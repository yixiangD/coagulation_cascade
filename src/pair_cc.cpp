/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include "pair_cc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "common.h"

using namespace LAMMPS_NS;
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON    1.0e-10
#define EPSILON_SQ 2.0e-20

/* ---------------------------------------------------------------------- */

PairCC::PairCC(LAMMPS *lmp) : Pair(lmp)
{
	random = NULL;
}

/* ---------------------------------------------------------------------- */

PairCC::~PairCC()
{
	if (allocated) {
		memory->destroy_2d_int_array(setflag);
		memory->destroy_2d_double_array(cutsq);

		memory->destroy_2d_double_array(cut);
		memory->destroy_2d_double_array(cutinv);
		memory->destroy_2d_double_array(cutc);

		memory->destroy_2d_double_array(a0);
		memory->destroy_2d_double_array(gamma);
		memory->destroy_2d_double_array(sigma);
		memory->destroy_2d_double_array(s1);
		memory->destroy_3d_double_array(kC);
		memory->destroy_3d_double_array(kappa);
		memory->destroy_3d_double_array(s2);
	}

	if (random) delete random;
}

/* ---------------------------------------------------------------------- */
void PairCC::compute(int eflag, int vflag)
{
	if (eflag || vflag) ev_setup(eflag,vflag);
	else evflag = vflag_fdotr = 0;

	double **x = atom->x;
	double **v = atom->v;
	double **f = atom->f;
	double **T = atom->T;
	double **Q = atom->Q;
	int  *type = atom->type;
	double *mass = atom->mass;
	int nlocal = atom->nlocal;
	int nall = nlocal + atom->nghost;
	double *special_lj = force->special_lj;
	int newton_pair = force->newton_pair;
	double dtinvsqrt = 1.0/sqrt(update->dt);

	int inum = list->inum;
	int *ilist = list->ilist;
	int *numneigh = list->numneigh;
	int **firstneigh = list->firstneigh;

	const double r_on = 0.01*cut_global;
	const double r_off= cut_global;

	// loop over neighbors of my atoms
	for (int p = 0; p < inum; p++)
	{
		int  i       = ilist     [p];
		int  itype   = type      [i];
		int  jnum    = numneigh  [i];
		int *jlist   = firstneigh[i];
		double xtmp  = x[i][0];
		double ytmp  = x[i][1];
		double ztmp  = x[i][2];
		double vxtmp = v[i][0];
		double vytmp = v[i][1];
		double vztmp = v[i][2];

		for (int q = 0; q < jnum; q++)
		{
			int j = jlist[q];
			double factor_dpd;
			if (j < nall) factor_dpd = 1.0;
			else
			{
				factor_dpd = special_lj[j/nall];
				j %= nall;
			}

			int jtype = type[j];
			double dx = xtmp - x[j][0];
			double dy = ytmp - x[j][1];
			double dz = ztmp - x[j][2];
			double rsq  = dx*dx + dy*dy + dz*dz;
			double SigmaIJ = sigma[itype][jtype];
			double GammaIJ = gamma[itype][jtype];
			double sf = s1[itype][jtype];

			if ( rsq < cutsq[itype][jtype])
			{
				double r    = sqrt(rsq);
				double rinv = 1.0 / r ;
				
				double wC, wR, wD;
		    if(r >= r_off ) { wC = wR = wD = 0.0;}
    		else if(r < r_on){ wC = wR = wD = 1.0;}
    		else
    		{
        	double r_rc = r * cutinv[itype][jtype];
        	wC = 1.0 - r_rc;
        	wC = MAX(0.0,wC);
        	wR = pow(wC, 0.5*sf);
        	wD = wR * wR;
    		}

				// conservative force
				double fC, eC;
				fC=a0[itype][jtype]*wC;
				eC=-0.5*a0[itype][jtype]*cut[itype][jtype] * wC *wC;

				// dissipative force
				double dvx = vxtmp - v[j][0];
				double dvy = vytmp - v[j][1];
				double dvz = vztmp - v[j][2];
				double dot = ( dx*dvx + dy*dvy + dz*dvz ) * rinv;
				double fD  = -GammaIJ * wD * dot ;

				// random force
				double RN = random->gaussian();
				RN = MAX(-4.0,MIN(RN,4.0));
				double fR = SigmaIJ * wR * dtinvsqrt * RN ;
				
				// sum it up
				double fpair = ( fC + fD + fR ) * factor_dpd * rinv ;

				f[i][0] += dx * fpair ;
				f[i][1] += dy * fpair ;
				f[i][2] += dz * fpair ;

				// chemical concentration transport
				if( r < cutc[itype][jtype])
    		{
					for(int k = 0; k < CTYPES; ++k)
					{
						double sC   = s2[itype][jtype][k];
			      double wCR = 1.0 - r/cutc[itype][jtype];
      			wCR = MAX(0.0,wCR);
			      wCR = pow(wCR, 0.5*sC);
					
						double dQc  = -kappa[itype][jtype][k] * wCR*wCR * ( T[i][k] - T[j][k] ); 	// q_cond

						if (atom->x[j][dir] >= outlet_loc || atom->x[j][dir] < domain->boxlo[dir]) //Alireza: zero-flux outflow BC 
							if (dir == 0) dQc -= dQc*(dx*dx)/rsq;
							else if (dir == 1) dQc -= dQc*(dy*dy)/rsq;
							else dQc -= dQc*(dz*dz)/rsq;

						Q[i][k] += ( dQc );
						if (newton_pair || j < nlocal)
							Q[j][k]  -= ( dQc );
					}
				}
				if (newton_pair || j < nlocal) {
					f[j][0] -= dx * fpair ;
					f[j][1] -= dy * fpair ;
					f[j][2] -= dz * fpair ;					
				}

				if (evflag) ev_tally(i,j,nlocal,newton_pair,eC,0.0,fpair,dx,dy,dz);
			}
		}
	}

	if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
	 allocate all arrays
------------------------------------------------------------------------- */

void PairCC::allocate()
{
	allocated = 1;
	int n = atom->ntypes;

	setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			setflag[i][j] = 0;

	cutsq   = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");
	cut     = memory->create_2d_double_array(n+1,n+1,"pair:cut");
	cutinv  = memory->create_2d_double_array(n+1,n+1,"pair:cutinv");
	cutc    = memory->create_2d_double_array(n+1,n+1,"pair:cutc");
	a0      = memory->create_2d_double_array(n+1,n+1,"pair:a0");
	gamma   = memory->create_2d_double_array(n+1,n+1,"pair:gamma");
	sigma   = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
	s1	= memory->create_2d_double_array(n+1,n+1,"pair:s1");
	kC	= memory->create_3d_double_array(n+1,n+1,CTYPES,"pair:kC");
	kappa   = memory->create_3d_double_array(n+1,n+1,CTYPES,"pair:kappa");
	s2	= memory->create_3d_double_array(n+1,n+1,CTYPES,"pair:s2");
}

/* ----------------------------------------------------------------------
	 global settings
------------------------------------------------------------------------- */

void PairCC::settings(int narg, char **arg)
{
	if (narg != 4) error->all("Illegal pair_style command");

	cut_global = atof(arg[0]);
	seed = atoi(arg[1]);
  dir = atoi(arg[2]);	//Alireza: orient the geometry such that the outlet is in x, y or z dir
  if (dir < 0 || dir > 2) error->all("PairCC: Outflow direction not valid");
	outlet_loc = atof(arg[3]);	//Alireza: exact location of outlet

	// initialize Marsaglia RNG with processor-unique seed

	if (seed <= 0) error->all("Illegal pair_style command");
	if (random) delete random;
	
	seed += int(time(0));
  random = new RanMars(lmp, (seed+comm->me) % 900000000 );
//	random = new RanMars(lmp,seed + comm->me);

	// reset cutoffs that have been explicitly set
	if (allocated)
	{
		int i,j;
		for (i = 1; i <= atom->ntypes; i++)
			for (j = i+1; j <= atom->ntypes; j++)
				if (setflag[i][j])
				{
					cut[i][j] = cut_global;
					cutinv[i][j] = 1.0 / cut[i][j] ;
				}
	}
}

/* ----------------------------------------------------------------------
	 set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCC::coeff(int narg, char **arg)
{
	if (narg != 7+1+3*CTYPES ) {
		error->all("Incorrect args for pair coefficients");
	}
	if (!allocated) allocate();

	int ilo,ihi,jlo,jhi;
	force->bounds(arg[0],atom->ntypes,ilo,ihi);
	force->bounds(arg[1],atom->ntypes,jlo,jhi);

	double a0_one 	   = atof(arg[2]);
	double gamma_one   = atof(arg[3]);
	double sigma_one   = atof(arg[4]);
	double s1_one      = atof(arg[5]);
	double cut_one     = atof(arg[6]);
	double cut_two     = atof(arg[7]);
	double kC_one[CTYPES], kappa_one[CTYPES], s2_one[CTYPES];
	for(int k=0; k<CTYPES; k++)
	{
    	kC_one[k]      = atof(arg[8+3*k]);
	    kappa_one[k]   = atof(arg[9+3*k]);
	    s2_one[k]      = atof(arg[10+3*k]);
	}

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
		for (int j = MAX(jlo,i); j <= jhi; j++) {
			a0     [i][j]    = a0_one;
			gamma  [i][j]    = gamma_one;
			sigma  [i][j]    = sigma_one;
			s1     [i][j]    = s1_one;
			cut    [i][j]    = cut_one;
			cutc   [i][j]    = cut_two;
			for(int k=0; k<CTYPES; k++)
			{
				kC [i][j][k] = kC_one[k];
			 	kappa [i][j][k] = kappa_one[k];
	    	s2 [i][j][k] = s2_one[k];
			}
	
			cutinv [i][j]    = 1.0 / cut_one;
			setflag[i][j]    = 1;
			count++;
		}
	}

	if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
	 init specific to this pair style
------------------------------------------------------------------------- */

void PairCC::init_style()
{
	if (atom->avec->ghost_velocity == 0)
		error->all("Pair CC requires ghost atoms store velocity");

	// if newton off, forces between atoms ij will be double computed
	// using different random numbers

	if (force->newton_pair == 0 && comm->me == 0) error->warning(
			"Pair CC needs newton pair on for momentum and heat conservation");

	int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
	 init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCC::init_one(int i, int j)
{
	if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

	cut[j][i]     = cut[i][j];
	cutinv[j][i]  = cutinv[i][j];
	cutc[j][i]    = cutc[i][j];
	a0[j][i]      = a0[i][j];
	gamma[j][i]   = gamma[i][j];
	sigma[j][i]   = sigma[i][j];
	s1[j][i]      = s1[i][j];
	for(int k=0; k<CTYPES; k++)
	{
		kC[j][i][k]      = kC[i][j][k];
		kappa[j][i][k]   = kappa[i][j][k];
		s2[j][i][k]      = s2[i][j][k];
	}

	return cut[i][j];
}

/* ----------------------------------------------------------------------
	 proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCC::write_restart(FILE *fp)
{
	write_restart_settings(fp);

	int i,j;
	for (i = 1; i <= atom->ntypes; i++)
	{
		for (j = i; j <= atom->ntypes; j++) {
			fwrite(&setflag[i][j],sizeof(int),1,fp);
			if (setflag[i][j]) {
				fwrite(&a0     [i][j],sizeof(double),1,fp);
				fwrite(&gamma  [i][j],sizeof(double),1,fp);
				fwrite(&sigma  [i][j],sizeof(double),1,fp);
				fwrite(&s1     [i][j],sizeof(double),1,fp);
				fwrite(&cut    [i][j],sizeof(double),1,fp);
				fwrite(&cutc   [i][j],sizeof(double),1,fp);
				for(int k=0; k<CTYPES; k++)
				{
					fwrite(&kC   [i][j][k],sizeof(double),1,fp);
					fwrite(&kappa[i][j][k],sizeof(double),1,fp);
					fwrite(&s2   [i][j][k],sizeof(double),1,fp);
				}
			}
		}
	}
}

/* ----------------------------------------------------------------------
	 proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCC::read_restart(FILE *fp)
{
	read_restart_settings(fp);

	allocate();

	int i,j;
	int me = comm->me;
	for (i = 1; i <= atom->ntypes; i++) {
		for (j = i; j <= atom->ntypes; j++) {
			if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
			MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
			if (setflag[i][j]) {
				if (me == 0) {
					fread(&a0     [i][j],sizeof(double),1,fp);
					fread(&gamma  [i][j],sizeof(double),1,fp);
					fread(&sigma  [i][j],sizeof(double),1,fp);
					fread(&s1     [i][j],sizeof(double),1,fp);
					fread(&cut    [i][j],sizeof(double),1,fp);
					fread(&cutc   [i][j],sizeof(double),1,fp);
					for(int k=0; k<CTYPES; k++)
					{
						fread(&kC   [i][j][k],sizeof(double),1,fp);
						fread(&kappa[i][j][k],sizeof(double),1,fp);
						fread(&s2   [i][j][k],sizeof(double),1,fp);
					}
				}
				MPI_Bcast(&a0     [i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&gamma  [i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&sigma  [i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&s1     [i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&cut    [i][j],1,MPI_DOUBLE,0,world);
				MPI_Bcast(&cutc   [i][j],1,MPI_DOUBLE,0,world);
				for(int k=0; k<CTYPES; k++)
				{
					MPI_Bcast(&kC   [i][j][k],1,MPI_DOUBLE,0,world);
					MPI_Bcast(&kappa[i][j][k],1,MPI_DOUBLE,0,world);
					MPI_Bcast(&s2   [i][j][k],1,MPI_DOUBLE,0,world);
				}
				cutinv[i][j] = 1.0 / cut[i][j];
			}
		}
	}
}

/* ----------------------------------------------------------------------
	 proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCC::write_restart_settings(FILE *fp)
{
	fwrite(&cut_global,sizeof(double),1,fp);
	fwrite(&seed,sizeof(int),1,fp);
	fwrite(&dir,sizeof(int),1,fp);
	fwrite(&outlet_loc,sizeof(double),1,fp);
	fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
	 proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCC::read_restart_settings(FILE *fp)
{
	if (comm->me == 0) {
		fread(&cut_global,sizeof(double),1,fp);
		fread(&seed,sizeof(int),1,fp);
		fread(&dir,sizeof(int),1,fp);
		fread(&outlet_loc,sizeof(double),1,fp);
		fread(&mix_flag,sizeof(int),1,fp);
	}
	MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
	MPI_Bcast(&seed,1,MPI_INT,0,world);
	MPI_Bcast(&dir,1,MPI_INT,0,world);
	MPI_Bcast(&outlet_loc,1,MPI_DOUBLE,0,world);
	MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

	// initialize Marsaglia RNG with processor-unique seed
	// same seed that pair_style command initially specified

	if (random) delete random;
  random = new RanMars(lmp, (seed+comm->me) % 900000000 );
}

/* ---------------------------------------------------------------------- */

double PairCC::single(int i, int j, int itype, int jtype, double rsq,
					 double factor_coul, double factor_dpd, double &fforce)
{
	double r,rinv,wr,phi;

	return 0.0;
}
