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

/* ----------------------------------------------------------------------
   Contributing author: Alireza Yazdani
------------------------------------------------------------------------- */
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_dpd_misc.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "common.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON 1.0e-10

/* ---------------------------------------------------------------------- */

PairDPDMisc::PairDPDMisc(LAMMPS *lmp) : Pair(lmp)
{
	single_enable = 0;
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPDMisc::~PairDPDMisc()
{
  if (allocated){
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);

    memory->destroy_2d_int_array(force_type);
    memory->destroy_2d_double_array(a0);
    memory->destroy_2d_double_array(gamma);
    memory->destroy_2d_double_array(sigma);
    memory->destroy_2d_double_array(cut_dpd);
    memory->destroy_2d_double_array(weight_exp);
    memory->destroy_2d_double_array(de);
    memory->destroy_2d_double_array(beta);
    memory->destroy_2d_double_array(r0m);
    memory->destroy_2d_double_array(morse1);
    memory->destroy_2d_double_array(cut_morse);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma_lj);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
    memory->destroy_2d_double_array(cut_lj);
    memory->destroy_2d_int_array(spec);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDMisc::compute(int eflag, int vflag)
{
  int i,j,k,l,m,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wr,randnum,factor_dpd,r2inv,r6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;
	int *neighs;
  double ff[6];
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  double shift = sqrt(3.0);
	double r_ce, aa, bb, dexp, dr;
	double dtv = update->dt;

	int *molecule = atom->molecule;
	int imol, jmol;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
		if (atom->molecule_flag) imol = molecule[i];
    
    if(stress_ind == 2)
      for(k = 0; k < n_stress_tot; ++k)
        modify->fix[stress_list[k]]->virial1(i);
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_dpd = 1.0;
      else {
				factor_dpd = special_lj[j/nall];
				j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
				r = sqrt(rsq);
				if (r < EPSILON) continue;     // r can be 0.0 in DPD systems

				if (force_type[itype][jtype] == 1){	// soft DPD conservative force + hard exponential force
					rinv = 1.0/r;
					delvx = vxtmp - v[j][0];
					delvy = vytmp - v[j][1];
					delvz = vztmp - v[j][2];
					dot = delx*delvx + dely*delvy + delz*delvz;
					wr = 1.0 - r/cut_dpd[itype][jtype];
      	  wd = pow(wr,weight_exp[itype][jtype]);
       		randnum = shift*(2.0*random->uniform()-1.0);

		// include hard/soft exponential potential between colloid-solvent, colloid-colloid
  	      if ((itype == 1 && jtype == 2 && atom->q[j] == 1.0) || (itype == 2 && jtype == 1 && atom->q[i] == 1.0)){
						aa = 1.0; bb = -60.0; r_ce = 1.0;
      	    if (r <= r_ce)
        	    wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
						else wr = 0.0;
        	}else if ((itype == 1 && jtype == 2 && atom->q[j] != 1.0) || (itype == 2 && jtype == 1 && atom->q[i] != 1.0)){
          	aa = 1.0; bb = -60.0; r_ce = 0.5;
          	if (r <= r_ce){
            	wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
							wd = pow((1.0-r/r_ce),weight_exp[itype][jtype]);
          	}else { wr = 0.0; wd = 0.0; }
        	}

        	if ((itype == 2 && jtype == 2) && (atom->q[i] >= 1.0 && atom->q[j] >= 1.0)){
						aa = 1.0; bb = -200.0; r_ce = 1.0;
          	if (r < r_ce)
            	wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
						else wr = 0.0;
        	}else if ((itype == 2 && jtype == 2 && atom->q[i] < 1.0) || (itype == 2 && jtype == 2 && atom->q[j] < 1.0)){
          	aa = 5.0; bb = -200.0; r_ce = 1.0;
          	if (r < r_ce)
            	wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
          	else wr = 0.0;
					}

#if 1
        if ((itype == 2 && jtype == 3 && atom->q[i] >= 1.0) || (itype == 3 && jtype == 2 && atom->q[j] >= 1.0)){
          aa = 1.0; bb = -200.0; r_ce = 0.5;
          if (r < r_ce)
            wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
          else wr = 0.0;
        }else if ((itype == 2 && jtype == 3 && atom->q[i] < 1.0) || (itype == 3 && jtype == 2 && atom->q[j] < 1.0))
          wr = 0.0;
#else
        if ((itype == 2 && jtype == 4 && atom->q[i] >= 1.0) || (itype == 4 && jtype == 2 && atom->q[j] >= 1.0)){
          aa = 1.0; bb = -200.0; r_ce = 0.5;
          if (r < r_ce)
            wr = aa/(1.0-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
          else wr = 0.0;
        }else if ((itype == 2 && jtype == 4 && atom->q[i] < 1.0) || (itype == 4 && jtype == 2 && atom->q[j] < 1.0))
          wr = 0.0;
#endif

		// conservative force = a0 * wr
		// drag force = -gamma * wd^2 * (delx dot delv) / r
		// random force = sigma * wd * rnd * dtinvsqrt;

					fpair = a0[itype][jtype]*wr;
					fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
					fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
					fpair *= factor_dpd*rinv;

					if (eflag) {
	  				evdwl = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype]*wr*wr;
	  				evdwl *= factor_dpd;
					}

		// Morse repulsive force to avoid cell overlap
				} else if (force_type[itype][jtype] == 2){
					if (atom->molecule_flag) jmol = molecule[j];
					if (imol && jmol && imol != jmol){
						double relax = 1.0;
						if ( ActState(imol, jmol) ) relax = 0.0;
					  dr = r - r0m[itype][jtype];
					  dexp = exp(-beta[itype][jtype] * dr);
					  fpair = relax * morse1[itype][jtype] * (dexp*dexp - dexp) / r;

					  if (eflag) evdwl = relax * de[itype][jtype] * (dexp*dexp - 2.0*dexp); 
					}

		// lj repulsive force to avoid cell overlap (Morse is prefered)
  			} else if (force_type[itype][jtype] == 3){
						if (atom->molecule_flag) jmol = molecule[j];
						if (imol && jmol && imol != jmol){
  	        	r2inv = 1.0/rsq;
    					r6inv = r2inv*r2inv*r2inv;
          		fpair = r6inv*(lj1[itype][jtype]*r6inv - lj2[itype][jtype])*r2inv;

					  if (eflag) evdwl = r6inv*(lj1[itype][jtype]*r6inv/12.0 - lj2[itype][jtype]/6.0);
						}

		// locally variable viscosity based on the local concentration
  			} else if (force_type[itype][jtype] == 4){
						int cc = spec[itype][jtype] - 1;
						double concen, factor_cc;
          	rinv = 1.0/r;
	          delvx = vxtmp - v[j][0];
  	        delvy = vytmp - v[j][1];
    	      delvz = vztmp - v[j][2];
      	    dot = delx*delvx + dely*delvy + delz*delvz;
        	  wr = 1.0 - r/cut_dpd[itype][jtype];
          	wd = pow(wr,weight_exp[itype][jtype]);
          	randnum = shift*(2.0*random->uniform()-1.0);

						concen = (atom->T[i][cc] + atom->T[j][cc]) / 2.0;
						if (concen == 0.0)
							factor_cc = 0.0;
						else
							factor_cc = fib_conv * 1.0 / (8e-12 * pow(concen/fib_thr, fib_pow) / (lscale*lscale));

		        fpair = a0[itype][jtype]*wr;
    		    fpair -= factor_cc*gamma[itype][jtype]*wd*wd*dot*rinv;
        		fpair += sqrt(factor_cc)*sigma[itype][jtype]*wd*randnum*dtinvsqrt;
        		fpair *= factor_dpd*rinv;

						if (eflag) {
		  				evdwl = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype]*wr*wr;
	  					evdwl *= factor_dpd;
						}

  			}	else
          	error->all("PairDPDMisc: Incorrect Force Type!!");

				f[i][0] += delx*fpair;
				f[i][1] += dely*fpair;
				f[i][2] += delz*fpair;
				if (newton_pair || j < nlocal) {
	  			f[j][0] -= delx*fpair;
	  			f[j][1] -= dely*fpair;
	  			f[j][2] -= delz*fpair;
				}

  			if(stress_ind == 2){
    			ff[0] = delx*delx*fpair;
	  			ff[1] = dely*dely*fpair;
	  			ff[2] = delz*delz*fpair;
	  			ff[3] = delx*dely*fpair;
	  			ff[4] = delx*delz*fpair;
	  			ff[5] = dely*delz*fpair;
	    		for(k = 0; k < n_stress_tot; ++k)
  	    		modify->fix[stress_list[k]]->virial2(j,ff);
  			}

				if (evflag) ev_tally(i,j,nlocal,newton_pair,
			     evdwl,0.0,fpair,delx,dely,delz);
      }	
    }
  }

  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairDPDMisc::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");
  force_type = memory->create_2d_int_array(n+1,n+1,"pair:force_type");
  a0 = memory->create_2d_double_array(n+1,n+1,"pair:a0");
  gamma = memory->create_2d_double_array(n+1,n+1,"pair:gamma");
  sigma = memory->create_2d_double_array(n+1,n+1,"pair:sigma");
  cut_dpd = memory->create_2d_double_array(n+1,n+1,"pair:cut_dpd");
  weight_exp = memory->create_2d_double_array(n+1,n+1,"pair:weight_exp");
  de = memory->create_2d_double_array(n+1,n+1,"pair:de");
  beta = memory->create_2d_double_array(n+1,n+1,"pair:beta");
  r0m = memory->create_2d_double_array(n+1,n+1,"pair:r0m");
  morse1 = memory->create_2d_double_array(n+1,n+1,"pair:morse1");
  cut_morse = memory->create_2d_double_array(n+1,n+1,"pair:cut_morse");
  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma_lj = memory->create_2d_double_array(n+1,n+1,"pair:sigma_lj");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
  cut_lj = memory->create_2d_double_array(n+1,n+1,"pair:cut_lj");
  spec = memory->create_2d_int_array(n+1,n+1,"pair:spec");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairDPDMisc::settings(int narg, char **arg)
{
  if (narg != 2) error->all("Illegal pair_style command");

  cut_global = atof(arg[0]);
  seed = atoi(arg[1]);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all("Illegal pair_style command");
  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDMisc::coeff(int narg, char **arg)
{
  int count, spec_one;
  double a0_one, gamma_one, sigma_one, cut_dpd_one, weight_exp_one;
  double de_one,r0m_one,beta_one,cut_morse_one,epsilon_one,sigma_lj_one,cut_lj_one;

  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  if (atoi(arg[2]) == 1) {
    if (narg != 8) error->all("Incorrect args in pair_coeff command");
    a0_one = atof(arg[3]);
    gamma_one = atof(arg[4]);
    sigma_one = atof(arg[5]);
    weight_exp_one = atof(arg[6]);
    cut_dpd_one = atof(arg[7]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 1;
        a0[i][j] = a0_one;
        gamma[i][j] = gamma_one;
        sigma[i][j] = sigma_one;
        weight_exp[i][j] = weight_exp_one;
        cut_dpd[i][j] = cut_dpd_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
  }
  else if (atoi(arg[2]) == 2) {
    if (narg != 7) error->all("Incorrect args in pair_coeff command");
    de_one = atof(arg[3]);
    beta_one = atof(arg[4]);
    r0m_one = atof(arg[5]);
    cut_morse_one = atof(arg[6]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 2;
        de[i][j] = de_one;
        beta[i][j] = beta_one;
        r0m[i][j] = r0m_one;
        cut_morse[i][j] = cut_morse_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
  }
  else if (atoi(arg[2]) == 3) {
    if (narg != 6) error->all("Incorrect args in pair_coeff command");
    epsilon_one = atof(arg[3]);
    sigma_lj_one = atof(arg[4]);
    cut_lj_one = atof(arg[5]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 3;
        epsilon[i][j] = epsilon_one;
        sigma_lj[i][j] = sigma_lj_one;
        cut_lj[i][j] = cut_lj_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
  }
  else if (atoi(arg[2]) == 4) {
    if (narg != 9) error->all("Incorrect args in pair_coeff command");
    a0_one = atof(arg[3]);
    gamma_one = atof(arg[4]);
    sigma_one = atof(arg[5]);
    weight_exp_one = atof(arg[6]);
    cut_dpd_one = atof(arg[7]);
    spec_one = atoi(arg[8]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 4;
        a0[i][j] = a0_one;
        gamma[i][j] = gamma_one;
        sigma[i][j] = sigma_one;
        weight_exp[i][j] = weight_exp_one;
        cut_dpd[i][j] = cut_dpd_one;
				spec[i][j] = spec_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
	}
  else
    error->all("Incorrect pair definition in pair_coeff command");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDMisc::init_style()
{
  if (atom->avec->ghost_velocity == 0)
    error->all("Pair dpd requires ghost atoms store velocity");

  // if newton off, forces between atoms ij will be double computed
  // using different random numbers

  if (force->newton_pair == 0 && comm->me == 0) error->warning(
      "Pair dpd needs newton pair on for momentum conservation");

  int irequest = neighbor->request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairDPDMisc::init_one(int i, int j)
{
	double value;
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

  a0[j][i] = a0[i][j];
  gamma[j][i] = gamma[i][j];
  sigma[j][i] = sigma[i][j];
  weight_exp[j][i] = weight_exp[i][j];
  cut_dpd[j][i] = cut_dpd[i][j];
  de[j][i] = de[i][j];
  beta[j][i] = beta[i][j];
  r0m[j][i] = r0m[i][j];
  morse1[i][j] = morse1[j][i] = 2.0*de[i][j]*beta[i][j];
  cut_morse[j][i] = cut_morse[i][j];
  epsilon[j][i] = epsilon[i][j];
  sigma_lj[j][i] = sigma_lj[i][j];
  lj1[i][j] = lj1[j][i] = 48.0 * epsilon[i][j] * pow(sigma_lj[i][j],12.0);
  lj2[i][j] = lj2[j][i] = 24.0 * epsilon[i][j] * pow(sigma_lj[i][j],6.0);
  cut_lj[j][i] = cut_lj[i][j];
  force_type[j][i] = force_type[i][j];
	spec[j][i] = spec[i][j];
  if (force_type[i][j] == 1)
		value = cut_dpd[i][j];
  else if (force_type[i][j] == 2)
		value = cut_morse[i][j];
	else if (force_type[i][j] == 3)
		value = cut_lj[i][j];
	else if (force_type[i][j] == 4)
		value = cut_dpd[i][j];
     
	return value;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDMisc::write_restart(FILE *fp)
{
  write_restart_settings(fp);
  int i,j;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&force_type[i][j],sizeof(int),1,fp);
        if (force_type[i][j] == 1){
    			fwrite(&a0[i][j],sizeof(double),1,fp);
          fwrite(&gamma[i][j],sizeof(double),1,fp);
          fwrite(&sigma[i][j],sizeof(double),1,fp);
          fwrite(&weight_exp[i][j],sizeof(double),1,fp);
          fwrite(&cut_dpd[i][j],sizeof(double),1,fp);
  			}
        else if (force_type[i][j] == 2){
          fwrite(&de[i][j],sizeof(double),1,fp);
          fwrite(&beta[i][j],sizeof(double),1,fp);
          fwrite(&r0m[i][j],sizeof(double),1,fp);
          fwrite(&cut_morse[i][j],sizeof(double),1,fp);
  			}
        else if (force_type[i][j] == 3){
          fwrite(&epsilon[i][j],sizeof(double),1,fp);
          fwrite(&sigma_lj[i][j],sizeof(double),1,fp);
          fwrite(&cut_lj[i][j],sizeof(double),1,fp);
  			}
				else{
	        fwrite(&a0[i][j],sizeof(double),1,fp);
          fwrite(&gamma[i][j],sizeof(double),1,fp);
          fwrite(&sigma[i][j],sizeof(double),1,fp);
          fwrite(&weight_exp[i][j],sizeof(double),1,fp);
          fwrite(&cut_dpd[i][j],sizeof(double),1,fp);
          fwrite(&spec[i][j],sizeof(int),1,fp);
				}
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDMisc::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
  			if (me == 0) {
          fread(&force_type[i][j],sizeof(int),1,fp);
          if (force_type[i][j] == 1){
      			fread(&a0[i][j],sizeof(double),1,fp);
            fread(&gamma[i][j],sizeof(double),1,fp);
            fread(&sigma[i][j],sizeof(double),1,fp);
            fread(&weight_exp[i][j],sizeof(double),1,fp);
            fread(&cut_dpd[i][j],sizeof(double),1,fp);
    			}
          else if (force_type[i][j] == 2){
            fread(&de[i][j],sizeof(double),1,fp);
            fread(&beta[i][j],sizeof(double),1,fp);
            fread(&r0m[i][j],sizeof(double),1,fp);
            fread(&cut_morse[i][j],sizeof(double),1,fp);
    			}
          else if (force_type[i][j] == 3){
            fread(&epsilon[i][j],sizeof(double),1,fp);
            fread(&sigma_lj[i][j],sizeof(double),1,fp);
            fread(&cut_lj[i][j],sizeof(double),1,fp);
    			}
					else{
            fread(&a0[i][j],sizeof(double),1,fp);
            fread(&gamma[i][j],sizeof(double),1,fp);
            fread(&sigma[i][j],sizeof(double),1,fp);
            fread(&weight_exp[i][j],sizeof(double),1,fp);
            fread(&cut_dpd[i][j],sizeof(double),1,fp);
            fread(&spec[i][j],sizeof(int),1,fp);
					}
  			}
        MPI_Bcast(&force_type[i][j],1,MPI_INT,0,world);
        if (force_type[i][j] == 1){
    			MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&weight_exp[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&cut_dpd[i][j],1,MPI_DOUBLE,0,world);
  			}
        else if (force_type[i][j] == 2){
          MPI_Bcast(&de[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&beta[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&r0m[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_morse[i][j],1,MPI_DOUBLE,0,world);
  			}
        else if (force_type[i][j] == 3){
          MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&sigma_lj[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
  			}
				else{
          MPI_Bcast(&a0[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&weight_exp[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_dpd[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&spec[i][j],1,MPI_INT,0,world);
				}
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDMisc::write_restart_settings(FILE *fp)
{
	fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDMisc::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
}

/*----------------------------------------------------------------------- */

int PairDPDMisc::ActState(int imol, int jmol)
{
  if (atom->molecular == 0) return 0;
	if (ifix == modify->nfix) return 0;
  if (molecule == NULL){
    for (ifix = 0; ifix < modify->nfix; ifix++)
      if (strcmp("molecule",modify->fix[ifix]->style) == 0) break;
    if (ifix != modify->nfix) molecule = dynamic_cast<FixMolecule *> (modify->fix[ifix]);
  }

	if (!molecule) return 0;

	if ( molecule->list_tact[imol] > 0.0 && molecule->list_tact[jmol] > 0.0 )
		return 1;
	else
		return 0;
}

/*----------------------------------------------------------------------- */
