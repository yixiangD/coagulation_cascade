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
   Contributing author: Alireza Yazdani/Dmitry Fedosov
------------------------------------------------------------------------- */
#include "mpi.h"
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_dpd_adh.h"
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

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EPSILON 1.0e-10
#define DELTA 10

/* ---------------------------------------------------------------------- */

PairDPDAdh::PairDPDAdh(LAMMPS *lmp) : Pair(lmp)
{
  random = NULL;
}

/* ---------------------------------------------------------------------- */

PairDPDAdh::~PairDPDAdh()
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
    memory->destroy_2d_double_array(kf0);
    memory->destroy_2d_double_array(sig);
    memory->destroy_2d_double_array(temp);
    memory->destroy_2d_double_array(r0);
    memory->destroy_2d_double_array(ks);
    memory->destroy_2d_double_array(cut_spring);
    memory->destroy_2d_double_array(de);
    memory->destroy_2d_double_array(r0m);
    memory->destroy_2d_double_array(beta);
    memory->destroy_2d_double_array(cut_morse);
    memory->destroy_2d_double_array(epsilon);
    memory->destroy_2d_double_array(sigma_lj);
    memory->destroy_2d_double_array(lj1);
    memory->destroy_2d_double_array(lj2);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void PairDPDAdh::compute(int eflag, int vflag)
{
  int i,j,k,l,m,ii,jj,inum,jnum,itype,jtype,offset;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  double rsq,r,rinv,dot,wd,wr,randnum,factor_dpd,kf,pp,r2inv,r6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;
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
	int *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  double shift = sqrt(3.0);
	double r_ce, aa, bb, dexp, dr;
  char warning[128];
	double dtv = update->dt;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms
	int createcount, ncreate = 0;
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

				if (force_type[itype][jtype] == 1){
					rinv = 1.0/r;
					delvx = vxtmp - v[j][0];
					delvy = vytmp - v[j][1];
					delvz = vztmp - v[j][2];
					dot = delx*delvx + dely*delvy + delz*delvz;
					wr = 1.0 - r/cut_dpd[itype][jtype];
      	  wd = pow(wr,weight_exp[itype][jtype]);
       		randnum = shift*(2.0*random->uniform()-1.0);

		// Alireza: include hard/soft exponential potential between colloid-solvent, colloid-colloid
  	      if ((itype == 1 && jtype == 2 && atom->q[j] == 1.0) || (itype == 2 && jtype == 1 && atom->q[i] == 1.0)){
						aa = 1.0; bb = -60.0; r_ce = 1.0;
      	    if (r <= r_ce)
        	    wr = aa/(1-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
						else wr = 0.0;
/*
          	if (r <= r_ce)
            	wd = pow((1.0-r/r_ce),weight_exp[itype][jtype]);
          	else wd = 0.0;
*/
        	} else if ((itype == 1 && jtype == 2 && atom->q[j] != 1.0) || (itype == 2 && jtype == 1 && atom->q[i] != 1.0)){
          	aa = 1.0; bb = -60.0; r_ce = 1.0;
          	if (r <= r_ce){
            	wr = aa/(1-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
							wd = pow((1.0-r/r_ce),weight_exp[itype][jtype]);
						} else {wr = 0.0; wd = 0.0;}
					}

        	if ((itype == 2 && jtype == 2) && (atom->q[i] >= 1.0 && atom->q[j] >= 1.0)){
						aa = 1.0; bb = -200.0; r_ce = 1.0;
          	if (r < r_ce)
            	wr = aa/(1-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
						else wr = 0.0;
        	}else if ((itype == 2 && jtype == 2 && atom->q[i] < 1.0) || (itype == 2 && jtype == 2 && atom->q[j] < 1.0)){
          	aa = 5.0; bb = -200.0; r_ce = 1.0;
          	if (r < r_ce)
            	wr = aa/(1-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
          	else wr = 0.0;
					}

					if ((itype == 2 && jtype == 3 && atom->q[i] >= 1.0) || (itype == 3 && jtype == 2 && atom->q[j] >= 1.0)){
          	aa = 1.0; bb = -200.0; r_ce = 0.5;
          	if (r < r_ce)
            	wr = aa/(1-exp(-r_ce*bb))*(exp(-r*bb)-exp(-r_ce*bb));
          	else wr = 0.0;
					}else if ((itype == 2 && jtype == 3 && atom->q[i] < 1.0) || (itype == 3 && jtype == 2 && atom->q[j] < 1.0))
						wr = 0.0;

		// conservative force = a0 * wr
		// drag force = -gamma * wd^2 * (delx dot delv) / r
		// random force = sigma * wd * rnd * dtinvsqrt;

					fpair = a0[itype][jtype]*wr;
					fpair -= gamma[itype][jtype]*wd*wd*dot*rinv;
					fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt;
					fpair *= factor_dpd*rinv;

		// Alireza: adhesive dynamics
				} else if (force_type[itype][jtype] == 2){
					if (atom->molecule[i] == atom->molecule[j]) continue;

          fpair = 0.0;
          delvx = vxtmp - v[j][0];
          delvy = vytmp - v[j][1];
          delvz = vztmp - v[j][2];
					double vs = sqrt(delvx*delvx+delvy*delvy+delvz*delvz);

          if (itype == ligand_type[itype][jtype]){
						if (atom->num_bond[i] > jmaxbond || atom->num_bond[j] > imaxbond) continue;
            kf = kf0[itype][jtype]*vs*exp((ks[itype][jtype]*fabs(r-r0[itype][jtype])*(sig[itype][jtype]-0.5*fabs(r-r0[itype][jtype])))/temp[itype][jtype]);	//Alireza: GPIb-vWF
            pp = 1.0 - exp(-kf*dtv);
            if (random->uniform() < pp){
		fprintf(stdout,"IM HERE (0), %d %d %d %d\n",atom->tag[i], atom->num_bond[i], atom->tag[j], atom->num_bond[j]);
							for (m = 0; m < Nreceptor; m++){
								if (atom->num_bond[i] > jmaxbond || atom->num_bond[j] > imaxbond) break;
                atom->bond_type[i][atom->num_bond[i]] = bond_type[itype][jtype];
                atom->bond_atom[i][atom->num_bond[i]] = atom->tag[j];
                atom->bond_length[i][atom->num_bond[i]] = r0[itype][jtype];
                atom->num_bond[i]++;

                atom->bond_type[j][atom->num_bond[j]] = bond_type[itype][jtype];
                atom->bond_atom[j][atom->num_bond[j]] = atom->tag[i];
                atom->bond_length[j][atom->num_bond[j]] = r0[itype][jtype];
                atom->num_bond[j]++;

								ncreate++;
                //fpair += -2.0*ks[itype][jtype]*(r-r0[itype][jtype])/r;
                neighbor->forced_reneigh = 1;
        			}
      			}
    			}
          else if (jtype == ligand_type[itype][jtype] && j < nlocal){
						if (atom->num_bond[i] > imaxbond || atom->num_bond[j] > jmaxbond) continue;
          	kf = kf0[itype][jtype]*vs*exp((ks[itype][jtype]*fabs(r-r0[itype][jtype])*(sig[itype][jtype]-0.5*fabs(r-r0[itype][jtype])))/temp[itype][jtype]);	//Alireza: GPIb-vWF
            pp = 1.0 - exp(-kf*dtv);
            if (random->uniform() < pp){
//		fprintf(stdout,"IM HERE (1), %d %d %d %d %d %d %d\n",Nreceptor, imaxbond, jmaxbond, atom->tag[i], atom->num_bond[i], atom->tag[j], atom->num_bond[j]);
							for (m = 0; m < Nreceptor; m++){
								if (atom->num_bond[i] > imaxbond || atom->num_bond[j] > jmaxbond) break;
                atom->bond_type[j][atom->num_bond[j]] = bond_type[itype][jtype];
                atom->bond_atom[j][atom->num_bond[j]] = atom->tag[i];
                atom->bond_length[j][atom->num_bond[j]] = r0[itype][jtype];
                atom->num_bond[j]++;

                atom->bond_type[i][atom->num_bond[i]] = bond_type[itype][jtype];
                atom->bond_atom[i][atom->num_bond[i]] = atom->tag[j];
                atom->bond_length[i][atom->num_bond[i]] = r0[itype][jtype];
                atom->num_bond[i]++;

								ncreate++;
                //fpair += -2.0*ks[itype][jtype]*(r-r0[itype][jtype])/r;
                neighbor->forced_reneigh = 1;
        			}
      			}
    			}
				} else if (force_type[itype][jtype] == 3){
          fpair = 2.0*de[itype][jtype]*beta[itype][jtype]*(exp(2.0*beta[itype][jtype]*(r0m[itype][jtype]-r))-exp(beta[itype][jtype]*(r0m[itype][jtype]-r)));
  			} else if (force_type[itype][jtype] == 4){
	          rinv = 1.0/r;
  	        r2inv = 1.0/rsq;
    				r6inv = r2inv*r2inv*r2inv;
      	    delvx = vxtmp - v[j][0];
        	  delvy = vytmp - v[j][1];
          	delvz = vztmp - v[j][2];
          	dot = delx*delvx + dely*delvy + delz*delvz;
          	wr = 1.0 - r/cut_dpd[itype][jtype];
          	wd = pow(wr,weight_exp[itype][jtype]);
          	randnum = shift*(2.0*random->uniform()-1.0);
          	fpair = r6inv*(lj1[itype][jtype]*r6inv - lj2[itype][jtype])*r2inv;
          	fpair -= gamma[itype][jtype]*wd*wd*dot*r2inv;
          	fpair += sigma[itype][jtype]*wd*randnum*dtinvsqrt*rinv;
          	fpair *= factor_dpd;
  			}	else
          	error->all("Incorrect Force Type!!");

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
					if (force_type[itype][jtype] == 2)
	          for(k = 0; k < n_stress_tot; ++k)
  	          modify->fix[stress_list[k]]->virial3(i,j,ff);
					else
	    			for(k = 0; k < n_stress_tot; ++k)
  	    			modify->fix[stress_list[k]]->virial2(j,ff);
  			}

				if (eflag) {
	  			evdwl = 0.5*a0[itype][jtype]*cut_dpd[itype][jtype]*wr*wr;
	  			evdwl *= factor_dpd;
				}

				if (evflag) ev_tally(i,j,nlocal,newton_pair,
			     evdwl,0.0,fpair,delx,dely,delz);
      }	
    }
  }

  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  atom->nbonds += createcount;

  l = 0;
  MPI_Allreduce(&neighbor->forced_reneigh,&l,1,MPI_INT,MPI_MAX,world);
  neighbor->forced_reneigh = l;

  if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairDPDAdh::allocate()
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
  kf0 = memory->create_2d_double_array(n+1,n+1,"pair:kf0");
  sig = memory->create_2d_double_array(n+1,n+1,"pair:sig");
  temp = memory->create_2d_double_array(n+1,n+1,"pair:temp");
  r0 = memory->create_2d_double_array(n+1,n+1,"pair:r0");
  ligand_type = memory->create_2d_int_array(n+1,n+1,"pair:ligand_type");
  bond_type = memory->create_2d_int_array(n+1,n+1,"pair:bond_type");
  ks = memory->create_2d_double_array(n+1,n+1,"pair:ks");
  cut_spring = memory->create_2d_double_array(n+1,n+1,"pair:cut_spring");
  de = memory->create_2d_double_array(n+1,n+1,"pair:de");
  r0m = memory->create_2d_double_array(n+1,n+1,"pair:r0m");
  beta = memory->create_2d_double_array(n+1,n+1,"pair:beta");
  cut_morse = memory->create_2d_double_array(n+1,n+1,"pair:cut_morse");
  epsilon = memory->create_2d_double_array(n+1,n+1,"pair:epsilon");
  sigma_lj = memory->create_2d_double_array(n+1,n+1,"pair:sigma_lj");
  lj1 = memory->create_2d_double_array(n+1,n+1,"pair:lj1");
  lj2 = memory->create_2d_double_array(n+1,n+1,"pair:lj2");
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairDPDAdh::settings(int narg, char **arg)
{
  short unsigned seed_h;
 
  if (narg != 1) error->all("Illegal pair_style command");

  seed = atoi(arg[0]);

  // initialize Marsaglia RNG with processor-unique seed

  if (seed <= 0) error->all("Illegal fix pair_style command");
  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
  seed_h = (comm->me+3)*seed;
  seed48(&seed_h);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairDPDAdh::coeff(int narg, char **arg)
{
  int count, ligand_one, bond_one, Nreceptor, imaxbond, jmaxbond;
  double a0_one, gamma_one, sigma_one, cut_dpd_one, weight_exp_one;
  double kf0_one, sig_one, r0_one, temp_one, ks_one, cut_spring_one;
  double de_one,r0m_one,beta_one,cut_morse_one,epsilon_one,sigma_lj_one;

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
    if (narg != 14) error->all("Incorrect args in pair_coeff command");
  	ligand_one = atoi(arg[3]);
  	bond_one = atoi(arg[4]);
    ks_one = atof(arg[5]);
    r0_one = atof(arg[6]);
    kf0_one = atof(arg[7]);
    sig_one = atof(arg[8]);
    temp_one = atof(arg[9]);
    cut_spring_one = atof(arg[10]);
    Nreceptor = atoi(arg[11]);
    imaxbond = atoi(arg[12]);
    jmaxbond = atoi(arg[13]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 2;
				ligand_type[i][j] = ligand_one;
				bond_type[i][j] = bond_one;
        ks[i][j] = ks_one;
        r0[i][j] = r0_one;
        kf0[i][j] = kf0_one;
        sig[i][j] = sig_one;
        temp[i][j] = temp_one;
        cut_spring[i][j] = cut_spring_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
  }
  else if (atoi(arg[2]) == 3) {
    if (narg != 7) error->all("Incorrect args in pair_coeff command");
    de_one = atof(arg[3]);
    r0m_one = atof(arg[4]);
    beta_one = atof(arg[5]);
    cut_morse_one = atof(arg[6]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 3;
        de[i][j] = de_one;
        r0m[i][j] = r0m_one;
        beta[i][j] = beta_one;
        cut_morse[i][j] = cut_morse_one;
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all("Incorrect args in pair_coeff command");
  }
  else if (atoi(arg[2]) == 4) {
    if (narg != 9) error->all("Incorrect args in pair_coeff command");
    epsilon_one = atof(arg[3]);
    sigma_lj_one = atof(arg[4]);
    gamma_one = atof(arg[5]);
    sigma_one = atof(arg[6]);
    weight_exp_one = atof(arg[7]);
    cut_dpd_one = atof(arg[8]);

    count = 0;
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo,i); j <= jhi; j++) {
        force_type[i][j] = 4;
        epsilon[i][j] = epsilon_one;
        sigma_lj[i][j] = sigma_lj_one;
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
  else
    error->all("Incorrect pair definition in pair_coeff command");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairDPDAdh::init_style()
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

double PairDPDAdh::init_one(int i, int j)
{
	double value;
  if (setflag[i][j] == 0) error->all("All pair coeffs are not set");

//  for (i = 1; i <= atom->ntypes; i++)
//    for (j = i+1; j <= atom->ntypes; j++)
//      if (setflag[i][j] == 0) error->all("All cross pair coeffs are not set");

  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      a0[j][i] = a0[i][j];
      gamma[j][i] = gamma[i][j];
      sigma[j][i] = sigma[i][j];
      weight_exp[j][i] = weight_exp[i][j];
      cut_dpd[j][i] = cut_dpd[i][j];
			ligand_type[j][i] = ligand_type[i][j];
			bond_type[j][i] = bond_type[i][j];
      ks[j][i] = ks[i][j];
      r0[j][i] = r0[i][j];
      kf0[j][i] = kf0[i][j];
      sig[j][i] = sig[i][j];
      temp[j][i] = temp[i][j];
      cut_spring[j][i] = cut_spring[i][j];
      de[j][i] = de[i][j];
      r0m[j][i] = r0m[i][j];
      beta[j][i] = beta[i][j];
      cut_morse[j][i] = cut_morse[i][j];
      force_type[j][i] = force_type[i][j];
      if (force_type[i][j] == 1 || force_type[i][j] == 4){
        cutsq[i][j] = cut_dpd[i][j]*cut_dpd[i][j];
				value = cut_dpd[i][j];
      } else if (force_type[i][j] == 2){
        cutsq[i][j] = cut_spring[i][j]*cut_spring[i][j];
				value = cut_spring[i][j];
      } else if (force_type[i][j] == 3){
        cutsq[i][j] = cut_morse[i][j]*cut_morse[i][j];
				value = cut_morse[i][j];
			}
      cutsq[j][i] = cutsq[i][j];
      epsilon[j][i] = epsilon[i][j];
      sigma_lj[j][i] = sigma_lj[i][j];
      lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma_lj[i][j],12.0);
      lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma_lj[i][j],6.0);
      lj1[j][i] = lj1[i][j];
      lj2[j][i] = lj2[i][j];
    }

	return value;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDAdh::write_restart(FILE *fp)
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
					fwrite(&ligand_type[i][j],sizeof(int),1,fp);
					fwrite(&bond_type[i][j],sizeof(int),1,fp);
          fwrite(&ks[i][j],sizeof(double),1,fp);
          fwrite(&r0[i][j],sizeof(double),1,fp);
          fwrite(&kf0[i][j],sizeof(double),1,fp);
          fwrite(&sig[i][j],sizeof(double),1,fp);
          fwrite(&temp[i][j],sizeof(double),1,fp);
          fwrite(&cut_spring[i][j],sizeof(double),1,fp);
					fwrite(&Nreceptor,sizeof(int),1,fp);
					fwrite(&imaxbond,sizeof(int),1,fp);
					fwrite(&jmaxbond,sizeof(int),1,fp);
  			}
        else if (force_type[i][j] == 3){
          fwrite(&de[i][j],sizeof(double),1,fp);
          fwrite(&r0m[i][j],sizeof(double),1,fp);
          fwrite(&beta[i][j],sizeof(double),1,fp);
          fwrite(&cut_morse[i][j],sizeof(double),1,fp);
  			}
        else{
          fwrite(&epsilon[i][j],sizeof(double),1,fp);
          fwrite(&sigma_lj[i][j],sizeof(double),1,fp);
          fwrite(&gamma[i][j],sizeof(double),1,fp);
          fwrite(&sigma[i][j],sizeof(double),1,fp);
          fwrite(&weight_exp[i][j],sizeof(double),1,fp);
          fwrite(&cut_dpd[i][j],sizeof(double),1,fp);
  			}
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDAdh::read_restart(FILE *fp)
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
						fread(&ligand_type[i][j],sizeof(int),1,fp);
						fread(&bond_type[i][j],sizeof(int),1,fp);
            fread(&ks[i][j],sizeof(double),1,fp);
            fread(&r0[i][j],sizeof(double),1,fp);
            fread(&kf0[i][j],sizeof(double),1,fp);
            fread(&sig[i][j],sizeof(double),1,fp);
            fread(&temp[i][j],sizeof(double),1,fp);
            fread(&cut_spring[i][j],sizeof(double),1,fp);
						fread(&Nreceptor,sizeof(int),1,fp);
						fread(&imaxbond,sizeof(int),1,fp);
						fread(&jmaxbond,sizeof(int),1,fp);
    			}
          else if (force_type[i][j] == 3){
            fread(&de[i][j],sizeof(double),1,fp);
            fread(&r0m[i][j],sizeof(double),1,fp);
            fread(&beta[i][j],sizeof(double),1,fp);
            fread(&cut_morse[i][j],sizeof(double),1,fp);
    			}
          else{
            fread(&epsilon[i][j],sizeof(double),1,fp);
            fread(&sigma_lj[i][j],sizeof(double),1,fp);
            fread(&gamma[i][j],sizeof(double),1,fp);
            fread(&sigma[i][j],sizeof(double),1,fp);
            fread(&weight_exp[i][j],sizeof(double),1,fp);
            fread(&cut_dpd[i][j],sizeof(double),1,fp);
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
        	MPI_Bcast(&ligand_type[i][j],1,MPI_INT,0,world);
        	MPI_Bcast(&bond_type[i][j],1,MPI_INT,0,world);
          MPI_Bcast(&ks[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&kf0[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&sig[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&temp[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_spring[i][j],1,MPI_DOUBLE,0,world);
        	MPI_Bcast(&Nreceptor,1,MPI_INT,0,world);
        	MPI_Bcast(&imaxbond,1,MPI_INT,0,world);
        	MPI_Bcast(&jmaxbond,1,MPI_INT,0,world);
  			}
        else if (force_type[i][j] == 3){
          MPI_Bcast(&de[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&r0m[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&beta[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&cut_morse[i][j],1,MPI_DOUBLE,0,world);
  			}
        else{
          MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
          MPI_Bcast(&sigma_lj[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&gamma[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
   			 	MPI_Bcast(&weight_exp[i][j],1,MPI_DOUBLE,0,world);
    			MPI_Bcast(&cut_dpd[i][j],1,MPI_DOUBLE,0,world);
  			}
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairDPDAdh::write_restart_settings(FILE *fp)
{
  fwrite(&seed,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairDPDAdh::read_restart_settings(FILE *fp)
{
  short unsigned seed_h;
 
  if (comm->me == 0) {
    fread(&seed,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&seed,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

  // initialize Marsaglia RNG with processor-unique seed
  // same seed that pair_style command initially specified

  if (random) delete random;
  random = new RanMars(lmp,seed + comm->me);
  seed_h = (comm->me+3)*seed;
  seed48(&seed_h);
}
