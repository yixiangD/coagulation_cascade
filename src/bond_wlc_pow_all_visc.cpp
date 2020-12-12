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
#include "stdlib.h"
#include "bond_wlc_pow_all_visc.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include "random_mars.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondWLC_POW_ALL_VISC::BondWLC_POW_ALL_VISC(LAMMPS *lmp) : Bond(lmp)
{
  random = new RanMars(lmp, 5684 + comm->me);
}

/* ---------------------------------------------------------------------- */

BondWLC_POW_ALL_VISC::~BondWLC_POW_ALL_VISC()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(temp);
    memory->sfree(r0);
    memory->sfree(mu_targ);
    memory->sfree(qp);
    memory->sfree(gamc);
    memory->sfree(gamt);
    memory->sfree(sigc);
    memory->sfree(sigt);
  }

	if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::compute(int eflag, int vflag)
{
  double ebond,fbond;
  int i1,i2,n,type,factor;
  double rr,rfactor,ra,rlogarg,kph,l0,lmax,mu,lambda;
  double dvx, dvy, dvz, vv;
  char warning[128];
  double ff[6];
  
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;  
  double **v = atom->v;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double *bondlist_length = neighbor->bondlist_length;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    l0 = bondlist_length[n];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    dvx = v[i1][0] - v[i2][0];
    dvy = v[i1][1] - v[i2][1];
    dvz = v[i1][2] - v[i2][2];
    domain->minimum_image(delx,dely,delz);

   	ra = sqrt(delx*delx + dely*dely + delz*delz);
   	lmax = l0*r0[type];
   	rr = 1.0/r0[type];  
   	kph = pow(l0,qp[type])*temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr);
   	mu = 0.25*sqrt(3.0)*(temp[type]*(-0.25/pow(1.0-rr,2) + 0.25 + 0.5*rr/pow(1.0-rr,3))/lmax/rr + kph*(qp[type]+1.0)/pow(l0,qp[type]+1.0));
   	lambda = mu/mu_targ[type];
   	kph = kph*mu_targ[type]/mu;
   	rr = ra/lmax; 
   	rlogarg = pow(ra,qp[type]+1.0);
   	vv = (delx*dvx + dely*dvy +  delz*dvz)/ra;
    
      if (rr >= 1.0) {
        sprintf(warning,"WLC bond too long: %d %d %d %d %d %d %g",
                update->ntimestep,atom->molecule[i1],atom->molecule[i2],atom->tag[i1],atom->tag[i2],type,rr);
        error->warning(warning);
      }

    generate_wrr();
    fbond = -temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr)/lambda/ra + kph/rlogarg + (sigc[type]*wrr[3] - gamc[type]*vv)/ra;

    // energy
    if (eflag) {
     	ebond += 0.25*temp[type]*lmax*(3.0*rr*rr-2.0*rr*rr*rr)/(1.0-rr)/lambda;
     	if (qp[type] == 1.0)
       	ebond -= kph*log(ra);
     	else
       	ebond += kph/(qp[type]-1.0)/pow(ra,qp[type]-1.0);
    }

   	// apply force to each of 2 atoms

   	if (newton_bond || i1 < nlocal) {
   	  f[i1][0] += delx*fbond - gamt[type]*dvx + sigt[type]*wrr[0]/ra;
      f[i1][1] += dely*fbond - gamt[type]*dvy + sigt[type]*wrr[1]/ra;
      f[i1][2] += delz*fbond - gamt[type]*dvz + sigt[type]*wrr[2]/ra;
   	}

   	if (newton_bond || i2 < nlocal) {
    	f[i2][0] -= delx*fbond - gamt[type]*dvx + sigt[type]*wrr[0]/ra;
      f[i2][1] -= dely*fbond - gamt[type]*dvy + sigt[type]*wrr[1]/ra;
      f[i2][2] -= delz*fbond - gamt[type]*dvz + sigt[type]*wrr[2]/ra;
    }
    
    // virial contribution

  	if (stress_ind == 2){
      ff[0] = delx*delx*fbond;
      ff[1] = dely*dely*fbond;
      ff[2] = delz*delz*fbond;
      ff[3] = delx*dely*fbond;
      ff[4] = delx*delz*fbond;
      ff[5] = dely*delz*fbond;
    
      for(int k = 0; k < n_stress_tot; ++k)
	      modify->fix[stress_list[k]]->virial3(i1, i2, ff);
  	}

  	if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
	}
}

/* ---------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  temp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:temp");
  r0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:r0");
  mu_targ = (double *) memory->smalloc((n+1)*sizeof(double),"bond:mu_targ");
  qp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:qp");
  gamc = (double *) memory->smalloc((n+1)*sizeof(double),"bond:gamc");
  gamt = (double *) memory->smalloc((n+1)*sizeof(double),"bond:gamt");
  sigc = (double *) memory->smalloc((n+1)*sizeof(double),"bond:sigc");
  sigt = (double *) memory->smalloc((n+1)*sizeof(double),"bond:sigt"); 
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::coeff(int narg, char **arg)
{
  if (narg != 7) error->all("Incorrect args in bond_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double temp_one = atof(arg[1]);
  double r0_one = atof(arg[2]);
  double mu_one = atof(arg[3]);
  double qp_one = atof(arg[4]);
  double gamc_one = atof(arg[5]);
  double gamt_one = atof(arg[6]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    temp[i] = temp_one;
    r0[i] = r0_one;
    mu_targ[i] = mu_one;
    qp[i] = qp_one;
    gamc[i] = gamc_one;
    gamt[i] = gamt_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args in bond_coeff command");
}
 
/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::init_style()
{
  double sdtt = sqrt(update->dt);

  for (int i = 1; i <= atom->nbondtypes; i++){
  //  if (setflag[i] == 0) error->all("All bond coeffs are not set");
    if (gamt[i] > 3.0*gamc[i]) error->one("Gamma_t > 3*Gamma_c");
    sigc[i] = sqrt(2.0*temp[i]*(3.0*gamc[i]-gamt[i]))/sdtt;
    sigt[i] = 2.0*sqrt(gamt[i]*temp[i])/sdtt;
  }
}

/* ---------------------------------------------------------------------- */

double BondWLC_POW_ALL_VISC::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::write_restart(FILE *fp)
{
  fwrite(&temp[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&mu_targ[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&qp[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamc[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&gamt[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&temp[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&mu_targ[1],sizeof(double),atom->nbondtypes,fp);
    fread(&qp[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gamc[1],sizeof(double),atom->nbondtypes,fp);
    fread(&gamt[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&temp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&mu_targ[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&qp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamc[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&gamt[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondWLC_POW_ALL_VISC::single(int type, double rsq, int i, int j)
{
  return 0.0;
}

/* ----------------------------------------------------------------------- */

void BondWLC_POW_ALL_VISC::generate_wrr()
{
  int i;
  double ww[3][3];
  double v1, v2, factor, ss;

  for (i=0; i<5; i++){
    ss = 100.0;
    while ( ss > 1.0 ){
      v1 = 2.0 * random->uniform() - 1.0;
      v2 = 2.0 * random->uniform() - 1.0;
      ss = v1*v1 + v2*v2;
    }
    factor = sqrt(-2.0 * log(ss)/ss);
    if (i < 3){
      ww[i][0] = factor*v1;
      ww[i][1] = factor*v2; 
    }
    else if (i == 3){
      ww[0][2] = factor*v1;
      ww[1][2] = factor*v2;
    }
    else
      ww[2][2] = factor*v1; 
  }
  wrr[3] = (ww[0][0]+ww[1][1]+ww[2][2])/3.0;
  wrr[0] = (ww[0][0]-wrr[3])*delx + 0.5*(ww[0][1]+ww[1][0])*dely + 0.5*(ww[0][2]+ww[2][0])*delz;
  wrr[1] = 0.5*(ww[1][0]+ww[0][1])*delx + (ww[1][1]-wrr[3])*dely + 0.5*(ww[1][2]+ww[2][1])*delz;
  wrr[2] = 0.5*(ww[2][0]+ww[0][2])*delx + 0.5*(ww[2][1]+ww[1][2])*dely + (ww[2][2]-wrr[3])*delz;
}

/* ----------------------------------------------------------------------- */
