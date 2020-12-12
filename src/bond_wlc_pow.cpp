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
#include "bond_wlc_pow.h"
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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondWLC_POW::BondWLC_POW(LAMMPS *lmp) : Bond(lmp)
{
}

/* ---------------------------------------------------------------------- */

BondWLC_POW::~BondWLC_POW()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(temp);
    memory->sfree(r0);
    memory->sfree(lambda);
    memory->sfree(qp);
    memory->sfree(kp);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLC_POW::compute(int eflag, int vflag)
{
  double ebond,fbond;
  
  int i1,i2,n,type,factor;
  double delx,dely,delz,rr,fforce,rfactor,ra,rlogarg,kph,l0;
  char warning[128];
  double ff[6];
  
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;  
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
    domain->minimum_image(delx,dely,delz);

    // force from log term
    
    ra = sqrt(delx*delx + dely*dely + delz*delz);
    rr = 1.0/r0[type]; 
    kph = pow(l0,qp[type])*temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr)/lambda[type];
    rr = ra/l0/r0[type]; 
    rlogarg = pow(ra,qp[type]+1.0);

    if (rr >= 1.0) {
      sprintf(warning,"WLC bond too long: %d %d %d %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],rr);
      error->warning(warning);
    }   

    /*if (rlogarg < 0.001) {
      sprintf(warning,"POW bond too short: %d %d %d %g",
              update->ntimestep,atom->tag[i1],atom->tag[i2],ra);
      error->warning(warning);
      rlogarg = 0.001;
    } */
    
    fforce = - temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr)/lambda[type]/ra + kph/rlogarg;


    // energy
    ebond = 0.0;
    if (eflag) {
      ebond += 0.25*temp[type]*l0*r0[type]*(3.0*rr*rr-2.0*rr*rr*rr)/(1.0-rr)/lambda[type];
      if (qp[type] == 1.0)
        ebond -= kph*log(ra);
      else
        ebond += kph/(qp[type]-1.0)/pow(ra,qp[type]-1.0);
    }

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fforce;
      f[i1][1] += dely*fforce;
      f[i1][2] += delz*fforce;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fforce;
      f[i2][1] -= dely*fforce;
      f[i2][2] -= delz*fforce;
    }
    
    // virial contribution

    if (stress_ind == 2){
      ff[0] = delx*delx*fforce;
      ff[1] = dely*dely*fforce;
      ff[2] = delz*delz*fforce;
      ff[3] = delx*dely*fforce;
      ff[4] = delx*delz*fforce;
      ff[5] = dely*delz*fforce;
    
      for(int k = 0; k < n_stress_tot; ++k)
        modify->fix[stress_list[k]]->virial3(i1, i2, ff);
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fforce,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondWLC_POW::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;
  temp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:temp");
  r0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:r0");
  lambda = (double *) memory->smalloc((n+1)*sizeof(double),"bond:lambda");
  kp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:kp");
  qp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:qp");
  
  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");

  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondWLC_POW::coeff(int narg, char **arg)
{
  if (narg != 6) error->all("Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double temp_one = atof(arg[1]);
  double r0_one = atof(arg[2]);
  double lambda_one = atof(arg[3]);
  double kp_one = atof(arg[4]);
  double qp_one = atof(arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    temp[i] = temp_one;
    r0[i] = r0_one;
    lambda[i] = lambda_one;
    kp[i] = kp_one;
    qp[i] = qp_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   check if special_bond settings are valid
------------------------------------------------------------------------- */

void BondWLC_POW::init_style()
{
  // special bonds should be 0 1 1
  /*
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) {
    if (comm->me == 0)
      error->warning("Use special bonds = 0,1,1 with bond style wlc");
  }
  */
}

/* ---------------------------------------------------------------------- */

double BondWLC_POW::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void BondWLC_POW::write_restart(FILE *fp)
{
  fwrite(&temp[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&lambda[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&kp[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&qp[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void BondWLC_POW::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&temp[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&lambda[1],sizeof(double),atom->nbondtypes,fp);
    fread(&kp[1],sizeof(double),atom->nbondtypes,fp);
    fread(&qp[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&temp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&lambda[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kp[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&qp[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondWLC_POW::single(int type, double rsq, int i, int j)
{
  double ra = sqrt(rsq);
  double rr = 1.0/r0[type]; 
  double l0 = 1.0;
  double kph = pow(l0,qp[type])*temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr)/lambda[type];
  rr = ra/l0/r0[type]; 
  double ebond = 0;
  ebond += 0.25*temp[type]*l0*r0[type]*(3.0*rr*rr-2.0*rr*rr*rr)/(1.0-rr)/lambda[type];
  if (qp[type] == 1.0)
    ebond -= kph*log(ra);
  else
    ebond += kph/(qp[type]-1.0)/pow(ra,qp[type]-1.0);

  return ebond;
}
