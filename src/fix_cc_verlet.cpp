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
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_cc_verlet.h"
#include "modify.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "iostream"
#include "fstream"
#include "math.h"
#include "common.h"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

FixCCVerlet::FixCCVerlet(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"cc/verlet") != 0 && narg != 6)
    error->all("Illegal fix CCVerlet command");
	time_integrate = 1;

	if (strcmp(arg[3],"x") == 0) dir = 0;
	else if (strcmp(arg[3],"y") == 0) dir = 1;
	else if (strcmp(arg[3],"z") == 0) dir = 2;
	else error->all("cc/verlet: Inflow direction not valid");	
	inlet_loc = atof(arg[4]);	// exact location of inlet

  char fname[FILENAME_MAX];
  sprintf(fname,arg[5]);

	// initial conditions for all species
  initial = (double *) memory->smalloc(CTYPES*sizeof(double),"cc/verlet:initial");

 	if (comm->me == 0){
  	fstream file_read;
    file_read.open(fname,ios::in);

    if(file_read.fail())
    	error->one("cc/verlet: Cannot read concentration initial condition file");

    for(int i = 0; i < CTYPES; i++)
    	file_read >> initial[i];

    file_read.close();
  }
  MPI_Bcast(initial,CTYPES,MPI_DOUBLE,0,world);
}

/* ---------------------------------------------------------------------- */

FixCCVerlet::~FixCCVerlet()
{
	memory->sfree(initial);
}

/* ---------------------------------------------------------------------- */

int FixCCVerlet::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCVerlet::init()
{
}

/* ---------------------------------------------------------------------- */

void FixCCVerlet::initial_integrate(int vflag)
{
	double dtTm;

	double **T = atom->T;
	double **Q = atom->Q;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

	for (int i = 0; i < nlocal; i++){
		if (mask[i] & groupbit){
			dtTm = 0.5*update->dt;

			for(int k = 0; k < CTYPES; k++)
			{	
				T[i][k] += Q[i][k] * dtTm ;
				T[i][k] = T[i][k] > 0 ? T[i][k] : 0.0;
			}
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixCCVerlet::final_integrate()
{
  double dtTm;

  double **T = atom->T;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			dtTm = 0.5*update->dt;                

			for(int k = 0; k < CTYPES; k++)
			{
				T[i][k] += dtTm * Q[i][k];
				T[i][k] = T[i][k] > 0 ? T[i][k] : 0.0;
			}
   	}
  }
}

/* ---------------------------------------------------------------------- */

void FixCCVerlet::pre_force(int vflag)
{
  double **x = atom->x;
  double **T = atom->T;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;

// zeroing concentrations for particles moving out of the boundary
  for (int i = 0; i < nall; i++)
  	if (mask[i] & groupbit)
			if (x[i][dir] <= inlet_loc || x[i][dir] > domain->boxhi[dir])
	    	for(int k = 0; k < CTYPES; k++)	T[i][k] = initial[k];
}

/* ---------------------------------------------------------------------- */
