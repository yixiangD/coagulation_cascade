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
#include "fix_cc_bcverlet.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "iostream"
#include "common.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixCCBCVerlet::FixCCBCVerlet(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (strcmp(style,"cc_bcverlet") != 0 && narg != 4)
    error->all("Illegal fix CCBCVerlet command");

  time_integrate = 1;
//  printf("caution:fix_cc_bcverlet, cpu is %d\n", comm->me);
//  MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

FixCCBCVerlet::~FixCCBCVerlet()
{
}

/* ---------------------------------------------------------------------- */

int FixCCBCVerlet::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCBCVerlet::init()
{
}

/* ---------------------------------------------------------------------- */

void FixCCBCVerlet::initial_integrate(int vflag)
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

			T[i][0] += Q[i][0] * dtTm ;
			T[i][0] = T[i][0] > 0 ? T[i][0] : 0.0;
			T[i][1] += Q[i][1] * dtTm ;
			T[i][1] = T[i][1] > 0 ? T[i][1] : 0.0;
			T[i][6] += Q[i][6] * dtTm ;
			T[i][6] = T[i][6] > 0 ? T[i][6] : 0.0;
			T[i][7] += Q[i][7] * dtTm ;
			T[i][7] = T[i][7] > 0 ? T[i][7] : 0.0;
			T[i][12]+= Q[i][12] * dtTm ;
			T[i][12] = T[i][12] > 0 ? T[i][12] : 0.0;
			T[i][13]+= Q[i][13] * dtTm ;
			T[i][13] = T[i][13] > 0 ? T[i][13] : 0.0;
			T[i][19]+= Q[i][19] * dtTm ;
			T[i][19] = T[i][19] > 0 ? T[i][19] : 0.0;
			
		}
	}
}

/* ---------------------------------------------------------------------- */

void FixCCBCVerlet::final_integrate()
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

		T[i][0] += Q[i][0] * dtTm ;
		T[i][0] = T[i][0] > 0 ? T[i][0] : 0.0;
		T[i][1] += Q[i][1] * dtTm ;
		T[i][1] = T[i][1] > 0 ? T[i][1] : 0.0;
		T[i][6] += Q[i][6] * dtTm ;
		T[i][6] = T[i][6] > 0 ? T[i][6] : 0.0;
		T[i][7] += Q[i][7] * dtTm ;
		T[i][7] = T[i][7] > 0 ? T[i][7] : 0.0;
		T[i][12]+= Q[i][12] * dtTm ;
		T[i][12] = T[i][12] > 0 ? T[i][12] : 0.0;
		T[i][13]+= Q[i][13] * dtTm ;
		T[i][13] = T[i][13] > 0 ? T[i][13] : 0.0;
		T[i][19]+= Q[i][19] * dtTm ;
		T[i][19] = T[i][19] > 0 ? T[i][19] : 0.0;
    	}
  }
}

/* ---------------------------------------------------------------------- */

