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
#include <cmath>
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_cc_source.h"
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

FixCCSource::FixCCSource(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int arg_index = 3;
  step = atoi( arg[arg_index++] );
  typec = atoi( arg[arg_index++] );
  wdim  = atoi( arg[arg_index++] );

  if(typec > CTYPES) 
	  error->all("Illegal fix cc_source command -- type option wrong.\n");   

  if (strcmp(arg[arg_index],"circle") == 0)       	index = 0;
  else if (strcmp(arg[arg_index],"rectangle") == 0)   index = 1;
  else error->all("Illegal fix cc_source command");
	arg_index++;

  if(index == 0){
  	if (narg != 11 ) error->all("Illegal fix cc_source command");
	
		center[0] = atof( arg[arg_index++] );
		center[1] = atof( arg[arg_index++] );
		radius = atof( arg[arg_index++] );
    value  = atof( arg[arg_index++] );
  }
  else if(index == 1){
  	if (narg != 12 ) error->all("Illegal fix cc_source command");
	
		center[0] = atof( arg[arg_index++] );
		center[1] = atof( arg[arg_index++] );
		length = atof( arg[arg_index++] );
		width  = atof( arg[arg_index++] );
    value  = atof( arg[arg_index++] );
  }
  else error->all("Illegal fix cc_source command");
}

/* ---------------------------------------------------------------------- */

FixCCSource::~FixCCSource()
{
}

/* ---------------------------------------------------------------------- */

int FixCCSource::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCSource::init()
{
}

/* ---------------------------------------------------------------------- */


void FixCCSource::post_force(int vflag)
{
  double **x = atom->x;
  double **T = atom->T;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double drx, dry, rsq;
  double radius_sq = radius*radius;

  if(update->ntimestep > step) 
  for (int i = 0; i < nlocal; i++) {

	if (mask[i] & groupbit) {
		if(index == 0){
			drx = x[i][(wdim+1)%3] - center[0];
		  dry = x[i][(wdim+2)%3] - center[1];
		  rsq = drx*drx + dry*dry;
		  if(rsq < radius_sq)			
				Q[i][typec-1] += value;
		}
		else if(index == 1){
			drx = x[i][(wdim+1)%3] - center[0];
		  dry = x[i][(wdim+2)%3] - center[1];
		  if(std::abs(drx) < length && std::abs(dry) < width)			
				Q[i][typec-1] += value;
		}
  }

  }
}

/* ---------------------------------------------------------------------- */

