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
   Modified by Zhen Li
   March 25 2013
   ----------------------------------------------------------------------*/
#include "string.h"
#include "stdlib.h"
#include "fix_cc_reflect.h"
#include "atom.h"
#include "modify.h"
#include "domain.h"
#include "error.h"
#include "update.h"
#include "comm.h"

using namespace LAMMPS_NS; 

/* ---------------------------------------------------------------------- */

FixCCReflect::FixCCReflect(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all("Illegal fix cc/reflect command");
  
  dim  = atoi(arg[3]);
  lvalue = atof(arg[4]);
  hvalue = atof(arg[5]);

  xflag = yflag = zflag = 0;
  if (dim == 0) xflag = 1;
  else if (dim == 1) yflag = 1;
  else if (dim == 2) zflag = 1;
  else error->all("Illegal fix cc/reflect command");

  if (xflag && domain->xperiodic )
    error->all("Cannot use wall in periodic dimension");
  if (yflag && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (zflag && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");

//  printf("caution:fix_cc_reflect, cpu is %d\n", comm->me);
//  MPI_Barrier(world);

}

/* ---------------------------------------------------------------------- */

int FixCCReflect::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCReflect::post_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if ( x[i][dim] < lvalue) {
	x[i][dim] = lvalue + (lvalue - x[i][dim]);
//	v[i][dim] = -v[i][dim];
        for(int j=0; j<3; j++)
	   v[i][j] = -v[i][j];	  
      }
      if ( x[i][dim] > hvalue) {
	x[i][dim] = hvalue - (x[i][dim] - hvalue);
//	v[i][dim] = -v[i][dim];
	for(int j=0; j<3; j++)
	   v[i][j] = -v[i][j];
      }
    }
  }
}
