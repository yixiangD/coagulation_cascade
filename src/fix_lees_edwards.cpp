/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "comm.h"
#include "fix_lees_edwards.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "stdio.h"
#include "string.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixLeesEdwards::FixLeesEdwards(LAMMPS *lmp, int narg, char **arg) : Fix(lmp,narg, arg)
{
  if (narg != 4) error->all("Illegal fix Lees/Edwards command");

  vel = atof(arg[3]);
	dir = 0;		// direction of flow: 0 along x; 1 along y
	adir = 1;
}

/* ---------------------------------------------------------------------- */

FixLeesEdwards::~FixLeesEdwards()
{
}

/* ---------------------------------------------------------------------- */

int FixLeesEdwards::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_EXCHANGE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLeesEdwards::initial_integrate(int vflag)
{
  int i;
  double time = update->dt*update->ntimestep;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for(i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if (x[i][adir] >= domain->boxhi[adir]){
        x[i][dir] -= time*vel;
        v[i][dir] -= vel;
      }
      if (x[i][adir] < domain->boxlo[adir]){
        x[i][dir] += time*vel;
        v[i][dir] += vel;
      }
  	}
}

/* ----------------------------------------------------------------------
   exchange atoms between processors 
------------------------------------------------------------------------- */

void FixLeesEdwards::pre_exchange()
{
  double **x = atom->x;
  int *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
			domain->remap(x[i],image[i]);

  comm->irregular();
}

/* ---------------------------------------------------------------------- 
   adjust the velocity of the ghost particles 
------------------------------------------------------------------------- */

void FixLeesEdwards::pre_force(int vflag)
{

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  for(int i = nlocal; i < nall; i++)
    if (mask[i] & groupbit){
    	if (x[i][adir] >= domain->boxhi[adir]) v[i][dir] += vel;
    	if (x[i][adir] < domain->boxlo[adir])	v[i][dir] -= vel;
  	}
}
