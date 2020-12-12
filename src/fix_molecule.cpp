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
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_molecule.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "iostream"
#include "fstream"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

FixMolecule::FixMolecule(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if(strcmp(style,"molecule") != 0)
  	error->all("Illegal fix Molecule command");
 
	nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixMolecule::~FixMolecule()
{
  memory->sfree(tempi);
	memory->sfree(tempd);
	memory->sfree(list_tact);
}

/* ---------------------------------------------------------------------- */

int FixMolecule::setmask()
{
  int mask = 0;
	mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMolecule::init()
{
	tempi = (int *) memory->smalloc(atom->n_mol*sizeof(int),"molecule:tempi");
	tempd = (double *) memory->smalloc(atom->n_mol*sizeof(double),"molecule:tempd");
	// time of activation of molecules
	list_tact = (double *) memory->smalloc(atom->n_mol*sizeof(double),"molecule:list_tact");
	for (int i = 0; i < atom->n_mol; i++) list_tact[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixMolecule::end_of_step()
{
  double **T = atom->T;
  double *tact = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  	if (mask[i] & groupbit){
			if ( tact[i] == 0.0 ) continue;
			if (list_tact[ atom->molecule[i] ] == 0.0)
				list_tact[ atom->molecule[i] ] = tact[i];
			else if (tact[i] < list_tact[ atom->molecule[i] ])
				list_tact[ atom->molecule[i] ] = tact[i];
//			printf("Molecule: activation time for molecule %d is: %f\n", atom->molecule[i], list_tact[atom->molecule[i]]);
		}

	MPI_Allreduce(&list_tact[0],&tempd[0],atom->n_mol,MPI_DOUBLE,MPI_MAX,world);

	for (int i = 0; i < atom->n_mol; i++) list_tact[i] = tempd[i];
}
 
/* ---------------------------------------------------------------------- */
