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

#ifndef DIHEDRAL_BEND_H
#define DIHEDRAL_BEND_H

#include "stdio.h"
#include "dihedral.h"
#include "fix_molecule.h"

namespace LAMMPS_NS {

class DihedralBend : public Dihedral {
 public:
  DihedralBend(class LAMMPS *);
  ~DihedralBend();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *k,*cos_shift,*sin_shift;
  int *sign,*multiplicity;
  double *theta0,*t_remod;
	class FixMolecule *molecule = NULL;

  void allocate();
  double cell_remodel(int, int, int, int, double);	// Alireza: cell membrane remodeling upon activation
};

}

#endif
