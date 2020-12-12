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

#ifndef ANGLE_AREA_VOLUME_H
#define ANGLE_AREA_VOLUME_H

#include "stdio.h"
#include "angle.h"
#include "fix_molecule.h"

namespace LAMMPS_NS {

class AngleAreaVolume : public Angle {
 public:
  AngleAreaVolume(class LAMMPS *);
  ~AngleAreaVolume();
  void compute(int, int);
  void coeff(int, int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 private:
  double *ka, *a0, *kv, *v0, *kl, *t_remod;
  int *ttyp;
  int init_on;
  class FixMolecule *molecule = NULL;

  void allocate();
  double cell_remodel(int, int, int, double);  // Alireza: cell membrane remodeling upon activation
};

}

#endif
