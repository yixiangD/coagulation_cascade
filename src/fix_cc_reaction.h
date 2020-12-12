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

#ifndef FIX_CC_REACTION_H
#define FIX_CC_REACTION_H

#include "fix.h"
#include "fix_molecule.h"

namespace LAMMPS_NS {

class FixCCReaction : public Fix {
 public:
  FixCCReaction(class LAMMPS *, int, char **);
  ~FixCCReaction();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void end_of_step();

 protected:
  int sys, Np, dir;
	double lo, hi;
  double *par, *Ri;
	class FixMolecule *molecule = NULL;

};

}

#endif
