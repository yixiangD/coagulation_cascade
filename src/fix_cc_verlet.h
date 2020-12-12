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

#ifndef FIX_CC_VERLET_H
#define FIX_CC_VERLET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCCVerlet : public Fix {
 public:
  FixCCVerlet(class LAMMPS *, int, char **);
  ~FixCCVerlet();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void pre_force(int);

 protected:
	int dir;
	double inlet_loc;
	double *initial;

};

}

#endif
