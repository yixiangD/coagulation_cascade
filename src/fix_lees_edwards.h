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

#ifndef FIX_LEES_EDWARDS_H
#define FIX_LEES_EDWARDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLeesEdwards : public Fix {
 public:
  FixLeesEdwards(class LAMMPS *, int, char **);
  ~FixLeesEdwards();
  int setmask();
  void initial_integrate(int);
  void pre_exchange();
  void pre_force(int);

 private:
	int dir, adir;
  double vel;
};

}

#endif
