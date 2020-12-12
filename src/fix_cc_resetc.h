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

#ifndef FIX_CC_RESETC_H
#define FIX_CC_RESETC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCCResetC : public Fix {
 public:
  FixCCResetC(class LAMMPS *, int, char **);
  ~FixCCResetC();
  int setmask();
  void init();
  void final_integrate();

 protected:
  int step;
  int typec;
  int index;
  double center[2], radius, length, width;
  double value;
};

}

#endif
