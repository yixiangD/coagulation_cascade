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

#ifndef FIX_CC_BCFLUX_H
#define FIX_CC_BCFLUX_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCCBcflux : public Fix {
 public:
  FixCCBcflux(class LAMMPS *, int, char **);
  ~FixCCBcflux();
  int setmask();
  void init();
  void post_force(int);

 protected:
  double Ri[7];
  double Diff[23];
  double *Diff_inv;
  double par[9];
  double cw[3];

};

}

#endif
