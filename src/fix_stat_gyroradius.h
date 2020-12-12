/* ----------------------------------------------------------------------
   Dmitry Fedosov 08/12/05 - accumulation of statistics

   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef FIX_STAT_GYRORADIUS_H
#define FIX_STAT_GYRORADIUS_H

#include "fix_stat.h"

namespace LAMMPS_NS{

class Fix_Stat_Gyroradius : public FixStat {
 public:
  Fix_Stat_Gyroradius(class LAMMPS *, int, char **);
  virtual ~Fix_Stat_Gyroradius();
  void end_of_step();
  void write_stat(int);
  int setmask();
  double *rad, *c_m, *c_mt;
  int nm, init_on;
};

}
#endif
