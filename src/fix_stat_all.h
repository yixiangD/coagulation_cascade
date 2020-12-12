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

#ifndef FIX_STAT_ALL_H
#define FIX_STAT_ALL_H

#include "fix_stat.h"

namespace LAMMPS_NS {

class FixStatAll : public FixStat {
 public:
  FixStatAll(class LAMMPS *, int, char **);
  ~FixStatAll();
  int setmask();
  void initial_integrate(int);
  void init();
  void end_of_step();
  void write_stat(int);
  void virial1(int);
  void virial2(int, double[]);
  void virial3(int, int, double[]);
  void virial4(int, int, int, double[]);
  void virial5(int, int, int, int, double[]);
  void virial_f(int, double[]);

  double ***num, ****vel, ****aT;
  double ****stress, ****pstress;

 protected:
	double ***tmp;
	double ***numtmp, ****veltmp, ****aTtmp;
	double ****v_avg, ****v_avg_tmp;
  double ****ss, ****ss_p, ****vv, ****vv_p;
  int poly_ind, is1, js1, ks1, map1;
};

}

#endif
