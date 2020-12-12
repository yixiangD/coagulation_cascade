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

#ifndef FIX_STAT_H
#define FIX_STAT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStat : public Fix {
 friend class FixSolidBound;
 friend class FixBondADCreate;

 public:
  FixStat(class LAMMPS *, int, char **);
  virtual ~FixStat();
  int setmask();
  virtual void end_of_step(){};
  virtual void init(){};
  int map_index(double, double, double);
  int map_index_cyl(double, double, double);
  virtual void write_stat(int){};

 protected:
  char fname[FILENAME_MAX], stype[FILENAME_MAX], dir[1]; 
  int nx, ny, nz;
  int st_start, dump_each, num_step;
  int is, js, ks;
  int xper, yper, zper;
  double xs, ys, zs;
  double xlo, ylo, zlo, xhi, yhi, zhi;
  double dxlo, dylo, dzlo, dxhi, dyhi, dzhi;
  double dxs, dys, dzs;
	double R, centx, centy, centz;
	double high, low, cent2, cent3;
	int n1, n2, n3;
	double *profile;
	int index;
};

}

#endif
