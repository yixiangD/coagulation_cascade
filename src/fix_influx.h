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

#ifndef FIX_INFLUX_H
#define FIX_INFLUX_H

#include "fix.h"
#include "fix_com_gyration.h"

namespace LAMMPS_NS {

class FixInflux : public Fix {
 public:
  FixInflux(class LAMMPS *, int, char **);
  ~FixInflux();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
	char fixcomid[128], coord;
	int frequency, type_rbc, type_plat, type_bank, platgroupbit, bankgroupbit, dir_1, dir_2, nmax, Nplat, chkpoint, setupbank, nbonds0;
  double flux, cutoff;
  double xlo, xhi, ylo, yhi, zlo, zhi, rlo, rhi, ycent, zcent, xprd, yprd, zprd;
  int *list_rbc, *list_plat, *list_bank, *list_mol, *loc_mol, *rls_mol;
	double *tempd, *Xo = NULL;

	class FixComGyration *comgyrat = NULL;

	int chkregion(int);
};

}

#endif
