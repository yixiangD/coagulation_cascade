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
/* --------------------------------------------------------------------
   changed by Huan Lei
   June 17 2009
   ---------------------------------------------------------------------*/


#ifndef FIX_ADD_FORCE_H
#define FIX_ADD_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAddForce : public Fix {
 public:
  FixAddForce(class LAMMPS *, int, char **);
  virtual ~FixAddForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);
  void interp_force(double *);

 private:
  double xvalue,yvalue,zvalue,ff;
  double foriginal[4],foriginal_all[4];
  double force, freq, dtv;
  int force_flag;
  int dir, index;         // 0 - force, 1 - torque, 2 - periodic force, 3 - periodic torque, 4 - RPF 
  int nlevels_respa;
  int dimx, dimy;
  double Xmin, Xmax, Ymin, Ymax;
  double **meshX, **meshY, **gradxP, **gradyP;
	double R0, alpha, len, X0;
};
  
}

#endif
