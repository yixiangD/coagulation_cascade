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

#ifndef PAIR_CC_H
#define PAIR_CC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCC : public Pair {
 public:
  PairCC(class LAMMPS *);
  ~PairCC();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 private:
  int seed, dir;
  double cut_global, outlet_loc;
  double **cut;
  double **cutinv;
  double **cutc;
  double **a0;
  double **gamma;
  double **sigma;
  double **s1;
  double ***kC;
  double ***kappa;
  double ***s2;
  class RanMars *random;

  void allocate();
};

}

#endif
