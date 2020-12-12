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

#ifndef FIX_BOND_AD_BREAK_H
#define FIX_BOND_AD_BREAK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondADBreak : public Fix {
 public:
  FixBondADBreak(class LAMMPS *, int, char **);
  ~FixBondADBreak();
  int setmask();
  void init();
	void setup(int);
  void post_integrate();

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double compute_vector(int);
  double memory_usage();

 private:
	char bmodel[128];
  int me;
  int btype,seed;
  double cutsq;
	double ks,r0,kr0,sig,temp;

  int breakcount,breakcounttotal;
  int nmax;
  int *partner;
  double *pp,*probability;

  class RanMars *random;
  int nlevels_respa;
};

}

#endif
