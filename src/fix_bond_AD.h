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

#ifndef FIX_BOND_AD_H
#define FIX_BOND_AD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBondAD : public Fix {
 public:
  FixBondAD(class LAMMPS *, int, char **);
  ~FixBondAD();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
  int me;
  int iatomtype,jatomtype;
  int btype,ligtype,seed;
  int Nlig,imaxbond,jmaxbond;
  double cutsq;
	double kf,ks,r0,kf0,sig,temp;
	double bondforce;

  int createcount,createcounttotal;
  int nmax;
  int *bondcount;

  class RanMars *random;
  class NeighList *list;
  int countflag;
};

}

#endif
