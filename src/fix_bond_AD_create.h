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

#ifndef FIX_BOND_AD_CREATE_H
#define FIX_BOND_AD_CREATE_H

#include "fix.h"
#include "fix_stat.h"
#include "fix_stat_all.h"

namespace LAMMPS_NS {

class FixBondADCreate : public Fix {
 public:
  FixBondADCreate(class LAMMPS *, int, char **);
  ~FixBondADCreate();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  double compute_vector(int);
  double memory_usage();

 private:
	char bmodel[128];
  int me;
  int iatomtype,jatomtype;
  int btype,seed;
  int imaxbond,jmaxbond;
	int Nbond;
  double cutsq;
	double ks,r0,kf0,sig,temp;

  int createcount,createcounttotal;
  int nmax;
  int *partner,*bondcount;
  double *pp,**probability;

  class RanMars *random;
  class NeighList *list;
  int countflag,reverseflag,partnerflag;

  char fixstatid[128];
  class FixStatAll *stat = NULL;
};

}

#endif
