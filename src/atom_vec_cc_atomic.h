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

#ifndef ATOM_VEC_CC_ATOMIC_H
#define ATOM_VEC_CC_ATOMIC_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecCCAtomic : public AtomVec {
 public:
  AtomVecCCAtomic(class LAMMPS *, int, char **);
  virtual ~AtomVecCCAtomic() {}
  virtual void grow(int);
  virtual void copy(int, int);
  virtual int pack_comm(int, int *, double *, int, int *);
  virtual void unpack_comm(int, int, double *);
  virtual int pack_reverse(int, int, double *);
  virtual void unpack_reverse(int, int *, double *);
  virtual int pack_border(int, int *, double *, int, int *);
  virtual void unpack_border(int, int, double *);
  virtual int pack_exchange(int, double *);
  virtual int unpack_exchange(double *);
  virtual int size_restart();
  virtual int pack_restart(int, double *);
  virtual int unpack_restart(double *);
  virtual void create_atom(int, double *);
  virtual void data_atom(double *, int, char **);
  virtual int data_atom_hybrid(int, char **);
  virtual double memory_usage();

 protected:
  int *tag,*type,*mask,*image;
  double **x,**v,**f;
	double *q;
  double **T;
  double **Q;
};

}

#endif
