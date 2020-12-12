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

#ifndef PAIR_DPD_ADH_H
#define PAIR_DPD_ADH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairDPDAdh : public Pair {
 public:
  PairDPDAdh(class LAMMPS *);
  ~PairDPDAdh();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);

 private:

  int seed, Nreceptor, imaxbond, jmaxbond;
  int **force_type, **ligand_type, **bond_type;
  double cut_dpd_global;
  double **a0,**gamma,**sigma,**cut_dpd,**weight_exp;
  double **ks,**r0,**kf0,**sig,**temp,**cut_spring;
  double **de, **r0m, **beta, **cut_morse;
  double **lj1, **lj2, **epsilon, **sigma_lj;
  class RanMars *random;

  void allocate();
};

}

#endif
