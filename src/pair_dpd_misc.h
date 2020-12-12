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

#ifndef PAIR_DPD_MISC_H
#define PAIR_DPD_MISC_H

#include "pair.h"
#include "fix_molecule.h"

namespace LAMMPS_NS {

class PairDPDMisc : public Pair {
 public:
  PairDPDMisc(class LAMMPS *);
  ~PairDPDMisc();
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
  double cut_global;
  int seed;
  int **force_type, **spec;
  double **a0,**gamma,**sigma,**cut_dpd,**weight_exp;
  double **de, **r0m, **beta, **cut_morse, **morse1;
  double **lj1, **lj2, **epsilon, **sigma_lj, **cut_lj;
  int ifix = -1;
  class RanMars *random;
	class FixMolecule *molecule = NULL;

  void allocate();
	int ActState(int, int);		// determine the activation state of either molecules
};

}

#endif
