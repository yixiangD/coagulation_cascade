/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef BOND_CELLULAR_H
#define BOND_CELLULAR_H

#include "stdio.h"
#include "bond.h"
#include "fix_molecule.h"

namespace LAMMPS_NS {

class BondCellular : public Bond {
 public:
  BondCellular(class LAMMPS *);
  ~BondCellular();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int);

 private:
  int *bond_index, onoff;
  double *ks,*r0,*temp,*dsig,*kr0,*rcs,*mu_targ,*qp,*t_remod;
  double dtv;
  double *gamc, *gamt, *sigc, *sigt;
	double bondforce;
  double wrr[4], delx, dely, delz;
	class RanMars *random;
  int fix_check;                   // # of fixes that induce reneigh
  int *fixchecklist;               // which fixes to check
	class FixMolecule *molecule = NULL;

  void allocate();
  void generate_wrr();
  double cell_remodel(int, int, double);	// Alireza: cell membrane remodeling upon activation
};

}
#endif
