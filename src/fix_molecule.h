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

#ifndef FIX_MOLECULE_H
#define FIX_MOLECULE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMolecule : public Fix {
 friend class FixCCReaction;
 friend class BondCellular;
 friend class DihedralBend;
 friend class AngleAreaVolume;
 friend class PairDPDMisc;

 public:
  FixMolecule(class LAMMPS *, int, char **);
  ~FixMolecule();
  int setmask();
  void init();
  void end_of_step();

 protected:
  int *tempi;
	double *tempd, *list_tact;
};

}

#endif
