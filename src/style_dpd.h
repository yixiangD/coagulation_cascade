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

#ifdef AtomInclude
#include "atom_vec_dpd_atomic.h"
#include "atom_vec_dpd_bond.h"
#include "atom_vec_dpd_angle.h"
#include "atom_vec_dpd_full.h"
#include "atom_vec_dpd_charge.h"
#include "atom_vec_dpd_molecular.h"
#endif

#ifdef AtomClass
AtomStyle(dpd/atomic,AtomVecDPDAtomic)
AtomStyle(dpd/bond,AtomVecDPDBond)
AtomStyle(dpd/angle,AtomVecDPDAngle)
AtomStyle(dpd/full,AtomVecDPDFull)
AtomStyle(dpd/charge,AtomVecDPDCharge)
AtomStyle(dpd/molecular,AtomVecDPDMolecular)
#endif

#ifdef PairInclude
#include "pair_dpd.h"
#include "pair_dpd_adh.h"
#include "pair_dpd_misc.h"
#endif

#ifdef PairClass
PairStyle(dpd,PairDPD)
PairStyle(dpd/adh,PairDPDAdh)
PairStyle(dpd/misc,PairDPDMisc)
#endif
