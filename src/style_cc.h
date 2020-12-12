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
#include "atom_vec_cc_atomic.h"
#include "atom_vec_cc_full.h"
#endif

#ifdef AtomClass
AtomStyle(cc/atomic,AtomVecCCAtomic)
AtomStyle(cc/full,AtomVecCCFull)
#endif

#ifdef PairInclude
#include "pair_cc.h"
#endif

#ifdef PairClass
PairStyle(cc,PairCC)
#endif

#ifdef FixInclude
#include "fix_cc_verlet.h"
#include "fix_cc_resetc.h"
#include "fix_dp_force.h"
#include "fix_cc_reflect.h"
#include "fix_cc_source.h"
#include "fix_cc_reaction.h"
#include "fix_cc_bcflux.h"
#include "fix_cc_bcverlet.h"
#endif

#ifdef FixClass
FixStyle(cc/verlet,FixCCVerlet)
FixStyle(cc/resetC,FixCCResetC)
FixStyle(dpforce,FixDPForce)
FixStyle(cc/reflect,FixCCReflect)
FixStyle(cc/source,FixCCSource)
FixStyle(cc/reaction,FixCCReaction)
FixStyle(cc/bcflux,FixCCBcflux)
FixStyle(cc/bcverlet,FixCCBCVerlet)
#endif


#ifdef DumpInclude
#include "dump_cc.h"
#endif

#ifdef DumpClass
DumpStyle(cc,DumpCC)
#endif
