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

#ifdef FixInclude
#include "fix_stat.h"
#include "fix_stat_com_dens.h"
#include "fix_stat_dens.h"
#include "fix_stat_dens_cyl.h"
#include "fix_stat_diff_displ.h"
#include "fix_stat_diff_poly.h"
#include "fix_stat_gyroradius.h"
#include "fix_stat_snap.h"
#include "fix_stat_stress.h"
#include "fix_stat_temperature.h"
#include "fix_stat_vel.h"
#include "fix_stat_vel_cyl.h"
#include "fix_stat_velT.h"
#include "fix_stat_velT_cyl.h"
#include "fix_stat_all.h"
#endif

#ifdef FixClass
FixStyle(stat,FixStat)
FixStyle(stat/com/dens,FixStatComDens)
FixStyle(stat/dens,FixStatDens)
FixStyle(stat/dens/cyl,FixStatDensCyl)
FixStyle(stat/diff_displ,Fix_Stat_Diff_Displ)
FixStyle(stat/diff_poly,Fix_Stat_Diff_Poly)
FixStyle(stat/gyroradius,Fix_Stat_Gyroradius)
FixStyle(stat/snap,FixStatSnap)
FixStyle(stat/stress,FixStatStress)
FixStyle(stat/temperature,FixStatTemp)
FixStyle(stat/vel,FixStatVel)
FixStyle(stat/vel/cyl,FixStatVelCyl)
FixStyle(stat/velT,FixStatVelT)
FixStyle(stat/velT/cyl,FixStatVelTCyl)
FixStyle(stat/all,FixStatAll)
#endif

