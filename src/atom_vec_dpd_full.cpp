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

#include "stdlib.h"
#include "atom_vec_dpd_full.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecDPDFull::AtomVecDPDFull(LAMMPS *lmp, int narg, char **arg) : 
  AtomVecFull(lmp, narg, arg)
{
  molecular = 1;
  bonds_allow = 1;
  angles_allow = 1;
  dihedrals_allow = 1;
  impropers_allow = 1;
  mass_type = 1;
  comm_x_only = 0;
  ghost_velocity = 1;
  size_comm = 6;
  size_reverse = 3;
  size_border = 12;
  size_data_atom = 7;
  size_data_vel = 4;
  xcol_data = 5;

  atom->molecule_flag = atom->q_flag = 1;
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::pack_comm(int n, int *list, double *buf,
			   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
      dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
      dz = pbc[2]*domain->zprd;
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::pack_comm_one(int i, double *buf)
{
  buf[0] = v[i][0];
  buf[1] = v[i][1];
  buf[2] = v[i][2];
  return 3;
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDFull::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++]; 
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::unpack_comm_one(int i, double *buf)
{
  v[i][0] = buf[0];
  v[i][1] = buf[1];
  v[i][2] = buf[2];
  return 3;
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::pack_border(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;
  int img;
  int pbc_index;
  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = static_cast<double> (tag[j]);
      buf[m++] = static_cast<double> (type[j]);
      buf[m++] = static_cast<double> (mask[j]);
      buf[m++] = q[j];
      buf[m++] = static_cast<double> (molecule[j]);
      buf[m++] = static_cast<double> (image[j]);
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = molecule[j];
      
      img = image[j]; 
      pbc_index = pbc[0];
      if(pbc_index > 0)
      {
        while(pbc_index != 0)
        {
          int idim = img & 1023;
          int otherdims = img ^ idim;
          idim--;
          idim &= 1023;
          img = otherdims | idim;
          --pbc_index;
        }
      }
      if(pbc_index < 0)
      {
        while(pbc_index != 0)
        {
          int idim = img & 1023;
          int otherdims = img ^ idim;
          idim++;
          idim &= 1023;
          img = otherdims | idim;
          ++pbc_index;
        }
      }
      
      pbc_index = pbc[1];
      if(pbc_index > 0)
      {
        while(pbc_index != 0)
        {
          int idim = (img >> 10) & 1023;
          int otherdims = img ^ (idim << 10);
          idim--;
          idim &= 1023;
          img = otherdims | (idim << 10);
          --pbc_index;
        }
      }
      if(pbc_index < 0)
      {
        while(pbc_index != 0)
        {
          int idim = (img >> 10) & 1023;
          int otherdims = img ^ (idim << 10);
          idim++;
          idim &= 1023;
          img = otherdims | (idim << 10);
          ++pbc_index;
        }
      }
      
      pbc_index = pbc[2];
      if(pbc_index > 0)
      {
        while(pbc_index != 0)
        {
          int idim = (img >> 20);
          int otherdims = img ^ (idim << 20);
          idim--;
          idim &= 1023;
          img = otherdims | (idim << 20);
          --pbc_index;
        }
      }
      if(pbc_index < 0)
      {
        while(pbc_index != 0)
        {
          int idim = (img >> 20) ;
          int otherdims = img ^ (idim << 20);
          idim++;
          idim &= 1023;
          img = otherdims | (idim << 20);
          ++pbc_index;
        }
      }
      buf[m++] = static_cast<double>(img);
      //buf[m++] = image[j] - (pbc[0]+pbc[1]*1000+pbc[2]*1000000);
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::pack_border_one(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = molecule[i];
  buf[2] = v[i][0];
  buf[3] = v[i][1];
  buf[4] = v[i][2];
  return 5;
}

/* ---------------------------------------------------------------------- */

void AtomVecDPDFull::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    q[i] = buf[m++];
    molecule[i] = static_cast<int> (buf[m++]);
    image[i] = static_cast<int> (buf[m++]);
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecDPDFull::unpack_border_one(int i, double *buf)
{
  q[i] = buf[0];
  molecule[i] = static_cast<int> (buf[1]);
  v[i][0] = buf[2];
  v[i][1] = buf[3];
  v[i][2] = buf[4];
  return 5;
}

/* ---------------------------------------------------------------------- */

