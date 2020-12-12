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

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <algorithm>
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "dihedral_bend.h"
#include "atom.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001

/* ---------------------------------------------------------------------- */

DihedralBend::DihedralBend(LAMMPS *lmp) : Dihedral(lmp) {}

/* ---------------------------------------------------------------------- */

DihedralBend::~DihedralBend()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(k);
    memory->sfree(theta0);
    memory->sfree(t_remod);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::compute(int eflag, int vflag)
{
  double edihedral,f1[3],f2[3],f3[3],f4[3];
  
  int n,i1,i2,i3,i4,type;
  double d21x,d21y,d21z,d31x,d31y,d31z,d32x,d32y,d32z;
  double d34x,d34y,d34z,d24x,d24y,d24z,d14x,d14y,d14z;
  double n1x,n1y,n1z,n2x,n2y,n2z,n1,n2,nn;
  double costheta,sintheta,mx;
  double alfa,a11,a12,a22;
  double s1x,s1y,s1z,s2x,s2y,s2z,s3x,s3y,s3z,s4x,s4y,s4z;
  double ff[6];
  
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

  edihedral = 0.0;
  energy = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 2-1 distance
    d21x = x[i2][0] - x[i1][0];
    d21y = x[i2][1] - x[i1][1];
    d21z = x[i2][2] - x[i1][2];
    domain->minimum_image(d21x,d21y,d21z);

    // 3-1 distance
    d31x = x[i3][0] - x[i1][0];
    d31y = x[i3][1] - x[i1][1];
    d31z = x[i3][2] - x[i1][2];
    domain->minimum_image(d31x,d31y,d31z);

    // 3-2 distance
    d32x = x[i3][0] - x[i2][0];
    d32y = x[i3][1] - x[i2][1];
    d32z = x[i3][2] - x[i2][2];
    domain->minimum_image(d32x,d32y,d32z);

    // 3-4 distance
    d34x = x[i3][0] - x[i4][0];
    d34y = x[i3][1] - x[i4][1];
    d34z = x[i3][2] - x[i4][2];
    domain->minimum_image(d34x,d34y,d34z);

    // 2-4 distance
    d24x = x[i2][0] - x[i4][0];
    d24y = x[i2][1] - x[i4][1];
    d24z = x[i2][2] - x[i4][2];
    domain->minimum_image(d24x,d24y,d24z);

    // 1-4 distance
    d14x = x[i1][0] - x[i4][0];
    d14y = x[i1][1] - x[i4][1];
    d14z = x[i1][2] - x[i4][2];
    domain->minimum_image(d14x,d14y,d14z);
    
    // calculate normals
    n1x = d21y*d31z - d31y*d21z;
    n1y = d31x*d21z - d21x*d31z;
    n1z = d21x*d31y - d31x*d21y;
    n2x = d34y*d24z - d24y*d34z;
    n2y = d24x*d34z - d34x*d24z;
    n2z = d34x*d24y - d24x*d34y;
    n1 = n1x*n1x + n1y*n1y + n1z*n1z;
    n2 = n2x*n2x + n2y*n2y + n2z*n2z;
    nn = sqrt(n1*n2);  

    // cos(theta) and sin(theta) calculation 
    costheta = (n1x*n2x + n1y*n2y + n1z*n2z)/nn; 
    if (costheta > 1.0) costheta = 1.0;
    if (costheta < -1.0) costheta = -1.0;
    sintheta = sqrt(1.0-costheta*costheta); 
    if (sintheta < SMALL) sintheta = SMALL;
    mx = (n1x-n2x)*d14x + (n1y-n2y)*d14y + (n1z-n2z)*d14z;
    if (mx < 0)
      sintheta = -sintheta;
     
    // coeffs calculation
    double relax = 1.0;
		if (t_remod[type] > 0.0) relax = cell_remodel(i1, i2, i3, i4, t_remod[type]);
    alfa = relax*k[type]*(cos(theta0[type])-costheta*sin(theta0[type])/sintheta);   
    a11 = -alfa*costheta/n1;
    a12 = alfa/nn;
    a22 = -alfa*costheta/n2;

    // forces calculation
    s1x = a11*(n1y*d32z - n1z*d32y) + a12*(n2y*d32z - n2z*d32y);
    s1y = a11*(n1z*d32x - n1x*d32z) + a12*(n2z*d32x - n2x*d32z);
    s1z = a11*(n1x*d32y - n1y*d32x) + a12*(n2x*d32y - n2y*d32x);
    s2x = a11*(n1z*d31y - n1y*d31z) + a22*(n2y*d34z - n2z*d34y) +  
          a12*(n2z*d31y - n2y*d31z + n1y*d34z - n1z*d34y);
    s2y = a11*(n1x*d31z - n1z*d31x) + a22*(n2z*d34x - n2x*d34z) +  
          a12*(n2x*d31z - n2z*d31x + n1z*d34x - n1x*d34z);
    s2z = a11*(n1y*d31x - n1x*d31y) + a22*(n2x*d34y - n2y*d34x) +  
          a12*(n2y*d31x - n2x*d31y + n1x*d34y - n1y*d34x);
    s3x = a11*(n1y*d21z - n1z*d21y) + a22*(n2z*d24y - n2y*d24z) +  
          a12*(n2y*d21z - n2z*d21y + n1z*d24y - n1y*d24z);     
    s3y = a11*(n1z*d21x - n1x*d21z) + a22*(n2x*d24z - n2z*d24x) +  
          a12*(n2z*d21x - n2x*d21z + n1x*d24z - n1z*d24x);
    s3z = a11*(n1x*d21y - n1y*d21x) + a22*(n2y*d24x - n2x*d24y) +  
          a12*(n2x*d21y - n2y*d21x + n1y*d24x - n1x*d24y);
    s4x = a22*(n2z*d32y - n2y*d32z) + a12*(n1z*d32y - n1y*d32z);
    s4y = a22*(n2x*d32z - n2z*d32x) + a12*(n1x*d32z - n1z*d32x);
    s4z = a22*(n2y*d32x - n2x*d32y) + a12*(n1y*d32x - n1x*d32y);

    if (eflag){
      mx = costheta*cos(theta0[type]) + sintheta*sin(theta0[type]);
      edihedral = relax*k[type]*(1.0-mx);
      energy += relax*k[type]*(1.0-mx);
    }
    
    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f1[0] = s1x;
      f1[1] = s1y;
      f1[2] = s1z;
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f2[0] = s2x;
      f2[1] = s2y;
      f2[2] = s2z;
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f3[0] = s3x;
      f3[1] = s3y;
      f3[2] = s3z;
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f4[0] = s4x;
      f4[1] = s4y;
      f4[2] = s4z;
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (stress_ind == 2){
      ff[0] = -d31x*s1x - d32x*s2x - d34x*s4x;
      ff[1] = -d31y*s1y - d32y*s2y - d34y*s4y;
      ff[2] = -d31z*s1z - d32z*s2z - d34z*s4z;
      ff[3] = -d31x*s1y - d32x*s2y - d34x*s4y;
      ff[4] = -d31x*s1z - d32x*s2z - d34x*s4z;
      ff[5] = -d31y*s1z - d32y*s2z - d34y*s4z;
      
      for(int k = 0; k < n_stress_tot; ++k)
        modify->fix[stress_list[k]]->virial5(i1, i2, i3, i4, ff);
    }

    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,f1,f3,f4,
	       -d21x,-d21y,-d21z,d32x,d32y,d32z,-d34x,-d34y,-d34z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBend::allocate()
{
  allocated = 1;

  int n = atom->ndihedraltypes;

  k = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:k");
  theta0 = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:theta0");
  t_remod = (double *) memory->smalloc((n+1)*sizeof(double),"dihedral:t_remod");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;

 }

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralBend::coeff(int which, int narg, char **arg)
{
  if (which != 0) error->all("Illegal coeffs for this dihedral style");
  if (narg != 4) error->all("Incorrect args in dihedral_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->ndihedraltypes,ilo,ihi);

  double k_one = atof(arg[1]);
  double theta0_one = atof(arg[2]);
  double t_remod_one = atof(arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one*M_PI/180.0;
    t_remod[i] = t_remod_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args in dihedral_coeff command");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void DihedralBend::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&theta0[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&t_remod[1],sizeof(double),atom->ndihedraltypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void DihedralBend::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&theta0[1],sizeof(double),atom->ndihedraltypes,fp);
    fread(&t_remod[1],sizeof(double),atom->ndihedraltypes,fp);
  }
  MPI_Bcast(&k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta0[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_remod[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------*/

double DihedralBend::cell_remodel(int i1, int i2, int i3, int i4, double tau)
{
  if (atom->molecular == 0) return 1.0;
  if (molecule == NULL){
    int ifix;
    for (ifix = 0; ifix < modify->nfix; ifix++)
      if (strcmp("molecule",modify->fix[ifix]->style) == 0) break;
  	if (ifix == modify->nfix) error->one("DihedralBend: FixMolecule style is not defined.");
 	 	molecule = dynamic_cast<FixMolecule *> (modify->fix[ifix]);
  }

  double tact = 0.0;
  if (molecule)
    tact = molecule->list_tact[ atom->molecule[i1] ];
  if (tact == 0.0) return 1.0;

  double coeff = exp( -(update->dt*update->ntimestep-tact) / tau );
  coeff = (coeff > 0.01) ? coeff : 0.01;

  return coeff;
}

/* ----------------------------------------------------------------------*/
