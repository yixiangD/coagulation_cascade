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
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "fix_stat_snap.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatSnap::FixStatSnap(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(nx,ny,nz,"fix_stat_snap:num");
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
        num[i][j][k] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixStatSnap::~FixStatSnap()
{
  memory->destroy_3d_double_array(num);
}

/* ---------------------------------------------------------------------- */

int FixStatSnap::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatSnap::end_of_step()
{
  int l;
  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;
	double *q = atom->q;
  int t_step = update->ntimestep;

  if (t_step > st_start){
    for (l=0; l<atom->nlocal; l++)
      if (mask[l] & groupbit)
        if (map_index(x[l][0],x[l][1],x[l][2]) && q[l] == 1.0)
          num[is][js][ks] = 1.0;
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatSnap:: write_stat(int step)
{
  int i, j, k, l;
  int total = nx*ny*nz;
  double x, y, z;
  double *ntmp, *tmp;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_snap:ntmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_snap:tmp");

  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=num[i][j][k];
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step);
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"Binary\" \n");
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);
    l = 0;
    for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          fprintf(out_stat,"%lf %lf %lf %lf \n",x, y, z, ntmp[l]);
          l++;
        }
    fclose(out_stat);
  }
  memory->sfree(ntmp);
  memory->sfree(tmp);
}

/* ---------------------------------------------------------------------- */
