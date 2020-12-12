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
#include "fix_stat_dens.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatDens::FixStatDens(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(nx,ny,nz,"fix_stat_dens:num");
  nmass = memory->create_3d_double_array(nx,ny,nz,"fix_stat_dens:nmass");
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++){
        num[i][j][k] = 0.0;
        nmass[i][j][k] = 0.0;
      }
}

/* ---------------------------------------------------------------------- */

FixStatDens::~FixStatDens()
{
  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(nmass);
}

/* ---------------------------------------------------------------------- */

int FixStatDens::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatDens::end_of_step()
{
  int l;
  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;
  double *mass = atom->mass;
  int t_step = update->ntimestep;

  if (t_step > st_start){
    for (l=0; l<atom->nlocal; l++)
      if (mask[l] & groupbit)
        if (map_index(x[l][0],x[l][1],x[l][2])){
          num[is][js][ks] += 1.0;
          nmass[is][js][ks] += mass[type[l]];
        }
    num_step++;
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatDens:: write_stat(int step)
{
  int i, j, k, l;
  int total = nx*ny*nz;
  double x, y, z;
  double *ntmp, *mtmp, *tmp;
  double vol = xs*ys*zs/total;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:ntmp");
  mtmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:mtmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:tmp");

  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=num[i][j][k]/vol/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        mtmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=nmass[i][j][k]/vol/num_step;
        nmass[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,mtmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  num_step = 0;

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step);
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"number density\",\"density\",\"particles number\" \n");
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);
    l = 0;
    for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          fprintf(out_stat,"%lf %lf %lf %lf %lf %lf \n",x, y, z, ntmp[l], mtmp[l], ntmp[l]*vol);
          l++;
        }
    fclose(out_stat);
  }
  memory->sfree(ntmp);
  memory->sfree(mtmp);
  memory->sfree(tmp);
}

/* ---------------------------------------------------------------------- */
