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
#include "fix_stat_vel.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatVel::FixStatVel(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(nx,ny,nz,"fix_stat_vel:num");
  vx = memory->create_3d_double_array(nx,ny,nz,"fix_stat_vel:vx");
  vy = memory->create_3d_double_array(nx,ny,nz,"fix_stat_vel:vy");
  vz = memory->create_3d_double_array(nx,ny,nz,"fix_stat_vel:vz");
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++){
        num[i][j][k] = 0.0;
        vx[i][j][k] = 0.0;
        vy[i][j][k] = 0.0;
        vz[i][j][k] = 0.0;
      }
}

/* ---------------------------------------------------------------------- */

FixStatVel::~FixStatVel()
{
  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(vx);
  memory->destroy_3d_double_array(vy);
  memory->destroy_3d_double_array(vz);
}

/* ---------------------------------------------------------------------- */

int FixStatVel::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatVel::end_of_step()
{
  int l;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  int t_step = update->ntimestep;
  if (t_step > st_start){
//    if((dump_each - (t_step%dump_each)) < 50000)
//    {
    for (l=0; l<atom->nlocal; l++)
      if (mask[l] & groupbit)
        if (map_index(x[l][0],x[l][1],x[l][2])){
          num[is][js][ks] += 1.0;
          vx[is][js][ks] += v[l][0];
          vy[is][js][ks] += v[l][1];
          vz[is][js][ks] += v[l][2];
        }
    num_step++;
//    }
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatVel:: write_stat(int step)
{
  int i, j, k, l;
  double x, y, z;
  int total = nx*ny*nz;
  double *ntmp, *vxtmp, *vytmp, *vztmp, *tmp;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_vel:ntmp");
  vxtmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_vel:vxtmp");
  vytmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_vel:vytmp");
  vztmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_vel:vztmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_vel:tmp");

  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=num[i][j][k]/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        vxtmp[l] = 0.0;
        vytmp[l] = 0.0;
        vztmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=vx[i][j][k]/num_step;
        vx[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vxtmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=vy[i][j][k]/num_step;
        vy[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vytmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=vz[i][j][k]/num_step;
        vz[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vztmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  num_step = 0;

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step);
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"v_x\",\"v_y\",\"v_z\" \n");
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);
    l = 0;
    for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          if (ntmp[l] > 0.0)
            fprintf(out_stat,"%lf %lf %lf %lf %lf %lf \n",x, y, z, vxtmp[l]/ntmp[l], vytmp[l]/ntmp[l], vztmp[l]/ntmp[l]);
          else
            fprintf(out_stat,"%lf %lf %lf 0.0 0.0 0.0 \n",x, y, z);
          l++;
        }
    fclose(out_stat);
  }
  memory->sfree(ntmp);
  memory->sfree(vxtmp);
  memory->sfree(vytmp);
  memory->sfree(vztmp);
  memory->sfree(tmp);
}
