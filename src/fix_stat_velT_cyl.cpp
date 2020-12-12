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
#include "fix_stat_velT_cyl.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatVelTCyl::FixStatVelTCyl(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(n1,n2,n3,"fix_stat_velT:num");
  vx = memory->create_3d_double_array(n1,n2,n3,"fix_stat_velT:vx");
  vy = memory->create_3d_double_array(n1,n2,n3,"fix_stat_velT:vy");
  vz = memory->create_3d_double_array(n1,n2,n3,"fix_stat_velT:vz");
  aT = memory->create_3d_double_array(n1,n2,n3,"fix_stat_velT:aT");
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      for (k=0; k<n3; k++){
        num[i][j][k] = 0.0;
        vx[i][j][k] = 0.0;
        vy[i][j][k] = 0.0;
        vz[i][j][k] = 0.0;
				aT[i][j][k] = 0.0;
      }
}

/* ---------------------------------------------------------------------- */

FixStatVelTCyl::~FixStatVelTCyl()
{
  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(vx);
  memory->destroy_3d_double_array(vy);
  memory->destroy_3d_double_array(vz);
  memory->destroy_3d_double_array(aT);
}

/* ---------------------------------------------------------------------- */

int FixStatVelTCyl::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatVelTCyl::end_of_step()
{
  int l;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **T = atom->T;
  int t_step = update->ntimestep;
  if (t_step > st_start){
//    if((dump_each - (t_step%dump_each)) < 60000)
//    {
    for (l=0; l<atom->nlocal; l++)
      if (mask[l] & groupbit)
        if (map_index_cyl(x[l][0],x[l][1],x[l][2])){
          num[is][js][ks] += 1.0;
          vx[is][js][ks] += v[l][0];
          vy[is][js][ks] += v[l][1];
          vz[is][js][ks] += v[l][2];
				  aT[is][js][ks] += T[l][index-1];
        }
    num_step++;
//    }
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatVelTCyl:: write_stat(int step)
{
  int i, j, k, l;
	double x, y, z, rr, theta;
  int total = n1*n2*n3;
  double *ntmp, *vxtmp, *vytmp, *vztmp, *aTtmp, *tmp;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:ntmp");
  vxtmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:vxtmp");
  vytmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:vytmp");
  vztmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:vztmp");
  aTtmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:aTtmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_velT:tmp");

  l = 0;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++){
        tmp[l]=num[i][j][k]/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        vxtmp[l] = 0.0;
        vytmp[l] = 0.0;
        vztmp[l] = 0.0;
				aTtmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++){
        tmp[l]=vx[i][j][k]/num_step;
        vx[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vxtmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++){
        tmp[l]=vy[i][j][k]/num_step;
        vy[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vytmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++){
        tmp[l]=vz[i][j][k]/num_step;
        vz[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,vztmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<n3; k++)
    for (j=0; j<n2; j++)
      for (i=0; i<n1; i++){
        tmp[l]=aT[i][j][k]/num_step;
        aT[i][j][k] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,aTtmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  num_step = 0;

  if (!(comm->me)){
    FILE* out_cc;
    sprintf(f_name,"%s_%02d.%d.plt",fname,index,step);
    out_cc=fopen(f_name,"w");
    fprintf(out_cc,"VARIABLES=\"x\",\"y\",\"z\",\"v_x\",\"v_y\",\"v_z\",\"T\" \n");
    fprintf(out_cc,"ZONE I=%d,J=%d,K=%d, F=POINT \n", n2, n3, n1);
    l = 0;
    for (i=0; i<n1; i++)
      for (k=0; k<n3; k++)
        for (j=0; j<n2; j++){
          x = low + (i+0.5)*(high-low)/n1;
          if (strcmp(stype,"cyl") == 0)
            rr = (j+0.5)*R/n2;
          else if (strcmp(stype,"sten") == 0)
            rr = (j+0.5)*profile[i]/n2;
          theta = (k+0.5)*2*M_PI/n3;
          y = cent2 + rr*cos(theta);
          z = cent3 + rr*sin(theta);
          if (ntmp[l] > 0.0)
            fprintf(out_cc,"%lf %lf %lf %lf %lf %lf %lf \n",x, y, z, vxtmp[l]/ntmp[l], vytmp[l]/ntmp[l], vztmp[l]/ntmp[l], aTtmp[l]/ntmp[l]);
          else
            fprintf(out_cc,"%lf %lf %lf 0.0 0.0 0.0 0.0\n",x, y, z);
          l++;
        }
    fclose(out_cc);
  }
  memory->sfree(ntmp);
  memory->sfree(vxtmp);
  memory->sfree(vytmp);
  memory->sfree(vztmp);
  memory->sfree(aTtmp);
  memory->sfree(tmp);
}
