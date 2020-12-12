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
#include "fix_stat_dens_cyl.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "math.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatDensCyl::FixStatDensCyl(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(n1,n2,n3,"fix_stat_dens:num");
  nmass = memory->create_3d_double_array(n1,n2,n3,"fix_stat_dens:nmass");
  for (i=0;i<n1;i++)
    for (j=0;j<n2;j++)
      for (k=0;k<n3;k++){
        num[i][j][k] = 0.0;
        nmass[i][j][k] = 0.0;
      }
}

/* ---------------------------------------------------------------------- */

FixStatDensCyl::~FixStatDensCyl()
{
  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(nmass);
}

/* ---------------------------------------------------------------------- */

int FixStatDensCyl::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatDensCyl::end_of_step()
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
        if (map_index_cyl(x[l][0],x[l][1],x[l][2])){
          num[is][js][ks] += 1.0;
          nmass[is][js][ks] += mass[type[l]];
        }
    num_step++;
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatDensCyl:: write_stat(int step)
{
  int i, j, k, l;
  int total = n1*n2*n3;
  double x, y, z, rr, theta;
  double *ntmp, *mtmp, *tmp;
  double vol;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:ntmp");
  mtmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:mtmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:tmp");

  l = 0;
  for (i=0; i<n1; i++)
    for (k=0; k<n3; k++)
      for (j=0; j<n2; j++){
        tmp[l]=num[i][j][k]/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        mtmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (i=0; i<n1; i++)
    for (k=0; k<n3; k++)
      for (j=0; j<n2; j++){
        tmp[l]=nmass[i][j][k]/num_step;
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
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", n2, n3, n1);
    l = 0;
    for (i=0; i<n1; i++)
      for (k=0; k<n3; k++)
        for (j=0; j<n2; j++){
          x = low + (i+0.5)*(high-low)/n1;
          if (strcmp(stype,"cyl") == 0){
	          vol = (high-low)/n1*M_PI/n3*(2.0*j+1)*R*R/n2/n2;  
            rr = (j+0.5)*R/n2;
          }else if (strcmp(stype,"sten") == 0){
          	vol = (high-low)/n1*M_PI/n3*(2.0*j+1)*profile[i]*profile[i]/n2/n2;  
            rr = (j+0.5)*profile[i]/n2;
					}
          theta = (k+0.5)*2*M_PI/n3; 
          y = cent2 + rr*cos(theta);
          z = cent3 + rr*sin(theta);
          fprintf(out_stat,"%lf %lf %lf %15.10lf %15.10lf %lf \n",x, y, z, ntmp[l]/vol, mtmp[l]/vol, ntmp[l]);
          l++;
	}
    fclose(out_stat);
  }

  memory->sfree(ntmp);
  memory->sfree(mtmp);
  memory->sfree(tmp);
}

/* ---------------------------------------------------------------------- */
