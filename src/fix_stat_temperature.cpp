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
#include "fix_stat_temperature.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatTemp::FixStatTemp(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k;
  num = memory->create_3d_double_array(nx,ny,nz,"stat_temperature:num");
  tx = memory->create_3d_double_array(nx,ny,nz,"stat_temperature:tx"); 
  ty = memory->create_3d_double_array(nx,ny,nz,"stat_temperature:ty"); 
  tz = memory->create_3d_double_array(nx,ny,nz,"stat_temperature:tz"); 
  for (i=0;i<nx;i++)
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++){
        num[i][j][k] = 0.0;
        tx[i][j][k] = 0.0;
        ty[i][j][k] = 0.0;
        tz[i][j][k] = 0.0;
      }

  num_step = 0;
}

/* ---------------------------------------------------------------------- */

FixStatTemp::~FixStatTemp()
{
  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(tx);
  memory->destroy_3d_double_array(ty);
  memory->destroy_3d_double_array(tz);
}

/* ---------------------------------------------------------------------- */

int FixStatTemp::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatTemp::end_of_step()
{
  int t_step = update->ntimestep;

  int i,j,k,l;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  double h_mass; 
  double v_avg[nx][ny][nz][4], v_avg_tmp[nx][ny][nz][4];  
  int nlocal = atom->nlocal;
  
  if(t_step > st_start) {
  
    for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
      for (k=0; k<nz; k++)
        for (l=0; l<4; l++){
          v_avg[i][j][k][l] = 0.0;
          v_avg_tmp[i][j][k][l] = 0.0;
	}

  for (l=0; l<nlocal; l++)
    if (mask[l] & groupbit)
      if (map_index(x[l][0],x[l][1],x[l][2])){
        v_avg_tmp[is][js][ks][0] += v[l][0];
        v_avg_tmp[is][js][ks][1] += v[l][1];
        v_avg_tmp[is][js][ks][2] += v[l][2];
        v_avg_tmp[is][js][ks][3] += 1.0;         
      }

  MPI_Allreduce(&v_avg_tmp,&v_avg,nx*ny*nz*4,MPI_DOUBLE,MPI_SUM,world);

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++)
        for (l=0; l<3; l++)
          v_avg[i][j][k][l] /= v_avg[i][j][k][3];


  for (int l=0;l<atom->nlocal;l++)
    if (mask[l] & groupbit)     
      if (map_index(x[l][0],x[l][1],x[l][2])){
        num[is][js][ks] += 1.0;
        h_mass = mass[type[l]];
        tx[is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][0]-v_avg[is][js][ks][0]);
        ty[is][js][ks] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][1]-v_avg[is][js][ks][1]);
        tz[is][js][ks] += h_mass*(v[l][2]-v_avg[is][js][ks][2])*(v[l][2]-v_avg[is][js][ks][2]);
      }
  num_step++;
  
  if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatTemp:: write_stat(int step)
{
  
  int i, j, k, l;
  double x, y, z;
  int total = nx*ny*nz;
  double *ntmp, *txtmp, *tytmp, *tztmp, *tmp;
  double tfactor = force->mvv2e/force->boltz;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"stat_temperature:ntmp");
  txtmp = (double *) memory->smalloc(total*sizeof(double),"stat_temperature:txtmp");
  tytmp = (double *) memory->smalloc(total*sizeof(double),"stat_temperature:tytmp");
  tztmp = (double *) memory->smalloc(total*sizeof(double),"stat_temperature:tztmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"stat_temperature:tmp");

  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=num[i][j][k]/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        txtmp[l] = 0.0;
        tytmp[l] = 0.0;
        tztmp[l] = 0.0;
	l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=tx[i][j][k]/num_step;
        tx[i][j][k] = 0.0;
	l++;
      }
  MPI_Reduce(tmp,txtmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=ty[i][j][k]/num_step;
        ty[i][j][k] = 0.0;
	l++;
      }
  MPI_Reduce(tmp,tytmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  l = 0;
  for (k=0; k<nz; k++)
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++){
        tmp[l]=tz[i][j][k]/num_step;
        tz[i][j][k] = 0.0;
	l++;
      }
  MPI_Reduce(tmp,tztmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  num_step = 0;

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step); 
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"T_x\",\"T_y\",\"T_z\",\"T_tot\"  \n");
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);
    l = 0;
    for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          if (ntmp[l] > 0.0)
            fprintf(out_stat,"%lf %lf %lf %lf %lf %lf %lf \n",x, y, z, tfactor*txtmp[l]/ntmp[l], tfactor*tytmp[l]/ntmp[l],tfactor*tztmp[l]/ntmp[l],tfactor*(txtmp[l]+tytmp[l]+tztmp[l])/ntmp[l]/3.0);
          else
            fprintf(out_stat,"%lf %lf %lf 0.0 0.0 0.0 0.0 \n",x, y, z);    
	  l++;
	}
    fclose(out_stat);
  }    
  memory->sfree(ntmp);
  memory->sfree(txtmp);
  memory->sfree(tytmp);
  memory->sfree(tztmp);
  memory->sfree(tmp);    

 }

/* ---------------------------------------------------------------------- */
