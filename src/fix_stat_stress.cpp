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
#include "fix_stat_stress.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "modify.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatStress::FixStatStress(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i,j,k,l;

  poly_ind = 0;
  if (force->bond || force->angle || force->dihedral) poly_ind = 1;  

  ss = memory->create_4d_double_array(nx,ny,nz,6,"stat_stress:ss");
  vv = memory->create_4d_double_array(nx,ny,nz,6,"stat_stress:vv");
  if (poly_ind){
    ss_p = memory->create_4d_double_array(nx,ny,nz,6,"stat_stress:ss_p");
    vv_p = memory->create_4d_double_array(nx,ny,nz,6,"stat_stress:vv_p");
  }

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++)
        for (l=0; l<6; l++){
          ss[i][j][k][l] = 0.0;
          vv[i][j][k][l] = 0.0;
          if (poly_ind){
            ss_p[i][j][k][l] = 0.0;
            vv_p[i][j][k][l] = 0.0;
	  }
	}
  num_step = 0;
  force->stress_ind = 1;
 }

/* ---------------------------------------------------------------------- */

FixStatStress::~FixStatStress()
{
  memory->destroy_4d_double_array(ss);
  memory->destroy_4d_double_array(vv);
  if (poly_ind){
    memory->destroy_4d_double_array(ss_p);
    memory->destroy_4d_double_array(vv_p);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::init()
{
  int i;
  for(i = 0; i < modify->nfix; ++i)
    if(strcmp(id, modify->fix[i]->id) == 0)
    {
      force->stress_list[force->n_stress_tot] = i;
      ++force->n_stress_tot;
    }
}

/* ---------------------------------------------------------------------- */

int FixStatStress::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatStress::initial_integrate(int vflag)
{
  int t_step = update->ntimestep;
  if(t_step > st_start)
    force->stress_ind = 2;  
}

/* ---------------------------------------------------------------------- */

void FixStatStress::end_of_step()
{
  int i,j,k,l;
  int *type = atom->type;
  int *mask = atom->mask;
  int *mol = atom->molecule;
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  double h_mass; 
  double v_avg[nx][ny][nz][4], v_avg_tmp[nx][ny][nz][4];  
  int nlocal = atom->nlocal;
  int t_step = update->ntimestep;

  if (t_step > st_start){
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

  for (l=0; l<nlocal; l++)
    if (mask[l] & groupbit)
      if (map_index(x[l][0],x[l][1],x[l][2])){
        h_mass = mass[type[l]];
        vv[is][js][ks][0] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][0]-v_avg[is][js][ks][0]);
        vv[is][js][ks][1] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][1]-v_avg[is][js][ks][1]);
        vv[is][js][ks][2] += h_mass*(v[l][2]-v_avg[is][js][ks][2])*(v[l][2]-v_avg[is][js][ks][2]);
        vv[is][js][ks][3] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][1]-v_avg[is][js][ks][1]);
        vv[is][js][ks][4] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][2]-v_avg[is][js][ks][2]);
        vv[is][js][ks][5] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][2]-v_avg[is][js][ks][2]);
        if (poly_ind && mol[l]){
          vv_p[is][js][ks][0] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][0]-v_avg[is][js][ks][0]);
          vv_p[is][js][ks][1] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][1]-v_avg[is][js][ks][1]);
          vv_p[is][js][ks][2] += h_mass*(v[l][2]-v_avg[is][js][ks][2])*(v[l][2]-v_avg[is][js][ks][2]);
          vv_p[is][js][ks][3] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][1]-v_avg[is][js][ks][1]);
          vv_p[is][js][ks][4] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][2]-v_avg[is][js][ks][2]);
          vv_p[is][js][ks][5] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][2]-v_avg[is][js][ks][2]);
	}
      }
  num_step++;

  if (t_step%dump_each == 0) write_stat(t_step);
 
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial1(int ii)
{
  map1 = map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]);

  if (map1 && (atom->mask[ii] & groupbit)){
    is1 = is; js1 = js; ks1 = ks;
  }
  else
    map1 = 0; 

}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial2(int ii, double ff[6])
{
  int i;  

  if (map1){
    for (i=0; i<6; i++) 
      ss[is1][js1][ks1][i] += 0.5*ff[i]; 
  }

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++) 
      ss[is][js][ks][i] += 0.5*ff[i];     
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial3(int ii, int jj, double ff[6])
{
  int i;  

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.5*ff[i];     
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.5*ff[i];     
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial4(int ii, int jj, int kk, double ff[6])
{
  int i;  

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += ff[i]/3.0;     
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += ff[i]/3.0;     
  }

  if (map_index(atom->x[kk][0],atom->x[kk][1],atom->x[kk][2]) && (atom->mask[kk] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += ff[i]/3.0;     
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial5(int ii, int jj, int kk, int ll, double ff[6])
{
  int i;  

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.25*ff[i];     
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.25*ff[i];     
  }

  if (map_index(atom->x[kk][0],atom->x[kk][1],atom->x[kk][2]) && (atom->mask[kk] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.25*ff[i];     
  }

  if (map_index(atom->x[ll][0],atom->x[ll][1],atom->x[ll][2]) && (atom->mask[ll] & groupbit)){
    for (i=0; i<6; i++) 
      ss_p[is][js][ks][i] += 0.25*ff[i];     
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress::virial_f(int ii, double ff[6])
{
  int i;

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss[is][js][ks][i] += ff[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixStatStress:: write_stat(int step)
{
  int i,j,k,l;
  double stmp[nx][ny][nz][6], sptmp[nx][ny][nz][6], tmp[nx][ny][nz][6]; 
  double vol = xs*ys*zs/nx/ny/nz;
  double x,y,z;  
  char f_name[FILENAME_MAX];

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++)
        for (l=0; l<6; l++){
          if (poly_ind)
            tmp[i][j][k][l] = -(ss[i][j][k][l] + ss_p[i][j][k][l] + vv[i][j][k][l])/num_step/vol;
          else 
            tmp[i][j][k][l] = -(ss[i][j][k][l] + vv[i][j][k][l])/num_step/vol; 
          ss[i][j][k][l] = 0.0;
          vv[i][j][k][l] = 0.0;
          stmp[i][j][k][l] = 0.0;
          sptmp[i][j][k][l] = 0.0; 
	}

  MPI_Reduce(&tmp,&stmp,nx*ny*nz*6,MPI_DOUBLE,MPI_SUM,0,world);

  if (poly_ind){
    for (i=0; i<nx; i++)
      for (j=0; j<ny; j++)
        for (k=0; k<nz; k++)
          for (l=0; l<6; l++){
            tmp[i][j][k][l] = -(ss_p[i][j][k][l] + vv_p[i][j][k][l])/num_step/vol;
            ss_p[i][j][k][l] = 0.0;
            vv_p[i][j][k][l] = 0.0;
	  }

    MPI_Reduce(&tmp,&sptmp,nx*ny*nz*6,MPI_DOUBLE,MPI_SUM,0,world);
  }    
  num_step = 0;

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step); 
    out_stat=fopen(f_name,"w");
    if (poly_ind)
      fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"Sxx\",\"Syy\",\"Szz\",\"Sxy\",\"Sxz\",\"Syz\",\"Press\",\"SPxx\",\"SPyy\",\"SPzz\",\"SPxy\",\"SPxz\",\"SPyz\",\"Press_p\"  \n");
    else
      fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"Sxx\",\"Syy\",\"Szz\",\"Sxy\",\"Sxz\",\"Syz\",\"Press\"  \n");

    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);

    for (k=0; k<nz; k++)
      for (j=0; j<ny; j++)
        for (i=0; i<nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          if (poly_ind)
            fprintf(out_stat,"%lf %lf %lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf \n",x, y, z,stmp[i][j][k][0],stmp[i][j][k][1],stmp[i][j][k][2],stmp[i][j][k][3],stmp[i][j][k][4],stmp[i][j][k][5],-(stmp[i][j][k][0] + stmp[i][j][k][1] + stmp[i][j][k][2])/3.0,sptmp[i][j][k][0],sptmp[i][j][k][1],sptmp[i][j][k][2],sptmp[i][j][k][3],sptmp[i][j][k][4],sptmp[i][j][k][5],-(sptmp[i][j][k][0] + sptmp[i][j][k][1] + sptmp[i][j][k][2])/3.0);
          else
            fprintf(out_stat,"%lf %lf %lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf %15.10lf \n",x, y, z,stmp[i][j][k][0],stmp[i][j][k][1],stmp[i][j][k][2],stmp[i][j][k][3],stmp[i][j][k][4],stmp[i][j][k][5],-(stmp[i][j][k][0] + stmp[i][j][k][1] + stmp[i][j][k][2])/3.0);            
	}
    fclose(out_stat);
    } 
}

/* ---------------------------------------------------------------------- */
