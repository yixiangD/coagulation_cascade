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
#include "fix_stat_all.h"
#include "atom.h"
#include "modify.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "common.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStatAll::FixStatAll(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  int i, j, k, l;
  num = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:num");
  numtmp = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:numtmp");
  tmp = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:tmp");
	vel = new double***[3];
	veltmp = new double***[3];
	for (i = 0; i < 3; i++){
  	vel[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:vel");
  	veltmp[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:veltmp");
	}
	aT = new double***[CTYPES];
	aTtmp = new double***[CTYPES];
	for (i = 0; i < CTYPES; i++){
  	aT[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:aT");
  	aTtmp[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:aTtmp");
	}
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++){
        num[i][j][k] = 0.0;
        numtmp[i][j][k] = 0.0;
				for (l=0; l<3; l++){
	        vel[l][i][j][k] = 0.0;
  	      veltmp[l][i][j][k] = 0.0;
				}
				for (l=0; l<CTYPES; l++){
					aT[l][i][j][k] = 0.0;
					aTtmp[l][i][j][k] = 0.0;
				}
      }

  poly_ind = 0;
  if (force->bond || force->angle || force->dihedral) poly_ind = 1;

	stress = new double***[6];
	pstress = new double***[6];
	ss = new double***[6];
	vv = new double***[6];
  for (i = 0; i < 6; i++){
  	stress[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:stress");
  	pstress[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:pstress");
  	ss[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:ss");
  	vv[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:vv");
	}
  if (poly_ind){
		ss_p = new double***[6];
		vv_p = new double***[6];
  	for (i = 0; i < 6; i++){
    	ss_p[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:ss_p");
    	vv_p[i] = memory->create_3d_double_array(nx,ny,nz,"fix_stat_all:vv_p");
	  }
	}

  for (l=0; l<6; l++)
  	for (i=0; i<nx; i++)
    	for (j=0; j<ny; j++)
      	for (k=0; k<nz; k++){
          ss[l][i][j][k] = 0.0;
          vv[l][i][j][k] = 0.0;
					stress[l][i][j][k] = 0.0;
					pstress[l][i][j][k] = 0.0;
          if (poly_ind){
            ss_p[l][i][j][k] = 0.0;
            vv_p[l][i][j][k] = 0.0;
    			}
  			}

	v_avg = memory->create_4d_double_array(nx,ny,nz,4,"fix_stat_all:v_avg");
	v_avg_tmp = memory->create_4d_double_array(nx,ny,nz,4,"fix_stat_all:v_avg_tmp");

  force->stress_ind = 1;
}

/* ---------------------------------------------------------------------- */

FixStatAll::~FixStatAll()
{
	int i;

  memory->destroy_3d_double_array(num);
  memory->destroy_3d_double_array(numtmp);
  memory->destroy_3d_double_array(tmp);
  for (i = 0; i < 3; i++){
	  memory->destroy_3d_double_array(vel[i]);
  	memory->destroy_3d_double_array(veltmp[i]);
	}
	delete[] vel;
	delete[] veltmp;
  for (i = 0; i < CTYPES; i++){
  	memory->destroy_3d_double_array(aT[i]);
  	memory->destroy_3d_double_array(aTtmp[i]);
	}
	delete[] aT;
	delete[] aTtmp;
  for (i = 0; i < 6; i++){
  	memory->destroy_3d_double_array(ss[i]);
  	memory->destroy_3d_double_array(vv[i]);
  	memory->destroy_3d_double_array(stress[i]);
  	memory->destroy_3d_double_array(pstress[i]);
  	if (poly_ind){
    	memory->destroy_3d_double_array(ss_p[i]);
    	memory->destroy_3d_double_array(vv_p[i]);
  	}
	}
	delete[] ss;
	delete[] vv;
	delete[] stress;
	delete[] pstress;
	if (poly_ind){
		delete[] ss_p;
		delete[] vv_p;
	}
  memory->destroy_4d_double_array(v_avg);
  memory->destroy_4d_double_array(v_avg_tmp);
}

/* ---------------------------------------------------------------------- */

void FixStatAll::init()
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

int FixStatAll::setmask()
{
  int mask = 0;
	mask |= INITIAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatAll::initial_integrate(int vflag)
{
  int t_step = update->ntimestep;
  if(t_step > st_start)
    force->stress_ind = 2;
}

/* ---------------------------------------------------------------------- */

void FixStatAll::end_of_step()
{
  int i, j, k, l;
  int *mask = atom->mask;
  int *type = atom->type;
  int *mol = atom->molecule;
  double *mass = atom->mass;
  double **x = atom->x;
  double **v = atom->v;
  double **T = atom->T;
  double h_mass;
  int t_step = update->ntimestep;
//	double v_avg[nx][ny][nz][4], v_avg_tmp[nx][ny][nz][4];

  if (t_step > st_start){
    for (i=0; i<nx; i++)
     for (j=0; j<ny; j++)
       for (k=0; k<nz; k++)
         for (l=0; l<4; l++){
           v_avg[i][j][k][l] = 0.0;
           v_avg_tmp[i][j][k][l] = 0.0;
         }

    for (l=0; l<atom->nlocal; l++)
      if (mask[l] & groupbit)
        if (map_index(x[l][0],x[l][1],x[l][2])){
          numtmp[is][js][ks] += 1.0;
					for (i=0; i<3; i++)
          	veltmp[i][is][js][ks] += v[l][i];
					for (i=0; i<CTYPES; i++)
					  aTtmp[i][is][js][ks] += T[l][i];
					for (i=0; i<3; i++)
        		v_avg_tmp[is][js][ks][i] += v[l][i];
        	v_avg_tmp[is][js][ks][3] += 1.0;
        }

  	MPI_Allreduce(&v_avg_tmp[0][0][0][0],&v_avg[0][0][0][0],nx*ny*nz*4,MPI_DOUBLE,MPI_SUM,world);

	  for (i=0; i<nx; i++)
  	  for (j=0; j<ny; j++)
    	  for (k=0; k<nz; k++)
      	  for (l=0; l<3; l++)
        	  v_avg[i][j][k][l] /= v_avg[i][j][k][3];

	  for (l=0; l<atom->nlocal; l++)
  	  if (mask[l] & groupbit)
    	  if (map_index(x[l][0],x[l][1],x[l][2])){
      	  h_mass = mass[type[l]];
        	vv[0][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][0]-v_avg[is][js][ks][0]);
        	vv[1][is][js][ks] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][1]-v_avg[is][js][ks][1]);
        	vv[2][is][js][ks] += h_mass*(v[l][2]-v_avg[is][js][ks][2])*(v[l][2]-v_avg[is][js][ks][2]);
        	vv[3][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][1]-v_avg[is][js][ks][1]);
        	vv[4][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][2]-v_avg[is][js][ks][2]);
        	vv[5][is][js][ks] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][2]-v_avg[is][js][ks][2]);
        	if (poly_ind && mol[l]){
          	vv_p[0][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][0]-v_avg[is][js][ks][0]);
          	vv_p[1][is][js][ks] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][1]-v_avg[is][js][ks][1]);
          	vv_p[2][is][js][ks] += h_mass*(v[l][2]-v_avg[is][js][ks][2])*(v[l][2]-v_avg[is][js][ks][2]);
          	vv_p[3][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][1]-v_avg[is][js][ks][1]);
          	vv_p[4][is][js][ks] += h_mass*(v[l][0]-v_avg[is][js][ks][0])*(v[l][2]-v_avg[is][js][ks][2]);
          	vv_p[5][is][js][ks] += h_mass*(v[l][1]-v_avg[is][js][ks][1])*(v[l][2]-v_avg[is][js][ks][2]);
  				}
      	}

    num_step++;
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatAll:: write_stat(int step)
{
  int i, j, k, l;
	int total = nx*ny*nz;
	double vol = xs*ys*zs/(double) total;

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      for (k=0; k<nz; k++){
        tmp[i][j][k] = numtmp[i][j][k]/(double) num_step;
        numtmp[i][j][k] = 0.0;
			}

  MPI_Allreduce(&tmp[0][0][0],&num[0][0][0],total,MPI_DOUBLE,MPI_SUM,world);

	for (l=0; l<3; l++){
  	for (i=0; i<nx; i++)
    	for (j=0; j<ny; j++)
      	for (k=0; k<nz; k++){
					if (num[i][j][k] != 0)
		      	tmp[i][j][k] = veltmp[l][i][j][k]/num[i][j][k]/(double) num_step;
					else
						tmp[i][j][k] = 0.0;
					veltmp[l][i][j][k] = 0.0;
				}

		MPI_Allreduce(&tmp[0][0][0],&vel[l][0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
	}

	for (l=0; l<CTYPES; l++){
	  for (i=0; i<nx; i++)
  	  for (j=0; j<ny; j++)
    	  for (k=0; k<nz; k++){
					if (num[i][j][k] != 0)
	      	  tmp[i][j][k] = aTtmp[l][i][j][k]/num[i][j][k]/(double) num_step;
					else
						tmp[i][j][k] = 0.0;
					aTtmp[l][i][j][k] = 0.0;
				}

		MPI_Allreduce(&tmp[0][0][0],&aT[l][0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
	}

  for (l=0; l<6; l++){
  	for (i=0; i<nx; i++)
    	for (j=0; j<ny; j++)
      	for (k=0; k<nz; k++){
          if (poly_ind)
            tmp[i][j][k] = -(ss[l][i][j][k] + ss_p[l][i][j][k] + vv[l][i][j][k])/vol/(double) num_step;
          else
            tmp[i][j][k] = -(ss[l][i][j][k] + vv[l][i][j][k])/vol/(double) num_step;
          ss[l][i][j][k] = 0.0;
          vv[l][i][j][k] = 0.0;
  			}

		MPI_Allreduce(&tmp[0][0][0],&stress[l][0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
	}

  if (poly_ind){
    for (l=0; l<6; l++){
    	for (i=0; i<nx; i++)
      	for (j=0; j<ny; j++)
        	for (k=0; k<nz; k++){
            tmp[i][j][k] = -(ss_p[l][i][j][k] + vv_p[l][i][j][k])/vol/(double) num_step;
            ss_p[l][i][j][k] = 0.0;
            vv_p[l][i][j][k] = 0.0;
    			}

			MPI_Allreduce(&tmp[0][0][0],&pstress[l][0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
		}
	}

	num_step = 0;
}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial1(int ii)
{
  map1 = map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]);

  if (map1 && (atom->mask[ii] & groupbit)){
    is1 = is; js1 = js; ks1 = ks;
  }
  else
    map1 = 0;

}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial2(int ii, double ff[6])
{
  int i;

  if (map1){
    for (i=0; i<6; i++)
      ss[i][is1][js1][ks1] += 0.5*ff[i];
  }

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss[i][is][js][ks] += 0.5*ff[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial3(int ii, int jj, double ff[6])
{
  int i;

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.5*ff[i];
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.5*ff[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial4(int ii, int jj, int kk, double ff[6])
{
  int i;

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += ff[i]/3.0;
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += ff[i]/3.0;
  }

  if (map_index(atom->x[kk][0],atom->x[kk][1],atom->x[kk][2]) && (atom->mask[kk] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += ff[i]/3.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial5(int ii, int jj, int kk, int ll, double ff[6])
{
  int i;

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.25*ff[i];
  }

  if (map_index(atom->x[jj][0],atom->x[jj][1],atom->x[jj][2]) && (atom->mask[jj] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.25*ff[i];
  }

  if (map_index(atom->x[kk][0],atom->x[kk][1],atom->x[kk][2]) && (atom->mask[kk] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.25*ff[i];
  }

  if (map_index(atom->x[ll][0],atom->x[ll][1],atom->x[ll][2]) && (atom->mask[ll] & groupbit)){
    for (i=0; i<6; i++)
      ss_p[i][is][js][ks] += 0.25*ff[i];
  }
}

/* ---------------------------------------------------------------------- */

void FixStatAll::virial_f(int ii, double ff[6])
{
  int i;

  if (map_index(atom->x[ii][0],atom->x[ii][1],atom->x[ii][2]) && (atom->mask[ii] & groupbit)){
    for (i=0; i<6; i++)
      ss[i][is][js][ks] += ff[i];
  }
}
