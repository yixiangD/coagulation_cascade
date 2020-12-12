/* ----------------------------------------------------------------------
  Dmitry Fedosov - 08/12/05  accumulation of statistics
 
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

------------------------------------------------------------------------- */
#include "mpi.h"
#include "fix_stat_gyroradius.h"
#include "atom.h"
#include "memory.h"
#include "comm.h"
#include "update.h"
#include "domain.h" 

using namespace LAMMPS_NS;
using namespace std;
/* ---------------------------------------------------------------------- */
Fix_Stat_Gyroradius::Fix_Stat_Gyroradius(LAMMPS *lmp, int narg, char **arg)
  :FixStat(lmp, narg,arg)
{
  init_on = 0;
  num_step = 0;
}

/* ---------------------------------------------------------------------- */

Fix_Stat_Gyroradius::~Fix_Stat_Gyroradius()
{
  if (init_on){
    memory->sfree(rad);
    memory->sfree(c_m);
    memory->sfree(c_mt);
  }
}

/* ---------------------------------------------------------------------- */

int Fix_Stat_Gyroradius::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


/* ---------------------------------------------------------------------- */

void Fix_Stat_Gyroradius::end_of_step()
{
  int i,j;
  double xx[3], mss;
  double **x = atom->x;
  double *mass = atom->mass;
  int *type = atom->type;
  int *n_atoms = atom->atoms_in_mol;
  int nlocal = atom->nlocal;
  int t_step = update->ntimestep;
  int mol;
  
  if (init_on == 0){
    nm = atom->n_mol;
    rad = (double *) memory->smalloc(nm*sizeof(double),"stat_gyroradius:rad");
    c_m = (double *) memory->smalloc(4*nm*sizeof(double),"stat_gyroradius:c_m");
    c_mt = (double *) memory->smalloc(4*nm*sizeof(double),"stat_gyroradius:c_mt");
    for(i=0; i<nm; i++)
      rad[i] = 0.0;   

    init_on = 1;   
  }
  if(t_step > st_start) { 
  for(i=0; i<4*nm; i++) {   
    c_m[i] = 0.0;
    c_mt[i] = 0.0;
  }

  for (i=0; i<nlocal; i++){
    mol = atom->molecule[i];     
    if (mol){
      mss = mass[type[i]];
      for (j=0; j<3; j++) xx[j] = x[i][j];
      domain->unmap(xx,atom->image[i]);
      for (j=0; j<3; j++) c_m[4*(mol-1)+j] += mss*xx[j];
      c_m[4*(mol-1)+3] += mss;
    }
  }  
  /*
  for (i=0; i<nm; i++)
    for (j=0; j<3; j++) 
      c_m[3*i+j] /= n_atoms[i];
   */
  MPI_Allreduce(c_m,c_mt,4*nm,MPI_DOUBLE,MPI_SUM,world);

  for (i=0; i<nm; i++)
    for (j=0; j<3; j++)
      c_mt[4*i+j] /= c_mt[4*i+3]; 
  for (i=0; i<nlocal; i++){
    mol = atom->molecule[i];   
    if (mol){
      for (j=0; j<3; j++) xx[j] = x[i][j];
      domain->unmap(xx,atom->image[i]);
      for (j=0; j<3; j++) 
        rad[mol-1] += (xx[j]-c_mt[4*(mol-1)+j])*(xx[j]-c_mt[4*(mol-1)+j]);
    }
  }
  num_step++;

  if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void Fix_Stat_Gyroradius:: write_stat(int step)
{
  int i;
  double *rtmp, *tmp;
  double radd = 0.0;
  char f_name[FILENAME_MAX];

  rtmp = (double *) memory->smalloc(nm*sizeof(double),"stat_gyroradius:rtmp");
  tmp = (double *) memory->smalloc(nm*sizeof(double),"stat_gyroradius:tmp");

  for (i=0; i<nm; i++){
    tmp[i] = rad[i]/num_step/atom->atoms_in_mol[i];
    rad[i] = 0.0;
    rtmp[i] = 0.0; 
  }
  MPI_Reduce(tmp,rtmp,nm,MPI_DOUBLE,MPI_SUM,0,world); 
  num_step = 0;

  if (!(comm->me)){
    for (i=0; i<nm; i++)
      radd += rtmp[i];
    radd /= nm;
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step);
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"Average radius of gyration squared = %lf \n",radd);
    fprintf(out_stat,"For each molecule = \n");
    for (i=0; i<nm; i++)
      fprintf(out_stat,"%d %lf \n",i+1,rtmp[i]);    
    fclose(out_stat);
  }    
  memory->sfree(rtmp);
  memory->sfree(tmp);        
}

/* ---------------------------------------------------------------------- */

