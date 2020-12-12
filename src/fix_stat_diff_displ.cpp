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
#include "fix_stat_diff_displ.h"
#include "atom.h"
#include "memory.h"
#include "comm.h"
#include "update.h"
#include "domain.h" 
#include <iostream>

using namespace LAMMPS_NS;
using namespace std;
/* ---------------------------------------------------------------------- */
Fix_Stat_Diff_Displ::Fix_Stat_Diff_Displ(LAMMPS *lmp, int narg, char **arg)
  :FixStat(lmp, narg,arg)
{
  
  int i,rcount;
  natoms = static_cast<int>(atom->natoms);
  index = 0;
  num_tot = 0;                             // change of a beginning cut
  //num_tot = static_cast<int> (10.0/update->dt + 1.0e-6);
  num_tot = dump_each - num_tot;

  displ = (double *) memory->smalloc(num_tot*sizeof(double),"fix_stat_diff_displ:displ");
  displ_tot = (double *) memory->smalloc(num_tot*sizeof(double),"fix_stat_diff_displ:displ_tot");
  x_init = (double *) memory->smalloc(3*natoms*sizeof(double), "fix_stat_diff_displ:x_init");
  for (i=0;i<num_tot;i++){
    displ[i] = 0.0;
    displ_tot[i] = 0.0;
  }
  num_c = 0;

  tot = 0;
  rcount = 0;
  for (i=0; i<atom->nlocal; i++)
    if (atom->mask[i] & groupbit)
      rcount++;

  MPI_Allreduce(&rcount,&tot,1,MPI_INT,MPI_SUM,world);  
  num_step = 0;
  n_skip = (st_start + dump_each - num_tot)%nevery;
  n_skip = nevery - n_skip - 1;
}

/* ---------------------------------------------------------------------- */

Fix_Stat_Diff_Displ::~Fix_Stat_Diff_Displ()
{
  memory->sfree(displ);
  memory->sfree(displ_tot);
  memory->sfree(x_init);
}

/* ---------------------------------------------------------------------- */

int Fix_Stat_Diff_Displ::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


/* ---------------------------------------------------------------------- */

void Fix_Stat_Diff_Displ::end_of_step()
{
  int l, count, k;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;
  double xx[3];
  int *tag = atom->tag;
  int t_step = update->ntimestep;
  int id;
  if (t_step > st_start){
    if (index == 0){
      double *rbuff;
      rbuff = (double *) memory->smalloc(3*natoms*sizeof(double), "fix_stat_diff_displ:rbuff");
      for(l = 0; l < 3*natoms; ++l)
        rbuff[l] = 0;

      for (l=0; l<nlocal; l++)
        if (mask[l] & groupbit){
          xx[0] = x[l][0];  xx[1] = x[l][1];  xx[2] = x[l][2];
          domain->unmap(xx,atom->image[l]);
          id = 3*(tag[l]-1);
          rbuff[id] = xx[0];
          rbuff[id+1] = xx[1];
          rbuff[id+2] = xx[2];
        }
      index = 1;
      MPI_Allreduce(rbuff,x_init,3*natoms,MPI_DOUBLE,MPI_SUM,world);
      memory->sfree(rbuff);
    }     
    count = (int)(t_step-st_start-1)%dump_each;
    if (count > dump_each - num_tot - 1){
      count = count - dump_each + num_tot;
      for (l=0; l<nlocal; l++) 
        if (mask[l] & groupbit){
          id = 3*(tag[l]-1);
          xx[0] = x[l][0];  xx[1] = x[l][1];  xx[2] = x[l][2]; 
          domain->unmap(xx,atom->image[l]);
          displ[count] += (x_init[id]-xx[0])*(x_init[id]-xx[0]) + (x_init[id+1]-xx[1])*(x_init[id+1]-xx[1]) + (x_init[id+2]-xx[2])*(x_init[id+2]-xx[2]);
        }
    }
    if ((t_step-st_start)%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void Fix_Stat_Diff_Displ:: write_stat(int step)
{
  int i;
  double *dtmp, *tmp;
  double time;
  double dtt = update->dt;
  char f_name[FILENAME_MAX];

  dtmp = (double *) memory->smalloc(num_tot*sizeof(double),"stat_diffusion_displ:dtmp");
  tmp = (double *) memory->smalloc(num_tot*sizeof(double),"stat_diffusion_displ:tmp");

  index = 0;
  num_c++;
  for (i=0; i<num_tot; ++i){
    tmp[i]=displ[i]/tot;
    displ[i] = 0.0;
    dtmp[i] = 0.0;
  }
  MPI_Reduce(tmp,dtmp,num_tot,MPI_DOUBLE,MPI_SUM,0,world);

  if (!(comm->me)){
    int step_tot = 1 + static_cast<int>((num_tot - 1 - n_skip)/nevery);
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step); 
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"t\",\"displ\"  \n");
    fprintf(out_stat,"ZONE I=%d, F=POINT \n", step_tot);

    for (i=n_skip; i<num_tot; i=i+nevery){
      displ_tot[i] += dtmp[i];
      time = dtt*(i + dump_each - num_tot - n_skip);
      fprintf(out_stat,"%lf %20.16lf \n",time, displ_tot[i]/num_c);
    }
    fclose(out_stat);
  }    
  memory->sfree(dtmp);
  memory->sfree(tmp);    
}

/* ---------------------------------------------------------------------- */

