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
#include "fix_stat_diff_poly.h"
#include "atom.h"
#include "memory.h"
#include "comm.h"
#include "update.h"
#include "domain.h" 

using namespace LAMMPS_NS;
using namespace std;
/* ---------------------------------------------------------------------- */
Fix_Stat_Diff_Poly::Fix_Stat_Diff_Poly(LAMMPS *lmp, int narg, char **arg)
  :FixStat(lmp, narg,arg)
{
  int i;
  
  index = 0;
  init_on = 0;
  num_tot = 0;
  //num_tot = static_cast<int> (10.0/update->dt + 1.0e-6);
  num_tot = dump_each - num_tot;

  displ = (double *) memory->smalloc(num_tot*sizeof(double),"stat_diffusion_poly:displ");
  displ_tot = (double *) memory->smalloc(num_tot*sizeof(double),"stat_diffusion_poly:displ_tot");
  for (i=0; i<num_tot; i++){
    displ[i] = 0.0;
    displ_tot[i] = 0.0;
  }
  num_c = 0;
  num_step = 0;
  n_skip = (st_start + dump_each - num_tot)%nevery;
  n_skip = nevery - n_skip - 1;
}

/* ---------------------------------------------------------------------- */

Fix_Stat_Diff_Poly::~Fix_Stat_Diff_Poly()
{
  memory->sfree(displ);
  memory->sfree(displ_tot);

  if (init_on){
    memory->sfree(c0);
    memory->sfree(c_m);
    memory->sfree(c_mt);
  }
}

/* ---------------------------------------------------------------------- */

int Fix_Stat_Diff_Poly::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}


/* ---------------------------------------------------------------------- */

void Fix_Stat_Diff_Poly::end_of_step()
{
  int i,j,mol, count;
  double **x = atom->x;
  double *mass = atom->mass;
  int *type = atom->type;
  int *n_atoms = atom->atoms_in_mol;
  int t_step = update->ntimestep;
  int nlocal = atom->nlocal;
  double xx[3],mss; 
 
  if (init_on == 0){
    nm = atom->n_mol;
    c0 = (double *) memory->smalloc(3*nm*sizeof(double),"stat_diffusion_poly:c0");
    c_m = (double *) memory->smalloc(4*nm*sizeof(double),"stat_diffusion_poly:c_m");
    c_mt = (double *) memory->smalloc(4*nm*sizeof(double),"stat_diffusion_poly:c_mt");
    init_on = 1;
  }

  if(t_step > st_start) {
  
  count = (int)(t_step-st_start-1)%dump_each;
  
  if (count > dump_each - num_tot-1 || count == 0){
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
        c_m[4*i+j] /= n_atoms[i];
   */
    MPI_Reduce(c_m,c_mt,4*nm,MPI_DOUBLE,MPI_SUM,0,world);


    if (comm->me == 0){
      for (i=0; i<nm; i++)
        for (j=0; j<3; j++)
          c_mt[4*i+j] /= c_mt[4*i+3]; 
      if (index == 0){ 
        for (i=0; i<nm; i++)
          for (j=0; j<3; j++)
            c0[3*i+j] = c_mt[4*i+j];
        index = 1;
      }     


      if (count > dump_each - num_tot-1){    
        count = count - dump_each + num_tot;
        for (i=0; i<nm; i++)
          for (j=0; j<3; j++)
            displ[count] += (c0[3*i+j]-c_mt[4*i+j])*(c0[3*i+j]-c_mt[4*i+j]);    
      }
    }
  }
      
  if ((t_step-st_start)%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void Fix_Stat_Diff_Poly:: write_stat(int step)
{
  int i;
  double time;
  double dtt = update->dt;
  char f_name[FILENAME_MAX];

  index = 0;
  num_c++;

  if (!(comm->me)){
    int step_tot = 1 + static_cast<int>((num_tot - 1 - n_skip)/nevery);
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step); 
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"t\",\"displ\"  \n");
    fprintf(out_stat,"ZONE I=%d, F=POINT \n", step_tot);

    for (i=n_skip; i<num_tot; i=i+nevery){
      displ_tot[i] += displ[i]/nm;
      time = dtt*(i + dump_each - num_tot - n_skip);
      fprintf(out_stat,"%lf %lf  \n",time, displ_tot[i]/num_c);
      displ[i] = 0.0;
    }
    fclose(out_stat);
  }        
}

/* ---------------------------------------------------------------------- */

