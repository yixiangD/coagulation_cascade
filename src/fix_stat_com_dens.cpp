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
#include "fix_stat_com_dens.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define BIG 2147483647

/* ---------------------------------------------------------------------- */

FixStatComDens::FixStatComDens(LAMMPS *lmp, int narg, char **arg) :
  FixStat(lmp, narg, arg)
{
  if (atom->molecular == 0)
    error->all("Fix stat/com/dens requires molecular atom style");
}

/* ---------------------------------------------------------------------- */

FixStatComDens::~FixStatComDens()
{
  memory->destroy_3d_double_array(num);
  memory->sfree(massproc);
  memory->sfree(masstotal);
  memory->destroy_2d_double_array(com);
  memory->destroy_2d_double_array(comall);
}

/* ---------------------------------------------------------------------- */

int FixStatComDens::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixStatComDens::init()
{
  int i, j, k;
  num = memory->create_3d_double_array(nx,ny,nz,"stat/com/dens:num");
  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
        num[i][j][k] = 0.0;

  // setup molecule-based data
  nmolecules = molecules_in_group(idlo,idhi);

  massproc = (double *) memory->smalloc(nmolecules*sizeof(double),"stat/com/dens:massproc");
  masstotal = (double *) memory->smalloc(nmolecules*sizeof(double),"stat/com/dens:masstotal");
  com = memory->create_2d_double_array(nmolecules,3,"stat/com/dens:com");
  comall = memory->create_2d_double_array(nmolecules,3,"stat/com/dens:comall");

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int imol;
  double massone;

  for (int i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixStatComDens::end_of_step()
{
  int l;
  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;
  int t_step = update->ntimestep;

	molcom();

  if (t_step > st_start){
    for (l = 0; l < nmolecules; l++)
      if (map_index(comall[l][0],comall[l][1],comall[l][2])){
        num[is][js][ks] += 1.0;
      }
    num_step++;
    if (t_step%dump_each == 0) write_stat(t_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixStatComDens:: write_stat(int step)
{
  int i, j, k, l;
  int total = nx*ny*nz;
  double x, y, z;
  double *ntmp, *mtmp, *tmp;
  double vol = xs*ys*zs/total;
  char f_name[FILENAME_MAX];

  ntmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:ntmp");
  tmp = (double *) memory->smalloc(total*sizeof(double),"fix_stat_dens:tmp");

  l = 0;
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++){
        tmp[l] = num[i][j][k]/vol/num_step;
        num[i][j][k] = 0.0;
        ntmp[l] = 0.0;
        l++;
      }
  MPI_Reduce(tmp,ntmp,total,MPI_DOUBLE,MPI_SUM,0,world);
  num_step = 0;

  if (!(comm->me)){
    FILE* out_stat;
    sprintf(f_name,"%s.%d.plt",fname,step);
    out_stat=fopen(f_name,"w");
    fprintf(out_stat,"VARIABLES=\"x\",\"y\",\"z\",\"COM-density\" \n");
    fprintf(out_stat,"ZONE I=%d,J=%d,K=%d, F=POINT \n", nx, ny, nz);
    l = 0;
    for (k = 0; k < nz; k++)
      for (j = 0; j < ny; j++)
        for (i = 0; i < nx; i++){
          x = xlo + (i+0.5)*xs/nx;
          y = ylo + (j+0.5)*ys/ny;
          z = zlo + (k+0.5)*zs/nz;
          fprintf(out_stat,"%lf %lf %lf %lf \n",x, y, z, ntmp[l]);
          l++;
        }
    fclose(out_stat);
  }
  memory->sfree(ntmp);
  memory->sfree(tmp);
}

/* ----------------------------------------------------------------------
   calculate per-molecule COM
------------------------------------------------------------------------- */

void FixStatComDens::molcom()
{
  int imol;
  double dx,dy,dz,massone;
  double unwrap[3];

  for (int i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      com[imol][0] += unwrap[0] * massone;
      com[imol][1] += unwrap[1] * massone;
      com[imol][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  for (int i = 0; i < nmolecules; i++){
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }
}

/* ---------------------------------------------------------------------- */
int FixStatComDens::molecules_in_group(int &idlo, int &idhi)
{
  int i;

  // find lo/hi molecule ID for any atom in group
  // warn if atom in group has ID = 0

  int *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int lo = BIG;
  int hi = -BIG;
  int flag = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if (molecule[i] == 0) flag = 1;
      lo = MIN(lo,molecule[i]);
      hi = MAX(hi,molecule[i]);
    }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning("Atom with molecule ID = 0 included in "
                   "compute molecule group");

  MPI_Allreduce(&lo,&idlo,1,MPI_INT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_INT,MPI_MAX,world);
  if (idlo == BIG) return 0;

  // molmap = vector of length nlen
  // set to 1 for IDs that appear in group across all procs, else 0

  int nlen_tag = idhi-idlo+1;
  if (nlen_tag > BIG)
    error->all("Too many molecules for compute");
  int nlen = (int) nlen_tag;

  molmap = (int *) memory->smalloc(nlen*sizeof(int),"stat/com/dens:molmap");
  for (i = 0; i < nlen; i++) molmap[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molmap[molecule[i]-idlo] = 1;

  int *molmapall;
  molmapall = (int *) memory->smalloc(nlen*sizeof(int),"stat/com/dens:molmapall");
  MPI_Allreduce(molmap,molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap[i] = nmolecules++;
    else molmap[i] = -1;
  memory->sfree(molmapall);

  // warn if any molecule has some atoms in group and some not in group

  flag = 0;
  for (i = 0; i < nlocal; i++){
    if (mask[i] & groupbit) continue;
    if (molecule[i] < idlo || molecule[i] > idhi) continue;
    if (molmap[molecule[i]-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning("One or more compute molecules has atoms not in group");

  // if molmap simply stores 1 to Nmolecules, then free it

  if (idlo == 1 && idhi == nmolecules && nlen == nmolecules){
    memory->sfree(molmap);
    molmap = NULL;
  }
  return nmolecules;
}
