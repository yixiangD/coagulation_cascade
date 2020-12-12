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
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "atom.h"
#include "domain.h"
#include "fix_com_gyration.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include <string>
#include <sys/stat.h>

using namespace LAMMPS_NS;
using namespace std;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define BIG 2147483647

/* ---------------------------------------------------------------------- */

FixComGyration::FixComGyration(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8) error->all("Illegal fix com/gyration command");

  int arg_index = 3;
  nevery = atoi(arg[arg_index++]);
  comflag = atoi(arg[arg_index++]);	// flag for computing center of mass and gyration tensor of the cells
	countflag = atoi(arg[arg_index++]);	// flag for computing number of transient bonds
	posflag = atoi(arg[arg_index++]);	// flag for computing minimum cell distance to the walls
	sprintf(fname,arg[arg_index++]);
	if (comflag)
	  tensorflag = atoi(arg[arg_index++]);
	if (countflag){
		btype = atoi(arg[arg_index++]);	// bondtype to monitor
		sprintf(grp,arg[arg_index++]);
		int igroup = group->find(grp);
    if (igroup == -1) error->one("Group ID does not exist");
    groupbitLig = group->bitmask[igroup];	// finding group for ligands
	}
	sten_50 = sten_75 = 0;	// vessels with or without stenosis
	if (posflag){
		if (posflag == 2) sten_50 = 1;
		if (posflag == 3) sten_75 = 1;
	}
  if (posflag)
    if (sten_75){
      R = 27.2, xl = 60., xh = 105., H = 30., Yc = -15.2;
      xc = (xl+xh)/2., yc = H/2.;
    } else if (sten_50){
      R = 34.8, xl = 60., xh = 105., H = 30., Yc = -26.6;
      xc = (xl+xh)/2., yc = H/2.;
    }

  if (atom->molecular == 0)
    error->all("Fix com/gyration requires molecular atom style");
}

/* ---------------------------------------------------------------------- */

FixComGyration::~FixComGyration()
{
	if (comflag){
	  memory->sfree(massproc);
	  memory->sfree(masstotal);
  	memory->destroy_2d_double_array(com);
  	memory->destroy_2d_double_array(comall);
		if (tensorflag){
  		memory->destroy_2d_double_array(rgt);
			memory->destroy_2d_double_array(array);
		} else{
			memory->sfree(rg);
			memory->sfree(vector);
		}
	}
	if (posflag){
  	memory->sfree(Ymin);
  	memory->sfree(Ymax);
  	memory->sfree(Yminmax);
	}
	if (countflag){
  	memory->sfree(bondcount);
  	memory->sfree(nbond);
  	memory->sfree(nbondtot);
	}
	for (int i = 0; i < nmolecules; i++)
		 if (!comm->me) fclose(fp[i]);
	if (!comm->me) fclose(fp0);
}

/* ---------------------------------------------------------------------- */

int FixComGyration::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixComGyration::init()
{

  // setup molecule-based data
  nmolecules = molecules_in_group(idlo,idhi);

		if (comflag){
  massproc = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:massproc");
  masstotal = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:masstotal");
  com = memory->create_2d_double_array(nmolecules,3,"com/gyration:com");
  comall = memory->create_2d_double_array(nmolecules,3,"com/gyration:comall");

  if (tensorflag){
    rgt = memory->create_2d_double_array(nmolecules,6,"com/gyration:rgt");
    array = memory->create_2d_double_array(nmolecules,6,"com/gyration:array");
    size_array_rows = nmolecules;
    size_array_cols = 6;
  } else{
    rg = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:rg");
    vector = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:vector");
    size_vector = nmolecules;
  }
		}

	if (posflag){
  	Ymin = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:Ymin");
  	Ymax = (double *) memory->smalloc(nmolecules*sizeof(double),"com/gyration:Ymax");
  	Yminmax = (double *) memory->smalloc(2*nmolecules*sizeof(double),"com/gyration:Yminmax");
	}

	if (countflag){
  	bondcount = (int *) memory->smalloc(atom->nmax*sizeof(int),"com/gyration:bondcount");
	  nbond = (int *) memory->smalloc(nmolecules*sizeof(int),"com/gyration:nbond");
  	nbondtot = (int *) memory->smalloc(nmolecules*sizeof(int),"com/gyration:nbondtot");
	}
  mollist = (int *) memory->smalloc(atom->n_mol*sizeof(int),"com/gyration:mollist");
  tempi = (int *) memory->smalloc(atom->n_mol*sizeof(int),"com/gyration:tempi");

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int imol;
  double massone;

		if (comflag){
  for (int i = 0; i < nmolecules; i++) massproc[i] = 0.0;
  for (int i = 0; i < atom->n_mol; i++) tempi[i] = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
			tempi[imol-1] = imol;
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(tempi,mollist,atom->n_mol,MPI_INT,MPI_MAX,world);

	int nfile = 0;
  for (int i = 0; i < atom->n_mol; i++){
		if (mollist[i] == 0) continue;
    sprintf(f_name,"%s_mol%04d.dat",fname,mollist[i]);
    if (!comm->me) fp[nfile++] = fopen(f_name,"a+");
  //  fprintf(fp[i],"VARIABLES=\"t\",\"xcm\",\"ycm\",\"zcm\",\"Rgxx\",\"Rgyy\",\"Rgzz\",\"Rgxy\",\"Rgxz\",\"Rgyz\"\n");
  }
		}

	if (countflag)
		if (!comm->me){
	    sprintf(f_name,"%s_bondcount.dat",fname);
    	fp0 = fopen(f_name,"a+");
		}
}

/* ---------------------------------------------------------------------- */

void FixComGyration::end_of_step()
{
  int i,j,k,m;
  int imol;
  double dx,dy,dz,massone;
  double unwrap[3];
  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

	updatemol();	// update number of cells in group and reallocate memory

	if (countflag) {
  	bondcount = (int *) memory->srealloc(bondcount,nall*sizeof(int),"com/gyration:bondcount");
  	for (i = 0; i < nall; i++) bondcount[i] = 0;

  	for (i = 0; i < nlocal; i++){
			if (mask[i] & groupbit){
    		for (j = 0; j < num_bond[i]; j++) {
      		if (bond_type[i][j] == btype) {
        bondcount[i]++;
  			if (newton_bond) {
    			m = atom->map(bond_atom[i][j]);
    			if (m < 0)
      			error->one("Could not count initial bonds in fix com/gyration");
    			bondcount[m]++;
  			}
      		}
    		}
			}
		}
  	// if newton_bond is set, need to communicate ghost counts
  	if (newton_bond) comm->reverse_comm_fix(this);
	}

	if (posflag)
		for (i = 0; i < nmolecules; i++){
			Ymin[i] = (double) BIG; 
			Ymax[i] = 0.0;
		}

	if (countflag)
		for (i = 0; i < nmolecules; i++)
			nbond[i] = 0;

	if (comflag){
		molcom();
		if (tensorflag == 0)
  		for (i = 0; i < nmolecules; i++) rg[i] = 0.0;
		else
    	for (i = 0; i < nmolecules; i++)
      	for (j = 0; j < 6; j++)
        	rgt[i][j] = 0.0;
	}

  for (i = 0; i < nlocal; i++)
   	if (mask[i] & groupbit){
     	imol = molecule[i];
     	if (molmap) imol = molmap[imol-idlo];
     	else imol--;
     	domain->unmap(x[i],image[i],unwrap);

				if (comflag){
			if (tensorflag == 0){
      	dx = unwrap[0] - comall[imol][0];
      	dy = unwrap[1] - comall[imol][1];
      	dz = unwrap[2] - comall[imol][2];
      	if (rmass) massone = rmass[i];
      	else massone = mass[type[i]];
      	rg[imol] += (dx*dx + dy*dy + dz*dz) * massone;
			} else{
        dx = unwrap[0] - comall[imol][0];
        dy = unwrap[1] - comall[imol][1];
        dz = unwrap[2] - comall[imol][2];
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        rgt[imol][0] += dx*dx * massone;
        rgt[imol][1] += dy*dy * massone;
        rgt[imol][2] += dz*dz * massone;
        rgt[imol][3] += dx*dy * massone;
        rgt[imol][4] += dx*dz * massone;
        rgt[imol][5] += dy*dz * massone;
			}
				}

			// finding minimum distance to the wall and number of created bonds for platelets
			if (posflag){
				if (sten_75 || sten_50){
					if (x[i][1] < yc && x[i][0] >= xl && x[i][0] <= xh)
						yw = Yc + sqrt(R*R - (x[i][0]-xc)*(x[i][0]-xc));
					else if (x[i][1] > yc && x[i][0] >= xl && x[i][0] <= xh)
						yw = fabs(Yc) + H - sqrt(R*R - (x[i][0]-xc)*(x[i][0]-xc));
					else
						yw = x[i][1];
					if( fabs(x[i][1]-yw) != 0.0 && fabs(x[i][1]-yw) < Ymin[imol] ) Ymin[imol] = fabs(x[i][1]-yw);
					if( fabs(x[i][1]-yw) != 0.0 && fabs(x[i][1]-yw) > Ymax[imol] ) Ymax[imol] = fabs(x[i][1]-yw);
				} else {
					if (unwrap[1] < Ymin[imol]) Ymin[imol] = unwrap[1];
					if (unwrap[1] > Ymax[imol]) Ymax[imol] = unwrap[1];
				}
			}

			if (countflag)
				nbond[imol] += bondcount[i];
   	}

	if (comflag)
		if (tensorflag == 0){
	  	MPI_Allreduce(rg,vector,nmolecules,MPI_DOUBLE,MPI_SUM,world);
 			for (i = 0; i < nmolecules; i++)
   			vector[i] = sqrt(vector[i]/masstotal[i]);
		} else{
   		MPI_Allreduce(&rgt[0][0],&array[0][0],nmolecules*6,
   	  		          MPI_DOUBLE,MPI_SUM,world);
  		for (i = 0; i < nmolecules; i++)
 	  		for (j = 0; j < 6; j++)
   	  		array[i][j] /= masstotal[i];
		}

	if (posflag){
   	MPI_Allreduce(&Ymin[0],&Yminmax[0],nmolecules,
   		            MPI_DOUBLE,MPI_MIN,world);
    MPI_Allreduce(&Ymax[0],&Yminmax[nmolecules],nmolecules,
 	   	            MPI_DOUBLE,MPI_MAX,world);
	}

	if (countflag)
  	MPI_Allreduce(nbond,nbondtot,nmolecules,MPI_INT,MPI_SUM,world);

	write_stat();	// writing the stat files
}

/* ----------------------------------------------------------------------
   calculate per-molecule COM
------------------------------------------------------------------------- */

void FixComGyration::molcom()
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

void FixComGyration:: write_stat()
{
  int i, j;
	double time = update->dt*update->ntimestep;

  if (!comm->me)
		for (i = 0; i < nmolecules; i++){
			if ( (comflag && tensorflag) && posflag && countflag )
	      fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %d\n",
								time, comall[i][0], comall[i][1], comall[i][2], array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5],
								Yminmax[i], Yminmax[nmolecules+i], nbondtot[i]);
			else if ( (comflag && tensorflag) && posflag )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n",
                time, comall[i][0], comall[i][1], comall[i][2], array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5],
                Yminmax[i], Yminmax[nmolecules+i]);
			else if ( (comflag && tensorflag) && countflag )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %d\n",
                time, comall[i][0], comall[i][1], comall[i][2], array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5],
                nbondtot[i]);
      else if ( (comflag && tensorflag) )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n",
                time, comall[i][0], comall[i][1], comall[i][2], array[i][0], array[i][1], array[i][2], array[i][3], array[i][4], array[i][5]);
      else if ( (comflag && !tensorflag) )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n", time, comall[i][0], comall[i][1], comall[i][2], vector[i]);
      else if ( !comflag && posflag && countflag )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf %d\n",time, Yminmax[i], Yminmax[nmolecules+i], nbondtot[i]);
      else if ( !comflag && posflag )
        fprintf(fp[i],"%9.4lf %9.4lf %9.4lf\n",time, Yminmax[i], Yminmax[nmolecules+i]);
      else if ( !comflag && countflag )
        fprintf(fp[i],"%9.4lf %d\n",time, nbondtot[i]);
			fflush(fp[i]);
		}
	if (!comm->me)
		if (countflag){
			fprintf(fp0,"%9.4lf %d\n", time, atom->nbonds);
			fflush(fp0);
		}
}

/* ----------------------------------------------------------------------
   identify molecule IDs with atoms in group
   warn if any atom in group has molecule ID = 0
   warn if any molecule has only some atoms in group
   return Ncount = # of molecules with atoms in group
   set molmap to NULL if molecule IDs include all in range from 1 to Ncount
   else: molecule IDs range from idlo to idhi
         set molmap to vector of length idhi-idlo+1
         molmap[id-idlo] = index from 0 to Ncount-1
         return idlo and idhi
------------------------------------------------------------------------- */

int FixComGyration::molecules_in_group(int &idlo, int &idhi)
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

  molmap = (int *) memory->smalloc(nlen*sizeof(int),"com/gyration:molmap");
  for (i = 0; i < nlen; i++) molmap[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molmap[molecule[i]-idlo] = 1;

  int *molmapall;
  molmapall = (int *) memory->smalloc(nlen*sizeof(int),"com/gyration:molmapall");
  MPI_Allreduce(molmap,molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmols = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmols = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap[i] = nmols++;
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

  if (idlo == 1 && idhi == nmols && nlen == nmols){
    memory->sfree(molmap);
    molmap = NULL;
  }
  return nmols;
}

/* ---------------------------------------------------------------------- */

int FixComGyration::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++)
    buf[m++] = bondcount[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixComGyration::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    bondcount[j] += static_cast<int> (buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

void FixComGyration::updatemol()
{

  // setup molecule-based data
  int nmolnew = molecules_in_group(idlo,idhi);

	if (nmolnew == nmolecules) return;

	nmolecules = nmolnew;

		if (comflag){
  massproc = (double *) memory->srealloc(massproc,nmolecules*sizeof(double),"com/gyration:massproc");
  masstotal = (double *) memory->srealloc(masstotal,nmolecules*sizeof(double),"com/gyration:masstotal");
  memory->destroy_2d_double_array(com);
  memory->destroy_2d_double_array(comall);
  com = memory->create_2d_double_array(nmolecules,3,"com/gyration:com");
  comall = memory->create_2d_double_array(nmolecules,3,"com/gyration:comall");

  if (tensorflag){
  	memory->destroy_2d_double_array(rgt);
  	memory->destroy_2d_double_array(array);
    rgt = memory->create_2d_double_array(nmolecules,6,"com/gyration:rgt");
    array = memory->create_2d_double_array(nmolecules,6,"com/gyration:array");
    size_array_rows = nmolecules;
    size_array_cols = 6;
  } else{
    rg = (double *) memory->srealloc(rg,nmolecules*sizeof(double),"com/gyration:rg");
    vector = (double *) memory->srealloc(vector,nmolecules*sizeof(double),"com/gyration:vector");
    size_vector = nmolecules;
  }
		}

  if (posflag){
    Ymin = (double *) memory->srealloc(Ymin,nmolecules*sizeof(double),"com/gyration:Ymin");
    Ymax = (double *) memory->srealloc(Ymax,nmolecules*sizeof(double),"com/gyration:Ymax");
    Yminmax = (double *) memory->srealloc(Yminmax,2*nmolecules*sizeof(double),"com/gyration:Yminmax");
  }

  if (countflag){
    nbond = (int *) memory->srealloc(nbond,nmolecules*sizeof(int),"com/gyration:nbond");
    nbondtot = (int *) memory->srealloc(nbondtot,nmolecules*sizeof(int),"com/gyration:nbondtot");
  }

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int imol;
  double massone;

		if (comflag){
  for (int i = 0; i < nmolecules; i++) massproc[i] = 0.0;
  for (int i = 0; i < atom->n_mol; i++) tempi[i] = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
			tempi[imol-1] = imol;
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(tempi,mollist,atom->n_mol,MPI_INT,MPI_MAX,world);

	int nfile = 0;
  for (int i = 0; i < atom->n_mol; i++){
		if (mollist[i] == 0) continue;
    sprintf(f_name,"%s_mol%04d.dat",fname,mollist[i]);
		//if ( exists(f_name) ) { nfile++; continue; }
    if (!comm->me) fp[nfile++] = fopen(f_name,"a+");
  }
		}
}

/* ---------------------------------------------------------------------- */

bool FixComGyration::exists(const std::string& name)
{
  struct stat buffer;   
  return ( stat(name.c_str(), &buffer) == 0 ); 
}

/* ---------------------------------------------------------------------- */
