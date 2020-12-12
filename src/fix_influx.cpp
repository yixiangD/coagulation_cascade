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
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "modify.h"
#include "fix_influx.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "group.h"
#include "iostream"
#include "fstream"

#define MAX_GROUP 32

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

FixInflux::FixInflux(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
  	error->all("Illegal fix Influx command");
 
	nevery = 1;
	force_reneighbor = 1;
	next_reneighbor = -1;
	setupbank = 0;
	nmax = 0;

  int arg_index = 3, igroup;
	flux = atof(arg[arg_index++]);	// particle flux number per DPD time unit
	frequency = atoi(arg[arg_index++]);	// adjustbale frequency of release
  type_rbc = atoi(arg[arg_index++]);
 	type_plat = atoi(arg[arg_index++]);
 	type_bank = atoi(arg[arg_index++]);
 	cutoff = atof(arg[arg_index++]);
  igroup = group->find(arg[arg_index++]);
  if (igroup == -1) error->all("FixInflux: Group ID does not exist");
  bankgroupbit = group->inversemask[igroup];
	sprintf(&fixcomid[0],arg[arg_index++]);
	sprintf(&coord,arg[arg_index++]);
	switch (toupper(coord)){
		case 'X':
			dir_1 = 0; break;
		case 'Y':
			dir_1 = 1; break;
		case 'Z':
			dir_1 = 2; break;
		case 'R':
			dir_1 = 1; dir_2 = 2; break;
		default:
 			error->all("Illegal fix Influx command");
	}
	if (toupper(coord) == 'R'){
  	xlo = atof( arg[arg_index++] );
 		xhi = atof( arg[arg_index++] );
  	rlo = atof( arg[arg_index++] );
 		rhi = atof( arg[arg_index++] );
  	ycent = atof( arg[arg_index++] );
 		zcent = atof( arg[arg_index++] );
	} else{
	  xlo = atof( arg[arg_index++] );
 		xhi = atof( arg[arg_index++] );
  	ylo = atof( arg[arg_index++] );
 		yhi = atof( arg[arg_index++] );
 		zlo = atof( arg[arg_index++] );
 		zhi = atof( arg[arg_index++] );
	}

	// assumes flow is in x-dir
	xprd = domain->boxhi[0] - domain->boxlo[0];
	yprd = domain->boxhi[1] - domain->boxlo[1];
	zprd = domain->boxhi[2] - domain->boxlo[2];
	
}

/* ---------------------------------------------------------------------- */

FixInflux::~FixInflux()
{
	memory->sfree(list_rbc);
 	memory->sfree(list_plat);
 	memory->sfree(list_bank);
 	memory->sfree(list_mol);
 	memory->sfree(loc_mol);
 	memory->sfree(rls_mol);
}

/* ---------------------------------------------------------------------- */

int FixInflux::setmask()
{
  int mask = 0;
	mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInflux::init()
{
	list_mol = (int *) memory->smalloc(atom->n_mol*sizeof(int),"influx:list_mol");
 	loc_mol = (int *) memory->smalloc(atom->n_mol*sizeof(int),"influx:loc_mol");
 	rls_mol = (int *) memory->smalloc(atom->n_mol*sizeof(int),"influx:rls_mol");
  tempd = (double *) memory->smalloc(atom->n_mol*sizeof(double),"influx:tempd");
	chkpoint = update->ntimestep + frequency;
	nbonds0 = atom->nbonds;
}

/* ---------------------------------------------------------------------- */

void FixInflux::setup(int vflag)
{
	list_rbc = (int *) memory->smalloc(atom->nmax*sizeof(int),"influx:list_rbc");
  list_plat = (int *) memory->smalloc(atom->nmax*sizeof(int),"influx:list_plat");
  list_bank = (int *) memory->smalloc(atom->nmax*sizeof(int),"influx:list_bank");

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(fixcomid,modify->fix[ifix]->id) == 0) break;
  if (ifix == modify->nfix) error->one("FixInflux: FixComGyration style is not defined.");
  comgyrat = dynamic_cast<FixComGyration *> (modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void FixInflux::end_of_step()
{
	int i, j, k;
  double **x = atom->x;
  double **v = atom->v;
  int *type = atom->type;
  int *mol = atom->molecule;
  int *mask = atom->mask;
	int nmol = atom->n_mol;
	int nrbc, nplat, nbank, ncount, totcount;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  if (update->ntimestep != chkpoint) return;
	chkpoint = update->ntimestep + frequency;
	
//	if (atom->nbonds == nbonds0) return;
//	nbonds0 = atom->nbonds;

	int ninflux, ncross = 0;
	if (Xo == NULL){
		Nplat = comgyrat->nmolecules;
  	Xo = (double *) memory->smalloc(Nplat*sizeof(double),"influx:Xo");
		for (int i = 0; i < Nplat; i++) Xo[i] = comgyrat->comall[i][0];
	}
	for (i = 0; i < Nplat; i++){
		int dX = (int)(comgyrat->comall[i][0]/xprd) - (int)(Xo[i]/xprd);
		if (dX > 0) ncross++;
		tempd[i] = comgyrat->comall[i][0];
	}
	if (Nplat < comgyrat->nmolecules){
		for (i = Nplat; i < comgyrat->nmolecules; i++) tempd[i] = comgyrat->comall[i][0];
		Nplat = comgyrat->nmolecules;
  	Xo = (double *) memory->srealloc(Xo, Nplat*sizeof(double),"influx:Xo");
		for (i = 0; i < Nplat; i++) Xo[i] = tempd[i];
	}
	ninflux = (int) (flux*(double)frequency*update->dt - (double)ncross);
	if (ninflux < 0) ninflux = 0;
	if (ncross == 0) {frequency += 1000; chkpoint = update->ntimestep + frequency;}
	if (comm->me == 0)
		fprintf(stdout, "FixInflux: Number of inserted platelets: %d; Frequency is adjusted to: %d timesteps\n", ninflux, frequency);

	if (ninflux == 0) return;

	if (setupbank == 0){
		for (i = 0; i < nmol; i++) list_mol[i] = 0;
    for (i = 0; i < nlocal; i++)
			if (mask[i] & groupbit)
				if (type[i] == type_bank)
					list_mol[ mol[i]-1 ] = mol[i];
		setupbank = 1;
	}
#if 0
  if (atom->nmax > nmax) {
    memory->sfree(list_rbc);
    memory->sfree(list_plat);
    memory->sfree(list_bank);
    nmax = atom->nmax;
		list_rbc = (int *) memory->smalloc(nmax*sizeof(int),"influx:list_rbc");
  	list_plat = (int *) memory->smalloc(nmax*sizeof(int),"influx:list_plat");
  	list_bank = (int *) memory->smalloc(nmax*sizeof(int),"influx:list_bank");
  }
#endif

	nrbc = nplat = nbank = 0;
  for (i = 0; i < nlocal; i++)
		if (mask[i] & groupbit)
	if ( chkregion(i) ){
    if (type[i] == type_rbc)
      list_rbc[nrbc++] = i;
    else if (type[i] == type_plat)
      list_plat[nplat++] = i;
    else if (type[i] == type_bank)
      list_bank[nbank++] = i;
	}

  ncount = totcount = 0;
	for (i = 0; i < nmol; i++) loc_mol[i] = rls_mol[i] = 0;

  while (totcount < ninflux){
  	for (k = 0; k < nmol; k++){
      int chk_mol = list_mol[k];
			if (chk_mol == 0) continue;
    	for (i = 0; i < nbank; i++){
      	int ip = list_bank[i];
				if (mol[ip] != chk_mol) continue;
      	for (j = 0; j < nrbc; j++){
        	int jp = list_rbc[j];
        	double rx = x[ip][0] - x[jp][0];
        	double ry = x[ip][1] - x[jp][1];
        	double rz = x[ip][2] - x[jp][2];
					domain->minimum_image(rx,ry,rz);
        	double tmp = rx*rx + ry*ry + rz*rz;
        	if (sqrt(tmp) < cutoff) break;
					if (toupper(coord) == 'R'){
						double rip = (x[ip][dir_1]-ycent)*(x[ip][dir_1]-ycent)+(x[ip][dir_2]-zcent)*(x[ip][dir_2]-zcent);
						double rjp = (x[jp][dir_1]-ycent)*(x[jp][dir_1]-ycent)+(x[jp][dir_2]-zcent)*(x[jp][dir_2]-zcent);
						if (rip < rjp) break;
					} else
						if ( fabs(x[ip][dir_1]-0.5*yprd) < fabs(x[jp][dir_1]-0.5*yprd) ) break; // implemented for y-dir only
      	}
				if (j < nrbc) break;
    	}
			if (i == nbank) { loc_mol[k] = chk_mol; ncount++; break; }
		}

    MPI_Allreduce(&ncount,&totcount,1,MPI_INT,MPI_SUM,world);

		if (k == nmol) break;
	}

	MPI_Allreduce(&loc_mol[0],&rls_mol[0],nmol,MPI_INT,MPI_MAX,world);

  for (i = 0; i < nlocal; i++)
		if ( type[i] == type_plat && (mask[i] & groupbit) ) break;

	int maskplat, masktemp = 0;
	if (i < nlocal) masktemp = mask[i];

	MPI_Allreduce(&masktemp,&maskplat,1,MPI_INT,MPI_MAX,world);

	for (i = 0; i < nmol, ninflux > 0; i++){
		if (rls_mol[i] == 0) continue;
		for (j = 0; j < nlocal; j++)
	 		if (mol[j] == rls_mol[i]){
				type[j] = type_plat;
				mask[j] &= bankgroupbit;
				for (int igroup = 0; igroup < MAX_GROUP; igroup++)
					if ( group->names[igroup] && (maskplat & group->bitmask[igroup]) )
						mask[j] |= group->bitmask[igroup];
			}
		list_mol[i] = 0;
		--ninflux;
	}
	next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

int FixInflux::chkregion(int ip){
	double **x = atom->x;

	if (toupper(coord) == 'R'){
		double rp = sqrt( (x[ip][dir_1]-ycent)*(x[ip][dir_1]-ycent)+(x[ip][dir_2]-zcent)*(x[ip][dir_2]-zcent) );
		if ( rp > rlo && rp < rhi && x[ip][0] > xlo && x[ip][0] < xhi ) return 1;
	}else
		if (x[ip][0] > xlo && x[ip][0] < xhi && x[ip][1] > ylo && x[ip][1] < yhi && x[ip][2] > zlo && x[ip][2] < zhi ) return 1; 

	return 0;
}

/* ---------------------------------------------------------------------- */
