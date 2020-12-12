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
#include "fix_cc_reaction.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "modify.h"
#include "iostream"
#include "fstream"
#include "common.h"

using namespace LAMMPS_NS;

#define DIM 100 // allowing maximum 100 reactions

/* ---------------------------------------------------------------------- */

FixCCReaction::FixCCReaction(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
	int Np;
	double dum;
  char fname[FILENAME_MAX];
	char buf[BUFSIZ];
	FILE *f_read;

  if(strcmp(style,"cc_reaction") != 0 && narg != 9)
  	error->all("Illegal fix CCReaction command");
    
  sys = atoi(arg[3]); // ODE system model No. for coagulation kinetics
												// 0: Generic 1: Full Anand (w/ plasmin) in the bulk 2: Partial Anand (wo/ plasmin) in the bulk
												// 3: On platelet surface + ADP + vWF
  Np = atoi(arg[4]);
	lo = atof(arg[5]);		// lower bound of the reactive region
	hi = atof(arg[6]);		// upper bound of the reactive region
	dir = atoi(arg[7]);		// flow direction x: 0, y: 1, z: 2
  sprintf(fname,arg[8]); // all the rates in physical units
	nevery = 1;
    
  par = (double *) memory->smalloc(Np*sizeof(double),"cc_reaction:par");

  if (comm->me == 0){
    f_read = fopen(fname,"r");
    if(f_read == (FILE*) NULL)
      error->one("cc_reaction: Cannot read reaction rate file");
    for (int i=0; i<Np; i++){
      fgets(buf,BUFSIZ,f_read);
      sscanf(buf,"%lf %lf",&par[i],&dum);
			if (static_cast<int> (dum) == 1) par[i] *= tscaleCC;
    }
    fclose(f_read);
	}
  MPI_Bcast(par,Np,MPI_DOUBLE,0,world);
	Ri = (double *) memory->smalloc(DIM*sizeof(double),"cc_reaction:Ri");
}

/* ---------------------------------------------------------------------- */

FixCCReaction::~FixCCReaction()
{
	memory->sfree(par);
	memory->sfree(Ri);
}

/* ---------------------------------------------------------------------- */

int FixCCReaction::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
	mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCReaction::init()
{
}

/* ---------------------------------------------------------------------- */

void FixCCReaction::setup(int vflag)
{
	if (sys != 3) return;

  int ifix;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp("molecule",modify->fix[ifix]->style) == 0) break;
  if (ifix == modify->nfix) error->one("FixCCReaction: FixMolecule style is not defined.");
  molecule = dynamic_cast<FixMolecule *> (modify->fix[ifix]);
}

/* ---------------------------------------------------------------------- */

void FixCCReaction::post_force(int vflag)
{
  double **x = atom->x;
  double **T = atom->T;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

	for (int i = 0; i < DIM; i++) Ri[i] = 0.0;

  for (int i = 0; i < nlocal; i++) {
  	for(int k = 0; k < CTYPES; k++)	T[i][k] = T[i][k] > 0 ? T[i][k] : 0.0;

	if ( (mask[i] & groupbit) && x[i][dir] >= lo && x[i][dir] <= hi ) {

		switch (sys){
	case 0 :
		break;
	case 1 :
		Ri[0] = par[0]*T[i][12]*T[i][1]/(par[1]+T[i][1]) - par[2]*T[i][0]*T[i][14];
		Ri[1] = - par[0]*T[i][12]*T[i][1]/(par[1]+T[i][1]);
		Ri[2] = par[3]*T[i][8]*T[i][3]/(par[4]+T[i][3])  - par[5]*T[i][2] - par[6]*T[i][16]*T[i][2]/(par[7]+T[i][2]);		
		Ri[3] = - par[3]*T[i][8]*T[i][3]/(par[4]+T[i][3]);		
		Ri[4] = par[8]*T[i][8]*T[i][5]/(par[9]+T[i][5]) - par[10]*T[i][4] - par[11]*T[i][16]*T[i][4]/(par[12]+T[i][4]);
		Ri[5] = - par[8]*T[i][8]*T[i][5]/(par[9]+T[i][5]);
		Ri[6] = par[13]*T[i][23]*T[i][7]/(par[14]+T[i][7]) - par[15]*T[i][6]*T[i][14] - par[16]*T[i][15]*T[i][6];
		Ri[7] = - par[13]*T[i][23]*T[i][7]/(par[14]+T[i][7]);    	
		Ri[8] = par[17]*T[i][24]*T[i][9]/(par[18]+T[i][9]) - par[19]*T[i][8]*T[i][14];
		Ri[9] = - par[17]*T[i][24]*T[i][9]/(par[18]+T[i][9]);
		Ri[10]= par[20]*T[i][8]*T[i][11]/(par[21]+T[i][11]) - par[22]*T[i][20]*T[i][10]/(par[23]+T[i][10]);		
		Ri[11]= - par[20]*T[i][8]*T[i][11]/(par[21]+T[i][11]);
		Ri[12]= par[24]*T[i][8]*T[i][13]/(par[25]+T[i][13]) - par[26]*T[i][12]*T[i][14] - par[27]*T[i][12]*T[i][18];
		Ri[13]= - par[24]*T[i][8]*T[i][13]/(par[25]+T[i][13]);
		Ri[14]= - ( par[2]*T[i][0] + par[15]*T[i][6] + par[19]*T[i][8] + par[26]*T[i][12] )*T[i][14];
		Ri[15]= - par[16]*T[i][15]*T[i][6];
		Ri[16]= par[28]*T[i][8]*T[i][17]/(par[29]+T[i][17]) - par[30]*T[i][16]*T[i][18];
		Ri[17]= - par[28]*T[i][8]*T[i][17]/(par[29]+T[i][17]);	
		Ri[18]= - par[30]*T[i][16]*T[i][18] - par[27]*T[i][12]*T[i][18];
		Ri[19]= 0.0;
		Ri[20]= par[31]*T[i][19]*T[i][21]/(par[32]+T[i][21]) - par[33]*T[i][20]*T[i][22];
		Ri[21]= - par[31]*T[i][19]*T[i][21]/(par[32]+T[i][21]);
		Ri[22]= - par[33]*T[i][20]*T[i][22];
		Ri[23]= T[i][2]*T[i][0]/par[36];
		Ri[24]= T[i][4]*T[i][6]/par[37];
		break;
	case 2 :
    Ri[0] = par[0]*T[i][12]*T[i][1]/(par[1]+T[i][1]) - par[2]*T[i][0]*T[i][14];
    Ri[1] = - par[0]*T[i][12]*T[i][1]/(par[1]+T[i][1]);
    Ri[2] = par[3]*T[i][8]*T[i][3]/(par[4]+T[i][3])  - par[5]*T[i][2] - par[6]*T[i][16]*T[i][2]/(par[7]+T[i][2]);
    Ri[3] = - par[3]*T[i][8]*T[i][3]/(par[4]+T[i][3]);
    Ri[4] = par[8]*T[i][8]*T[i][5]/(par[9]+T[i][5]) - par[10]*T[i][4] - par[11]*T[i][16]*T[i][4]/(par[12]+T[i][4]);
    Ri[5] = - par[8]*T[i][8]*T[i][5]/(par[9]+T[i][5]);
    Ri[6] = par[13]*T[i][21]*T[i][7]/(par[14]+T[i][7]) - par[15]*T[i][6]*T[i][14] - par[16]*T[i][15]*T[i][6];
    Ri[7] = - par[13]*T[i][21]*T[i][7]/(par[14]+T[i][7]);
    Ri[8] = par[17]*T[i][22]*T[i][9]/(par[18]+T[i][9]) - par[19]*T[i][8]*T[i][14];
    Ri[9] = - par[17]*T[i][22]*T[i][9]/(par[18]+T[i][9]);
    Ri[10]= par[20]*T[i][8]*T[i][11]/(par[21]+T[i][11]);
    Ri[11]= - par[20]*T[i][8]*T[i][11]/(par[21]+T[i][11]);
    Ri[12]= par[24]*T[i][8]*T[i][13]/(par[25]+T[i][13]) - par[26]*T[i][12]*T[i][14] - par[27]*T[i][12]*T[i][18];
    Ri[13]= - par[24]*T[i][8]*T[i][13]/(par[25]+T[i][13]);
    Ri[14]= - ( par[2]*T[i][0] + par[15]*T[i][6] + par[19]*T[i][8] + par[26]*T[i][12] )*T[i][14];
    Ri[15]= - par[16]*T[i][15]*T[i][6];
    Ri[16]= par[28]*T[i][8]*T[i][17]/(par[29]+T[i][17]) - par[30]*T[i][16]*T[i][18];
    Ri[17]= - par[28]*T[i][8]*T[i][17]/(par[29]+T[i][17]);
    Ri[18]= - par[30]*T[i][16]*T[i][18] - par[27]*T[i][12]*T[i][18];
    Ri[19] = 0.0;
    Ri[20] = 0.0;
    Ri[21]= T[i][2]*T[i][0]/par[31];
    Ri[22]= T[i][4]*T[i][6]/par[32];
		break;
	case 3 :
		double tact = 0.0;
		if (molecule) tact = molecule->list_tact[ atom->molecule[i] ];
		if (tact > 0.0)
			Ri[20] = relcon * exp( -(update->dt*update->ntimestep-tact) / reltau ) / reltau;
		break;
		}

		for(int j=0; j<CTYPES; j++)
			Q[i][j] += Ri[j];

	}

  }
}

/* ---------------------------------------------------------------------- */

void FixCCReaction::end_of_step()
{
	if (sys != 3) return;

  double **x = atom->x;
  double **T = atom->T;
  double *tact = atom->q; // using variable q to set the time particle i becomes activated
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
  	if (mask[i] & groupbit){
			if ( x[i][dir] < lo || x[i][dir] > hi ){
				if (molecule) molecule->list_tact[ atom->molecule[i] ] = 0.0;
				tact[i] = 0.0;
				continue;
			}
			if (ADP_th > 0.0){
				if ( (T[i][8]/IIa_th + T[i][20]/ADP_th) > 1.0 )	// setting a threshold value for platelet activation based on thrombin & ADP concentrations
					if (tact[i] == 0.0) tact[i] = update->dt*update->ntimestep;
			}
			else{
				if ( T[i][8]/IIa_th > 1.0 ) // setting a threshold value for platelet activation based on thrombin concentration
					if (tact[i] == 0.0) tact[i] = update->dt*update->ntimestep;
			}
		}
}

/* ---------------------------------------------------------------------- */
