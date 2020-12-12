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
#include "math.h"
#include "mpi.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "fix_cc_bcflux.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "error.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include <iostream>
#include <fstream>
#include "common.h"

using namespace LAMMPS_NS;
using namespace std;

/* ---------------------------------------------------------------------- */

FixCCBcflux::FixCCBcflux(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    if(strcmp(style,"cc_bcflux") != 0 && narg != 4)
    	error->all("Illegal fix CCBcflux command");

    if(CTYPES != 25)
			error->all("fix CCBcflux expects 23+2 species");

    char fname[FILENAME_MAX];
    
    sprintf(fname,arg[3]);

    int Np=23;

    Diff_inv = (double *) memory->smalloc(Np*sizeof(double),"cc_bcflux:Diff_inv");
    if (comm->me == 0){
    	fstream file_read; 
  		file_read.open(fname,ios::in);

    	if(file_read.fail()) 
				error->all("cc_bcflux: Diffusion coeff. file should be prepared in this folder.");

			for(int i=0; i<Np; i++){
				file_read >> Diff[i];
				Diff[i] *= CC_FACTOR;
				Diff_inv[i]=1/Diff[i];
			}
			file_read.close();
    }    
    MPI_Bcast(Diff_inv,Np,MPI_DOUBLE,0,world);

    par[0] = 0.034;	//phi_11
    par[1] = 2000;	//Phi_11M
    par[2] = 32.4;	//k_7,9
    par[3] = 24.0;	//K_7,9M
    par[4] = 103;	//k_7,10
    par[5] = 240;	//K_7,10M
    par[6] = 6.52E-13;  //k_tPA^C
    par[7] = 9.27E-12;  //k_tPA^IIa
    par[8] = 5.059E-18; //k_tPA^Ia

    cw[0]  = 0.25;	//[TF-VIIa]
    cw[1]  = 2.0E+09;	//[ENDO]
    cw[2]  = 375;	//[XIIa]

//    printf("caution:fix_cc_bcflux, cpu is %d\n", comm->me);
//    MPI_Barrier(world);
}

/* ---------------------------------------------------------------------- */

FixCCBcflux::~FixCCBcflux()
{
  memory->sfree(Diff_inv); 
}

/* ---------------------------------------------------------------------- */

int FixCCBcflux::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCCBcflux::init()
{

}

/* ---------------------------------------------------------------------- */

void FixCCBcflux::post_force(int vflag)
{
  double **T = atom->T;
  double **Q = atom->Q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double tnow = update->ntimestep*update->dt;
  double L = domain->boxhi[0]-domain->boxlo[0];
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double t_par = par[7]*exp(-134.8*tnow);

  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
		double temp = par[2]*T[i][1]*cw[0]/(par[3]+T[i][1])*L;
		Ri[0] =  temp*Diff_inv[0];
		Ri[1] = -temp*Diff_inv[1];
		
		temp = par[3]*T[i][7]*cw[0]/(par[5]+T[i][7])*L;
		Ri[2] =  temp*Diff_inv[6];		
		Ri[3] = -temp*Diff_inv[7];

		temp = par[0]*T[i][13]*cw[2]/(par[1]+T[i][13])*L;		
		Ri[4] =  temp*Diff_inv[12];
		Ri[5] = -temp*Diff_inv[13];
		Ri[6] = (par[6]+par[7]*T[i][8]+par[8]*T[i][10])*cw[1]*L*Diff_inv[19];

		Q[i][0] += Ri[0];
		Q[i][1] += Ri[1];
		Q[i][6] += Ri[2];
		Q[i][7] += Ri[3];
		Q[i][12]+= Ri[4];
		Q[i][13]+= Ri[5];
		Q[i][19]+= Ri[6];
	}

  }
}

/* ---------------------------------------------------------------------- */

