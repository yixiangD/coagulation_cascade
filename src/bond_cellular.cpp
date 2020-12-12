/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Alireza Yazdani
------------------------------------------------------------------------- */

#include <algorithm>
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "bond_cellular.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "integrate.h"
#include "force.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "random_mars.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondCellular::BondCellular(LAMMPS *lmp) : Bond(lmp)
{
  random = new RanMars(lmp, 5684 + comm->me);
}

/* ----------------------------------------------------------------------
   free all arrays 
------------------------------------------------------------------------- */

BondCellular::~BondCellular()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(bond_index);
    memory->sfree(ks);
    memory->sfree(r0);
    memory->sfree(temp);
    memory->sfree(dsig);
    memory->sfree(kr0);
    memory->sfree(rcs);
    memory->sfree(mu_targ);
    memory->sfree(qp); 
    memory->sfree(gamc);
    memory->sfree(gamt);
    memory->sfree(sigc);
    memory->sfree(sigt);
    memory->sfree(t_remod);
  }

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void BondCellular::compute(int eflag, int vflag)
{
  int i1,i2,n,m,l,type,atom1,cc;
  double rsq,r,dr,fbond,kr,pp,lmax,rr;
  double rlogarg,kph,l0,mu,lambda;
  double dvx, dvy, dvz, vv;
  char warning[128];
  double ff[6], ebond;

  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

	ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  int nt = update->ntimestep;
  for (int i = 0; i < fix_check; i++)
  	if (nt == modify->fix[fixchecklist[i]]->next_reneighbor){
			comm->exchange();
			comm->borders();
			neighbor->build();
			break;
		}

  double **x = atom->x;
	double **f = atom->f;
  double **v = atom->v;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  double *bondlist_length = neighbor->bondlist_length;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

	int k,n1,n3,*slist;
  int **nspecial = atom->nspecial;
  int **special = atom->special;

	int breakcount, nbreak = 0;

  for (n = 0; n < nbondlist; n++) {
    
    cc = 0;
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    l0 = bondlist_length[n];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    dvx = v[i1][0] - v[i2][0];
    dvy = v[i1][1] - v[i2][1];
    dvz = v[i1][2] - v[i2][2];
    domain->minimum_image(delx,dely,delz);

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);

    if (bond_index[type] == 1){ // WLC bonds
			double relax = 1.0;
			if (t_remod[type] > 0.0) relax = cell_remodel(i1, i2, t_remod[type]);
      lmax = l0*r0[type];
      rr = 1.0/r0[type];  
      kph = pow(l0,qp[type])*temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr);
      mu = 0.25*sqrt(3.0)*(temp[type]*(-0.25/pow(1.0-rr,2) + 0.25 + 0.5*rr/pow(1.0-rr,3))/lmax/rr + kph*(qp[type]+1.0)/pow(l0,qp[type]+1.0));
      lambda = mu / (mu_targ[type] * relax);
      kph = kph * mu_targ[type] * relax / mu;
      rr = r/lmax;
      rlogarg = pow(r,qp[type]+1.0); 
      vv = (delx*dvx + dely*dvy +  delz*dvz)/r;
     
      if (rr >= 1.0) {
        sprintf(warning,"WLC bond too long: %d %d %d %d %d %d %g",
                update->ntimestep,atom->molecule[i1],atom->molecule[i2],atom->tag[i1],atom->tag[i2],type,rr);
        error->warning(warning);
      }  
       
      generate_wrr();
      fbond = -temp[type]*(0.25/(1.0-rr)/(1.0-rr)-0.25+rr)/lambda/r + kph/rlogarg + (sigc[type]*wrr[3] - gamc[type]*vv)/r;

      // energy
      if (eflag) {
          ebond += 0.25*temp[type]*lmax*(3.0*rr*rr-2.0*rr*rr*rr)/(1.0-rr)/lambda;
        if (qp[type] == 1.0)
          ebond -= kph*log(r);
        else
          ebond += kph/(qp[type]-1.0)/pow(r,qp[type]-1.0);
      }
    }

    else if (bond_index[type] == 2) { // Harmonic bond break
			fbond = 0.0;
			if (r > rcs[type]) r = rcs[type];
					if (!onoff) {
			dr = r - r0[type];
	   	//if (dr > 0.0) fbond = -2.0 * ks[type] * dr / r;
	   	fbond = -2.0 * ks[type] * dr / r;
  	 	if (eflag) ebond = ks[type] * dr * dr;
		//fprintf(stdout,"i=%d j=%d r= %e fbond = %e\n",atom->tag[i1],atom->tag[i2],r,fbond);
					} else {
			dr = r - r0[type];
			bondforce = ks[type]*fabs(dr);
      kr = kr0[type] * exp(dsig[type]*bondforce/temp[type]);

      pp = 1.0 - exp(-kr*dtv);
      m = 0;
      while (m < atom->num_bond[i1]){
  	  	atom1 = atom->map(atom->bond_atom[i1][m]);
	  		if (atom1 == i2)
			    if ( (dr >= 0.0) && (random->uniform() < pp || r > rcs[type]) ){
		fprintf(stdout,"bond broken [dr], [kr], [prob], [mol], [i], [j]: %e %e %e %d %d %d\n",dr, kr, pp, atom->molecule[atom1], atom->tag[i1], atom->tag[i2]);
            cc++; nbreak++;
	    			l = atom->num_bond[i1];
	    			atom->bond_atom[i1][m] = atom->bond_atom[i1][l-1];
	    			atom->bond_type[i1][m] = atom->bond_type[i1][l-1];
                    atom->bond_length[i1][m] = atom->bond_length[i1][l-1];
	    			atom->num_bond[i1]--;

    // remove J from special bond list for atom I
    // atom J will also do this

    slist = special[i1];
    n1 = nspecial[i1][0];
    n3 = nspecial[i1][2];
    for (k = 0; k < n1; k++)
      if (slist[k] == i2) break;
    for (; k < n3-1; k++) slist[k] = slist[k+1];
    nspecial[i1][0]--;
    nspecial[i1][1]--;
    nspecial[i1][2]--;

    slist = special[i2];
    n1 = nspecial[i2][0];
    n3 = nspecial[i2][2];
    for (k = 0; k < n1; k++)
      if (slist[k] == i1) break;
    for (; k < n3-1; k++) slist[k] = slist[k+1];
    nspecial[i2][0]--;
    nspecial[i2][1]--;
    nspecial[i2][2]--;

	  			}
					else{
						// force & energy
	        	if (dr > 0.0) fbond += -2.0 * ks[type] * dr / r;
  	 				if (eflag) ebond += ks[type] * dr * dr;
						m++;
					}
				else m++;
    	}
     	if (cc != 0)	neighbor->forced_reneigh = 1;

				}
    } else fbond = 0.0;

    if (bond_index[type] == 1){
      if (newton_bond || i1 < nlocal){
        f[i1][0] += delx*fbond - gamt[type]*dvx + sigt[type]*wrr[0]/r;
        f[i1][1] += dely*fbond - gamt[type]*dvy + sigt[type]*wrr[1]/r;
        f[i1][2] += delz*fbond - gamt[type]*dvz + sigt[type]*wrr[2]/r;
      }

      if (newton_bond || i2 < nlocal){
        f[i2][0] -= delx*fbond - gamt[type]*dvx + sigt[type]*wrr[0]/r;
        f[i2][1] -= dely*fbond - gamt[type]*dvy + sigt[type]*wrr[1]/r;
        f[i2][2] -= delz*fbond - gamt[type]*dvz + sigt[type]*wrr[2]/r;
      }
    } 
		else if (bond_index[type] == 2){
      if (newton_bond || i1 < nlocal){
        f[i1][0] += delx*fbond;
        f[i1][1] += dely*fbond;
        f[i1][2] += delz*fbond;
      }

      if (newton_bond || i2 < nlocal){
        f[i2][0] -= delx*fbond;
        f[i2][1] -= dely*fbond;
        f[i2][2] -= delz*fbond;
      }
    }

    // virial contribution
    if (stress_ind == 2){
      ff[0] = delx*delx*fbond;
      ff[1] = dely*dely*fbond;
      ff[2] = delz*delz*fbond;
      ff[3] = delx*dely*fbond;
      ff[4] = delx*delz*fbond;
      ff[5] = dely*delz*fbond;

      for(int k = 0; k < n_stress_tot; ++k)
        modify->fix[stress_list[k]]->virial3(i1, i2, ff);
    }

		if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }

  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);
  atom->nbonds -= breakcount;

  l = 0;
  MPI_Allreduce(&neighbor->forced_reneigh,&l,1,MPI_INT,MPI_MAX,world);
  neighbor->forced_reneigh = l;
}

/* ---------------------------------------------------------------------- */

void BondCellular::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  ks = (double *) memory->smalloc((n+1)*sizeof(double),"bond:ks");
  r0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:r0");
  temp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:temp");
  dsig = (double *) memory->smalloc((n+1)*sizeof(double),"bond:dsig");
  kr0 = (double *) memory->smalloc((n+1)*sizeof(double),"bond:kr0");
  rcs = (double *) memory->smalloc((n+1)*sizeof(double),"bond:rcs");
  mu_targ = (double *) memory->smalloc((n+1)*sizeof(double),"bond:lambda");
  qp = (double *) memory->smalloc((n+1)*sizeof(double),"bond:qp");
  bond_index = (int *) memory->smalloc((n+1)*sizeof(int),"bond:bond_index");
  gamc = (double *) memory->smalloc((n+1)*sizeof(double),"bond:gamc");
  gamt = (double *) memory->smalloc((n+1)*sizeof(double),"bond:gamt");
  sigc = (double *) memory->smalloc((n+1)*sizeof(double),"bond:sigc");
  sigt = (double *) memory->smalloc((n+1)*sizeof(double),"bond:sigt");
  t_remod = (double *) memory->smalloc((n+1)*sizeof(double),"bond:t_remod");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs from one line in input script
------------------------------------------------------------------------- */

void BondCellular::coeff(int narg, char **arg)
{
  int count, i;
  double temp_one, rcs_one, mu_targ_one, qp_one, ks_one, r0_one, dsig_one, kr0_one, gamc_one, gamt_one, t_remod_one;

  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  if (strcmp(arg[1],"rbc") == 0 || strcmp(arg[1],"plat") == 0) {
    if (narg != 9) error->all("Incorrect args in bond_coeff command");
    temp_one = atof(arg[2]);
    r0_one = atof(arg[3]);
    mu_targ_one = atof(arg[4]);
    qp_one = atof(arg[5]);
    gamc_one = atof(arg[6]);
    gamt_one = atof(arg[7]);
    t_remod_one = atof(arg[8]);

    count = 0;
    for (i = ilo; i <= ihi; i++) {
      bond_index[i] = 1;
      temp[i] = temp_one;
      r0[i] = r0_one;
      mu_targ[i] = mu_targ_one;
      qp[i] = qp_one;
      gamc[i] = gamc_one;
      gamt[i] = gamt_one;
      t_remod[i] = t_remod_one;
      setflag[i] = 1;
      count++;
    }
    if (count == 0) error->all("Incorrect args in bond_coeff command");
  }
  else if (strcmp(arg[1],"bond") == 0) {
    if (narg != 9) error->all("Incorrect args in bond_coeff command");
    ks_one = atof(arg[2]);
    r0_one = atof(arg[3]);
    kr0_one = atof(arg[4]);
    dsig_one = atof(arg[5]);
    temp_one = atof(arg[6]);
    rcs_one = atof(arg[7]);
		onoff = atoi(arg[8]);	// bond formation/dissociation flag is set to 0 (fix_bond_AD_create/break does the job)

    count = 0;
    for (i = ilo; i <= ihi; i++) {
      bond_index[i] = 2;
      ks[i] = ks_one;
      r0[i] = r0_one;
      kr0[i] = kr0_one; 
      dsig[i] = dsig_one;
      temp[i] = temp_one;
      rcs[i] = rcs_one;
      setflag[i] = 1;
      count++;
    }
    if (count == 0) error->all("Incorrect args in bond_coeff command");
  } else
    error->all("Incorrect bond definition in bond_coeff command");  
}

/* ----------------------------------------------------------------------
   error check and initialize all values needed for force computation
------------------------------------------------------------------------- */

void BondCellular::init_style()
{
  double sdtt = sqrt(update->dt);

  if (!allocated) error->all("Bond coeffs are not set");
  for (int i = 1; i <= atom->nbondtypes; i++){
//    if (setflag[i] == 0) error->all("All bond coeffs are not set");
    if (bond_index[i] == 1){
      if (gamt[i] > 3.0*gamc[i]) error->all("Gamma_t > 3*Gamma_c");
      sigc[i] = sqrt(2.0*temp[i]*(3.0*gamc[i]-gamt[i]))/sdtt;
      sigt[i] = 2.0*sqrt(gamt[i]*temp[i])/sdtt;
    }
  }

  dtv = update->dt;

  fixchecklist = NULL;
  fixchecklist = new int[modify->nfix];

  fix_check = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->force_reneighbor)
      fixchecklist[fix_check++] = i;

}

/* ----------------------------------------------------------------------
   return an equilbrium bond length 
------------------------------------------------------------------------- */

double BondCellular::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file 
------------------------------------------------------------------------- */

void BondCellular::write_restart(FILE *fp)
{
  int i;

  for (i = 1; i <= atom->nbondtypes; i++) {
    fwrite(&bond_index[i],sizeof(int),1,fp);
    if (bond_index[i] == 1) {
      fwrite(&temp[i],sizeof(double),1,fp);
      fwrite(&r0[i],sizeof(double),1,fp);
      fwrite(&mu_targ[i],sizeof(double),1,fp);
      fwrite(&qp[i],sizeof(double),1,fp);
      fwrite(&gamc[i],sizeof(double),1,fp);
      fwrite(&gamt[i],sizeof(double),1,fp);
      fwrite(&t_remod[i],sizeof(double),1,fp);
    } else {
      fwrite(&ks[i],sizeof(double),1,fp);
      fwrite(&r0[i],sizeof(double),1,fp);
      fwrite(&kr0[i],sizeof(double),1,fp);
      fwrite(&dsig[i],sizeof(double),1,fp);
      fwrite(&temp[i],sizeof(double),1,fp);
      fwrite(&rcs[i],sizeof(double),1,fp);
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void BondCellular::read_restart(FILE *fp)
{
  int i;

  allocate();

	for (i = 1; i <= atom->nbondtypes; i++){
		if (comm->me == 0) {
       fread(&bond_index[i],sizeof(int),1,fp);  
       if (bond_index[i] == 1){
         fread(&temp[i],sizeof(double),1,fp);
         fread(&r0[i],sizeof(double),1,fp);
         fread(&mu_targ[i],sizeof(double),1,fp);
         fread(&qp[i],sizeof(double),1,fp);
         fread(&gamc[i],sizeof(double),1,fp);
         fread(&gamt[i],sizeof(double),1,fp);
         fread(&t_remod[i],sizeof(double),1,fp);
       } else {
         fread(&ks[i],sizeof(double),1,fp);
         fread(&r0[i],sizeof(double),1,fp);
         fread(&kr0[i],sizeof(double),1,fp);
         fread(&dsig[i],sizeof(double),1,fp);
         fread(&temp[i],sizeof(double),1,fp);
         fread(&rcs[i],sizeof(double),1,fp);
       }
  	}
  	MPI_Bcast(&bond_index[i],1,MPI_INT,0,world);
		if (bond_index[i] == 1){
  		MPI_Bcast(&temp[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&r0[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&mu_targ[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&qp[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&gamc[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&gamt[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&t_remod[i],1,MPI_DOUBLE,0,world);
		} else {
	  	MPI_Bcast(&ks[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&r0[i],1,MPI_DOUBLE,0,world);
	  	MPI_Bcast(&kr0[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&dsig[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&temp[i],1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&rcs[i],1,MPI_DOUBLE,0,world);
		}
	}

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double BondCellular::single(int type, double rsq, int i, int j)
{
  return 0.0; 
}

/* ---------------------------------------------------------------------- */

void BondCellular::generate_wrr()
{
  int i;
  double ww[3][3];
  double v1, v2, factor, ss;

  for (i=0; i<5; i++){
    ss = 100.0;
    while ( ss > 1.0 ){
      v1 = 2.0 * random->uniform() - 1.0;
      v2 = 2.0 * random->uniform() - 1.0;
      ss = v1*v1 + v2*v2;
    }
    factor = sqrt(-2.0 * log(ss)/ss);
    if (i < 3){
      ww[i][0] = factor*v1;
      ww[i][1] = factor*v2;
    }
    else if (i == 3){
      ww[0][2] = factor*v1;
      ww[1][2] = factor*v2;
    }
    else
      ww[2][2] = factor*v1;
  }
  wrr[3] = (ww[0][0]+ww[1][1]+ww[2][2])/3.0;
  wrr[0] = (ww[0][0]-wrr[3])*delx + 0.5*(ww[0][1]+ww[1][0])*dely + 0.5*(ww[0][2]+ww[2][0])*delz;
  wrr[1] = 0.5*(ww[1][0]+ww[0][1])*delx + (ww[1][1]-wrr[3])*dely + 0.5*(ww[1][2]+ww[2][1])*delz;
  wrr[2] = 0.5*(ww[2][0]+ww[0][2])*delx + 0.5*(ww[2][1]+ww[1][2])*dely + (ww[2][2]-wrr[3])*delz;
}

/* ----------------------------------------------------------------------*/

double BondCellular::cell_remodel(int i1, int i2, double tau)
{
	if (atom->molecular == 0) return 1.0;
	if (molecule == NULL){
  	int ifix;
  	for (ifix = 0; ifix < modify->nfix; ifix++)
    	if (strcmp("molecule",modify->fix[ifix]->style) == 0) break;
  	if (ifix == modify->nfix) error->one("BondCellular: FixMolecule style is not defined.");
 	 	molecule = dynamic_cast<FixMolecule *> (modify->fix[ifix]);
	}
	if (atom->molecule[i1] != atom->molecule[i2]) return 1.0;
	
	double tact = 0.0;
	if (molecule)
		tact = molecule->list_tact[ atom->molecule[i1] ];
	if (tact == 0.0) return 1.0;

	double coeff = exp( -(update->dt*update->ntimestep-tact) / tau );
	coeff = (coeff > 0.01) ? coeff : 0.01;

	//fprintf(stdout,"bond relaxation is done in molecule %d relax is = %f\n",atom->molecule[i1],coeff);
  return coeff;
}

/* ----------------------------------------------------------------------*/
