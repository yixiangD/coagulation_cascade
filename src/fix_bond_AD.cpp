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
#include "string.h"
#include "stdlib.h"
#include "fix_bond_AD.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixBondAD::FixBondAD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 17) error->all("Illegal fix bond/adh command");

	me = comm->me;

  nevery = 1;

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  scalar_vector_freq = 1;
  extvector = 1;

  iatomtype = atoi(arg[3]);
  jatomtype = atoi(arg[4]);
  double cutoff = atof(arg[5]);
  ligtype = atoi(arg[6]);
  btype = atoi(arg[7]);
	ks = atof(arg[8]);
	r0 = atof(arg[9]);
	kf0 = atof(arg[10]);
	sig = atof(arg[11]);
	temp = atof(arg[12]);
  Nlig = atoi(arg[13]);
  imaxbond = atoi(arg[14]);
  jmaxbond = atoi(arg[15]);
  int seed = atoi(arg[16]);

  if (iatomtype < 1 || iatomtype > atom->ntypes || 
      jatomtype < 1 || jatomtype > atom->ntypes)
    error->all("Invalid atom type in fix bond/adh command");
  if (cutoff < 0.0) error->all("Illegal fix bond/adh command");
  if (btype < 1 || btype > atom->nbondtypes)
    error->all("Invalid bond type in fix bond/adh command");

  cutsq = cutoff*cutoff;

  // error check

  if (atom->molecular == 0)
    error->all("Cannot use fix bond/adh with non-molecular systems");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // perform initial allocation of atom-based arrays
  // register with Atom class
  // bondcount values will be initialized in setup()

  bondcount = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  countflag = 0;

  // set comm sizes needed by this fix

  comm_forward = 2;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;

  // zero out stats

  createcount = 0;
  createcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondAD::~FixBondAD()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete random;

  // delete locally stored arrays

  memory->sfree(bondcount);
}

/* ---------------------------------------------------------------------- */

int FixBondAD::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondAD::init()
{
  // check cutoff for iatomtype,jatomtype

  if (force->pair == NULL || cutsq > force->pair->cutsq[iatomtype][jatomtype]) 
    error->all("Fix bond/adh cutoff is longer than pairwise cutoff");
  
  if (force->angle || force->dihedral || force->improper) {
    if (me == 0) 
      error->warning("Created bonds will not create angles, "
		     "dihedrals, or impropers");
  }

  // need a half neighbor list, built when ever re-neighboring occurs

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
}

/* ---------------------------------------------------------------------- */

void FixBondAD::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondAD::setup(int vflag)
{
  int i,j,m;

  // compute initial bondcount if this is first run
  // can't do this earlier, like in constructor or init, b/c need ghost info

  if (countflag) return;
  countflag = 1;

  // count bonds stored with each bond I own
  // if newton bond is not set, just increment count on atom I
  // if newton bond is set, also increment count on atom J even if ghost
  // bondcount is long enough to tally ghost atom counts
  
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int nall = nlocal + nghost;
  int newton_bond = force->newton_bond;

  for (i = 0; i < nall; i++) bondcount[i] = 0;

  for (i = 0; i < nlocal; i++)
    for (j = 0; j < num_bond[i]; j++) {
      if (bond_type[i][j] == btype) {
        bondcount[i]++;
	if (newton_bond) {
	  m = atom->map(bond_atom[i][j]);
	  if (m < 0)
	    error->one("Could not count initial bonds in fix bond/adh");
	  bondcount[m]++;
	}
      }
    }

  // if newton_bond is set, need to communicate ghost counts

  if (newton_bond) comm->reverse_comm_fix(this);

	neighbor->build();
}

/* ---------------------------------------------------------------------- */

void FixBondAD::post_integrate()
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,r;
  double vxtmp,vytmp,vztmp,delvx,delvy,delvz;
  int *ilist,*jlist,*numneigh,**firstneigh;

  // need updated ghost atom positions

  comm->communicate();

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  // loop over neighbors of my atoms
  // setup possible partner list of bonds to create

  double **x = atom->x;
  double **v = atom->v;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int *type = atom->type;

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  double **bond_length = atom->bond_length;
  int *num_bond = atom->num_bond;
  int newton_bond = force->newton_bond;

	int k,n1,n3,*slist;
  int **nspecial = atom->nspecial;
  int **special = atom->special;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  int ncreate = 0;
	double dtv = update->dt;

	//Alireza: GPIba-vWF
	//const double lscale = 1.0e-6, tscale = 1.802e-4, fscale = 4.5e-14, escale = 4.5e-20;
	//const double Ftrans = 1.0e-11/fscale, alpha = 1.0+0.05/2.25, kf01 = 5.8e0*tscale, kf02 = 52.0e0*tscale, 
	//			 sigf1 = 1.58e-9/lscale, sigf2 = 1.31e-9/lscale, dgf1 = 6.17e-21/escale, dgf2 = 6.17e-21/escale;
	const double IIathresh = 2.5;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
//    if (atom->molecule[i] > 0)
//      if (atom->T[i][8]/IIathresh < 1.0) continue; //Alireza: setting a threshold value for platelet activation based on thrombin concentration
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    vxtmp = v[i][0];
    vytmp = v[i][1];
    vztmp = v[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (!(mask[j] & groupbit)) continue;
      if (atom->molecule[i] == atom->molecule[j]) continue;
//      if (atom->molecule[j] > 0)
//        if (atom->T[j][8]/IIathresh < 1.0) continue; //Alireza: setting a threshold value for platelet activation based on thrombin concentration
      jtype = type[j];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;
			double r = sqrt(rsq);

      delvx = vxtmp - v[j][0];
      delvy = vytmp - v[j][1];
      delvz = vztmp - v[j][2];
      double vs = sqrt(delvx*delvx+delvy*delvy+delvz*delvz);

      if (itype == ligtype){
	    	if (bondcount[i] > jmaxbond || bondcount[j] > imaxbond) continue;
				bondforce = ks*fabs(r-r0);	//Alireza: GPIba-vWF
//				if (bondforce < Ftrans)
//					kf = kf01*pow(alpha,1.5)*(1.0+sigf1*bondforce/2.0/dgf1)*
//							exp((dgf1/temp)*(1.0-alpha*(sigf1*bondforce/2.0/dgf1)*(sigf1*bondforce/2.0/dgf1)));
//				else
//					kf = kf02*pow(alpha,1.5)*(1.0+sigf2*bondforce/2.0/dgf2)*
//							exp((dgf2/temp)*(1.0-alpha*(sigf2*bondforce/2.0/dgf2)*(sigf2*bondforce/2.0/dgf2)));
       	kf = kf0 * exp((bondforce*(sig-0.5*fabs(r-r0)))/temp);

        double pp = 1.0 - exp(-kf*dtv);
        for (m = 0; m < Nlig; m++){
        	if (random->uniform() < pp){
//	if (bondforce < Ftrans)
//		fprintf(stdout,"bond formed itype: %e %e %d %d kf01\n",r-r0,pp,atom->tag[i],atom->tag[j]);
//	else
//		fprintf(stdout,"bond formed itype: %e %e %d %d kf02\n",r-r0,pp,atom->tag[i],atom->tag[j]);
	fprintf(stdout,"bond formed [dr], [kf], [prob], [mol], [i], [j]: %e %e %e %d %d %d\n",r-r0, kf, pp, atom->molecule[j], atom->tag[i], atom->tag[j]);
						if (bondcount[i] > jmaxbond || bondcount[j] > imaxbond) break;
           	bond_type[i][num_bond[i]] = btype;
            bond_atom[i][num_bond[i]] = atom->tag[j];
            bond_length[i][num_bond[i]] = r0;
            num_bond[i]++;

						bondcount[i]++;
						bondcount[j]++;
            ncreate++;

    // add a 1-2 neighbor to special bond list for atom I
    // atom J will also do this

    slist = special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    if (n3 == atom->maxspecial)
      error->warning("New bond exceeded special list size in fix bond/adh");
    for (k = n3; k > n1; k--) slist[k+1] = slist[k];
    slist[n1] = atom->tag[j];
    nspecial[i][0]++;
    nspecial[i][1]++;
    nspecial[i][2]++;

    slist = special[j];
    n1 = nspecial[j][0];
    n3 = nspecial[j][2];
    if (n3 == atom->maxspecial)
      error->warning("New bond exceeded special list size in fix bond/adh");
    for (k = n3; k > n1; k--) slist[k+1] = slist[k];
    slist[n1] = atom->tag[i];
    nspecial[j][0]++;
    nspecial[j][1]++;
    nspecial[j][2]++;

          }
        }
			}

      if (jtype == ligtype){
        if (bondcount[i] > imaxbond || bondcount[j] > jmaxbond) continue;
				bondforce = ks*fabs(r-r0);	//Alireza: GPIba-vWF
//        if (bondforce < Ftrans)
//          kf = kf01*pow(alpha,1.5)*(1.0+sigf1*bondforce/2.0/dgf1)*
//              exp((dgf1/temp)*(1.0-alpha*(sigf1*bondforce/2.0/dgf1)*(sigf1*bondforce/2.0/dgf1)));
//        else
//          kf = kf02*pow(alpha,1.5)*(1.0+sigf2*bondforce/2.0/dgf2)*
//              exp((dgf2/temp)*(1.0-alpha*(sigf2*bondforce/2.0/dgf2)*(sigf2*bondforce/2.0/dgf2)));
       	kf = kf0 * exp((bondforce*(sig-0.5*fabs(r-r0)))/temp);
        
				double pp = 1.0 - exp(-kf*dtv);
        for (m = 0; m < Nlig; m++){
        	if (random->uniform() < pp){
//	if (bondforce < Ftrans)
//		fprintf(stdout,"bond formed jtype: %e %e %d %d kf01\n",r-r0,pp,atom->tag[i],atom->tag[j]);
//	else
//		fprintf(stdout,"bond formed jtype: %e %e %d %d kf02\n",r-r0,pp,atom->tag[i],atom->tag[j]);
		fprintf(stdout,"bond formed [dr], [kf], [prob], [mol], [i], [j]: %e %e %e %d %d %d\n",r-r0, kf, pp, atom->molecule[i], atom->tag[i], atom->tag[j]);
      	    if (bondcount[i] > imaxbond || bondcount[j] > jmaxbond) break;
            bond_type[j][num_bond[j]] = btype;
            bond_atom[j][num_bond[j]] = atom->tag[i];
            bond_length[j][num_bond[j]] = r0;
            num_bond[j]++;
                
            bondcount[j]++;
            bondcount[i]++;
            ncreate++;

    // add a 1-2 neighbor to special bond list for atom I
    // atom J will also do this

    slist = special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    if (n3 == atom->maxspecial)
      error->warning("New bond exceeded special list size in fix bond/adh");
    for (k = n3; k > n1; k--) slist[k+1] = slist[k];
    slist[n1] = atom->tag[j];
    nspecial[i][0]++;
    nspecial[i][1]++;
    nspecial[i][2]++;

    slist = special[j];
    n1 = nspecial[j][0];
    n3 = nspecial[j][2];
    if (n3 == atom->maxspecial)
      error->warning("New bond exceeded special list size in fix bond/adh");
    for (k = n3; k > n1; k--) slist[k+1] = slist[k];
    slist[n1] = atom->tag[i];
    nspecial[j][0]++;
    nspecial[j][1]++;
    nspecial[j][2]++;

        	}
      	}
     	}
		}
	}

  // tally stats

  MPI_Allreduce(&ncreate,&createcount,1,MPI_INT,MPI_SUM,world);
  createcounttotal += createcount;
  atom->nbonds += createcount;

  // trigger reneighboring if any bonds were formed

  if (createcount) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixBondAD::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++) 
    buf[m++] = bondcount[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void FixBondAD::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    bondcount[j] += static_cast<int> (buf[m++]);
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixBondAD::grow_arrays(int nmax)
{
  bondcount = (int *)
    memory->srealloc(bondcount,nmax*sizeof(int),"bond/adh:bondcount");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixBondAD::copy_arrays(int i, int j)
{
  bondcount[j] = bondcount[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixBondAD::pack_exchange(int i, double *buf)
{
  buf[0] = bondcount[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixBondAD::unpack_exchange(int nlocal, double *buf)
{
  bondcount[nlocal] = static_cast<int> (buf[0]);
  return 1;
}

/* ---------------------------------------------------------------------- */

double FixBondAD::compute_vector(int n)
{
  if (n == 1) return (double) createcount;
  return (double) createcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixBondAD::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*2 * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
