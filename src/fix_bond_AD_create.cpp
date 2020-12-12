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
#include "fix_bond_AD_create.h"
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
#include "common.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixBondADCreate::FixBondADCreate(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 17) error->all("Illegal fix bond/AD/create command");

  MPI_Comm_rank(world,&me);

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
  btype = atoi(arg[6]);
  ks = atof(arg[7]);
  r0 = atof(arg[8]);
  kf0 = atof(arg[9]);
  sig = atof(arg[10]);
  temp = atof(arg[11]);
  Nbond = atoi(arg[12]);
  imaxbond = atoi(arg[13]);
  jmaxbond = atoi(arg[14]);
	int seed = atoi(arg[15]);
  sprintf(bmodel,arg[16]);

  if (iatomtype < 1 || iatomtype > atom->ntypes || 
      jatomtype < 1 || jatomtype > atom->ntypes)
    error->all("Invalid atom type in fix bond/AD/create command");
  if (cutoff < 0.0) error->all("Illegal fix bond/AD/create command");
  if (btype < 1 || btype > atom->nbondtypes)
    error->all("Invalid bond type in fix bond/AD/create command");

  cutsq = cutoff*cutoff;

  // error check

  if (atom->molecular == 0)
    error->all("Cannot use fix bond/create with non-molecular systems");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // perform initial allocation of atom-based arrays
  // register with Atom class
  // bondcount values will be initialized in setup()

  bondcount = NULL;
//  grow_arrays(atom->nmax);
  atom->add_callback(0);
  countflag = 0;

  // set comm sizes needed by this fix

  comm_forward = Nbond;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = NULL;
  pp = NULL;
  probability = NULL;

  // zero out stats

  createcount = 0;
  createcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondADCreate::~FixBondADCreate()
{
  // unregister callbacks to this fix from Atom class

  atom->delete_callback(id,0);

  delete random;

  // delete locally stored arrays

  memory->sfree(bondcount);
  memory->sfree(partner);
  memory->sfree(pp);
  memory->destroy_2d_double_array(probability);
}

/* ---------------------------------------------------------------------- */

int FixBondADCreate::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondADCreate::init()
{
  // check cutoff for iatomtype,jatomtype

  if (force->pair == NULL || cutsq > force->pair->cutsq[iatomtype][jatomtype]) 
    error->all("Fix bond/AD/create cutoff is longer than pairwise cutoff");
  
  // warn if angles, dihedrals, impropers are being used

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

void FixBondADCreate::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixBondADCreate::setup(int vflag)
{
  int i,j,m;

  nmax = atom->nmax;
  bondcount = (int *)
    memory->smalloc(nmax*sizeof(int),"bond/AD/create:bondcount");
  partner = (int *)
    memory->smalloc(nmax*sizeof(int),"bond/AD/create:partner");
  pp = (double *)
    memory->smalloc(nmax*sizeof(double),"bond/AD/create:pp");
  probability = memory->create_2d_double_array(nmax,Nbond,"bond/AD/create:probability");

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
	    error->one("Could not count initial bonds in fix bond/AD/create");
	  bondcount[m]++;
	}
      }
    }

  // if newton_bond is set, need to communicate ghost counts
  // use reverseflag to toggle operations inside pack/unpack methods

  reverseflag = 0;
  if (newton_bond) comm->reverse_comm_fix(this);
  reverseflag = 1;

	neighbor->build();
}

/* ---------------------------------------------------------------------- */

void FixBondADCreate::post_integrate()
{
  int i,j,k,m,ii,jj,inum,jnum,itype,jtype,n1,n3,possible;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,min,max;
  int *ilist,*jlist,*numneigh,**firstneigh,*slist;

  if (update->ntimestep % nevery) return;

  // need updated ghost atom positions

  comm->communicate();

  // resize bond partner list and initialize it
  // probability array overlays pp array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    partner = (int *)
      memory->srealloc(partner,nmax*sizeof(int),"bond/AD/create:partner");
    pp = (double *)
      memory->srealloc(pp,nmax*sizeof(double),"bond/AD/create:pp");
  	probability = memory->grow_2d_double_array(probability,nmax,Nbond,"bond/AD/create:probability");
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
		pp[i] = 0.0;
	}

	double dtv = update->dt;

  // loop over neighbors of my atoms
  // setup possible partner list of bonds to create

  double **x = atom->x;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
		if (atom->molecule[i] > 0)
			if (atom->q[i] == 0.0) continue; // platelet must be activated
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      if (!(mask[j] & groupbit)) continue;
			if (atom->molecule[i] == atom->molecule[j]) continue;
			if (atom->molecule[j] > 0)
				if (atom->q[j] == 0.0) continue; // platelet must be activated
      jtype = type[j];

      possible = 0;
      if (itype == iatomtype && jtype == jatomtype) {
	if ((imaxbond == 0 || bondcount[i] < imaxbond) &&
           (jmaxbond == 0 || bondcount[j] < jmaxbond)) 
           possible = 1;
      } else if (itype == jatomtype && jtype == iatomtype) {
	if ((jmaxbond == 0 || bondcount[i] < jmaxbond)&& 
          (imaxbond == 0 || bondcount[j] < imaxbond))
        possible = 1;
      }
      if (!possible) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq >= cutsq) continue;

      double r = sqrt(rsq);
      double bondforce = ks*fabs(r-r0);
      double kf, pbond;
      if (!strcmp(bmodel,"flex"))
	      if (bondforce < Ftrans){
  	    	kf = kf01*pow(alpha,1.5)*(1.0+sigf1*bondforce/2.0/dgf1)*
    	         exp((dgf1/temp)*(1.0-alpha*(sigf1*bondforce/2.0/dgf1)*(sigf1*bondforce/2.0/dgf1)));
					kf = (kf < kf01) ? kf01:kf;
          }	else{
        	kf = kf02*pow(alpha,1.5)*(1.0+sigf2*bondforce/2.0/dgf2)*
          	   exp((dgf2/temp)*(1.0-alpha*(sigf2*bondforce/2.0/dgf2)*(sigf2*bondforce/2.0/dgf2)));
					kf = (kf < kf02) ? kf02:kf;
		  }
       else if (!strcmp(bmodel,"catch")){
				kf = kf0 * exp((bondforce*(sig-0.5*fabs(r-r0)))/temp);
				kf = (kf < kf0) ? kf0:kf;
       }

      pbond = 1.0 - exp(-kf*dtv);

      if (pbond > pp[i]) {
          partner[i] = tag[j];
          pp[i] = pbond;
      }
      if (pbond > pp[j]) {
          partner[j] = tag[i];
          pp[j] = pbond;
      }
    }
  }

  // reverse comm of partner info

  if (force->newton_pair) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  for (i = 0; i < nlocal; i++)
  	if (partner[i]){
			for (j = 0; j < Nbond; j++)
				probability[i][j] = random->uniform();
		}

  partnerflag = 1;
	comm->comm_fix(this);
  partnerflag = 0;
	comm->comm_fix(this);

  // create bonds
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  int *num_bond = atom->num_bond;
  double **bond_length = atom->bond_length;
  int **nspecial = atom->nspecial;
  int **special = atom->special;
  int newton_bond = force->newton_bond;

  int ncreate = 0;
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    if (partner[j] != tag[i]) continue;

		for (k = 0; k < Nbond; k++) {
	    // apply probability constraint
  	  // MIN,MAX insures values are added in same order on different procs

    	min = MIN(probability[i][k],probability[j][k]);
    	max = MAX(probability[i][k],probability[j][k]);
    	//if (0.5*(min+max) >= pp[i]) continue;
    	if (min >= pp[i]) continue;

	    // store bond with atom I
  	  // if newton_bond is set, only store with I or J

    	if (!newton_bond || tag[i] < tag[j]) {
      	bond_type[i][num_bond[i]] = btype;
      	bond_atom[i][num_bond[i]] = tag[j];
				bond_length[i][num_bond[i]] = r0;
      	num_bond[i]++;
				fprintf(stdout,"bond formed [prob], [mol(i)], [mol(j)], [btype], [i], [j]: %e %d %d %d %d %d\n", pp[i], atom->molecule[i], atom->molecule[j], btype, tag[i], tag[j]);
    	}

	    // add a 1-2 neighbor to special bond list for atom I
	    // atom J will also do this

	    slist = atom->special[i];
  	  n1 = nspecial[i][0];
    	n3 = nspecial[i][2];
    	if (n3 == atom->maxspecial)
      	error->one("New bond exceeded special list size in fix bond/AD/create");
    	for (m = n3; m > n1; m--) slist[m+1] = slist[m];
    	slist[n1] = tag[j];
    	nspecial[i][0]++;
    	nspecial[i][1]++;
    	nspecial[i][2]++;

    	// increment bondcount

    	bondcount[i]++;

	    // count the created bond once

  	  if (tag[i] < tag[j]) ncreate++;
			break;
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

int FixBondADCreate::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
  int i,j,k,m;

  m = 0;
	if (partnerflag) {
	  for (i = 0; i < n; i++) {
  	  j = list[i];
    	buf[m++] = partner[j];
		}
  	return 1;

	}	else {
		for (i = 0; i < n; i++) {
			j = list[i];
			for (k = 0; k < Nbond; k++)
    		buf[m++] = probability[j][k];
  	}
		return Nbond;
	}
}

/* ---------------------------------------------------------------------- */

void FixBondADCreate::unpack_comm(int n, int first, double *buf)
{
  int i,j,m,last;

  m = 0;
  last = first + n;
	if (partnerflag)
	  for (i = first; i < last; i++)
  	  partner[i] = static_cast<int> (buf[m++]);
	else
		for (i = first; i < last; i++)
			for (j = 0; j < Nbond; j++)
    		probability[i][j] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixBondADCreate::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (reverseflag) {
    for (i = first; i < last; i++) {
      buf[m++] = partner[i];
      buf[m++] = pp[i];
    }
    return 2;

  } else {
    for (i = first; i < last; i++) 
      buf[m++] = bondcount[i];
    return 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixBondADCreate::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;

  if (reverseflag) {
    for (i = 0; i < n; i++) {
      j = list[i];
      if (buf[m+1] > pp[j]) {
	partner[j] = static_cast<int> (buf[m++]);
	pp[j] = buf[m++];
      } else m += 2;
    }

  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      bondcount[j] += static_cast<int> (buf[m++]);
    }
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays 
------------------------------------------------------------------------- */

void FixBondADCreate::grow_arrays(int nmax)
{
  bondcount = (int *)
 	  memory->srealloc(bondcount,nmax*sizeof(int),"bond/AD/create:bondcount");
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays 
------------------------------------------------------------------------- */

void FixBondADCreate::copy_arrays(int i, int j)
{
  bondcount[j] = bondcount[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc 
------------------------------------------------------------------------- */

int FixBondADCreate::pack_exchange(int i, double *buf)
{
  buf[0] = bondcount[i];
  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc 
------------------------------------------------------------------------- */

int FixBondADCreate::unpack_exchange(int nlocal, double *buf)
{
  bondcount[nlocal] = static_cast<int> (buf[0]);
  return 1;
}

/* ---------------------------------------------------------------------- */

double FixBondADCreate::compute_vector(int n)
{
  if (n == 1) return (double) createcount;
  return (double) createcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixBondADCreate::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax*2 * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
