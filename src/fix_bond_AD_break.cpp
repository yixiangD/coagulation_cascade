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
#include "fix_bond_AD_break.h"
#include "update.h"
#include "respa.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "common.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixBondADBreak::FixBondADBreak(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 12) error->all("Illegal fix bond/AD/break command");

  MPI_Comm_rank(world,&me);

  nevery = 1;

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = 2;
  scalar_vector_freq = 1;
  extvector = 1;

  double cutoff = atof(arg[3]);
  btype = atoi(arg[4]);
  ks = atof(arg[5]);
  r0 = atof(arg[6]);
  kr0 = atof(arg[7]);
  sig = atof(arg[8]);
  temp = atof(arg[9]);
	int seed = atoi(arg[10]);
  sprintf(bmodel,arg[11]);

  if (btype < 1 || btype > atom->nbondtypes)
    error->all("Invalid bond type in fix bond/AD/break command");
  if (cutoff < 0.0) error->all("Illegal fix bond/AD/break command");

  cutsq = cutoff*cutoff;

  // error check

  if (atom->molecular == 0)
    error->all("Cannot use fix bond/break with non-molecular systems");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + me);

  // set comm sizes needed by this fix

  comm_forward = 2;
  comm_reverse = 2;

  // allocate arrays local to this fix

  nmax = 0;
  partner = NULL;
  pp = NULL;
  probability = NULL;

  // zero out stats

  breakcount = 0;
  breakcounttotal = 0;
}

/* ---------------------------------------------------------------------- */

FixBondADBreak::~FixBondADBreak()
{
  delete random;

  // delete locally stored arrays

  memory->sfree(partner);
  memory->sfree(pp);
  memory->sfree(probability);
}

/* ---------------------------------------------------------------------- */

int FixBondADBreak::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBondADBreak::init()
{
  // warn if angles, dihedrals, impropers are being used

  if (force->angle || force->dihedral || force->improper) {
    if (me == 0) 
      error->warning("Broken bonds will not alter angles, "
		     "dihedrals, or impropers");
  }
}

/* ---------------------------------------------------------------------- */

void FixBondADBreak::setup(int vflag)
{
    nmax = atom->nmax;
    partner = (int *)
      memory->smalloc(nmax*sizeof(int),"bond/AD/break:partner");
    pp = (double *)
      memory->smalloc(nmax*sizeof(double),"bond/AD/break:pp");
    probability = (double *)
      memory->smalloc(nmax*sizeof(double),"bond/AD/break:probability");
}

/* ---------------------------------------------------------------------- */

void FixBondADBreak::post_integrate()
{
  int i,j,k,m,n,i1,i2,n1,n3,possible,type;
  double delx,dely,delz,rsq,min,max;
  int *slist;

  if (update->ntimestep % nevery) return;

  // need updated ghost atom positions

  comm->communicate();

  // resize bond partner list and initialize it
  // probability array overlays pp array
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    partner = (int *)
      memory->smalloc(nmax*sizeof(int),"bond/AD/break:partner");
    pp = (double *)
      memory->smalloc(nmax*sizeof(double),"bond/AD/break:pp");
    probability = (double *)
      memory->smalloc(nmax*sizeof(double),"bond/AD/break:probability");
  }

  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    partner[i] = 0;
    pp[i] = 0.0;
  }

	double dtv = update->dt;

  // loop over bond list
  // setup possible partner list of bonds to break

  double **x = atom->x;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];
    if (!(mask[i1] & groupbit)) continue;
    if (!(mask[i2] & groupbit)) continue;
    if (type != btype) continue;

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];
    domain->minimum_image(delx,dely,delz);
    rsq = delx*delx + dely*dely + delz*delz;

		double r = sqrt(rsq);
		if (r < r0) continue;
    double bondforce = ks*fabs(r-r0);
		double kr, pbond;
		if (!strcmp(bmodel,"flex"))
	    if (bondforce < Ftrans)
  	    kr = kr01*exp(sigr1*bondforce/temp);
    	else
      	kr = kr02*exp(sigr2*bondforce/temp);
		else if (!strcmp(bmodel,"slip"))
	    kr = kr0*exp(sig*bondforce/temp);

    pbond = 1.0 - exp(-kr*dtv);
		if (rsq > cutsq) pbond = 1.0;

    if (pbond > pp[i1]) {
      partner[i1] = tag[i2];
      pp[i1] = pbond;
    }
    if (pbond > pp[i2]) {
      partner[i2] = tag[i1];
      pp[i2] = pbond;
    }
  }

  // reverse comm of partner info

  if (force->newton_bond) comm->reverse_comm_fix(this);

  // each atom now knows its winning partner
  // for prob check, generate random value for each atom with a bond partner
  // forward comm of partner and random value, so ghosts have it

  for (i = 0; i < nlocal; i++)
    if (partner[i]) probability[i] = random->uniform();

  comm->comm_fix(this);

  // break bonds
  // if both atoms list each other as winning bond partner
  // and probability constraint is satisfied

  int **bond_type = atom->bond_type;
  int **bond_atom = atom->bond_atom;
  double **bond_length = atom->bond_length;
  int *num_bond = atom->num_bond;
  int **nspecial = atom->nspecial;
  int **special = atom->special;

  int nbreak = 0;
  for (i = 0; i < nlocal; i++) {
    if (partner[i] == 0) continue;
    j = atom->map(partner[i]);
    if (partner[j] != tag[i]) continue;

    // apply probability constraint
    // MIN,MAX insures values are added in same order on different procs

    min = MIN(probability[i],probability[j]);
    max = MAX(probability[i],probability[j]);
    if (0.5*(min+max) >= pp[i]) continue;
    //if (max >= pp[i]) continue;

    // delete bond from atom I if I stores it
    // atom J will also do this

    for (m = 0; m < num_bond[i]; m++) {
      if (bond_atom[i][m] == partner[i]) {
      	k = num_bond[i];
        bond_atom[i][m] = bond_atom[i][k-1];
        bond_type[i][m] = bond_type[i][k-1];
        bond_length[i][m] = bond_length[i][k-1];
        num_bond[i]--;
				fprintf(stdout,"bond broken [prob], [mol(i)], [mol(j)], [btype], [i], [j]: %e %d %d %d %d %d\n", pp[i], atom->molecule[i], atom->molecule[j], btype, tag[i], tag[j]);
				break;
      }
    }

    // remove J from special bond list for atom I
    // atom J will also do this

    slist = atom->special[i];
    n1 = nspecial[i][0];
    n3 = nspecial[i][2];
    for (m = 0; m < n1; m++)
      if (slist[m] == partner[i]) break;
    for (; m < n3-1; m++) slist[m] = slist[m+1];
    nspecial[i][0]--;
    nspecial[i][1]--;
    nspecial[i][2]--;
    
    // count the broken bond once

    if (tag[i] < tag[j]) nbreak++;
  }

  // tally stats

  MPI_Allreduce(&nbreak,&breakcount,1,MPI_INT,MPI_SUM,world);
  breakcounttotal += breakcount;
  atom->nbonds -= breakcount;

  // trigger reneighboring if any bonds were formed

  if (breakcount) next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixBondADBreak::pack_comm(int n, int *list, double *buf,
			     int pbc_flag, int *pbc)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = partner[j];
    buf[m++] = probability[j];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondADBreak::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    partner[i] = static_cast<int> (buf[m++]);
    probability[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixBondADBreak::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = partner[i];
    buf[m++] = pp[i];
  }
  return 2;
}

/* ---------------------------------------------------------------------- */

void FixBondADBreak::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    if (buf[m+1] > pp[j]) {
      partner[j] = static_cast<int> (buf[m++]);
      pp[j] = buf[m++];
    } else m += 2;
  }
}

/* ---------------------------------------------------------------------- */

double FixBondADBreak::compute_vector(int n)
{
  if (n == 1) return (double) breakcount;
  return (double) breakcounttotal;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double FixBondADBreak::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax * sizeof(double);
  return bytes;
}
