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
#include "stdlib.h"
#include "string.h"
#include "fix_print.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixPrint::FixPrint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix print command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix print command");

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  string = new char[n];
  strcpy(string,arg[4]);

  // parse optional args

  fp = NULL;
  screenflag = 1;

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"file") == 0 || strcmp(arg[iarg],"append") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix print command");
      if (me == 0) {
	if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
	else fp = fopen(arg[iarg+1],"a");
	if (fp == NULL) {
	  char str[128];
	  sprintf(str,"Cannot open fix print file %s",arg[iarg+1]);
	  error->one(str);
	}
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"screen") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix print command");
      if (strcmp(arg[iarg+1],"yes") == 0) screenflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) screenflag = 0;
      else error->all("Illegal fix print command");
      iarg += 2;
    } else error->all("Illegal fix print command");
  }

  // print header into file

  if (fp && me == 0) fprintf(fp,"# Fix print output for fix %s\n",id);

  copy = new char[MAXLINE];
  work = new char[MAXLINE];
}

/* ---------------------------------------------------------------------- */

FixPrint::~FixPrint()
{
  delete [] string;
  delete [] copy;
  delete [] work;

  if (fp && me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

int FixPrint::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrint::end_of_step()
{
  // make a copy of string to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  strcpy(copy,string);
  input->substitute(copy,0);
  strcat(copy,"\n");

  modify->addstep_compute(update->ntimestep + nevery);

  if (me == 0) {
    if (screenflag && screen) fprintf(screen,copy);
    if (screenflag && logfile) fprintf(logfile,copy);
    if (fp) {
      fprintf(fp,copy);
      fflush(fp);
    }
  }
}
