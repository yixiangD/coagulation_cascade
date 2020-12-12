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

#ifndef FIX_COM_GYRATION_H
#define FIX_COM_GYRATION_H

#include "fix.h"
#include <string>

namespace LAMMPS_NS {

class FixComGyration : public Fix {
 friend class FixInflux;

 public:
  FixComGyration(class LAMMPS *, int, char **);
  ~FixComGyration();
  int setmask();
  void init();
  void end_of_step();
  void write_stat();

  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 private:
  FILE *fp0, *fp[9999];
	char fname[FILENAME_MAX],f_name[FILENAME_MAX],grp[FILENAME_MAX];
  int comflag,tensorflag,posflag,sten_50,sten_75,countflag;
  int nmolecules;
  int idlo,idhi;
  int size_array_rows;      // rows in global array
  int size_array_cols;      // columns in global array
	int btype,groupbitLig;
  double R, xl, xh, H, Yc, xc, yc, yw;

  double *massproc,*masstotal;
  double **com,**comall;
  double *rg;
  double **rgt;
	double *vector;
	double **array;
	double *Ymin,*Ymax,*Yminmax;
	int		 *bondcount,*nbond,*nbondtot;
	int		 *tempi,*mollist;

  void molcom();

  int *molmap;                 // convert molecule ID to local index

  int molecules_in_group(int &, int &);
	void updatemol();
	bool exists(const std::string&);
};

}

#endif
