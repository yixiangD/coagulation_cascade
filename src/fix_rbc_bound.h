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

#ifndef FIX_BOUND_RBC_H
#define FIX_BOUND_RBC_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixBoundRBC : public Fix {
 public:
  FixBoundRBC(class LAMMPS *, int, char **);
  ~FixBoundRBC();
  int setmask();
  void setup(int);
  void initial_integrate(int);
  void final_integrate();

 private:
  int min_ind, max_ind, ind_tot, ntri, n_tri_check, n_point_check;
  int *ind_ntri, *ind_h, **ind_tri, **tri_ind, *tri_mol; 
  double **ind_x, **ind_v;
  int groupbit_comm, groupbit_inner, max_count;

  int nnmax, max_faces, ind_bounce;
  double d_cut, d_cut_sq, comm_cut;
  double **norm, **edge1, **edge2, **dif1, **dif2, **dd12;
  int *ind_norm, *tri_check, *point_check;
  int *num_faces, **ind_faces;
  double prd[3];
  
  int size_border, size_comm;
  int maxsend, maxrecv, maxswap;
  int *sendnum, *recvnum, *size_comm_send, *size_comm_recv; 
  int *maxsendlist, **sendlist;
  double *buf_send, *buf_recv; 
  double *slablo,*slabhi;

  void borders();
  void communicate();
  void local_faces(int);
  void grow_send(int);
  void grow_list(int, int);
  void grow_recv(int);
  int pack_border(int, double *, int, int *);
  int unpack_border(double *);
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *); 
  void face_decide();  
  void calc_norm(int);
  double find_root(double *);
  double newton(double, double *, double *);
  double solve_quadratic(double *);
};
 

}

#endif
