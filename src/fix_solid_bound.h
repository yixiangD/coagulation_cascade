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

#ifndef FIX_SOLID_BOUND_H
#define FIX_SOLID_BOUND_H

#include "fix.h"
#include "fix_stat.h"
#include "fix_stat_all.h"

namespace LAMMPS_NS {

class FixSolidBound : public Fix {
 public:
  FixSolidBound(class LAMMPS *, int, char **);
  virtual ~FixSolidBound();
  virtual int setmask();
	virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void pre_force(int);
  virtual void post_force(int);
	virtual void end_of_step();

 protected:
  int ind_shear, ind_out, ind_press, ind_drift,num_shapes, nmax, mirror, 
			n_per, n_press, n_shear, n_normal, n_drift, iter, iter_coupled, num_points, inject, num_at_types;  
  int max_faces, s_apply, p_apply, numt, nprob, prob_num, groupbit_s, groupbit_p;
  int *ptype, *refl, **ndiv, **period, *ind_coupled, *ind_bc_cc, **info_coupled, *at_type, *at_groupbit;
  int *num_faces, **ind_faces, **prob_ind, *num_faces_wal, *num_faces_out;
  int *num_faces_flux, **ind_faces_flux;
  double *prob, ***area, *at_dens;
  double **x0, **aa, ****vel, ***rot, *f_press, *f_shear_t, *f_shear_n, *bc_cc;
	double *r_d, *v_d, min_, max_;
/*	related to concentrations */
/*******************************************/
	int ind_cc, ind_param, ind_space, n_ccD, n_ccN, groupbit_c;
	double *concen_D, *concen_N, r_ccD, r_ccN;
	double *kappa, *diff;
/*	 problem-specific 25 species */
	int div[3], nstat[4], num_step;
	double X_bc[6];
	double ***tmp, ***ntmp, ***num_avg, ****aT, ****aTtmp, ****Ri;
  int is, js, ks;
  int xper, yper, zper;
  double xs, ys, zs;
  double xlo, ylo, zlo, xhi, yhi, zhi;
  double dxlo, dylo, dzlo, dxhi, dyhi, dzhi;
  double dxs, dys, dzs;
	int map_index(double, double, double);
	void bc_flux();
  void post_force_cc(int);
/*******************************************/
  double coeff, power, r_cut, r_cutshear, r_cutflux, r_shear, r_out, r_press, rr_shear, coeff_press, coeff_flux, power_out, num_dens, kbt;
  double ****velx, ****vely, ****velz, ****num, ****dens;
  double ****fsx, ****fsy, ****fsz;
  double *vcoupled, *xcoupled;
  class RanMars *ranmars0;
  double ****ncount;
	double **cent_loc, **cent_all;
	int *counter, *angle_per_mol;
	///////////////////////////////////////////////////
  //outflow adaptive force test
  double num_out, num_in, natoms_0; 
  double num_out_total, num_in_total;
  double num_out_desire;
  double factor;
  double *fac_ave;
  double **adap_ave;
  int dump_each;
  double ***velboundx, ***velboundy, ***velboundz, ****velz_image;
  int numdiff;
  int adaptive;
  double ****accum, ****beta, num_accum;
  int Nloop;
  // Womersley flow test
  double frequency, kvis,wom, hy, xx1, xx2, dpressure;
  double vfield(double &, double &);
  //////////////////////////////////////////////////
  void post_force_vel(int);
  void setup_rot(int, double[], double[]);
  void rot_forward(double &, double &, double &, int);
  void rot_back(double &, double &, double &, int);
  void face_decide();
  void recalc_force();
  void couple_run();
	void cent_mass();
	int in_cell(double *, int *, int, double);

	char fixsbid[128], fixstatid[128];
	class FixSolidBound *sbound = NULL;
	class FixStatAll *stat = NULL;
};

}
#endif
