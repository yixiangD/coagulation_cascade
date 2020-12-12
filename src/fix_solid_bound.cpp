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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "fix_solid_bound.h"
#include "random_mars.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "group.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "neighbor.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "common.h"

#define  NUM 1000
#define  NUM1 10000
#define  EPS 1e-6
#define  pressure_grid 1.0

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixSolidBound::FixSolidBound(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  int i, j, k, l, ttyp, igroup, n_dpd_fluid;
  double rdata[NUM1];
  double dummy,a1[3],a2[3];
  char grp[50];
  char buf[BUFSIZ];
  char fname[FILENAME_MAX];
  char fname1[FILENAME_MAX];
  char fname2[FILENAME_MAX];
  char fname3[FILENAME_MAX];
  char fname4[FILENAME_MAX];
  char fname5[FILENAME_MAX];
  char fname6[FILENAME_MAX];
  char fname7[FILENAME_MAX];
  FILE *f_read, *f_read1, *f_read2, *f_read3, *f_read4, *f_read5, *f_read6, *f_read7;

  ranmars0 = new RanMars(lmp, 4366*(2+comm->me));
  if (narg > 6) error->all("Illegal fix solid/boundary command");
	ind_cc = atoi(arg[3]);
    sprintf(fname,arg[4]);
  if (narg == 6) sprintf(fixstatid,arg[5]);
	nevery = 1;
  if (ind_cc == 0) ind_param = ind_space = n_ccD = n_ccN = 0;

  if (comm->me == 0){
    l = 0;
    f_read = fopen(fname,"r");
    if(f_read == (FILE*) NULL)
      error->one("Could not open input boundary file");
    fgets(buf,BUFSIZ,f_read);
    sscanf(buf,"%d %d",&num_shapes,&num_at_types);
    for (i=0; i<num_at_types; i++){
      fgets(buf,BUFSIZ,f_read);
      sscanf(buf,"%d %lf %s",&ttyp,&dummy,grp);
      rdata[l++] = static_cast<double> (ttyp);
      rdata[l++] = dummy;
      igroup = group->find(grp);
      if (igroup == -1) error->one("Group ID does not exist");
      rdata[l++] = static_cast<double> (group->bitmask[igroup]);
    }
    fgets(buf,BUFSIZ,f_read);
    sscanf(buf,"%lf %lf %lf %lf %d %d %d",&r_cut, &r_cutshear, &r_cutflux, &kbt, &mirror, &inject, &n_dpd_fluid);
    fgets(buf,BUFSIZ,f_read);
    sscanf(buf,"%d",&ind_shear);
    if (ind_shear){	// imposing no-slip condition (old method)
      sscanf(buf,"%d %lf %lf %lf %d %d %d %d %s",&ind_shear,&r_shear,&coeff,&power,&n_per,&iter,&iter_coupled,&s_apply,grp);
      igroup = group->find(grp);
      groupbit_s = group->bitmask[igroup];
    }
    fgets(buf,BUFSIZ,f_read);
    sscanf(buf,"%d",&ind_out);
    if (ind_out){	// imposing Neumann condition at the outlet (not very robust)
      sscanf(buf,"%d %lf %lf %lf %d %d %d",&ind_out,&r_out,&coeff,&power_out,&n_per,&iter,&iter_coupled);
    }
    fgets(buf,BUFSIZ,f_read);
    sscanf(buf,"%d",&ind_press);
    if (ind_press){	// wall no-slip condition (preferable)
      sscanf(buf,"%d %d %lf %d %lf %lf %d %s %s %s %s",&ind_press,&n_press,&r_press,&n_shear,&rr_shear,&coeff_press,&p_apply,fname1,fname2,fname3,grp);
      igroup = group->find(grp);
      groupbit_p = group->bitmask[igroup];
    }
	if (ind_cc){	// flux and Dirichlet BCs for concentration at the wall
    	fgets(buf,BUFSIZ,f_read);
		sscanf(buf,"%d %d %lf %d %lf %s %s %s %s",&ind_param,&n_ccD,&r_ccD,&n_ccN,&r_ccN,fname5,fname6,fname7,grp);
		igroup = group->find(grp);
		groupbit_c = group->bitmask[igroup];
	}
  }
  MPI_Bcast(&num_shapes,1,MPI_INT,0,world);
  MPI_Bcast(&num_at_types,1,MPI_INT,0,world);
  MPI_Bcast(&l,1,MPI_INT,0,world);
  MPI_Bcast(&rdata[0],l,MPI_DOUBLE,0,world);
  at_type = (int *) memory->smalloc(num_at_types*sizeof(int),"fix_solid_bound:at_type");
  at_groupbit = (int *) memory->smalloc(num_at_types*sizeof(int),"fix_solid_bound:at_groupbit");
  at_dens = (double *) memory->smalloc(num_at_types*sizeof(double),"fix_solid_bound:at_dens");
  l = 0;
  for (i=0; i<num_at_types; i++){
    at_type[i] = static_cast<int> (rdata[l++]);
    at_dens[i] = rdata[l++];
    at_groupbit[i] = static_cast<int> (rdata[l++]);
    if (at_type[i] == 1)
      num_dens = at_dens[i];
  }
  MPI_Bcast(&r_cut,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_cutshear,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&r_cutflux,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kbt,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mirror,1,MPI_INT,0,world);
  MPI_Bcast(&inject,1,MPI_INT,0,world);
  MPI_Bcast(&n_dpd_fluid,1,MPI_INT,0,world);
  natoms_0 = static_cast<double> (n_dpd_fluid);
  MPI_Bcast(&ind_shear,1,MPI_INT,0,world);
  if (ind_shear){
    MPI_Bcast(&r_shear,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&coeff,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&power,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&n_per,1,MPI_INT,0,world);
    MPI_Bcast(&iter,1,MPI_INT,0,world);
    MPI_Bcast(&iter_coupled,1,MPI_INT,0,world);
    MPI_Bcast(&s_apply,1,MPI_INT,0,world);
    MPI_Bcast(&groupbit_s,1,MPI_INT,0,world);
  }
  MPI_Bcast(&ind_out,1,MPI_INT,0,world);
  if (ind_out){
    MPI_Bcast(&r_out,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&coeff,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&power_out,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&n_per,1,MPI_INT,0,world);
    MPI_Bcast(&iter,1,MPI_INT,0,world);
    MPI_Bcast(&iter_coupled,1,MPI_INT,0,world);
  }
  MPI_Bcast(&ind_press,1,MPI_INT,0,world);
  if (ind_press){
    MPI_Bcast(&r_press,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&n_press,1,MPI_INT,0,world);
    MPI_Bcast(&rr_shear,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&n_shear,1,MPI_INT,0,world);
    MPI_Bcast(&coeff_press,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&p_apply,1,MPI_INT,0,world);
    MPI_Bcast(&groupbit_p,1,MPI_INT,0,world);
    f_press = (double *) memory->smalloc(n_press*sizeof(double),"fix_solid_bound:f_press");
    f_shear_t = (double *) memory->smalloc(n_shear*sizeof(double),"fix_solid_bound:f_shear_t");
    f_shear_n = (double *) memory->smalloc(n_shear*sizeof(double),"fix_solid_bound:f_shear_n");

    if (comm->me == 0){
      f_read1 = fopen(fname1,"r");
      if(f_read1 == (FILE*) NULL)
        error->one("Could not open input pressure file");
      for (j=0; j<n_press; j++){
        fgets(buf,BUFSIZ,f_read1);
        sscanf(buf,"%lf %lf",&dummy,&rdata[j]);
      }
      fclose(f_read1);
    }
    MPI_Bcast(&rdata,n_press,MPI_DOUBLE,0,world);
    for (j=0; j<n_press; j++)
      f_press[j] = rdata[j];
   // read adaptive shear force file  
    if (comm->me == 0){
      f_read2 = fopen(fname2,"r");
      if(f_read2 == (FILE*) NULL)
        error->one("Could not open input tangential shear file");
      for (j=0; j<n_shear; j++){
        fgets(buf,BUFSIZ,f_read2);
        sscanf(buf,"%lf %lf",&dummy,&rdata[j]);
      }
      fclose(f_read2);
    }
    MPI_Bcast(&rdata,n_shear,MPI_DOUBLE,0,world);
    for (j=0; j<n_shear; j++)
      f_shear_t[j] = rdata[j];
    
    if (comm->me == 0){
      f_read3 = fopen(fname3,"r");
      if(f_read3 == (FILE*) NULL)
        error->one("Could not open input normal shear file");
      for (j=0; j<n_press; j++){
        fgets(buf,BUFSIZ,f_read3);
        sscanf(buf,"%lf %lf",&dummy,&rdata[j]);
      }
      fclose(f_read3);
    }
    MPI_Bcast(&rdata,n_shear,MPI_DOUBLE,0,world);
    for (j=0; j<n_shear; j++)
      f_shear_n[j] = rdata[j];
  }  

  if (ind_cc){
    MPI_Bcast(&ind_param,1,MPI_INT,0,world);
    MPI_Bcast(&n_ccD,1,MPI_INT,0,world);
    MPI_Bcast(&r_ccD,1,MPI_DOUBLE,0,world);
    MPI_Bcast(&n_ccN,1,MPI_INT,0,world);
    MPI_Bcast(&r_ccN,1,MPI_DOUBLE,0,world);
		MPI_Bcast(&groupbit_c,1,MPI_INT,0,world);
		if (n_ccD) concen_D = (double *) memory->smalloc(n_ccD*sizeof(double),"fix_solid_bound:concen_D");
		if (n_ccN) concen_N = (double *) memory->smalloc(n_ccN*sizeof(double),"fix_solid_bound:concen_N");
    if (comm->me == 0 && n_ccD){
      f_read5 = fopen(fname5,"r");
      if(f_read5 == (FILE*) NULL)
        error->one("Could not open input Dirichlet concentration BC file");
      for (j=0; j<n_ccD; j++){
        fgets(buf,BUFSIZ,f_read5);
        sscanf(buf,"%lf %lf",&dummy,&rdata[j]);
      }
      fclose(f_read5);
    }
		if (n_ccD){
	  	MPI_Bcast(&rdata,n_ccD,MPI_DOUBLE,0,world);
  	  for (j=0; j<n_ccD; j++)
    	  concen_D[j] = rdata[j];
		}
    if (comm->me == 0 && n_ccN){
      f_read6 = fopen(fname6,"r");
      if(f_read6 == (FILE*) NULL)
        error->one("Could not open input Neumann concentration BC file");
      for (j=0; j<n_ccN; j++){
        fgets(buf,BUFSIZ,f_read6);
        sscanf(buf,"%lf %lf",&dummy,&rdata[j]);
      }
      fclose(f_read6);
    }
		if (n_ccN){
	  	MPI_Bcast(&rdata,n_ccN,MPI_DOUBLE,0,world);
	    for (j=0; j<n_ccN; j++)
  	    concen_N[j] = rdata[j];
		}
    kappa = (double *) memory->smalloc(NPDE*sizeof(double),"fix_solid_bound:kappa");
    diff = (double *) memory->smalloc(NPDE*sizeof(double),"fix_solid_bound:diff");
		if (comm->me == 0 && ind_param){
      f_read7 = fopen(fname7,"r");
      if(f_read7 == (FILE*) NULL)
        error->one("Could not open input concentration parameter file");
      for (j=0; j<NPDE+1; j++){
        fgets(buf,BUFSIZ,f_read7);
        sscanf(buf,"%lf",&rdata[j]);
      }
			ind_space = static_cast<int> (rdata[NPDE]);
			if (ind_space == 2){
				fgets(buf,BUFSIZ,f_read7);
				sscanf(buf,"%lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d",&X_bc[0],&X_bc[1],&X_bc[2],&X_bc[3],&X_bc[4],&X_bc[5],&div[0],&div[1],&div[2],&nstat[0],&nstat[1],&nstat[2],&nstat[3]);
			}
      fclose(f_read7);
		}
    if (ind_param){	// parameters and diffusion coefficients related to Dirichlet BCs are given here
    	MPI_Bcast(&ind_space,1,MPI_INT,0,world);
    	MPI_Bcast(&rdata,NPDE,MPI_DOUBLE,0,world);
        for(j=0; j<NPDE; j++){
            kappa[j] = CC_FACTOR*rdata[j];
			diff[j] = D_ksi+D_Fick*kappa[j];
		}
    } else ind_space = 0;
		if (ind_space == 2){	// time-average of concentrations adjacent to the wall is computed to estimate surface reactions
										// problem-specific (Anand's Model)
		  MPI_Bcast(&X_bc,6,MPI_DOUBLE,0,world);
		  MPI_Bcast(&div,3,MPI_INT,0,world);
		  MPI_Bcast(&nstat,4,MPI_INT,0,world);
			xlo = X_bc[0]; xhi = X_bc[1];
			ylo = X_bc[2]; yhi = X_bc[3];
			zlo = X_bc[4]; zhi = X_bc[5];
			if ((div[0] < 1)||(div[1] < 1)||(div[2] < 1)) error->one("fix_solid_bound:Illegal division for statistics");
		  if ((xlo >= xhi)||(ylo >= yhi)||(zlo >= zhi)) error->one("fix_solid_bound:Illegal coordinates for statistics");
			if (nstat[3] < 1) error->one("fix_solid_bound:number of species for flux BC must be at least 1");
		  xs = xhi - xlo;
		  ys = yhi - ylo;
		  zs = zhi - zlo;
		  xper = domain->xperiodic;
		  yper = domain->yperiodic;
		  zper = domain->zperiodic;
		  dxlo = domain->boxlo[0];
		  dxhi = domain->boxhi[0];
		  dylo = domain->boxlo[1];
		  dyhi = domain->boxhi[1];
		  dzlo = domain->boxlo[2];
		  dzhi = domain->boxhi[2];
		  dxs = domain->xprd;
		  dys = domain->yprd;
		  dzs = domain->zprd;
		  num_step = 0;
		  nevery = nstat[1];
		  num_avg = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:num_avg");
      tmp = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:tmp");
      ntmp = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:ntmp");
			aT = new double***[nstat[3]];
      aTtmp = new double***[nstat[3]];
      Ri = new double***[n_sur_reac];
			for (i = 0; i < nstat[3]; i++){
	  		aT[i] = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:aT[i]");
        aTtmp[i] = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:aTtmp[i]");
			}
			for (i = 0; i < n_sur_reac; i++)
				Ri[i] = memory->create_3d_double_array(div[0],div[1],div[2],"fix_solid_bound:Ri[i]");
      for (i = 0; i < div[0]; i++)
     		for (j = 0; j < div[1]; j++)
        	for (k = 0; k < div[2]; k++){
						tmp[i][j][k] = 0.0;
						ntmp[i][j][k] = 0.0;
						num_avg[i][j][k] = 0.0;
			    	for (l = 0; l < nstat[3]; l++){
              aT[l][i][j][k] = 0.0;
              aTtmp[l][i][j][k] = 0.0;
						}
						for (l = 0; l < n_sur_reac; l++)
							Ri[l][i][j][k] = 0.0;
					}
		}
 }

  numt = 0;
  nprob = 0;
  num_points = 0;
  max_faces = 5;
  x0 = memory->create_2d_double_array(num_shapes,3,"fix_solid_bound:x0");
  aa = memory->create_2d_double_array(num_shapes,4,"fix_solid_bound:aa");
  rot = memory->create_3d_double_array(num_shapes,3,3,"fix_solid_bound:rot");
  ndiv = memory->create_2d_int_array(num_shapes,2,"fix_solid_bound:ndiv");
  period = memory->create_2d_int_array(num_shapes,2,"fix_solid_bound:period");
  ptype = (int *) memory->smalloc(num_shapes*sizeof(int),"fix_solid_bound:ptype");
  refl = (int *) memory->smalloc(num_shapes*sizeof(int),"fix_solid_bound:refl");
  ind_coupled = (int *) memory->smalloc(num_shapes*sizeof(int),"fix_solid_bound:ind_coupled");
  ind_bc_cc = (int *) memory->smalloc(num_shapes*sizeof(int),"fix_solid_bound:ind_bc_cc");
  bc_cc = (double *) memory->smalloc(num_shapes*sizeof(double),"fix_solid_bound:bc_cc");

  vel = new double***[num_shapes];
  area = new double**[num_shapes];
  ncount = new double***[num_shapes];

  if (ind_out || ind_shear){
    velx = new double***[num_shapes];
    vely = new double***[num_shapes];
    velz = new double***[num_shapes];
    num = new double***[num_shapes];
  }
  if(ind_shear){
    fsx = new double***[num_shapes];
    fsy = new double***[num_shapes];
    fsz = new double***[num_shapes];
  } else if (ind_out){
    dens = new double***[num_shapes];
    beta = new double***[num_shapes];
    accum = new double***[num_shapes];
  }

  for (i=0; i<num_shapes; i++){
    if (comm->me == 0){
      fgets(buf,BUFSIZ,f_read);
      sscanf(buf,"%d",&ttyp);
    }
    MPI_Bcast(&ttyp,1,MPI_INT,0,world);
    ptype[i] = ttyp;
    switch (ptype[i]){
      case 1 :
        if (comm->me == 0)
          sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&rdata[0],&rdata[1],&rdata[2],&rdata[3],&rdata[4],&rdata[5],&rdata[6],&rdata[7],&rdata[8],&rdata[9],&rdata[10],&rdata[11],&rdata[12],&rdata[13]); 
        MPI_Bcast(&rdata,14,MPI_DOUBLE,0,world);
        for (j=0; j<3; j++){
          x0[i][j] = rdata[j+1];
          a1[j] = rdata[j+4] - rdata[j+1];
          a2[j] = rdata[j+7] - rdata[j+1];
		}
        ndiv[i][0] = 1;
        ndiv[i][1] = 1;
        refl[i] = static_cast<int> (rdata[10]);
        ind_coupled[i] = static_cast<int> (rdata[11]);
        ind_bc_cc[i] = static_cast<int> (rdata[12]);	// '1' for Dirichlet; '2' for Neumann
        bc_cc[i] = rdata[13];		// fixed Neumann or Dirichlet BC value
        area[i] = memory->create_2d_double_array(ndiv[i][0], ndiv[i][1], "fix_solid_bound:area[i]");
        setup_rot(i,a1,a2);
        if(ind_coupled[i] == 1)
          num_points += 1;
        break;
      case 2 :
        if (comm->me == 0)
          sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&rdata[0],&rdata[1],&rdata[2],&rdata[3],&rdata[4],&rdata[5],&rdata[6],&rdata[7],&rdata[8],&rdata[9],&rdata[10],&rdata[11],&rdata[12],&rdata[13],&rdata[14],&rdata[15],&rdata[16],&rdata[17]); 
        MPI_Bcast(&rdata,18,MPI_DOUBLE,0,world);
        for (j=0; j<3; j++){
          x0[i][j] = rdata[j+1];
          a1[j] = rdata[j+4] - rdata[j+1];
          a2[j] = rdata[j+7] - rdata[j+1];
		}
        ndiv[i][0] = static_cast<int> (rdata[10]);
        ndiv[i][1] = static_cast<int> (rdata[11]);
        refl[i] = static_cast<int> (rdata[12]);
        period[i][0] = static_cast<int> (rdata[13]);
        period[i][1] = static_cast<int> (rdata[14]);
        ind_coupled[i] = static_cast<int> (rdata[15]);
        ind_bc_cc[i] = static_cast<int> (rdata[16]);
        bc_cc[i] = rdata[17];
        area[i] = memory->create_2d_double_array(ndiv[i][0], ndiv[i][1], "fix_solid_bound:area[i]");
        setup_rot(i,a1,a2);
        if(ind_coupled[i] == 1)
          num_points += ndiv[i][0] * ndiv[i][1];
        break;
      case 3 :
        if (comm->me == 0)
          sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&rdata[0],&rdata[1],&rdata[2],&rdata[3],&rdata[4],&rdata[5],&rdata[6],&rdata[7],&rdata[8],&rdata[9],&rdata[10],&rdata[11],&rdata[12],&rdata[13],&rdata[14]); 
        MPI_Bcast(&rdata,15,MPI_DOUBLE,0,world);
        for (j=0; j<3; j++){
          x0[i][j] = rdata[j+1];
          a1[j] = rdata[j+4] - rdata[j+1];
		}
        aa[i][1] = rdata[7];
        ndiv[i][0] = static_cast<int> (rdata[8]);
        ndiv[i][1] = static_cast<int> (rdata[9]);
        refl[i] = static_cast<int> (rdata[10]);
        period[i][0] = static_cast<int> (rdata[11]);
        ind_coupled[i] = static_cast<int> (rdata[12]);
        ind_bc_cc[i] = static_cast<int> (rdata[13]);
        bc_cc[i] = rdata[14];
        area[i] = memory->create_2d_double_array(ndiv[i][0], ndiv[i][1], "fix_solid_bound:area[i]");
        setup_rot(i,a1,a2);
        if(ind_coupled[i] == 1)
          num_points += ndiv[i][0] * ndiv[i][1];
        break;
      case 4 :
        if (comm->me == 0)
          sscanf(buf,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&rdata[0],&rdata[1],&rdata[2],&rdata[3],&rdata[4],&rdata[5],&rdata[6],&rdata[7],&rdata[8],&rdata[9],&rdata[10]);
        MPI_Bcast(&rdata,11,MPI_DOUBLE,0,world);
        for (j=0; j<3; j++){
          x0[i][j] = rdata[j+1];
		}
        aa[i][0] = rdata[4];
        ndiv[i][0] = static_cast<int> (rdata[5]);
        ndiv[i][1] = static_cast<int> (rdata[6]);
        refl[i] = static_cast<int> (rdata[7]);
        ind_coupled[i] = static_cast<int> (rdata[8]);
        ind_bc_cc[i] = static_cast<int> (rdata[9]);
        bc_cc[i] = rdata[10];
        area[i] = memory->create_2d_double_array(ndiv[i][0], ndiv[i][1], "fix_solid_bound:area[i]");
        setup_rot(i,a1,a2);
        if(ind_coupled[i] == 1)
          num_points += ndiv[i][0] * ndiv[i][1];
        break;
    }
    vel[i] = memory->create_3d_double_array(ndiv[i][0],ndiv[i][1],3,"fix_solid_bound:vel[i]");
    ncount[i] = memory->create_3d_double_array(ndiv[i][0],ndiv[i][1],num_at_types,"fix_solid_bound:ncount[i]");
    for(j = 0; j < ndiv[i][0]; ++j)
      for(k = 0; k < ndiv[i][1]; ++k)
        for(l = 0; l < num_at_types; ++l)
          ncount[i][j][k][l] = 0.0;

    if (ind_out || ind_shear){
      velx[i] = memory->create_3d_double_array(2*n_per,ndiv[i][0],ndiv[i][1],"fix_solid_bound:velx[i]");
      vely[i] = memory->create_3d_double_array(2*n_per,ndiv[i][0],ndiv[i][1],"fix_solid_bound:vely[i]");
      velz[i] = memory->create_3d_double_array(2*n_per,ndiv[i][0],ndiv[i][1],"fix_solid_bound:velz[i]");
      num[i] = memory->create_3d_double_array(2*n_per,ndiv[i][0],ndiv[i][1],"fix_solid_bound:num[i]");
    }
    if (ind_shear){
      fsx[i] = memory->create_3d_double_array(2,ndiv[i][0],ndiv[i][1],"fix_solid_bound:fsx[i]");
      fsy[i] = memory->create_3d_double_array(2,ndiv[i][0],ndiv[i][1],"fix_solid_bound:fsy[i]");
      fsz[i] = memory->create_3d_double_array(2,ndiv[i][0],ndiv[i][1],"fix_solid_bound:fsz[i]");
    }
    else if (ind_out){
      dens[i] = memory->create_3d_double_array(2,ndiv[i][0],ndiv[i][1],"fix_solid_bound:dens[i]");
      beta[i] = memory->create_3d_double_array(2,ndiv[i][0],ndiv[i][1],"fix_solid_bound:beta[i]");
      accum[i] = memory->create_3d_double_array(2000,ndiv[i][0],ndiv[i][1],"fix_solid_bound:accum[i]");
    }

    numt += 4*2*n_per*ndiv[i][0]*ndiv[i][1]; 
    nprob += ndiv[i][0] * ndiv[i][1];
    for(j = 0; j < ndiv[i][0]; ++j)
      for(k = 0; k < ndiv[i][1]; ++k)
      {
        if(comm->me == 0){
          fgets(buf,BUFSIZ,f_read);
          sscanf(buf, "%lf %lf %lf", &rdata[0], &rdata[1], &rdata[2]);
        }
        MPI_Bcast(&rdata, 3, MPI_DOUBLE, 0, world);
        rot_forward(rdata[0], rdata[1], rdata[2], i);
        vel[i][j][k][0] = rdata[0];
        vel[i][j][k][1] = rdata[1];
        vel[i][j][k][2] = rdata[2];
      }
  }
  
  vcoupled = (double *)memory->smalloc(3*num_points*sizeof(double),"fix_solid_bound:vcoupled");
  xcoupled = (double *)memory->smalloc(3*num_points*sizeof(double),"fix_solid_bound:xcoupled");
  info_coupled = memory->create_2d_int_array(num_points, 3, "fix_solid_bound:couple");
 
  if (inject){
    prob = (double *) memory->smalloc(nprob*sizeof(double),"fix_solid_bound:prob");
    prob_ind = memory->create_2d_int_array(nprob,3,"fix_solid_bound:prob_ind"); 
    fac_ave = (double *) memory->smalloc(100*sizeof(double),"fix_solid_bound:fac_ave");
    adap_ave = memory->create_2d_double_array(1000,20,"fix_solid_bound:adap_ave");
  }

  if (comm->me == 0)
    fclose(f_read);

	cent_loc = cent_all = NULL; 
	counter = angle_per_mol = NULL;

  adaptive = ind_out;
  num_accum = 0.0;
  num_out = 0.0;
  num_in = 0.0;
  factor = 0.0;
  coeff_flux = 1.0;
  numdiff = 0;
  Nloop = 20;
  dump_each = 200000;
  // test for womersley flow /////////////////
  frequency = M_PI / 160.0;
  hy = 10.0;
  kvis = 0.54;
  wom = hy * sqrt(frequency/kvis/2.0);
  xx1 = cos(wom)*cosh(wom);
  xx2 = sin(wom)*sinh(wom);
  //dpressure = -0.1296;
  dpressure = 0.1296;
  /////////////////////////////////////////////
}

/* ---------------------------------------------------------------------- */

FixSolidBound::~FixSolidBound()
{
  int i;

  if (ind_press)
  {
    memory->sfree(f_press);
    memory->sfree(f_shear_t);
    memory->sfree(f_shear_n);
  }
  if (ind_cc)
	{
  	if (n_ccD) memory->sfree(concen_D);
  	if (n_ccN) memory->sfree(concen_N);
		memory->sfree(kappa);
		if (ind_space == 2)
		{
			memory->destroy_3d_double_array(tmp);
			memory->destroy_3d_double_array(ntmp);
			memory->destroy_3d_double_array(num_avg);
			for(i = 0; i < nstat[3]; ++i)
			{
				memory->destroy_3d_double_array(aT[i]);
				memory->destroy_3d_double_array(aTtmp[i]);
			}
			for(i = 0; i < 7; ++i)
				memory->destroy_3d_double_array(Ri[i]);
			delete[] aT;
			delete[] aTtmp;
			delete[] Ri;
		}
		if (ind_space == 1)
			memory->destroy_4d_double_array(Ri);
	}

  for(i = 0; i < num_shapes; ++i)
  {
    memory->destroy_3d_double_array(vel[i]);
    memory->destroy_2d_double_array(area[i]);
    memory->destroy_3d_double_array(ncount[i]);
  }
  delete[] vel;
  delete[] area;
  delete[] ncount;

  if (ind_out || ind_shear){
    for (i=0; i<num_shapes; i++){
      memory->destroy_3d_double_array(velx[i]);
      memory->destroy_3d_double_array(vely[i]);
      memory->destroy_3d_double_array(velz[i]);
      memory->destroy_3d_double_array(num[i]);
    }
    delete[] velx;
    delete[] vely;
    delete[] velz;
    delete[] num;
  }
  if (ind_shear){
    for (i=0; i<num_shapes; i++){
      memory->destroy_3d_double_array(fsx[i]);
      memory->destroy_3d_double_array(fsy[i]);
      memory->destroy_3d_double_array(fsz[i]);
    }
    delete[] fsx;
    delete[] fsy;
    delete[] fsz;
  } else if (ind_out){
      for (i=0; i<num_shapes; i++){
        memory->destroy_3d_double_array(dens[i]);
        memory->destroy_3d_double_array(beta[i]);
        memory->destroy_3d_double_array(accum[i]);
      }
    delete[] dens;
    delete[] beta;
    delete[] accum;
  }
    
  if (inject){
    memory->sfree(prob);
    memory->destroy_2d_int_array(prob_ind);
  }

	if (atom->molecule_flag){
		memory->destroy_2d_double_array(cent_loc);
		memory->destroy_2d_double_array(cent_all);
		memory->sfree(counter);
		memory->sfree(angle_per_mol);
	}
  
  memory->destroy_2d_double_array(x0);
  memory->destroy_2d_double_array(aa);
  memory->destroy_3d_double_array(rot);
  memory->destroy_2d_int_array(ndiv); 
  memory->destroy_2d_int_array(period);
  memory->destroy_2d_int_array(info_coupled);
  memory->destroy_2d_int_array(ind_faces);
  memory->destroy_2d_int_array(ind_faces_flux);
  memory->sfree(ptype);
  memory->sfree(refl);
  memory->sfree(ind_coupled);
  memory->sfree(ind_bc_cc);
  memory->sfree(bc_cc);
  memory->sfree(xcoupled);
  memory->sfree(vcoupled);
  memory->sfree(at_type);
  memory->sfree(at_dens);
  memory->sfree(at_groupbit);
  memory->sfree(num_faces);
  memory->sfree(num_faces_flux);
}

/* ---------------------------------------------------------------------- */

int FixSolidBound::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::setup(int vflag)
{

	int ifix;
/*
	if (strcmp(fixsbid,"NULL") == 0 || strcmp(fixsbid,"null") == 0)
		return;
  for (ifix = 0; ifix < modify->nfix; ifix++)
    if (strcmp(fixsbid,modify->fix[ifix]->id) == 0) break;
  if (ifix == modify->nfix) error->one("FixSolidBound: FixSolidBound id is not defined.");
  sbound = dynamic_cast<FixSolidBound *> (modify->fix[ifix]);
*/
	if (ind_space == 1){
  	for (ifix = 0; ifix < modify->nfix; ifix++)
    	if (strcmp(fixstatid,modify->fix[ifix]->id) == 0) break;
  	if (ifix == modify->nfix) error->one("FixSolidBound: FixStatAll id is not defined.");
  	stat = dynamic_cast<FixStatAll *> (modify->fix[ifix]);

		Ri = memory->create_4d_double_array(n_sur_reac,stat->nx,stat->ny,stat->nz,"fix_solid_bound:Ri");
	}

  int i,j,k,l;
  double ssum, theta, veln, phi, xtemp, ytemp, ztemp;

	nmax = atom->nmax;
  ind_faces = memory->create_2d_int_array(nmax,max_faces,"fix_solid_bound:ind_faces");
  num_faces = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces");
  ind_faces_flux = memory->create_2d_int_array(nmax,max_faces,"fix_solid_bound:ind_faces_flux");
  num_faces_flux = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_flux");
  num_faces_wal = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_wal");
  num_faces_out = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_out");

  if (ind_shear)
    for(i = 0; i < num_shapes; i++)
      for(k = 0; k < ndiv[i][0]; k++)
        for(l = 0; l < ndiv[i][1]; l++){
          for(j = 0; j < 2*n_per; j++){
            velx[i][j][k][l] = 0.0;
            vely[i][j][k][l] = 0.0;
            velz[i][j][k][l] = 0.0;
            num[i][j][k][l] = 0.0;
	  			}
          for(j = 0; j < 2; j++){
            fsx[i][j][k][l] = 0.0;
            fsy[i][j][k][l] = 0.0;
            fsz[i][j][k][l] = 0.0;
	  			}
				}

  if (ind_out)
    for(i = 0; i < num_shapes; i++)
      for(k = 0; k < ndiv[i][0]; k++)
        for(l = 0; l < ndiv[i][1]; l++)
          for(j = 0; j < 2*n_per; j++){
            velx[i][j][k][l] = 0.0;
            vely[i][j][k][l] = 0.0;
            velz[i][j][k][l] = 0.0;
            num[i][j][k][l] = 0.0;
          }

  l = 0;
  for(i = 0; i < num_shapes; ++i)
  {
    if(ind_coupled[i] == 1)
    {
      if(ptype[i] == 1)
      {
        xtemp = (aa[i][0] + aa[i][1]) / 3.0;
        ytemp = aa[i][2] / 3.0;
        ztemp = 0;
        rot_back(xtemp, ytemp, ztemp, i);
        xcoupled[l] = xtemp + x0[i][0];
        xcoupled[l+1] = ytemp + x0[i][1];
        xcoupled[l+2] = ztemp + x0[i][2];
        info_coupled[l/3][0] = i;
        info_coupled[l/3][1] = 0;
        info_coupled[l/3][2] = 0;
        l = l + 3;
      }
      else if(ptype[i] == 2)
      {
        for(j = 0; j < ndiv[i][0]; ++j)
          for(k = 0; k < ndiv[i][1]; ++k)
          {
            xtemp = aa[i][0]*(j+0.5)/ndiv[i][0] + aa[i][1]*(k+0.5)/ndiv[i][1];
            ytemp = aa[i][2]*(k+0.5)/ndiv[i][1];
            ztemp = 0;
            rot_back(xtemp, ytemp, ztemp, i);
            xcoupled[l] = xtemp + x0[i][0];
            xcoupled[l+1] = ytemp + x0[i][1];
            xcoupled[l+2] = ztemp + x0[i][2];
            info_coupled[l/3][0] = i;
            info_coupled[l/3][1] = j;
            info_coupled[l/3][2] = k;
            l = l + 3;
          }
      }
      else if(ptype[i] == 3)
      {
        for(j  = 0; j < ndiv[i][0]; ++j)
          for(k = 0; k < ndiv[i][1]; ++k)
          {
            theta = 2.0*M_PI*static_cast<double>(k+0.5)/static_cast<double>(ndiv[i][1]);
            xtemp = aa[i][1] * cos(theta);
            ytemp = aa[i][1] * sin(theta);
            ztemp = (j+0.5)*aa[i][0] / static_cast<double>(ndiv[i][0]);
            rot_back(xtemp, ytemp, ztemp, i);
            xcoupled[l] = xtemp + x0[i][0];
            xcoupled[l+1] = ytemp + x0[i][1];
            xcoupled[l+2] = ztemp + x0[i][2];
            info_coupled[l/3][0] = i;
            info_coupled[l/3][1] = j;
            info_coupled[l/3][2] = k;
            l = l + 3;
          }
      }
      else if(ptype[i] == 4)
      {
        for(j = 0; j < ndiv[i][0]; ++j)
          for(k = 0; k < ndiv[i][1]; ++k)
          {
            theta = M_PI*static_cast<double>(j+0.5)/static_cast<double>(ndiv[i][0]);
            phi = 2.0*M_PI*static_cast<double>(k+0.5)/static_cast<double>(ndiv[i][1]);
            xtemp = aa[i][0] * sin(theta) * cos(phi);
            ytemp = aa[i][0] * sin(theta) * sin(phi);
            ztemp = aa[i][0] * cos(theta);
            rot_back(xtemp, ytemp, ztemp, i);
            xcoupled[l] = xtemp + x0[i][0];
            xcoupled[l+1] = ytemp + x0[i][1];
            xcoupled[l+2] = ztemp + x0[i][2];
            info_coupled[l/3][0] = i;
            info_coupled[l/3][1] = j;
            info_coupled[l/3][2] = k;
            l = l + 3;
          }
      }
    }
  }
  if(l != 3*num_points)
    error->all("error in the couple initialization!\n");

  if (inject){
    prob_num = 0;
    ssum = 0.0;
    for(i = 0; i < num_shapes; i++)
      for(j = 0; j < ndiv[i][0]; j++)
        for(k = 0; k < ndiv[i][1]; ++k) 
          {
            if(ptype[i] < 3)
            {
              veln = vel[i][j][k][2];
            }
            else if(ptype[i] == 3)
            {
              theta = 2.0*M_PI*static_cast<double>(k+0.5)/static_cast<double>(ndiv[i][1]);
              veln = vel[i][j][k][0]*cos(theta) + vel[i][j][k][1]*sin(theta);
            }
            else
            {
              theta = M_PI*static_cast<double>(j+0.5)/static_cast<double>(ndiv[i][0]);
              phi = 2.0*M_PI*static_cast<double>(k+0.5)/static_cast<double>(ndiv[i][1]);
              veln = vel[i][j][k][0]*sin(theta)*cos(phi) + vel[i][j][k][1]*sin(theta)*sin(phi)
                + vel[i][j][k][2]*cos(theta);
            }
            if(veln > 0)
              {
                if(prob_num == 0)
                  prob[prob_num] = veln*area[i][j][k];
                else
                  prob[prob_num] = prob[prob_num-1] + veln*area[i][j][k];
                ssum += veln*area[i][j][k];
                prob_ind[prob_num][0] = i;
                prob_ind[prob_num][1] = j;
                prob_ind[prob_num][2] = k;
                prob_num++;
              }
          }
    if (prob_num == 0) printf("probability problem in init() in step %d\n", update->ntimestep);
    for(i = 0; i < prob_num; i++)
      prob[i] /= ssum;

		if (adaptive){
	    num_in_total = num_out_desire = ssum * num_dens;
  	  num_out_total = ssum * num_dens;
		}
  }

	face_decide();
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::initial_integrate(int vflag)
{
  int i,j,k,kk,ll,ind,cond,i_x[NUM],cc,bounce,m,cplane;
  int ix, iy, l, n, pnum, pnum1;
  double dl[3], t_x[NUM], xd[2], dd[3],xp[3],xh[3],vh[3],vh1[3],vv[3],norm[3], fs[6];
  double tt,dot,d1,d2,dtt,u_inv,u0,u1,rr,rr1,ntemp,nfluid,nplat;
  double theta,theta1,phi,hh,randnum1,randnum2,randnum3,xtemp,ytemp,ztemp,vxtemp,vytemp,vztemp,vtemp,xx[3]; 
  double part_dat[NUM1], part_comm[8*NUM1];
  int rcounts[comm->nprocs], displs[comm->nprocs], offset;
  bool inside;
  double dtv = update->dt; 
  int step_t = update->ntimestep;
  double **x = atom->x;
  double **v = atom->v;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag;

  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

	if (ind_cc == 1) return;

  for(i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      cond = 1;
      dtt = dtv;
      cplane = -1;
      while (cond){
        ind = 0;
        for(m = 0; m < num_faces[i]; m++){
          kk = ind_faces[i][m];
          if (kk != cplane){
            for(j = 0; j < 3; j++){
              xh[j] = x[i][j] - x0[kk][j];
              vh[j] = v[i][j];
	    			}  
            rot_forward(xh[0],xh[1],xh[2],kk);
            rot_forward(vh[0],vh[1],vh[2],kk);  
            for(j = 0; j < 3; j++){
              dl[j] = vh[j]*dtt;
              xp[j] = xh[j] - dl[j]; 
	    			}
           
            tt = -1.0;
            if (ptype[kk]<3){
              u1 = xp[1]/aa[kk][2]; 
              u0 = (xp[0]-u1*aa[kk][1])/aa[kk][0];
              ix = static_cast<int> (u0*ndiv[kk][0]);
              iy = static_cast<int> (u1*ndiv[kk][1]);
              if(ix < 0)
                ix = 0;
              if(ix > (ndiv[kk][0] - 1))
                ix = ndiv[kk][0] - 1;
              if(iy < 0)
                iy = 0;
              if(iy > (ndiv[kk][1] - 1))
                iy = ndiv[kk][1] - 1;
              for(j = 0; j < 3; ++j)
                xh[j] -= vel[kk][ix][iy][j] * dtt;
	      
              if (xh[2]*xp[2] < 0.0)
                if (refl[kk] == 3 || (xp[2]>0.0 && refl[kk] == 1) || (xp[2]<0.0 && refl[kk] == 2))
                  tt = xp[2]/(xp[2]-xh[2]);  
	      
            } else if (ptype[kk] == 3){  
              rr = sqrt(xp[0]*xp[0] + xp[1]*xp[1]);
              ix = static_cast<int> (xp[2]*ndiv[kk][0]/aa[kk][0]);
              if(ix < 0)
                ix = 0;
              if(ix > (ndiv[kk][0] - 1))
                ix = ndiv[kk][0] - 1;
              theta = acos(xp[0]/rr);
              theta1 = asin(xp[1]/rr);
              if (theta1 < 0.0)
                theta = 2.0*M_PI - theta;
              iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);
              if(iy < 0)
		          { 
		            std::cout<<"in initial iy is "<<iy<<" theta is "<<theta<<" time step is "<<update->ntimestep<<" cpu is "<<comm->me<<"\n";
		            iy = 0;
		          }
              if(iy > ndiv[kk][1] -1)
		          {
                  std::cout<<"in initial iy is "<<iy<<" theta is "<<theta<<" time step is "<<update->ntimestep<<" cpu is "<<comm->me<<"\n";
                  iy = ndiv[kk][1] - 1;
		          }
              for(j = 0; j < 3; ++j)
		          {
		            xh[j] -= vel[kk][ix][iy][j] * dtt;
		            dl[j] -= vel[kk][ix][iy][j] * dtt;
		          }
              d1 = sqrt(xh[0]*xh[0] + xh[1]*xh[1]);
              d2 = sqrt(xp[0]*xp[0] + xp[1]*xp[1]);
              if ((d1-aa[kk][1])*(d2-aa[kk][1]) < 0.0)
                if (refl[kk] == 3 || (d2>aa[kk][1] && refl[kk] == 1) || (d2<aa[kk][1] && refl[kk] == 2)){
                  vv[0] = dl[0]*dl[0] + dl[1]*dl[1];
                  vv[1] = (xp[0]*dl[0] + xp[1]*dl[1])/vv[0];
                  vv[2] = (d2*d2-aa[kk][1]*aa[kk][1])/vv[0];
                  dot =  vv[1]*vv[1] - vv[2];
                  if (dot > 0.0){
                    if (vv[2] > 0) 
	              			tt = -vv[1] - sqrt(dot);
	            			else
	              			tt = -vv[1] + sqrt(dot);
		    
                    if (tt<1.0e-6) tt = -1.0;  
		  						}
	        			}
	    			} else {
              rr = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);
              rr1 = sqrt(xp[0]*xp[0] + xp[1]*xp[1]);
              theta = acos(xp[2]/rr);
              ix = static_cast<int> (theta*ndiv[kk][0]/M_PI);
              theta = acos(xp[0]/rr1);
              theta1 = asin(xp[1]/rr1);
              if (theta1 < 0.0)
                theta = 2.0*M_PI - theta;
              iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);
              for(j = 0; j < 3; ++j)
              {
                xh[j] -= vel[kk][ix][iy][j] * dtt;
                dl[j] -= vel[kk][ix][iy][j] * dtt;
              }
              d1 = sqrt(xh[0]*xh[0] + xh[1]*xh[1] + xh[2]*xh[2]);
              d2 = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2]);
              if ((d1-aa[kk][0])*(d2-aa[kk][0]) < 0.0)
                if (refl[kk] == 3 || (d2>aa[kk][0] && refl[kk] == 1) || (d2<aa[kk][0] && refl[kk] == 2)){
                  vv[0] = dl[0]*dl[0] + dl[1]*dl[1] + dl[2]*dl[2];
                  vv[1] = (xp[0]*dl[0] + xp[1]*dl[1] + xp[2]*dl[2])/vv[0];
                  vv[2] = (d2*d2-aa[kk][0]*aa[kk][0])/vv[0];
                  dot =  vv[1]*vv[1] - vv[2];
                  if (dot > 0.0){
                    if (vv[2] > 0) 
	              			tt = -vv[1] - sqrt(dot);
	            			else
	              			tt = -vv[1] + sqrt(dot);

                    if (tt<1.0e-6) tt = -1.0;  
	  	  					}
	        			}
	    			}
          
	    			if (tt>0.0 && tt<=1.0){ 
            	t_x[ind] = tt;
	        		i_x[ind] = kk;
	        		ind++;
	    			}
	  			}
				}
	
        if (ind > 0){
          for(j = 1; j < ind; j++){          
            tt = 2.0; 
            for(k = 0; k < ind-j+1; k++)
              if (t_x[k]<tt) {
                tt = t_x[k];
                kk = k;      
	      			}
            cc = i_x[kk];
            t_x[kk] = t_x[ind-j];
            i_x[kk] = i_x[ind-j];
            t_x[ind-j] = tt;
            i_x[ind-j] = cc;   
          }
	  			while(ind){
            bounce = 1;
            tt = t_x[ind - 1];
	    			kk = i_x[ind - 1];
            for(k = 0; k < 3; k++){
              xp[k] = x[i][k] - v[i][k]*dtt - x0[kk][k];
              vh[k] = v[i][k];
	    		}  
            rot_forward(xp[0],xp[1],xp[2],kk);
            rot_forward(vh[0],vh[1],vh[2],kk);
            for (k = 0; k < 3; k++)
	      			dd[k] = xp[k] + tt*vh[k]*dtt;
	    
            switch (ptype[kk]){
              case 1 : 
                u1 = dd[1]/aa[kk][2]; 
                u0 = (dd[0]-u1*aa[kk][1])/aa[kk][0];
                ix = 0;
                iy = 0;
                if (u0 < 0.0 || u1 < 0.0 || u0+u1 > 1.0)
                  bounce = 0;
                break;
              case 2 : 
                u1 = dd[1]/aa[kk][2]; 
                u0 = (dd[0]-u1*aa[kk][1])/aa[kk][0];
                ix = static_cast<int> (u0*ndiv[kk][0]);
                iy = static_cast<int> (u1*ndiv[kk][1]);
                if(ix < 0)
                ix = 0;
                if(ix > (ndiv[kk][0] - 1))
                ix = ndiv[kk][0] - 1;
                if(iy < 0)
                iy = 0;
                if(iy > (ndiv[kk][1] - 1))
                iy = ndiv[kk][1] - 1;
                
                if((period[kk][0] == 0) && ((u0 < 0.0)||(u0 > 1.0)))
                  bounce = 0;
                if((period[kk][0] == 1) && (u0 > 1.0))
                  bounce = 0;
                if((period[kk][0] == 2) && (u0 < 0.0))
                  bounce = 0;
                
                if((period[kk][1] == 0) && ((u1 < 0.0)||(u1 > 1.0)))
                  bounce = 0;
                if((period[kk][1] == 1) && (u1 > 1.0))
                  bounce = 0;
                if((period[kk][1] == 2) && (u1 < 0.0))
                  bounce = 0;
		
                break;
              case 3 :
                rr = sqrt(dd[0]*dd[0] + dd[1]*dd[1]);
                ix = static_cast<int> (dd[2]*ndiv[kk][0]/aa[kk][0]);
                
                if(ix < 0)
		  					ix = 0;
                if(ix > (ndiv[kk][0] - 1))
                ix = ndiv[kk][0] - 1;
		
                theta = acos(dd[0]/rr);
                theta1 = asin(dd[1]/rr);
                if (theta1 < 0.0)
                  theta = 2.0*M_PI - theta;
                iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);
                
                if(iy < 0)
                {
                  std::cout<<"in initial iy is "<<iy<<" theta is "<<theta<<" time step is "<<update->ntimestep<<" cpu is "<<comm->me<<"\n";
                  iy = 0;
                }
                if(iy > ndiv[kk][1] -1)
                {
                  std::cout<<"in initial iy is "<<iy<<" theta is "<<theta<<" time step is "<<update->ntimestep<<" cpu is "<<comm->me<<"\n";
                  iy = ndiv[kk][1] - 1;
                }
                
                if((period[kk][0] == 0) && ((dd[2] < 0)||(dd[2] > aa[kk][0])))
                  bounce = 0;
                if((period[kk][0] == 1) && (dd[2] > aa[kk][0]))
                  bounce = 0;
                if((period[kk][0] == 2) && (dd[2] < aa[kk][0]))
                  bounce = 0;

                norm[0] = dd[0]/rr;
                norm[1] = dd[1]/rr; 
                norm[2] = 0.0;
                break;

              case 4 :
                rr = sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);
                rr1 = sqrt(dd[0]*dd[0] + dd[1]*dd[1]);
                theta = acos(dd[2]/rr);
                ix = static_cast<int> (theta*ndiv[kk][0]/M_PI);
                theta = acos(dd[0]/rr1);
                theta1 = asin(dd[1]/rr1);
                if (theta1 < 0.0)
                  theta = 2.0*M_PI - theta;
                iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);
                for(k = 0; k < 3; k++)
                  norm[k] = dd[k]/rr;          
                break;
	    				}
            
	    			if (bounce){
              cplane = kk; 
              if (mirror){
                if (ptype[kk]<3){
                  // fix for the pressure
                  if(stress_ind == 2)
                  {
                    fs[0] = fabs(2.0*vel[kk][ix][iy][2]-2.0*vh[2]);
                    fs[0] = 0.5*fs[0]*pressure_grid/dtv;
                    fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=0;
                  }
                  //
                  vh[2] = 2.0 * vel[kk][ix][iy][2]-vh[2];
                } 
                else
                {
                  d1 = 2.0 * ((vh[0]-vel[kk][ix][iy][0])*norm[0] + (vh[1]-vel[kk][ix][iy][1])*norm[1] 
                      + (vh[2]-vel[kk][ix][iy][2])*norm[2]);
                  d2 = 2.0 * ((vh1[0]-vel[kk][ix][iy][0])*norm[0] + (vh1[1]-vel[kk][ix][iy][1])*norm[1] 
                      + (vh1[2]-vel[kk][ix][iy][2])*norm[2]);
                  for (k = 0; k < 3; k++){ 
                    vh[k] += -d1*norm[k];     
		  						}
								}
	      			}
              
              else {
                for (k = 0; k < 3; k++){
                  vh[k] = -vh[k] + 2.0*vel[kk][ix][iy][k];     
								}
	      			}
                
              for (k = 0; k < 3; k++){
                xp[k] = dd[k];      
                xh[k] = xp[k] + (1.0-tt)*dtt*vh[k];      // verlet
	      			}
              rot_back(xh[0],xh[1],xh[2],kk);
              rot_back(vh[0],vh[1],vh[2],kk);
              for (k = 0; k < 3; k++){
                x[i][k] = xh[k] + x0[kk][k];
                v[i][k] = vh[k];
              } 
              // fix for the pressure
              if(stress_ind == 2)
              {
                for(int k = 0; k < n_stress_tot; ++k)
                  modify->fix[stress_list[k]]->virial_f(i,fs);
              }
              
              dtt = (1.0-tt)*dtt;  
              ind = 0; 
	    			}
	    			else{
              ind--;
              if (ind == 0) 
                cond = 0;
	    			}
	  			}
				}
        else{
          cond = 0;
        }
      }
    }
  
 // Alireza: Periodic Deletion & Insertion 

#if 0
    if(inject){

    i = 0;
    vtemp = 1.0;
    pnum = pnum1 = 0;
    while(i < nlocal){
      if(mask[i] & groupbit) {
        inside = true;
        for(kk = 0; kk < num_shapes; ++kk)
        { 
          for(j = 0; j < 3; ++j)
            xh[j] = x[i][j] - x0[kk][j];
          rot_forward(xh[0],xh[1],xh[2],kk);
          if(ptype[kk] < 3)
            hh = xh[2];
          else if(ptype[kk] == 3)
          {
            rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1]);
            hh = aa[kk][1] - rr;
          }
          else
          {
            rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1] + xh[2]*xh[2]);
            hh = rr - aa[kk][0];
          }
          
          if(hh < 0)
          {
            inside = false;
            break;
          }
        }

        if(inside == false)
        {
          num_out = num_out + 1.0;
          randnum1 = ranmars0->uniform();
          k = 0;
          while(prob[k] < randnum1)
            k++;
          if(k >= prob_num){
            printf("some problem here, randnum1 is %f, prob[prob_num-1] is %f\n", randnum1, prob[prob_num-1]);
            k = prob_num - 1;
          }
          l = prob_ind[k][0];
          m = prob_ind[k][1];
          n = prob_ind[k][2];
          
          randnum1 = ranmars0->uniform();
          randnum2 = ranmars0->uniform();
          randnum3 = ranmars0->uniform();
          if(ptype[l] == 1)
          { 
            while((randnum1 + randnum2) > 1.0)
            {
              randnum1 = ranmars0->uniform();
              randnum2 = ranmars0->uniform();
            }
            xtemp = aa[l][0]*(m + randnum1) / static_cast<double>(ndiv[l][0]) + aa[l][1]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
            ytemp = aa[l][2]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
            ztemp = vel[l][m][n][2]*dtv*randnum3; 
          }
          else if(ptype[l] == 2)
          {
            xtemp = aa[l][0]*(m + randnum1) / static_cast<double>(ndiv[l][0]) + aa[l][1]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
            ytemp = aa[l][2]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
            ztemp = vel[l][m][n][2]*dtv*randnum3;
          }
          else if(ptype[l] == 3)
          {
            theta = (n+0.5)*2.0*M_PI/static_cast<double>(ndiv[l][1]);
            rr = vel[l][m][n][0]*cos(theta) + vel[l][m][n][1]*sin(theta);
            rr = rr * dtv * randnum1 + aa[l][1];
            theta = (n + randnum2)*2.0*M_PI/static_cast<double>(ndiv[l][1]);
            xtemp = rr * cos(theta);
            ytemp = rr * sin(theta);
            ztemp = (m + randnum3) * aa[l][0] / static_cast<double>(ndiv[l][0]); 
          }
          else 
          {
            theta = (m+0.5)*M_PI/static_cast<double>(ndiv[i][0]);
            phi = (n+0.5)*2.0*M_PI/static_cast<double>(ndiv[i][1]);
            rr = vel[l][m][n][0]*sin(theta)*cos(phi) + vel[l][m][n][1]*sin(theta)*sin(phi) + vel[l][m][n][2]*cos(theta);
            rr = rr * dtv * randnum1 + aa[l][0];
            theta = (m + randnum2) * M_PI/static_cast<double>(ndiv[l][0]);
            phi = (n + randnum3)*2.0*M_PI/static_cast<double>(ndiv[i][1]);
            xtemp = rr * sin(theta)*cos(phi);
            ytemp = rr * sin(theta)*sin(phi);
            ztemp = rr * cos(theta);
          }
            vxtemp = ranmars0->gaussian() + vel[l][m][n][0];
            vytemp = ranmars0->gaussian() + vel[l][m][n][1];
            vztemp = ranmars0->gaussian() + vel[l][m][n][2]; 
            //vxtemp = vel[l][m][n][0];
            //vytemp = vel[l][m][n][1];
            //vztemp = vel[l][m][n][2]; 
            rot_back(xtemp,ytemp,ztemp,l);
            rot_back(vxtemp,vytemp,vztemp,l);
            xtemp = x0[l][0] + xtemp;
            ytemp = x0[l][1] + ytemp;
            ztemp = x0[l][2] + ztemp;
            part_dat[pnum*8] = xtemp;
            part_dat[pnum*8 + 1] = ytemp;
            part_dat[pnum*8 + 2] = ztemp;
            part_dat[pnum*8 + 3] = vxtemp;
            part_dat[pnum*8 + 4] = vytemp;
            part_dat[pnum*8 + 5] = vztemp;
            part_dat[pnum*8 + 6] = static_cast<double>(tag[i]);
            part_dat[pnum*8 + 7] = static_cast<double>(atom->type[i]);
            atom->avec->copy(nlocal-1,i);
            nlocal--;
            pnum++;
        }//inside == false loop
        else
          i++;
      }//mask loop
      else
        i++;
    }//nlocal loop
    atom->nlocal = nlocal;
    MPI_Allgather(&pnum, 1, MPI_INT, rcounts, 1, MPI_INT, world);
    offset = 0;
    for( i = 0; i < comm->nprocs; ++i)
    {
      pnum1 += rcounts[i];
      displs[i] = offset;
      rcounts[i] *= 8;
      offset += rcounts[i];
    }
    if(pnum1 > NUM1) error->all("problem with number of particles, change NUM1 !!!!\n");
    MPI_Allgatherv(&part_dat, 8*pnum, MPI_DOUBLE, &part_comm, rcounts, displs, MPI_DOUBLE, world);
    for( i = 0; i < pnum1; ++i)
    {
      xx[0] = part_comm[i*8];
      xx[1] = part_comm[i*8+1];
      xx[2] = part_comm[i*8+2];
      if (xx[0] >= domain->sublo[0] && xx[0] < domain->subhi[0] && xx[1] >= domain->sublo[1] && xx[1] < domain->subhi[1] && xx[2] >= domain->sublo[2] && xx[2] < domain->subhi[2]) {
        l = static_cast<int>(part_comm[i*8 + 7]);
      
      atom->avec->create_atom(l,xx);
      k = atom->nlocal - 1;
      mask[k] = 1 | groupbit;
      tag[k] = static_cast<int>(part_comm[i*8 + 6]);
      v[k][0] = part_comm[i*8 + 3];
      v[k][1] = part_comm[i*8 + 4];
      v[k][2] = part_comm[i*8 + 5];
      }
    }
    
  }//inject loop

#else
    if(inject){
    i = 0;
		k = 0;
    while(i < nlocal){
			if(!adaptive)
				if(atom->type[i] != 2){
					i++;
					k = i;
					continue;
				}
      if(mask[i] & groupbit) {
        inside = true;
        for(m = 0; m < num_faces[k]; m++){
          kk = ind_faces[k][m];
					if (ind_coupled[kk] == 2){
          	for(j = 0; j < 3; ++j)
            	xh[j] = x[i][j] - x0[kk][j];
          	rot_forward(xh[0],xh[1],xh[2],kk);
          	if(ptype[kk] < 3)
            	hh = xh[2];
          	else if(ptype[kk] == 3)
          	{
            	rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1]);
            	hh = aa[kk][1] - rr;
          	}
          	else
          	{
            	rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1] + xh[2]*xh[2]);
            	hh = rr - aa[kk][0];
          	}
          	if(hh < 0){
            	inside = false;
            	break;
          	}
					}
        }
        if(inside == false){
          if(atom->type[i] == 1) num_out += 1.0; //Alireza: exclude other particles than plasma
          atom->avec->copy(nlocal-1,i);
					k = nlocal - 1;
          nlocal--;
        }//inside == false loop
        else{
          i++;
					k = i;
				}
      }//mask loop
      else{
        i++;
				k = i;
			}
    }//nlocal loop
    atom->nlocal = nlocal;

  	int *cell, count = 0, inout = 0;
		double in_plane[3] = {5.0, 0.0, 0.0}, a_S = 3.2, dist;
		cell = NULL;
		if (atom->molecule_flag){
  		int nm = atom->n_mol;
  		if (!cell) cell = new int[nm];
			cent_mass();	//Alireza: finding center of mass of cells before inserting platelets
  		for (i = 0; i < nm; i++){
//	if (comm->me==0) printf("cell no. = %d cell_x = %f\n",i,cent_all[i][0]);
    		if ( fabs(cent_all[i][0] - in_plane[0]) < 3.0*a_S) 
					{ cell[count] = i; count++; }
    		else	continue;
  		}
		}

#if 1
    //insert particle
    for(kk = 0; kk < prob_num; kk++){
      l = prob_ind[kk][0];
      m = prob_ind[kk][1];
      n = prob_ind[kk][2];
//	Alireza: hard-coded to insert platelets in the anulus
//			double radsq = x0[l][1]*x0[l][1]+x0[l][2]*x0[l][2];
//			double rad = 7.5;
//			if (sqrt(radsq) < 0.8*rad) continue;

      for(ll = 0; ll < num_at_types; ++ll){
      	if(!adaptive)
        	if(at_type[ll] == 1) continue;
        ncount[l][m][n][ll] += coeff_flux*at_dens[ll]*dtv*area[l][m][n]*vel[l][m][n][2];
//	if (comm->me == 0) printf("l m n ncount vel: %d %d %d %f %f\n", l ,m, n, ncount[l][m][n][ll], vel[l][m][n][2]);
        if (ncount[l][m][n][ll] < 0.0)
          ncount[l][m][n][ll] = 0.0;
        while (ncount[l][m][n][ll] > 0.99999999){
          ncount[l][m][n][ll] -= 1.0;
          if (comm->me == 0) {
            randnum1 = ranmars0->uniform();
            randnum2 = ranmars0->uniform();
            randnum3 = ranmars0->uniform();
            if(ptype[l] == 1)
            {
              while((randnum1 + randnum2) > 1.0)
              {
                randnum1 = ranmars0->uniform();
                randnum2 = ranmars0->uniform();
              }
              xtemp = aa[l][0]*randnum1 + aa[l][1]*randnum2;
              ytemp = aa[l][2]*randnum2;
              ztemp = vel[l][m][n][2]*dtv*randnum3;
            }
            else if(ptype[l] == 2)
            {
              xtemp = aa[l][0]*(m + randnum1) / static_cast<double>(ndiv[l][0]) + aa[l][1]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
              ytemp = aa[l][2]*(n + randnum2) / static_cast<double>(ndiv[l][1]);
              ztemp = vel[l][m][n][2]*dtv*randnum3;
            }
            else if(ptype[l] == 3)
            {
              theta = (n+0.5)*2.0*M_PI/static_cast<double>(ndiv[l][1]);
              rr = vel[l][m][n][0]*cos(theta) + vel[l][m][n][1]*sin(theta);
              rr = rr * dtv * randnum1 + aa[l][1];
              theta = (n + randnum2)*2.0*M_PI/static_cast<double>(ndiv[l][1]);
              xtemp = rr * cos(theta);
              ytemp = rr * sin(theta);
              ztemp = (m + randnum3) * aa[l][0] / static_cast<double>(ndiv[l][0]);
            }
            else
            {
              theta = (m+0.5)*M_PI/static_cast<double>(ndiv[i][0]);
              phi = (n+0.5)*2.0*M_PI/static_cast<double>(ndiv[i][1]);
              rr = vel[l][m][n][0]*sin(theta)*cos(phi) + vel[l][m][n][1]*sin(theta)*sin(phi) + vel[l][m][n][2]*cos(theta);
              rr = rr * dtv * randnum1 + aa[l][0];
              theta = (m + randnum2) * M_PI/static_cast<double>(ndiv[l][0]);
              phi = (n + randnum3)*2.0*M_PI/static_cast<double>(ndiv[i][1]);
              xtemp = rr * sin(theta)*cos(phi);
              ytemp = rr * sin(theta)*sin(phi);
              ztemp = rr * cos(theta);
            }
            rot_back(xtemp, ytemp, ztemp, l);
            xx[0] = xtemp + x0[l][0];
            xx[1] = ytemp + x0[l][1];
            xx[2] = ztemp + x0[l][2];

						if (atom->molecule_flag)
							inout = in_cell(xx,cell,count-1,a_S);
							
/*
					// Alireza: hardcoded for micro-channel / tube
						double height = (max_-min_)/2.0, cent_ = (min_+max_)/2.0;
						if (fabs(xx[2]-cent_)/height < 0.8){
							if ((xx[2]-cent_) > 0)
								xx[2] += 0.75*height;
							else
								xx[2] -= 0.75*height;
						}
*/
//					////////////////////////////////////////
					}
					if (atom->molecule_flag){
						MPI_Bcast(&inout, 1, MPI_INT, 0, world);
						if (inout){	ncount[l][m][n][ll] += 1.0; break; }
					}

          MPI_Bcast(&xx[0], 3, MPI_DOUBLE, 0, world);
      	  if (xx[0] >= domain->sublo[0] && xx[0] < domain->subhi[0] && xx[1] >= domain->sublo[1]
            && xx[1] < domain->subhi[1] && xx[2] >= domain->sublo[2] && xx[2] < domain->subhi[2]) {

            atom->avec->create_atom(at_type[ll], xx);
            k = atom->nlocal - 1;
            mask[k] = 1 | groupbit;
            if (groupbit != at_groupbit[ll])
              mask[k] |= at_groupbit[ll];

            rr = sqrt(kbt);
            vxtemp = rr*ranmars0->gaussian() + vel[l][m][n][0];
            vytemp = rr*ranmars0->gaussian() + vel[l][m][n][1];
            vztemp = rr*ranmars0->gaussian() + vel[l][m][n][2]; 

            rot_back(vxtemp, vytemp, vztemp, l);
            v[k][0] = vxtemp;
            v[k][1] = vytemp;
            v[k][2] = vztemp;
            if (at_type[ll] == 1) num_in += 1.0;
          }
        }
      }
    }

//*************** Problem Specific (Plane Channel Flow) ********************//    
#else
    //insert particle
    for(i = 0; i < num_shapes; i++){
      for (j=0; j<ndiv[i][1]; j++){
      ncount[i][j] += 3.0*dtv*0.5*10*vel[i][0][j][2];
      if (ncount[i][j]<0) ncount[i][j]=0;
      while (ncount[i][j]>0.99999999){
        ncount[i][j] -= 1.0;
        if(comm->me == 0) {
        randnum1 = ranmars0->uniform();
        randnum2 = ranmars0->uniform();
        randnum3 = ranmars0->uniform();
        xtemp = aa[i][0]*(0 + randnum1) / static_cast<double>(ndiv[i][0]) 
          + aa[i][1]*(j + randnum2) / static_cast<double>(ndiv[i][1]);
        ytemp = aa[i][2]*(j + randnum2) / static_cast<double>(ndiv[i][1]);
        ztemp = vel[i][0][j][2]*dtv*randnum3; 
        rot_back(xtemp, ytemp, ztemp, i);
        
        xx[0] = xtemp + x0[i][0];
        xx[1] = ytemp + x0[i][1];
        xx[2] = ztemp + x0[i][2];
        }

        MPI_Bcast(xx, 3, MPI_DOUBLE, 0, world);
        if (xx[0] >= domain->sublo[0] && xx[0] < domain->subhi[0] && xx[1] >= domain->sublo[1] 
          && xx[1] < domain->subhi[1] && xx[2] >= domain->sublo[2] && xx[2] < domain->subhi[2]) {

        atom->avec->create_atom(1,xx);
        k = atom->nlocal - 1;
        mask[k] = 1 | groupbit;
          
        //vxtemp = ranmars0->gaussian() + vel[i][0][j][0];
        //vytemp = ranmars0->gaussian() + vel[i][0][j][1];
        //vztemp = ranmars0->gaussian() + vel[i][0][j][2]; 

        vxtemp = vel[i][0][j][0];
        vytemp = vel[i][0][j][1];
        vztemp = vel[i][0][j][2]; 

        rot_back(vxtemp, vytemp, vztemp, i);
        v[k][0] = vxtemp;
        v[k][1] = vytemp;
        v[k][2] = vztemp;
        
      }
      }   
      }
    }  
#endif
  }//inject loop

	if (inject){
		if (!adaptive){
			ntemp = 0.0;
			for (i = 0; i < atom->nlocal; i++)
      	if (atom->type[i] == 2) ntemp += 1.0;
			MPI_Allreduce(&ntemp,&nplat,1,MPI_DOUBLE,MPI_SUM,world);
			if(!comm->me && update->ntimestep % 5000 == 0)
				printf("Number of platelets: %d in step %d...\n", static_cast<int>(nplat), update->ntimestep);
		}else{
	  ntemp = static_cast<double>(atom->nlocal);
  	MPI_Allreduce(&ntemp,&atom->natoms,1,MPI_DOUBLE,MPI_SUM,world);
  	ntemp = 0.0;
  	for (i = 0; i < atom->nlocal; i++)
    	if (atom->type[i] == 1)	ntemp += 1.0;
  	MPI_Allreduce(&ntemp,&nfluid,1,MPI_DOUBLE,MPI_SUM,world);
  	if(!comm->me && update->ntimestep % iter == 0)
    	printf("Total number of atoms: %d & fluid particles: %d in step %d...\n", static_cast<int>(atom->natoms), static_cast<int>(nfluid), update->ntimestep);
  	num_accum += nfluid;
		}
#endif 
	}
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::pre_force(int vflag)
{
#if 1 
  if (atom->nlocal > nmax){
	  nmax = atom->nlocal;
		memory->destroy_2d_int_array(ind_faces);
    ind_faces = memory->create_2d_int_array(nmax,max_faces,"fix_solid_bound:ind_faces");
		memory->sfree(num_faces);
   	num_faces = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces");
		memory->destroy_2d_int_array(ind_faces_flux);
   	ind_faces_flux = memory->create_2d_int_array(nmax,max_faces,"fix_solid_bound:ind_faces_flux");
		memory->sfree(num_faces_flux);
   	num_faces_flux = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_flux");
		memory->sfree(num_faces_wal);
   	num_faces_wal = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_wal");
		memory->sfree(num_faces_out);
  	num_faces_out = (int *) memory->smalloc(nmax*sizeof(int),"fix_solid_bound:num_faces_out");
	}
#endif

  if (neighbor->ago == 0)
		face_decide();
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::post_force(int vflag)
{
	if (ind_press) post_force_vel(vflag);

	if (ind_cc) post_force_cc(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::end_of_step()
{
	int i, j, k, l, iz, index;
  double **x = atom->x;
  double **v = atom->v;
  double **T = atom->T;
	double *q = atom->q;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;
  int step_t = update->ntimestep;

	if (ind_space == 2){

  int total = div[0]*div[1]*div[2];

  if (step_t > nstat[0]){
    for (i = 0; i < atom->nlocal; i++)
      if (mask[i] & groupbit)
        if (map_index(x[i][0],x[i][1],x[i][2])){
          num_avg[is][js][ks] += 1.0;
					for (l = 0; l < nstat[3]; l++){
						index = spec[l];
	          aT[l][is][js][ks] += T[i][index];
					}
        }
		num_step++;
    if (step_t%nstat[2] == 0){
  		for (i = 0; i < div[0]; i++)
		    for (j = 0; j < div[1]; j++)
		      for (k = 0; k < div[2]; k++){
     				tmp[i][j][k] = num_avg[i][j][k]/num_step;
		        num_avg[i][j][k] = 0.0;
     			  ntmp[i][j][k] = 0.0;
						for (l = 0; l < nstat[3]; l++)
				       aTtmp[l][i][j][k] = 0.0;
      		}
  		MPI_Allreduce(&tmp[0][0][0],&ntmp[0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
			for (l = 0; l < nstat[3]; l++){
	      for (i = 0; i < div[0]; i++)
  	      for (j = 0; j < div[1]; j++)
    	      for (k = 0; k < div[2]; k++){
            	tmp[i][j][k] = aT[l][i][j][k]/num_step;
							aT[l][i][j][k] = 0.0;
						}
				MPI_Allreduce(&tmp[0][0][0],&aTtmp[l][0][0][0],total,MPI_DOUBLE,MPI_SUM,world);
			}
  		num_step = 0;
      for (l = 0; l < nstat[3]; l++)
        for (i = 0; i < div[0]; i++)
          for (j = 0; j < div[1]; j++)
            for (k = 0; k < div[2]; k++)
              aTtmp[l][i][j][k] = aTtmp[l][i][j][k]/ntmp[i][j][k];
		}
	}

	}
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::setup_rot(int id, double a1[], double a2[])
{
  int i,j;
  double norm[3],mr[3][3],a3[2];
  double nr,rx,theta;

  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++){
      rot[id][i][j] = 0.0;
      mr[i][j] = 0.0;  
    }

  if (ptype[id] < 3){  
    norm[0] = a1[1]*a2[2] - a1[2]*a2[1];
    norm[1] = a1[2]*a2[0] - a1[0]*a2[2];
    norm[2] = a1[0]*a2[1] - a1[1]*a2[0]; 
    nr = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    rx = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    if (rx > EPS){ 
      mr[0][0] = norm[0]*norm[2]/nr/rx;
      mr[0][1] = norm[1]*norm[2]/nr/rx;
      mr[0][2] = -rx/nr;
      mr[1][0] = -norm[1]/rx;
      mr[1][1] = norm[0]/rx;
      mr[1][2] = 0.0;
      mr[2][0] = norm[0]/nr;
      mr[2][1] = norm[1]/nr;
      mr[2][2] = norm[2]/nr;
    }
    else {
      mr[0][0] = 1.0;
      mr[1][1] = 1.0;
      mr[2][2] = 1.0;
      if (norm[2] < 0.0)
        mr[2][2] = -1.0;
    }
    a3[0] = mr[0][0]*a1[0] + mr[0][1]*a1[1] + mr[0][2]*a1[2];
    a3[1] = mr[1][0]*a1[0] + mr[1][1]*a1[1] + mr[1][2]*a1[2];
    rx = sqrt(a3[0]*a3[0] + a3[1]*a3[1]);
    rot[id][0][0] = (a3[0]*mr[0][0] + a3[1]*mr[1][0])/rx;
    rot[id][0][1] = (a3[0]*mr[0][1] + a3[1]*mr[1][1])/rx;
    rot[id][0][2] = (a3[0]*mr[0][2] + a3[1]*mr[1][2])/rx;
    rot[id][1][0] = (-a3[1]*mr[0][0] + a3[0]*mr[1][0])/rx;
    rot[id][1][1] = (-a3[1]*mr[0][1] + a3[0]*mr[1][1])/rx;
    rot[id][1][2] = (-a3[1]*mr[0][2] + a3[0]*mr[1][2])/rx;
    rot[id][2][0] = mr[2][0];
    rot[id][2][1] = mr[2][1];
    rot[id][2][2] = mr[2][2];
    aa[id][0] = rx;
    aa[id][1] = rot[id][0][0]*a2[0] + rot[id][0][1]*a2[1] + rot[id][0][2]*a2[2];
    aa[id][2] = rot[id][1][0]*a2[0] + rot[id][1][1]*a2[1] + rot[id][1][2]*a2[2];
    aa[id][3] = sqrt(aa[id][1]*aa[id][1] + aa[id][2]*aa[id][2]);
    if(ptype[id] == 1)
      area[id][0][0] = nr / 2.0;
    else
      for(i = 0; i < ndiv[id][0]; ++i)
        for(j = 0; j < ndiv[id][1]; ++j)
          area[id][i][j] = nr / ndiv[id][0] / ndiv[id][1];
  } else if (ptype[id] == 3){
    norm[0] = a1[0];
    norm[1] = a1[1];
    norm[2] = a1[2];
    nr = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    rx = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    if (rx > EPS){ 
      rot[id][0][0] = norm[0]*norm[2]/nr/rx;
      rot[id][0][1] = norm[1]*norm[2]/nr/rx;
      rot[id][0][2] = -rx/nr;
      rot[id][1][0] = -norm[1]/rx;
      rot[id][1][1] = norm[0]/rx;
      rot[id][1][2] = 0.0;
      rot[id][2][0] = norm[0]/nr;
      rot[id][2][1] = norm[1]/nr;
      rot[id][2][2] = norm[2]/nr;
    } else {
      rot[id][0][0] = 1.0;
      rot[id][1][1] = 1.0;
      rot[id][2][2] = 1.0;
      if (norm[2] < 0.0)
        rot[id][2][2] = -1.0;
    }
    aa[id][0] = nr;
    for(i = 0; i < ndiv[id][0]; ++i)
      for(j = 0; j < ndiv[id][1]; ++j)
        area[id][i][j] = 2.0*M_PI*aa[id][1]*aa[id][0]/ndiv[id][0]/ndiv[id][1]; 
  } else {
    for(i = 0; i < 3; i++)
      rot[id][i][i] = 1.0;
    for(i = 0; i < ndiv[id][0]; ++i)
    {
      theta = M_PI*static_cast<double>(i)/static_cast<double>(ndiv[id][0]);
      if(theta < 0.5*M_PI)
        area[id][i][0] = 2.0*M_PI*aa[id][0]*(1.0 - cos(theta))*aa[id][0];
      else
        area[id][i][0] = 2.0*M_PI*aa[id][0]*(1.0 + cos(theta))*aa[id][0];
    }
    for(i = 0; i < ndiv[id][0]; ++i)
      {
        if(i < ndiv[id][0] - 1)
          area[id][i][0] = area[id][i+1][0] - area[id][i][0];
        else
          area[id][i][0] = 4.0*M_PI*aa[id][0]*aa[id][0] - area[id][i][0];
        area[id][i][0] /= ndiv[id][1];
      }
    for(i = 0; i < ndiv[id][0]; ++i)
      for(j = 1; j < ndiv[id][1]; ++j)
        area[id][i][j] = area[id][i][0];
  }
} 

/* ---------------------------------------------------------------------- */

void FixSolidBound::rot_forward(double &x, double &y, double &z, int id)
{
  double x1,y1,z1;

  x1 = x; y1 = y; z1 = z;
  x = rot[id][0][0]*x1 + rot[id][0][1]*y1 + rot[id][0][2]*z1;
  y = rot[id][1][0]*x1 + rot[id][1][1]*y1 + rot[id][1][2]*z1;
  z = rot[id][2][0]*x1 + rot[id][2][1]*y1 + rot[id][2][2]*z1;
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::rot_back(double &x, double &y, double &z, int id)
{
  double x1,y1,z1;

  x1 = x; y1 = y; z1 = z;
  x = rot[id][0][0]*x1 + rot[id][1][0]*y1 + rot[id][2][0]*z1;
  y = rot[id][0][1]*x1 + rot[id][1][1]*y1 + rot[id][2][1]*z1;
  z = rot[id][0][2]*x1 + rot[id][1][2]*y1 + rot[id][2][2]*z1;
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::face_decide()
{
  int i,j,k,l,m;
  double dd[3],xd[2],vv[3];
  double u_inv,u0,u1;

  double **x = atom->x;
  int *mask = atom->mask;
  int **ind_tmp1, **ind_tmp2;

  for (i=0; i<atom->nlocal; i++){
    num_faces[i] = num_faces_flux[i] = 0; 
    if (mask[i] & groupbit){
//        if (x[i][0] < 110. || x[i][0] > 150.)
//					if ((x[i][1]*x[i][1]+x[i][2]*x[i][2]) < (RAD-r_cut)*(RAD-r_cut)) continue; // Alireza: hard-coded to expedite the search process; geometry is a tube w/ x-axis of symmetry
//        if (x[i][0] < 60.0 || x[i][0] > 105.0)
// 	        if ( fabs(x[i][1]-DY/2.0) < (DY/2.0-r_cut) && fabs(x[i][2]-DZ/2.0) < (DZ/2.0-r_cut) ) continue;
      if ( fabs(x[i][1]-DY/2.0) < (DY/2.0-r_cut) ) continue;
      for (j = 0; j < num_shapes; j++){
        for (k = 0; k < 3; k++)
          dd[k] = x[i][k] - x0[j][k];
        rot_forward(dd[0],dd[1],dd[2],j);
        switch (ptype[j]){
          case 1 : 
            u1 = dd[1]/aa[j][2]; 
            u0 = (dd[0]-u1*aa[j][1])/aa[j][0]; 
            if (u0 > -r_cutshear/aa[j][0] && u1 > -r_cutshear/aa[j][3] && u0+u1 < 1.0 + r_cutshear/aa[j][0] + r_cutshear/aa[j][3] && dd[2] < r_cut && dd[2] > -r_cut){
              ind_faces[i][num_faces[i]] = j;
              num_faces[i]++;
	    			}
            if (u0 > -r_cutflux/aa[j][0] && u1 > -r_cutflux/aa[j][3] && u0+u1 < 1.0 + r_cutflux/aa[j][0] + r_cutflux/aa[j][3] && dd[2] < r_cut && dd[2] > -r_cut){
              ind_faces_flux[i][num_faces_flux[i]] = j;
              num_faces_flux[i]++;
	    			}
            break;
          case 2 : 
            u1 = dd[1]/aa[j][2]; 
            u0 = (dd[0]-u1*aa[j][1])/aa[j][0];   
            if (u0 > -r_cutshear/aa[j][0] && u1 > -r_cutshear/aa[j][3] && u0 < 1.0 + r_cutshear/aa[j][0] && u1 < 1.0 + r_cutshear/aa[j][3] && dd[2] < r_cut && dd[2] > -r_cut){
              ind_faces[i][num_faces[i]] = j;
              num_faces[i]++;
	    			}           
            if (u0 > -r_cutflux/aa[j][0] && u1 > -r_cutflux/aa[j][3] && u0 < 1.0 + r_cutflux/aa[j][0] && u1 < 1.0 + r_cutflux/aa[j][3] && dd[2] < r_cut && dd[2] > -r_cut){
              ind_faces_flux[i][num_faces_flux[i]] = j;
              num_faces_flux[i]++;
	    			}           
            break;
          case 3 : 
            xd[0] = sqrt(dd[0]*dd[0] + dd[1]*dd[1]);
            if (xd[0] > aa[j][1]-r_cut && xd[0] < aa[j][1]+r_cut && dd[2] > -r_cut && dd[2] < aa[j][0]+r_cut){
              ind_faces[i][num_faces[i]] = j;
              num_faces[i]++;
	    			}             
            break;
          case 4 :
            xd[0] = sqrt(dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2]);
            if (xd[0] > aa[j][0]-r_cut && xd[0] < aa[j][0]+r_cut){
              ind_faces[i][num_faces[i]] = j;
              num_faces[i]++;
	    			}             
            break;
				}  
        if (num_faces[i] == max_faces || num_faces_flux[i] == max_faces) {
          m = max_faces;
          ind_tmp1 = memory->create_2d_int_array(nmax,m,"fix_solid_bound:ind_tmp1");
          ind_tmp2 = memory->create_2d_int_array(nmax,m,"fix_solid_bound:ind_tmp2");
          for (k=0; k<nmax; k++)
            for (l=0; l<m; l++){
              ind_tmp1[k][l] = ind_faces[k][l];
              ind_tmp2[k][l] = ind_faces_flux[k][l];
						}
          max_faces += 5;
          ind_faces = memory->grow_2d_int_array(ind_faces,nmax,max_faces,"fix_solid_bound:ind_faces");
          ind_faces_flux = memory->grow_2d_int_array(ind_faces_flux,nmax,max_faces,"fix_solid_bound:ind_faces_flux");
          for (k=0; k<nmax; k++)
            for (l=0; l<m; l++){
              ind_faces[k][l] = ind_tmp1[k][l];
              ind_faces_flux[k][l] = ind_tmp2[k][l];
						}
          memory->destroy_2d_int_array(ind_tmp1);
          memory->destroy_2d_int_array(ind_tmp2);
				}
      }
    }
  }

#if 0
    if (comm->me == 0)
    {
      ofstream f_stream;
      f_stream.open("faces.dat");
      if(f_stream.is_open())
      {
        for(i = 0; i < atom->nlocal; ++i){
          for(j = 0; j < num_faces[i]; ++j)
          {
            if (num_faces[i] != 0) f_stream << "atom: " << i << " face: " << ind_faces[i][j] << " index: " << ind_coupled[ind_faces[i][j]] << endl;
          }
        }
      }
      f_stream.close();
    }
#endif

#if 1
  for (i=0; i < atom->nlocal; i++){
    l = m = 0;
    if (mask[i] & groupbit){
      for (j = 0; j < num_faces[i]; j++){
        k = ind_faces[i][j];
        if (ind_coupled[k] == -1) l++;
        if (ind_coupled[k] == 0 || ind_coupled[k] == 2) m++;
      }
    }
    num_faces_wal[i] = l;
    num_faces_out[i] = m;
  }
#endif
} 

/* ---------------------------------------------------------------------- */

void FixSolidBound::recalc_force()
{
  int i,j,k,l,m,kk,shift;
  double tmp[numt],ntmp[numt],vv;
  double ssum, theta, phi, veln;
  int index = update->ntimestep % dump_each;
	int step_t = update->ntimestep;
  
  index = index / iter;
  //fac_ave[index] = factor;

 if (ind_shear){
  l = 0;
  for(i = 0; i < num_shapes; i++)
    if (ind_coupled[i] == -1){
      for (j=0; j<ndiv[i][0]; j++)
        for (k=0; k<ndiv[i][1]; k++)
          for (m=0; m<2*n_per; m++){
            tmp[l] = num[i][m][j][k];
            tmp[l+1] = velx[i][m][j][k];
            tmp[l+2] = vely[i][m][j][k];
            tmp[l+3] = velz[i][m][j][k];
            num[i][m][j][k] = 0.0;
            velx[i][m][j][k] = 0.0;
            vely[i][m][j][k] = 0.0;
            velz[i][m][j][k] = 0.0;
            ntmp[l] = 0.0;
            ntmp[l+1] = 0.0;
            ntmp[l+2] = 0.0;
            ntmp[l+3] = 0.0;
            l += 4;
          }
    }

  if (l != numt) error->all("Wrong count in fix_solid_bound command");
  MPI_Allreduce(&tmp,&ntmp,l,MPI_DOUBLE,MPI_SUM,world);

  l = 0;
  for(i = 0; i < num_shapes; i++)
    if (ind_coupled[i] == -1){
      for (j=0; j<ndiv[i][0]; j++)
        for (k=0; k<ndiv[i][1]; k++)
          for (m=0; m<2*n_per; m++){
            if (m == 0 && ntmp[l] > 0.0 && ntmp[l+4] > 0.0){
              vv=0.5*(3.0*ntmp[l+1]/ntmp[l]-ntmp[l+5]/ntmp[l+4]);
              fsx[i][0][j][k] += coeff*(vv-vel[i][j][k][0]);
              vv=0.5*(3.0*ntmp[l+2]/ntmp[l]-ntmp[l+6]/ntmp[l+4]);
              fsy[i][0][j][k] += coeff*(vv-vel[i][j][k][1]);
              vv=0.5*(3.0*ntmp[l+3]/ntmp[l]-ntmp[l+7]/ntmp[l+4]);
              fsz[i][0][j][k] += coeff*(vv-vel[i][j][k][2]);
            }
            if (m == n_per && ntmp[l] > 0.0 && ntmp[l+4] > 0.0){
              vv=0.5*(3.0*ntmp[l+1]/ntmp[l]-ntmp[l+5]/ntmp[l+4]);
              fsx[i][1][j][k] += coeff*(vv-vel[i][j][k][0]);
              vv=0.5*(3.0*ntmp[l+2]/ntmp[l]-ntmp[l+6]/ntmp[l+4]);
              fsy[i][1][j][k] += coeff*(vv-vel[i][j][k][1]);
              vv=0.5*(3.0*ntmp[l+3]/ntmp[l]-ntmp[l+7]/ntmp[l+4]);
              fsz[i][1][j][k] += coeff*(vv-vel[i][j][k][2]);
            }
            l += 4;
          }
    }

 } else if (ind_out){
  l = 0;
  for(i = 0; i < num_shapes; i++)
      for (j=0; j<ndiv[i][0]; j++)
        for (k=0; k<ndiv[i][1]; k++)
          for (m=0; m<2*n_per; m++){
            tmp[l] = num[i][m][j][k];
            tmp[l+1] = velx[i][m][j][k];
            tmp[l+2] = vely[i][m][j][k];
            tmp[l+3] = velz[i][m][j][k];
            num[i][m][j][k] = 0.0;
            velx[i][m][j][k] = 0.0;
            vely[i][m][j][k] = 0.0;
            velz[i][m][j][k] = 0.0;
            ntmp[l] = 0.0;
            ntmp[l+1] = 0.0;
            ntmp[l+2] = 0.0;
            ntmp[l+3] = 0.0;
            l += 4;
	  			}

  MPI_Allreduce(&tmp,&ntmp,l,MPI_DOUBLE,MPI_SUM,world);
  l = 0;
  for(i = 0; i < num_shapes; i++)
      for (j=0; j<ndiv[i][0]; j++)
        for (k=0; k<ndiv[i][1]; k++)
          for (m=0; m<2*n_per; m++){
            if (m == 0 && ntmp[l] > 0.0 && ntmp[l+4] > 0.0){
              vv = ntmp[l+3]/ntmp[l]-ntmp[l+7]/ntmp[l+4];
              beta[i][0][j][k] = coeff*vv;
              dens[i][0][j][k] = (ntmp[l]+ntmp[l+4])/(area[i][j][k]*r_out)/iter;
            }
            l += 4;
          }

  MPI_Allreduce(&num_out, &num_out_total, 1, MPI_DOUBLE, MPI_SUM,  world);
  MPI_Allreduce(&num_in, &num_in_total, 1, MPI_DOUBLE, MPI_SUM,  world);
  num_out = num_in = 0.0;
  num_out_total = num_out_total / (iter * update->dt);
  num_in_total = num_in_total / (iter * update->dt);

  if (num_accum/iter < natoms_0*0.9995)
    coeff_flux += fabs(num_accum/iter - natoms_0)/natoms_0;
  else if (num_accum/iter > natoms_0*1.0005)
    coeff_flux -= fabs(num_accum/iter - natoms_0)/natoms_0;

  coeff_flux = 1.0;

/*
	if (step_t < 100000){
		if (step_t%10000 == 0 && step_t != 0)
			coeff_flux += 0.1;
	}else 
		coeff_flux = 1.0;
	num_out_desire = num_in_total;
*/

  if(fabs(num_out_total - num_out_desire) > 1.0)
    factor += coeff_press * (num_out_total - num_out_desire);

  if(numdiff < (Nloop + 1))
    ++numdiff;

  //determine num_accum by total num of particle of the system  
  if((num_accum/iter < natoms_0*0.9995) && numdiff == (Nloop + 1))
    ++Nloop;
  if((num_accum/iter > natoms_0*1.0005) && numdiff == (Nloop + 1) && (Nloop > 1))
    --Nloop;

  shift = numdiff - Nloop;
  for (kk = 0; kk < num_shapes; kk++){
    if (ind_coupled[kk] == 2)
      if(shift < 1){
        for(i = 0; i < ndiv[kk][0]; ++i)
          for(j = 0; j < ndiv[kk][1]; ++j){
            accum[kk][numdiff-1][i][j] = beta[kk][0][i][j];
            beta[kk][0][i][j] = 0.0;
          }
  
        for(k = 0; k < numdiff; ++k)
          for(i = 0; i < ndiv[kk][0]; ++i)
            for(j = 0; j < ndiv[kk][1]; ++j)
              beta[kk][0][i][j] += accum[kk][k][i][j];
            
      } else{
        for(k = 0; k < Nloop-1; ++k)
          for(i = 0; i < ndiv[kk][0]; ++i)
            for(j = 0; j < ndiv[kk][1]; ++j)
              accum[kk][k][i][j] = accum[kk][k+shift][i][j];
            
        for(i = 0; i < ndiv[kk][0]; ++i)
          for(j = 0; j < ndiv[kk][1]; ++j){
            accum[kk][Nloop-1][i][j] = beta[kk][0][i][j];
	    			beta[kk][0][i][j] = 0.0;
          }

        for(k = 0; k < Nloop; ++k)
          for(i = 0; i < ndiv[kk][0]; ++i)
            for(j = 0; j < ndiv[kk][1]; ++j)
              beta[kk][0][i][j] += accum[kk][k][i][j];
    }
  }  
  if(!comm->me){
  	printf("num_in_total: %d, num_out_total: %d, Nloop: %d, shift: %d, factor: %lf\n", static_cast<int>(num_in_total), static_cast<int>(num_out_total), Nloop, shift, factor);
  }

  //reset numdiff to Nloop + 1
  if (shift >= 1)
    numdiff = Nloop + 1;
  //reset num_accum to zero
  num_accum = 0.0;
 }
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::couple_run()
{
  int i, ishape, j, k;
  double vx, vy, vz;
  for(i = 0; i < 3*num_points; i=i+3)
  {
    ishape = info_coupled[i/3][0];  
    j = info_coupled[i/3][1];
    k = info_coupled[i/3][2];
    vx = vcoupled[i];
    vy = vcoupled[i+1];
    vz = vcoupled[i+2];
    rot_forward(vx, vy, vz, ishape);
    vel[ishape][j][k][0] = vx;
    vel[ishape][j][k][1] = vy;
    vel[ishape][j][k][2] = vz;
  }
}

/*---------------------------------------------------------------------- */
double FixSolidBound::vfield(double &y, double &t)
{
  double u1 = dpressure * sin(frequency*t) / frequency / 3.0;
  double u2 = (xx1*cos(wom*y)*cosh(wom*y) + xx2*sin(wom*y)*sinh(wom*y))*sin(frequency*t);
  double u3 = (xx2*cos(wom*y)*cosh(wom*y) - xx1*sin(wom*y)*sinh(wom*y))*cos(frequency*t);
  double u4 = -dpressure * (u2 - u3)/(xx1*xx1 + xx2*xx2) / frequency / 3.0;

  double u = u1 + u4;
  return u;
}

/*---------------------------------------------------------------------- */
void FixSolidBound::cent_mass()
{
	int i1, i2, i3, type, i, j, m, n;
	double xx1[3], xx2[3], xx3[3];
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nm = atom->n_mol;
	double **x = atom->x;

	if (!counter){
  	cent_loc = memory->create_2d_double_array(nm,3,"fix_solid_bound:cent_loc");
  	cent_all = memory->create_2d_double_array(nm,3,"fix_solid_bound:cent_all");
    counter = (int *) memory->smalloc(nm*sizeof(int),"fix_solid_bound:counter");
    angle_per_mol = (int *) memory->smalloc(nm*sizeof(int),"fix_solid_bound:angle_per_mol");
	}

 	for (i = 0; i < nm; i++){
		counter[i] = 0;
		angle_per_mol[i] = 0;
  	for (j = 0; j < 3; j++){
    	cent_loc[i][j] = 0.0;
    	cent_all[i][j] = 0.0;
		}
	}

  for (n = 0; n < nanglelist; n++){
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    m = atom->molecule[i1]-1;

   	if (type == 1){
      for (j = 0; j < 3; j++){
        xx1[j] = x[i1][j];
        xx2[j] = x[i2][j];
        xx3[j] = x[i3][j];
      }
      cent_loc[m][0] +=  xx1[0] + xx2[0] + xx3[0];
      cent_loc[m][1] +=  xx1[1] + xx2[1] + xx3[1];
      cent_loc[m][2] +=  xx1[2] + xx2[2] + xx3[2];
			counter[m]++;
		}
	}
	MPI_Allreduce(&cent_loc[0][0],&cent_all[0][0],nm*3,MPI_DOUBLE,MPI_SUM,world);
	MPI_Allreduce(&counter[0],&angle_per_mol[0],nm,MPI_INT,MPI_SUM,world);
	
	for (i = 0; i < nm; i++)
		for (j = 0; j < 3; j++)
			cent_all[i][j] /= 3.0*angle_per_mol[i];
}

/*---------------------------------------------------------------------- */
int FixSolidBound::in_cell(double *in_loc, int *cell, int count, double a_S)
{
	if (count == 0) return 0;

	int id, m, inside = 0;
	double dist;

	for (m = 0; m < count; m++){
		id = cell[m]; 
    dist = sqrt( (cent_all[id][0]-in_loc[0])*(cent_all[id][0]-in_loc[0])
          + (cent_all[id][1]-in_loc[1])*(cent_all[id][1]-in_loc[1]) + (cent_all[id][2]-in_loc[2])*(cent_all[id][2]-in_loc[2]));
		if (dist < 1.25*a_S){ inside = 1; break; }
	}

	return inside;
}

/* ---------------------------------------------------------------------- */

int FixSolidBound::map_index(double x, double y, double z)
{
  int ind = 0;

  if (x<dxlo || x>=dxhi || y<dylo || y>=dyhi || z<dzlo || z>=dzhi){
    if (xper) {
      while (x >= dxhi)
        x -= dxs;
      while (x < dxlo)
        x += dxs;
    }
    if (yper) {
      while (y >= dyhi)
        y -= dys;
      while (y < dylo)
        y += dys;
    }
    if (zper) {
      while (z >= dzhi)
        z -= dzs;
      while (z < dzlo)
        z += dzs;
    }
  }
  if (x>=xlo && x<xhi && y>=ylo && y<yhi && z>=zlo && z<zhi){
    is = static_cast<int> ((x - xlo)*div[0]/xs);
    js = static_cast<int> ((y - ylo)*div[1]/ys);
    ks = static_cast<int> ((z - zlo)*div[2]/zs);
    ind = 1;
    if (is > div[0]-1) is--;
    if (js > div[1]-1) js--;
    if (ks > div[2]-1) ks--;
  }
  return ind;
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::bc_flux()
{
  double **T = atom->T;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double tnow = update->ntimestep*update->dt;

	if (ind_space == 1){
    double t_vals = vals[7]*exp(-vals[9]*(tnow-stat->st_start*update->dt));
    for (int i = 0; i < stat->nx; i++)
      for (int j = 0; j < stat->ny; j++)
        for (int k = 0; k < stat->nz; k++){
          double temp = vals[2]*stat->aT[1][i][j][k]*cw[0]/(vals[3]+stat->aT[1][i][j][k])*Linj;
          Ri[0][i][j][k] =  temp;
          Ri[1][i][j][k] = -temp;

          temp = vals[4]*stat->aT[7][i][j][k]*cw[0]/(vals[5]+stat->aT[7][i][j][k])*Linj;
          Ri[2][i][j][k] =  temp;
          Ri[3][i][j][k] = -temp;

          temp = vals[0]*stat->aT[13][i][j][k]*cw[2]/(vals[1]+stat->aT[13][i][j][k])*Linj;
          Ri[4][i][j][k] =  temp;
          Ri[5][i][j][k] = -temp;
          Ri[6][i][j][k] = (vals[6]+t_vals*stat->aT[8][i][j][k]+vals[8]*stat->aT[10][i][j][k])*cw[1]*Linj;
        }
	} else if (ind_space == 2){
		double t_vals = vals[7]*exp(-vals[9]*(tnow-nstat[0]*update->dt));
  	for (int i = 0; i < div[0]; i++)
  		for (int j = 0; j < div[1]; j++)
    		for (int k = 0; k < div[2]; k++){
	  			double temp = vals[2]*aTtmp[0][i][j][k]*cw[0]/(vals[3]+aTtmp[0][i][j][k])*Linj;
	  	  	Ri[0][i][j][k] =  temp;
   	 			Ri[1][i][j][k] = -temp;

	    		temp = vals[4]*aTtmp[1][i][j][k]*cw[0]/(vals[5]+aTtmp[1][i][j][k])*Linj;
		  		Ri[2][i][j][k] =  temp;
  		  	Ri[3][i][j][k] = -temp;

	  			temp = vals[0]*aTtmp[4][i][j][k]*cw[2]/(vals[1]+aTtmp[4][i][j][k])*Linj;
	  	  	Ri[4][i][j][k] =  temp;
		    	Ri[5][i][j][k] = -temp;
					Ri[6][i][j][k] = (vals[6]+t_vals*aTtmp[2][i][j][k]+vals[8]*aTtmp[3][i][j][k])*cw[1]*Linj;
				}
	}
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::post_force_vel(int vflag)
{
  int i,m,kk,ll,j,k,l,ix,iy,iz,iz1;
  double xh[3], xh1[3], vh[3],vd[3],ff[3],fs[6],u0,u1,weight,rr,ac,as,theta,theta1,hh,rr1,sigma_t,sigma_n,frand;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int step_t = update->ntimestep;
  double dtv = update->dt;
  double dtinvsqrt = 1.0/sqrt(update->dt);
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;

  if (ind_shear || ind_press)
    for(i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
      l = num_faces_wal[i];
      ll = num_faces_out[i];
      for(m = 0; m < num_faces_flux[i]; m++){
        kk = ind_faces_flux[i][m];
        for(j = 0; j < 3; j++){
          xh[j] = x[i][j] - x0[kk][j];
          vh[j] = v[i][j];
        }
        rot_forward(xh[0],xh[1],xh[2],kk);
        rot_forward(vh[0],vh[1],vh[2],kk);
        switch (ptype[kk]){
          case 1 :
            u1 = xh[1]/aa[kk][2];
            u0 = (xh[0]-u1*aa[kk][1])/aa[kk][0];
            if (u0 >= 0.0 && u1 >= 0.0 && u0+u1 <= 1.0){
              if (ind_shear && (mask[i] & groupbit_s) && ind_coupled[kk] == -1){
                if (xh[2]>0.0 && xh[2]<r_shear && (s_apply == 1 || s_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_per/r_shear);

                  if(iz > n_per -1)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = n_per - 1;
                  }
                  if(iz < 0)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = 0;
                  }

                  velx[kk][iz][0][0] += vh[0];
                  vely[kk][iz][0][0] += vh[1];
                  velz[kk][iz][0][0] += vh[2];
                  num[kk][iz][0][0] += 1.0;
                  weight = 1.0-xh[2]/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][0][0][0]*weight;
                  ff[1] = fsy[kk][0][0][0]*weight;
                  ff[2] = fsz[kk][0][0][0]*weight;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
                if (xh[2]>-r_shear && xh[2]<0.0 && (s_apply == 2 || s_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_per/r_shear);

                  if(iz > n_per -1)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = n_per - 1;
                  }
                  if(iz < 0)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = 0;
                  }

                  velx[kk][iz+n_per][0][0] += vh[0];
                  vely[kk][iz+n_per][0][0] += vh[1];
                  velz[kk][iz+n_per][0][0] += vh[2];
                  num[kk][iz+n_per][0][0] += 1.0;
                  weight = 1.0+xh[2]/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][1][0][0]*weight;
                  ff[1] = fsy[kk][1][0][0]*weight;
                  ff[2] = fsz[kk][1][0][0]*weight;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
              }
              if (ind_out && (mask[i] & groupbit_p)){
                if (xh[2]>0.0 && xh[2]<r_out && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (fabs(xh[2])*n_per/r_out);
                  velx[kk][iz][0][0] += vh[0];
                  vely[kk][iz][0][0] += vh[1];
                  velz[kk][iz][0][0] += vh[2];
                  num[kk][iz][0][0] += 1.0;

                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  ff[2] = 0.0;

                  if (ind_coupled[kk] == 2){
                    weight = 1.0-fabs(xh[2])/r_out;
                    weight = pow(weight,power_out);
                    ff[2] = beta[kk][0][0][0]*weight;
                    iz = static_cast<int> (fabs(xh[2])*n_press/r_press);
//                    weight = f_press[iz]/f_press[0];
                    ff[2] += factor*f_press[iz];

                    //fix for the pressure
                    if(stress_ind == 2)
                    {
                      fs[0] = 0.5*fabs(ff[2])*pressure_grid;
                      fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=0;
                      for(k = 0; k < n_stress_tot; ++k)
                        modify->fix[stress_list[k]]->virial_f(i,fs);
                    }
                  }
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
              }
              if (ind_press && !ind_shear && (mask[i] & groupbit_p) && ind_coupled[kk] == -1){
                if (xh[2]>0.0 && xh[2]<rr_shear && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_shear/rr_shear);
									iz = iz > 10 ? iz : 10;
                  sigma_t = sqrt(-2.0*kbt*f_shear_t[iz]);
                  sigma_n = sqrt(-2.0*kbt*f_shear_n[iz]);
                  for(j = 0; j < 3; j++){
                    vd[j] = vh[j]-vel[kk][0][0][j];
                    if (j < 2)
                      ff[j] = f_shear_t[iz]*vd[j] + sigma_t*ranmars0->gaussian()*dtinvsqrt;
                    else
                      ff[j] = f_shear_n[iz]*vd[j] + sigma_n*ranmars0->gaussian()*dtinvsqrt;
                  }
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
                if (xh[2]>-rr_shear && xh[2]<0.0 && (p_apply == 2 || p_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_shear/rr_shear);
									iz = iz > 10 ? iz : 10;
                  sigma_t = sqrt(-2.0*kbt*f_shear_t[iz]);
                  sigma_n = sqrt(-2.0*kbt*f_shear_n[iz]);
                  for(j = 0; j < 3; j++){
                    vd[j] = vh[j]-vel[kk][0][0][j];
                    if (j < 2)
                      ff[j] = f_shear_t[iz]*vd[j] + sigma_t*ranmars0->gaussian()*dtinvsqrt;
                    else
                      ff[j] = f_shear_n[iz]*vd[j] + sigma_n*ranmars0->gaussian()*dtinvsqrt;
                  }
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
              }
              if (ind_press && (mask[i] & groupbit_p)){
                if (xh[2]>0.0 && xh[2]<r_press && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_press/r_press);
                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  if(ind_coupled[kk] == 0 || ind_coupled[kk] == 2){
                    if(ind_out)
                      ff[2] = f_press[iz] * (dens[kk][0][0][0]/num_dens);
                    else
                      ff[2] = 0.0;
                  }else
                      ff[2] = f_press[iz];
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];

                  //fix for the pressure
                  if(stress_ind == 2)
                  {
                    fs[0] = 0.5*f_press[iz]*pressure_grid;
                    fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=0;
                    for(k = 0; k < n_stress_tot; ++k)
                      modify->fix[stress_list[k]]->virial_f(i,fs);
                  }
                }
                if (xh[2]>-r_press && xh[2]<0.0 && (p_apply == 2 || p_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_press/r_press);
                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  if(ind_coupled[kk] == 0 || ind_coupled[kk] == 2){
                    if(ind_out)
                      ff[2] = -f_press[iz] * (dens[kk][0][0][0]/num_dens);
                    else
                      ff[2] = 0.0;
                  }else
                      ff[2] = -f_press[iz];
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
              }
            }
            break;
          case 2 :
            u1 = xh[1]/aa[kk][2];
            u0 = (xh[0]-u1*aa[kk][1])/aa[kk][0];
            if (u0 >= 0.0 && u1 >= 0.0 && u0 < 1.0 && u1 < 1.0){
              ix = static_cast<int> (u0*ndiv[kk][0]);
              iy = static_cast<int> (u1*ndiv[kk][1]);
              if (ind_shear && (mask[i] & groupbit_s) && ind_coupled[kk] == -1){
                if (xh[2]>0.0 && xh[2]<r_shear && (s_apply == 1 || s_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_per/r_shear);
                  velx[kk][iz][ix][iy] += vh[0];
                  vely[kk][iz][ix][iy] += vh[1];
                  velz[kk][iz][ix][iy] += vh[2];
                  num[kk][iz][ix][iy] += 1.0;
                  weight = 1.0-xh[2]/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][0][ix][iy]*weight;
                  ff[1] = fsy[kk][0][ix][iy]*weight;
                  ff[2] = fsz[kk][0][ix][iy]*weight;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
                if (xh[2]>-r_shear && xh[2]<0.0 && (s_apply == 2 || s_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_per/r_shear);
                  velx[kk][iz+n_per][ix][iy] += vh[0];
                  vely[kk][iz+n_per][ix][iy] += vh[1];
                  velz[kk][iz+n_per][ix][iy] += vh[2];
                  num[kk][iz+n_per][ix][iy] += 1.0;
                  weight = 1.0+xh[2]/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][1][ix][iy]*weight;
                  ff[1] = fsy[kk][1][ix][iy]*weight;
                  ff[2] = fsz[kk][1][ix][iy]*weight;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
              }
              if (ind_out && (mask[i] & groupbit_p)){
                if (xh[2]>0.0 && xh[2]<r_out && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (fabs(xh[2])*n_per/r_out);
                  velx[kk][iz][ix][iy] += vh[0];
                  vely[kk][iz][ix][iy] += vh[1];
                  velz[kk][iz][ix][iy] += vh[2];
                  num[kk][iz][ix][iy] += 1.0;

                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  ff[2] = 0.0;

                  if (ind_coupled[kk] == 2){
                    weight = 1.0-fabs(xh[2])/r_out;
                    weight = pow(weight,power_out);
                    ff[2] = beta[kk][0][ix][iy]*weight;
                    iz = static_cast<int> (fabs(xh[2])*n_press/r_press);
//                    weight = f_press[iz]/f_press[0];
                    ff[2] += factor*f_press[iz];

                    //fix for the pressure
                    if(stress_ind == 2)
                    {
                      fs[0] = 0.5*fabs(ff[2])*pressure_grid;
                      fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=0;
                      for(k = 0; k < n_stress_tot; ++k)
                        modify->fix[stress_list[k]]->virial_f(i,fs);
                    }
                  }

                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
              }
              if (ind_press && !ind_shear && (mask[i] & groupbit_p) && ind_coupled[kk] == -1){
                if (xh[2]>0.0 && xh[2]<rr_shear && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_shear/rr_shear);
									iz = iz > 10 ? iz : 10;
                  sigma_t = sqrt(-2.0*kbt*f_shear_t[iz]);
                  sigma_n = sqrt(-2.0*kbt*f_shear_n[iz]);
                  for(j = 0; j < 3; j++){
                    vd[j] = vh[j]-vel[kk][ix][iy][j];
                    if (j < 2)
                      ff[j] = f_shear_t[iz]*vd[j] + sigma_t*ranmars0->gaussian()*dtinvsqrt;
                    else
                      ff[j] = f_shear_n[iz]*vd[j] + sigma_n*ranmars0->gaussian()*dtinvsqrt;
                  }
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
                if (xh[2]>-rr_shear && xh[2]<0.0 && (p_apply == 2 || p_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_shear/rr_shear);
									iz = iz > 10 ? iz : 10;
                  sigma_t = sqrt(-2.0*kbt*f_shear_t[iz]);
                  sigma_n = sqrt(-2.0*kbt*f_shear_n[iz]);
                  for(j = 0; j < 3; j++){
                    vd[j] = vh[j]-vel[kk][ix][iy][j];
                    if (j < 2)
                      ff[j] = f_shear_t[iz]*vd[j] + sigma_t*ranmars0->gaussian()*dtinvsqrt;
                    else
                      ff[j] = f_shear_n[iz]*vd[j] + sigma_n*ranmars0->gaussian()*dtinvsqrt;
                  }
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
              }
              if (ind_press && (mask[i] & groupbit_p)){
                if (xh[2]>0.0 && xh[2]<r_press && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (xh[2]*n_press/r_press);
                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  if(ind_coupled[kk] == 0 || ind_coupled[kk] == 2){
                    if(ind_out)
                      ff[2] = f_press[iz] * (dens[kk][0][ix][iy]/num_dens);
                    else
                      ff[2] = 0.0;
                  }else
                      ff[2] = f_press[iz];
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];

                  //fix for the pressure
                  if(stress_ind == 2)
                  {
                    fs[0] = 0.5*f_press[iz]*pressure_grid;
                    fs[1]=fs[2]=fs[3]=fs[4]=fs[5]=0;
                    for(k = 0; k < n_stress_tot; ++k)
                      modify->fix[stress_list[k]]->virial_f(i,fs);
                  }
                }
                if (xh[2]>-r_press && xh[2]<0.0 && (p_apply == 2 || p_apply == 0)){
                  iz = static_cast<int> (-xh[2]*n_press/r_press);
                  ff[0] = 0.0;
                  ff[1] = 0.0;
                  if(ind_coupled[kk] == 0 || ind_coupled[kk] == 2){
                    if(ind_out)
                      ff[2] = -f_press[iz] * (dens[kk][0][ix][iy]/num_dens);
                    else
                      ff[2] = 0.0;
                  }else
                      ff[2] = -f_press[iz];
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
              }
            }
            break;
          case 3 :
            if (xh[2]>=0.0 && xh[2] < aa[kk][0]){
              rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1]);
              ix = static_cast<int> (xh[2]*ndiv[kk][0]/aa[kk][0]);

              if(ix > ndiv[kk][0] -1)
              {
                std::cout<<"strange ix "<<ix<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                ix = ndiv[kk][0] - 1;
              }
              if(ix < 0)
              {
                std::cout<<"strange ix "<<ix<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                ix = 0;
              }

              theta = acos(xh[0]/rr);
              theta1 = asin(xh[1]/rr);
              if (theta1 < 0.0)
                theta = 2.0*M_PI - theta;
              iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);

              if(iy > ndiv[kk][1] -1)
              {
                std::cout<<"strange iy "<<iy<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                iy = ndiv[kk][1] - 1;
              }
              if(iy < 0)
              {
                std::cout<<"strange iy "<<iy<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                iy = 0;
              }

              hh = rr - aa[kk][1];
              if (ind_shear && (mask[i] & groupbit_s)){
                if (hh>0.0 && hh<r_shear && (s_apply == 1 || s_apply == 0)){
                  iz = static_cast<int> (hh*n_per/r_shear);

                  if(iz > n_per -1)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = n_per - 1;
                  }
                  if(iz < 0)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = 0;
                  }

                  velx[kk][iz][ix][iy] += vh[0];
                  vely[kk][iz][ix][iy] += vh[1];
                  velz[kk][iz][ix][iy] += vh[2];
                  num[kk][iz][ix][iy] += 1.0;
                  weight = 1.0-hh/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][0][ix][iy]*weight;
                  ff[1] = fsy[kk][0][ix][iy]*weight;
                  ff[2] = fsz[kk][0][ix][iy]*weight;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
                }
                if (hh>-r_shear && hh<0.0 && (s_apply == 2 || s_apply == 0)){
                  iz = static_cast<int> (-hh*n_per/r_shear);

                  if(iz > n_per -1)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = n_per - 1;
                  }
                  if(iz < 0)
                  {
                    std::cout<<"strange iz "<<iz<<" cpu is "<<comm->me<<" time step "<<update->ntimestep<<"\n";
                    iz = 0;
                  }

                  velx[kk][iz+n_per][ix][iy] += vh[0];
                  vely[kk][iz+n_per][ix][iy] += vh[1];
                  velz[kk][iz+n_per][ix][iy] += vh[2];
                  num[kk][iz+n_per][ix][iy] += 1.0;
                  weight = 1.0+hh/r_shear;
                  weight = pow(weight,power);
                  ff[0] = fsx[kk][1][ix][iy]*weight;
                  ff[1] = fsy[kk][1][ix][iy]*weight;
                  ff[2] = fsz[kk][1][ix][iy]*weight;

                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] -= ff[j];
              }
            }
              if (ind_press && (mask[i] & groupbit_p)){
                if (hh>0.0 && hh<r_press && (p_apply == 1 || p_apply == 0)){
                  iz = static_cast<int> (hh*n_press/r_press);
                  ff[0] = f_press[iz]*xh[0]/rr;
                  ff[1] = f_press[iz]*xh[1]/rr;
                  ff[2] = 0.0;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
                if (hh>-r_press && hh<0.0 && (p_apply == 2 || p_apply == 0)){
                  iz = static_cast<int> (-hh*n_press/r_press);
                  ff[0] = -f_press[iz]*xh[0]/rr;
                  ff[1] = -f_press[iz]*xh[1]/rr;
                  ff[2] = 0.0;
                  rot_back(ff[0],ff[1],ff[2],kk);
                  for(j = 0; j < 3; j++)
                    f[i][j] += ff[j];
                }
              }
            }
            break;
          case 4 :
            rr = sqrt(xh[0]*xh[0] + xh[1]*xh[1] + xh[2]*xh[2]);
            rr1 = sqrt(xh[0]*xh[0] + xh[1]*xh[1]);
            theta = acos(xh[2]/rr);
            ix = static_cast<int> (theta*ndiv[kk][0]/M_PI);
            theta = acos(xh[0]/rr1);
            theta1 = asin(xh[1]/rr1);
            if (theta1 < 0.0)
              theta = 2.0*M_PI - theta;
            iy = static_cast<int> (0.5*theta*ndiv[kk][1]/M_PI);
            hh = rr - aa[kk][0];
            if (ind_shear && (mask[i] & groupbit_s)){
              if (hh>0.0 && hh<r_shear && (s_apply == 1 || s_apply == 0)){
                iz = static_cast<int> (hh*n_per/r_shear);
                velx[kk][iz][ix][iy] += vh[0];
                vely[kk][iz][ix][iy] += vh[1];
                velz[kk][iz][ix][iy] += vh[2];
                num[kk][iz][ix][iy] += 1.0;
                weight = 1.0-hh/r_shear;
                weight = pow(weight,power);
                ff[0] = fsx[kk][0][ix][iy]*weight;
                ff[1] = fsy[kk][0][ix][iy]*weight;
                ff[2] = fsz[kk][0][ix][iy]*weight;
                rot_back(ff[0],ff[1],ff[2],kk);
                for(j = 0; j < 3; j++)
                  f[i][j] -= ff[j];
              }
              if (hh>-r_shear && hh<0.0 && (s_apply == 2 || s_apply == 0)){
                iz = static_cast<int> (-hh*n_per/r_shear);
                velx[kk][iz+n_per][ix][iy] += vh[0];
                vely[kk][iz+n_per][ix][iy] += vh[1];
                velz[kk][iz+n_per][ix][iy] += vh[2];
                num[kk][iz+n_per][ix][iy] += 1.0;
                weight = 1.0+hh/r_shear;
                weight = pow(weight,power);
                ff[0] = fsx[kk][1][ix][iy]*weight;
                ff[1] = fsy[kk][1][ix][iy]*weight;
                ff[2] = fsz[kk][1][ix][iy]*weight;
                rot_back(ff[0],ff[1],ff[2],kk);
                for(j = 0; j < 3; j++)
                  f[i][j] -= ff[j];
              }
            }
            if (ind_press && (mask[i] & groupbit_p)){
              if (hh>0.0 && hh<r_press && (p_apply == 1 || p_apply == 0)){
                iz = static_cast<int> (hh*n_press/r_press);
                ff[0] = f_press[iz]*xh[0]/rr;
                ff[1] = f_press[iz]*xh[1]/rr;
                ff[2] = f_press[iz]*xh[2]/rr;
                rot_back(ff[0],ff[1],ff[2],kk);
                for(j = 0; j < 3; j++)
                  f[i][j] += ff[j];
              }
              if (hh>-r_press && hh<0.0 && (p_apply == 2 || p_apply == 0)){
                iz = static_cast<int> (-hh*n_press/r_press);
                ff[0] = -f_press[iz]*xh[0]/rr;
                ff[1] = -f_press[iz]*xh[1]/rr;
                ff[2] = -f_press[iz]*xh[2]/rr;
                rot_back(ff[0],ff[1],ff[2],kk);
                for(j = 0; j < 3; j++)
                  f[i][j] += ff[j];
              }
            }
            break;
        }
        if (ind_coupled[kk] == -1)
          if (l > 1)
            l--;
      }
    		}
		}

  if (ind_shear)
    if(step_t%iter == 0)
      recalc_force();

  if (ind_out){
    if(step_t%iter == 0 && step_t != 0)
      recalc_force();
  }
}

/* ---------------------------------------------------------------------- */

void FixSolidBound::post_force_cc(int vflag)
{
  int i,m,kk,j,k,ix,iy,iz,iz1;
  double xh[3], xh1[3], vh[3],vd[3],u0,u1;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int step_t = update->ntimestep;
	int reac = 0;

  double **T = atom->T;
  double **Q = atom->Q;

  if (ind_space == 1 && step_t>stat->st_start && step_t%stat->dump_each == 0)
		bc_flux();
  if (ind_space == 2 && step_t>nstat[0] && step_t%nstat[2] == 0)
    bc_flux();

  if (ind_cc)
    for(i = 0; i < nlocal; i++){
				if (mask[i] & groupbit){
      for(m = 0; m < num_faces_flux[i]; m++){
        kk = ind_faces_flux[i][m];
        for(j = 0; j < 3; j++){
          xh[j] = x[i][j] - x0[kk][j];
          vh[j] = v[i][j];
        }
        rot_forward(xh[0],xh[1],xh[2],kk);
        rot_forward(vh[0],vh[1],vh[2],kk);
        switch (ptype[kk]){
          case 1 :
            u1 = xh[1]/aa[kk][2];
            u0 = (xh[0]-u1*aa[kk][1])/aa[kk][0];
            if (u0 >= 0.0 && u1 >= 0.0 && u0+u1 <= 1.0){
              if (n_ccD && (mask[i] & groupbit_c) && (ind_bc_cc[kk] == 1))
                if (fabs(xh[2]) < r_ccD){
                  iz = static_cast<int> (xh[2]*n_ccD/r_ccD);
                  for (j = 0; j < NPDE; j++)
                    Q[i][j] += kappa[j] * (bc_cc[kk] - T[i][j]) * concen_D[iz];
                }
              if (n_ccN && (mask[i] & groupbit_c) && (ind_bc_cc[kk] == 2))
                if (fabs(xh[2]) < r_ccN){
                  iz = static_cast<int> (xh[2]*n_ccN/r_ccN);
                  for (j = 0; j < NPDE; j++)
                    Q[i][j] += diff[j] * bc_cc[kk] * concen_N[iz];
									if (ind_space == 1 && (stat->map_index(x[i][0],x[i][1],x[i][2])))
										{is = stat->is; js = stat->js; ks = stat->ks; reac = 1;}
                  if (ind_space == 2 && (map_index(x[i][0],x[i][1],x[i][2]))) reac = 1;
									if (reac == 1){   //Alireza: hard-coded for 25-species model in a channel
										for (j = 0; j < n_sur_reac; j++)
                    	Q[i][sur_reac[j]] += Ri[j][is][js][ks] * concen_N[iz];
										reac = 0;
									}
//	fprintf(stdout, "%e %e %e %e\n", stat->aT[1][is][js][ks], stat->aT[7][is][js][ks], stat->aT[13][is][js][ks], Ri[0][is][js][ks]);
                }
            }
            break;
          case 2 :
            u1 = xh[1]/aa[kk][2];
            u0 = (xh[0]-u1*aa[kk][1])/aa[kk][0];
            if (u0 >= 0.0 && u1 >= 0.0 && u0 < 1.0 && u1 < 1.0){
              ix = static_cast<int> (u0*ndiv[kk][0]);
              iy = static_cast<int> (u1*ndiv[kk][1]);
              if (n_ccD && (mask[i] & groupbit_c) && (ind_bc_cc[kk] == 1))
                if (fabs(xh[2]) < r_ccD){
                  iz = static_cast<int> (fabs(xh[2])*n_ccD/r_ccD);
                  for (j = 0; j < NPDE; j++)
                    Q[i][j] += kappa[j] * (bc_cc[kk] - T[i][j]) * concen_D[iz];
                }
              if (n_ccN && (mask[i] & groupbit_c) && (ind_bc_cc[kk] == 2))
                if (fabs(xh[2]) < r_ccN){
                  iz = static_cast<int> (fabs(xh[2])*n_ccN/r_ccN);
                  for (j = 0; j < NPDE; j++)
                    Q[i][j] += diff[j] * bc_cc[kk] * concen_N[iz];
                  if (ind_space == 1 && (stat->map_index(x[i][0],x[i][1],x[i][2])))
                    {is = stat->is; js = stat->js; ks = stat->ks; reac = 1;}
                  if (ind_space == 2 && (map_index(x[i][0],x[i][1],x[i][2]))) reac = 1;
                  if (reac == 1){   //Alireza: hard-coded for 25-species model in a channel
                    for (j = 0; j < n_sur_reac; j++)
                      Q[i][sur_reac[j]] += Ri[j][is][js][ks] * concen_N[iz];
                    reac = 0;
                  }
                }
            }
            break;
          case 3 :
            error->one("ADR boundary conditions are not set up for cylinderical wall boundaries...\n");
            break;
          case 4 :
            error->one("ADR boundary conditions are not set up for spherical wall boundaries...\n");
            break;
        }
      }
    		}
		}
}
