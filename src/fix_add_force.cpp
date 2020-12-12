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
/* ----------------------------------------------------------------------
   changed by Huan Lei
   June 17 2009
   ----------------------------------------------------------------------*/
#include "mpi.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_add_force.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "comm.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixAddForce::FixAddForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  char buf[BUFSIZ];
  char fname[FILENAME_MAX];
  FILE *f_read;
  int i, j, k, m;

  if (narg < 7) error->all("Illegal fix addforce command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  scalar_vector_freq = 1;
  extscalar = 1;
  extvector = 1;

  index = atoi(arg[3]);

	switch (index){
		case 0: case 1: case 2: case 3: case 4:
  		xvalue = atof(arg[4]);
  		yvalue = atof(arg[5]);
  		zvalue = atof(arg[6]);
  		dir = static_cast<int> (xvalue);
  		force = yvalue;
  		freq = zvalue; 
		break;
		case 5 :		// Alireza: pressure gradient as the body force
      xvalue = atof(arg[4]);
      yvalue = atof(arg[5]);
      zvalue = atof(arg[6]);
      ff = atof(arg[7]);
    	sprintf(fname,arg[8]);
  		if (comm->me == 0){
    		f_read = fopen(fname,"r");
    		if(f_read == (FILE*) NULL)
      		error->one("Could not open input forcing file");
    		fgets(buf,BUFSIZ,f_read);
    		sscanf(buf,"%d %d %lf %lf %lf %lf",&dimx,&dimy,&Xmin,&Xmax,&Ymin,&Ymax);
  		}
  		MPI_Bcast(&dimx,1,MPI_INT,0,world);
  		MPI_Bcast(&dimy,1,MPI_INT,0,world);
  		MPI_Bcast(&Xmin,1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&Xmax,1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&Ymin,1,MPI_DOUBLE,0,world);
  		MPI_Bcast(&Ymax,1,MPI_DOUBLE,0,world);

  		meshX = memory->create_2d_double_array(dimy,dimx,"fix_add_force:meshX");
  		meshY = memory->create_2d_double_array(dimy,dimx,"fix_add_force:meshY");
  		gradxP = memory->create_2d_double_array(dimy,dimx,"fix_add_force:gradxP");
  		gradyP = memory->create_2d_double_array(dimy,dimx,"fix_add_force:gradyP");

  		if (comm->me == 0){
    		for (i = 0; i < dimy; i++){
      		for (j = 0; j < dimx; j++){
        		fgets(buf,BUFSIZ,f_read);
        		sscanf(buf,"%lf %lf %lf %lf",&meshX[i][j],&meshY[i][j],&gradxP[i][j],&gradyP[i][j]);
      		}
    		}
				
  			fclose(f_read);
  		}
  		MPI_Bcast(&meshX[0][0],dimx*dimy,MPI_DOUBLE,0,world);
  		MPI_Bcast(&meshY[0][0],dimx*dimy,MPI_DOUBLE,0,world);
  		MPI_Bcast(&gradxP[0][0],dimx*dimy,MPI_DOUBLE,0,world);
  		MPI_Bcast(&gradyP[0][0],dimx*dimy,MPI_DOUBLE,0,world);

  		if (comm->me == 0) printf("READ PRESSURE GRADIENT DATA FILE................... DONE\n");
		break;
		case 6 :		//Alireza: added for stenosis
      xvalue = atof(arg[4]);
      yvalue = atof(arg[5]);
      zvalue = atof(arg[6]);
      ff = atof(arg[7]);
			R0 = atof(arg[8]);
			alpha = atof(arg[9]);
			len	= atof(arg[10]);
			X0 = atof(arg[11]);
		break;
	}
  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixAddForce::~FixAddForce()
{
	if (index == 5){
	  memory->destroy_2d_double_array(gradxP);
	  memory->destroy_2d_double_array(gradyP);
	  memory->sfree(meshX);
	  memory->sfree(meshY);
	}
}

/* ---------------------------------------------------------------------- */

int FixAddForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddForce::init()
{
  
  dtv = update->dt;
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAddForce::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force(int vflag)
{
  int i, j;
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double xh[3],cm[4],cmt[4];
  double time;
  double z_half = 0.5*domain->zprd;
  double z_min = domain->boxlo[2];
	double theta, R, u;

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  if (index == 1 || index == 3){
    for(i=0; i<4; i++){
      cm[i] = 0.0;
      cmt[i] = 0.0;
    }
    for (i=0; i<nlocal; i++)
      if (mask[i] & groupbit) {
        for (j=0; j<3; j++) xh[j] = x[i][j];
        domain->unmap(xh,atom->image[i]);
        for (j=0; j<3; j++) cmt[j] += xh[j];
        cmt[3] += 1.0;
      }
    MPI_Allreduce(&cmt[0],&cm[0],4,MPI_DOUBLE,MPI_SUM,world);
    for(i=0; i<3; i++)
      cm[i] /= cm[3];
  }
  if (index == 2 || index == 3){
    time = dtv*update->ntimestep;
    ff = force; //*cos(2.0*M_PI*time*freq); //commented by LG to get steady flow
  }
  // foriginal[0] = - x dot f = "potential" for added force
  // foriginal[123] = force on atoms before extra force added
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit){
      foriginal[1] += f[i][0];
      foriginal[2] += f[i][1];
      foriginal[3] += f[i][2];
      switch (index){
        case 0 :
          foriginal[0] -= xvalue*x[i][0] + yvalue*x[i][1] + zvalue*x[i][2];
          f[i][0] += xvalue;
          f[i][1] += yvalue;
          f[i][2] += zvalue;
	  		break;
        case 1 :
          for (j=0; j<3; j++) xh[j] = x[i][j];
          domain->unmap(xh,atom->image[i]);
          for(j=0; j<3; j++) cmt[j] = xh[j] - cm[j];
          if (dir == 1){
            foriginal[0] -= force*(cmt[2]*x[i][1] - cmt[1]*x[i][2]);
            f[i][1] += cmt[2]*force;
            f[i][2] -= cmt[1]*force;
          }
          else if (dir == 2){
            foriginal[0] -= force*(-cmt[2]*x[i][0] + cmt[0]*x[i][2]);
            f[i][0] -= cmt[2]*force;
            f[i][2] += cmt[0]*force;
          }
          else{
            foriginal[0] -= force*(cmt[1]*x[i][0] - cmt[0]*x[i][1]);
            f[i][0] += cmt[1]*force;
            f[i][1] -= cmt[0]*force;
          }
	  		break;
        case 2 :
          if (dir == 1){
            foriginal[0] -= ff*x[i][0];
            f[i][0] += ff;
          }
          else if (dir == 2){
            foriginal[0] -= ff*x[i][1];
            f[i][1] += ff;
          }
          else{
            foriginal[0] -= ff*x[i][2];
            f[i][2] += ff;
          }
	  		break;
        case 3 :
          for (j=0; j<3; j++) xh[j] = x[i][j];
          domain->unmap(xh,atom->image[i]);
          for(j=0; j<3; j++) cmt[j] = xh[j] - cm[j];
          if (dir == 1){
            foriginal[0] -= ff*(cmt[2]*x[i][1] - cmt[1]*x[i][2]);
            f[i][1] += cmt[2]*ff;
            f[i][2] -= cmt[1]*ff;
          }
          else if (dir == 2){
            foriginal[0] -= ff*(-cmt[2]*x[i][0] + cmt[0]*x[i][2]);
            f[i][0] -= cmt[2]*ff;
            f[i][2] += cmt[0]*ff;
          }
          else{
            foriginal[0] -= ff*(cmt[1]*x[i][0] - cmt[0]*x[i][1]);
            f[i][0] += cmt[1]*ff;
            f[i][1] -= cmt[0]*ff;
          }
	  		break;
        case 4 :
        	j = static_cast<int>((x[i][2]-z_min)/z_half + 2.0);
          if (j%2 == 0) {
            foriginal[0] -= xvalue*x[i][0];
            f[i][0] += xvalue;
	  			}
          else {
            foriginal[0] += xvalue*x[i][0];
            f[i][0] -= xvalue;
          }
	  		break;
        case 5 :
          for (j=0; j<3; j++) xh[j] = x[i][j];
          interp_force(xh);

          foriginal[0] -= ff*xvalue*x[i][0] + ff*yvalue*x[i][1] + zvalue*x[i][2];
          f[i][0] += ff*xvalue;
          f[i][1] += ff*yvalue;
          f[i][2] += zvalue;
        break;
				case 6 :		// Alireza: flow along x-axis
					if (x[i][0] <= X0 || x[i][0] >= (X0+2.0*len*R0)){
	          foriginal[0] -= xvalue*x[i][0] + yvalue*x[i][1] + zvalue*x[i][2];
  	        f[i][0] += xvalue;
    	      f[i][1] += yvalue;
      	    f[i][2] += zvalue;
					}else{
      			//u = (x[i][0] - X0 - len*R0) / (len*R0);
						//theta = M_PI*u;
      			//R = R0 * (1.0 - alpha*(1.0+cos(theta))/2.0);
						//ff = xvalue*pow(R0/R,4.0);

            foriginal[0] -= ff*xvalue*x[i][0] + yvalue*x[i][1] + zvalue*x[i][2];
            f[i][0] += ff*xvalue;
            f[i][1] += yvalue;
            f[i][2] += zvalue;

					}
				break;
      }
    }
    
}

/* ---------------------------------------------------------------------- */

void FixAddForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixAddForce::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixAddForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

#if 0
/* ----------------------------------------------------------------------
   return components of added force interpolated at the locations of atoms
------------------------------------------------------------------------- */

void FixAddForce::interp_force(double *xh)
{
  int i, j, cnt;
  int ix, iy, iz = 0;
  double fx, fy, fz;

  ix = static_cast<int> (fabs(xh[0]-Xmin)*dimx/fabs(Xmax-Xmin));
  iy = static_cast<int> (fabs(xh[1]-Ymin)*dimy/fabs(Ymax-Ymin));

  if (ix == dimx-1) ix--;
  if (iy == dimy-1) iy--;

  fx = fy = fz = 0.0;
  cnt = 0;
  if (gradxP[iy][ix] != 999.0){
    fx += gradxP[iy][ix];
    fy += gradyP[iy][ix];
    cnt++;
  }
  if (gradxP[iy+1][ix] != 999.0){
    fx += gradxP[iy+1][ix];
    fy += gradyP[iy+1][ix];
    cnt++;
		}
  }
  if (gradxP[iy][ix+1] != 999.0){
    fx += gradxP[iy][ix+1];
    fy += gradyP[iy][ix+1];
    cnt++;
  }
  if (gradxP[iy+1][ix+1] != 999.0){
    fx += gradxP[iy+1][ix+1];
    fy += gradyP[iy+1][ix+1];
    cnt++;
  }
  if (cnt != 0){
   fx /= static_cast<double> (cnt);
   fy /= static_cast<double> (cnt);
  }
  xvalue = fx; yvalue = fy; zvalue = fz;
}

#else
void FixAddForce::interp_force(double *xh)
{
  int i, j, cnt;
  int ix, iy, iz = 0;
  double fx, fy, fz;
  double min_dist = 1000.0, dist;

  ix = static_cast<int> (fabs(xh[0]-Xmin)*dimx/(Xmax-Xmin));
  if (ix >= dimx-1) ix = dimx-2;
  iy = static_cast<int> (fabs(xh[1]-meshY[0][ix])*dimy/(meshY[dimy-1][ix]-meshY[0][ix]));
  if (iy >= dimy-1) iy = dimy-2;

  fx = fy = fz = 0.0;
  if (gradxP[iy][ix] != 999.0){
    dist = sqrt((xh[0]-meshX[iy][ix])*(xh[0]-meshX[iy][ix])+(xh[1]-meshY[iy][ix])*(xh[1]-meshY[iy][ix]));
    if (dist < min_dist){
    	min_dist = dist;
    	fx = gradxP[iy][ix];
    	fy = gradyP[iy][ix];
    }
  }
  if (gradxP[iy+1][ix] != 999.0){
    dist = sqrt((xh[0]-meshX[iy+1][ix])*(xh[0]-meshX[iy+1][ix])+(xh[1]-meshY[iy+1][ix])*(xh[1]-meshY[iy+1][ix]));
    if (dist < min_dist){
      min_dist = dist;
    	fx = gradxP[iy+1][ix];
    	fy = gradyP[iy+1][ix];
    }
  }
  if (gradxP[iy][ix+1] != 999.0){
    dist = sqrt((xh[0]-meshX[iy][ix+1])*(xh[0]-meshX[iy][ix+1])+(xh[1]-meshY[iy][ix+1])*(xh[1]-meshY[iy][ix+1]));
    if (dist < min_dist){
      min_dist = dist;
	    fx = gradxP[iy][ix+1];
  	  fy = gradyP[iy][ix+1];
		}
  }
  if (gradxP[iy+1][ix+1] != 999.0){
    dist = sqrt((xh[0]-meshX[iy+1][ix+1])*(xh[0]-meshX[iy+1][ix+1])+(xh[1]-meshY[iy+1][ix+1])*(xh[1]-meshY[iy+1][ix+1]));
    if (dist < min_dist){
      min_dist = dist;
	    fx = gradxP[iy+1][ix+1];
  	  fy = gradyP[iy+1][ix+1];
		}
  }
  xvalue = fx; yvalue = fy; zvalue = fz;
}
#endif
