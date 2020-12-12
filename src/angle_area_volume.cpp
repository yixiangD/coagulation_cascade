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

#include <algorithm>
#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "angle_area_volume.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include "update.h"

using namespace LAMMPS_NS;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleAreaVolume::AngleAreaVolume(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleAreaVolume::~AngleAreaVolume()
{
  if (allocated) {
    memory->sfree(setflag);
    memory->sfree(ka);
    memory->sfree(a0);
    memory->sfree(kv);
    memory->sfree(v0);
    memory->sfree(kl);
    memory->sfree(t_remod);   
  }
  if (init_on){
    memory->sfree(h_area);
    memory->sfree(h_volume);
    memory->sfree(ttyp);
  }

}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::compute(int eflag, int vflag)
{
  double eangle,f1[3],f3[3];

  int i1,i2,i3,n,m,j,type;
  double d21x,d21y,d21z,d31x,d31y,d31z,d32x,d32y,d32z;
  double nx,ny,nz,nn,mx,my,mz,aa,vv,ar0;
  double coefl, coefa, coefv, coefca;
  double s1x,s1y,s1z,s2x,s2y,s2z,s3x,s3y,s3z;
  double s1xv,s1yv,s1zv,s2xv,s2yv,s2zv,s3xv,s3yv,s3zv;
  int nm = atom->n_mol;
  int ttyp1[nm];
  double dath[2*nm],datt[2*nm],xx1[3],xx2[3],xx3[3]; 
  double ff[6];

  double delx1,dely1,delz1,delx2,dely2,delz2;
  double dtheta,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
  
  int stress_ind = force->stress_ind;
  int n_stress_tot = force->n_stress_tot;
  int *stress_list = force->stress_list;
  
  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;


  if (init_on == 0){
    h_area = (double *) memory->smalloc(nm*sizeof(double),"angle_area_volume:h_area");
    h_volume = (double *) memory->smalloc(nm*sizeof(double),"angle_area_volume:h_volume");
    ttyp = (int *) memory->smalloc(nm*sizeof(int),"angle_area_volume:ttyp");
    for (n = 0; n < nm; n++){
      ttyp[n] = 0;
      ttyp1[n] = 0;
    }
    for (n = 0; n < neighbor->nanglelist; n++) {
      i1 = neighbor->anglelist[n][0];
      m = atom->molecule[i1]-1;
      ttyp1[m] = neighbor->anglelist[n][3]; 
    }
    init_on = 1;
    MPI_Allreduce(&ttyp1,ttyp,nm,MPI_INT,MPI_MAX,world);
  }
  for (n = 0; n < nm; n++) {
    h_area[n] = 0.0;
    h_volume[n] = 0.0;
  }

  energy = 0.0;
  energy_a = 0.0;
  energy_v = 0.0;
  energy_l = 0.0;
    
  for (n = 0; n < 2*nm; n++){
    dath[n] = 0.0;
    datt[n] = 0.0;   
  }

  int *n_atoms = atom->atoms_in_mol;
  double *anglelist_area = neighbor->anglelist_area;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    m = atom->molecule[i1]-1;

		if(type == 1 || type == 2 || type == 3){
    	// 2-1 distance
    	d21x = x[i2][0] - x[i1][0];
    	d21y = x[i2][1] - x[i1][1];
    	d21z = x[i2][2] - x[i1][2];
    	domain->minimum_image(d21x,d21y,d21z);

    	// 3-1 distance
    	d31x = x[i3][0] - x[i1][0];
    	d31y = x[i3][1] - x[i1][1];
    	d31z = x[i3][2] - x[i1][2];
    	domain->minimum_image(d31x,d31y,d31z);

    	// calculate normal
    	nx = d21y*d31z - d31y*d21z;
    	ny = d31x*d21z - d21x*d31z;
    	nz = d21x*d31y - d31x*d21y;
    	nn = sqrt(nx*nx + ny*ny + nz*nz);
    
    	// calculate center
    	for (j = 0; j < 3; j++){
      	xx1[j] = x[i1][j];
      	xx2[j] = x[i2][j];
      	xx3[j] = x[i3][j]; 
    	}
    
    	domain->unmap(xx1,atom->image[i1]);
    	domain->unmap(xx2,atom->image[i2]);
    	domain->unmap(xx3,atom->image[i3]);
    
    	mx =  xx1[0] + xx2[0] + xx3[0];
    	my =  xx1[1] + xx2[1] + xx3[1];
    	mz =  xx1[2] + xx2[2] + xx3[2];
    
    	// calculate area and volume
    	aa = 0.5*nn;
    	vv = (nx*mx + ny*my + nz*mz)/18.0;
    	dath[m] += aa;
    	dath[m+nm] += vv;
    	h_area[m] += aa;
    	h_volume[m] += vv;
  	}
	}

  MPI_Allreduce(&dath[0],&datt[0],2*nm,MPI_DOUBLE,MPI_SUM,world);
  if (eflag) {
    for (m = 0; m < nm; m++){
      type = ttyp[m];
      energy_a += 0.5*ka[type]*(a0[type]-datt[m])*(a0[type]-datt[m])/a0[type];
      energy_v += 0.5*kv[type]*(v0[type]-datt[m+nm])*(v0[type]-datt[m+nm])/v0[type];
      if(m == 1 && update->ntimestep % 500 == 0)
        if(comm->me == 0) printf("total area is %f, total volume is %f, desired is %f %f, in step %d\n", datt[m], datt[m+nm], a0[type], v0[type], update->ntimestep); 
    }
  } 

  for (n = 0; n < nanglelist; n++) {
    i1 = anglelist[n][0];
    i2 = anglelist[n][1];
    i3 = anglelist[n][2];
    type = anglelist[n][3];
    if(type == 1 || type == 2 || type == 3) {
    	ar0 = anglelist_area[n];
    	m = atom->molecule[i1]-1;

    	// 2-1 distance
    	d21x = x[i2][0] - x[i1][0];
    	d21y = x[i2][1] - x[i1][1];
    	d21z = x[i2][2] - x[i1][2];
    	domain->minimum_image(d21x,d21y,d21z);

    	// 3-1 distance
    	d31x = x[i3][0] - x[i1][0];
    	d31y = x[i3][1] - x[i1][1];
    	d31z = x[i3][2] - x[i1][2];
    	domain->minimum_image(d31x,d31y,d31z);

    	// 3-2 distance
    	d32x = x[i3][0] - x[i2][0];
    	d32y = x[i3][1] - x[i2][1];
    	d32z = x[i3][2] - x[i2][2];
    	domain->minimum_image(d32x,d32y,d32z);
    
    	// calculate normal
    	nx = d21y*d31z - d31y*d21z;
    	ny = d31x*d21z - d21x*d31z;
    	nz = d21x*d31y - d31x*d21y;
    	nn = sqrt(nx*nx + ny*ny + nz*nz);
    
    	// calculate center
    	for (j = 0; j < 3; j++){
      	xx1[j] = x[i1][j];
       	xx2[j] = x[i2][j];
        xx3[j] = x[i3][j]; 
    	}
    	domain->unmap(xx1,atom->image[i1]);
    	domain->unmap(xx2,atom->image[i2]);
    	domain->unmap(xx3,atom->image[i3]);
    	mx =  xx1[0] + xx2[0] + xx3[0];
    	my =  xx1[1] + xx2[1] + xx3[1];
    	mz =  xx1[2] + xx2[2] + xx3[2];
        
    	// calculate coeffs
			double relax = 1.0;
		  if (t_remod[type] > 0.0) relax = cell_remodel(i1, i2, i3, t_remod[type]);
    	coefl = 0.5*kl[type]*(ar0-0.5*nn)/ar0/nn; 
    	coefa = 0.5*ka[type]*(a0[type]-datt[m])/a0[type]/nn;
    	coefca = coefl + coefa;        
    	coefv = relax*kv[type]*(v0[type]-datt[m+nm])/v0[type]/18.0;  

    	if (eflag) {
      	energy_l += 0.5*kl[type]*(ar0-0.5*nn)*(ar0-0.5*nn)/ar0;
        eangle = 0.5*kl[type]*(ar0-0.5*nn)*(ar0-0.5*nn)/ar0;
    	}

    	s1x = coefca*(ny*d32z - nz*d32y);
    	s1y = coefca*(nz*d32x - nx*d32z);    
    	s1z = coefca*(nx*d32y - ny*d32x);
    	s2x = coefca*(nz*d31y - ny*d31z);
    	s2y = coefca*(nx*d31z - nz*d31x);
    	s2z = coefca*(ny*d31x - nx*d31y);
    	s3x = coefca*(ny*d21z - nz*d21y);
    	s3y = coefca*(nz*d21x - nx*d21z);
    	s3z = coefca*(nx*d21y - ny*d21x);
    
    	s1xv = coefv*(nx + d32z*my - d32y*mz);
    	s1yv = coefv*(ny - d32z*mx + d32x*mz);    
    	s1zv = coefv*(nz + d32y*mx - d32x*my);
    	s2xv = coefv*(nx - d31z*my + d31y*mz);
    	s2yv = coefv*(ny + d31z*mx - d31x*mz);
    	s2zv = coefv*(nz - d31y*mx + d31x*my);
    	s3xv = coefv*(nx + d21z*my - d21y*mz);
    	s3yv = coefv*(ny - d21z*mx + d21x*mz);
    	s3zv = coefv*(nz + d21y*mx - d21x*my);
    
    	// apply force to each of 3 atoms
    	if (newton_bond || i1 < nlocal) {
      	    f1[0] = s1x + s1xv;
      	    f1[1] = s1y + s1yv;
            f1[2] = s1z + s1zv;
            f[i1][0] += f1[0];
            f[i1][1] += f1[1];
            f[i1][2] += f1[2];
      }

      if (newton_bond || i2 < nlocal) {
            f[i2][0] += s2x+s2xv;
            f[i2][1] += s2y+s2yv;
            f[i2][2] += s2z+s2zv;
      }

      if (newton_bond || i3 < nlocal) {
            f3[0] = s3x+s3xv;
            f3[1] = s3y+s3yv;
            f3[2] = s3z+s3zv;
            f[i3][0] += f3[0];
            f[i3][1] += f3[1];
            f[i3][2] += f3[2];
      }

    	if (stress_ind == 2){
      	    vv = 2.0*datt[m+nm]*coefv/n_atoms[m];
            ff[0] = d21x*s2x + d31x*s3x + (d21x*(s2xv-s1xv)+d31x*(s3xv-s1xv)+d32x*(s3xv-s2xv))/3.0 + vv;
            ff[1] = d21y*s2y + d31y*s3y + (d21y*(s2yv-s1yv)+d31y*(s3yv-s1yv)+d32y*(s3yv-s2yv))/3.0 + vv;
            ff[2] = d21z*s2z + d31z*s3z + (d21z*(s2zv-s1zv)+d31z*(s3zv-s1zv)+d32z*(s3zv-s2zv))/3.0 + vv;
            ff[3] = d21x*s2y + d31x*s3y + (d21x*(s2yv-s1yv)+d31x*(s3yv-s1yv)+d32x*(s3yv-s2yv))/3.0;
            ff[4] = d21x*s2z + d31x*s3z + (d21x*(s2zv-s1zv)+d31x*(s3zv-s1zv)+d32x*(s3zv-s2zv))/3.0;
            ff[5] = d21y*s2z + d31y*s3z + (d21y*(s2zv-s1zv)+d31y*(s3zv-s1zv)+d32y*(s3zv-s2zv))/3.0;
      
            for(int k = 0; k < n_stress_tot; ++k)
                modify->fix[stress_list[k]]->virial4(i1, i2, i3, ff);
      }

      if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
			 -d21x,-d21y,-d21z,d32x,d32y,d32z);
    }
    
    else { ///// this is for normal angle potential
        
				// 1st bond
    	delx1 = x[i1][0] - x[i2][0];
    	dely1 = x[i1][1] - x[i2][1];
    	delz1 = x[i1][2] - x[i2][2];
    	domain->minimum_image(delx1,dely1,delz1);
      rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
      r1 = sqrt(rsq1);

        // 2nd bond
			delx2 = x[i3][0] - x[i2][0];
    	dely2 = x[i3][1] - x[i2][1];
    	delz2 = x[i3][2] - x[i2][2];
    	domain->minimum_image(delx2,dely2,delz2);
    	rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
    	r2 = sqrt(rsq2);

    	// angle (cos and sin)
			c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    	c /= r1*r2;

    	if (c > 1.0) c = 1.0;
    	if (c < -1.0) c = -1.0;

    	s = sqrt(1.0 - c*c);
    	if (s < SMALL) s = SMALL;
    	s = 1.0/s;

    	// force & energy

    	dtheta = acos(c) - a0[type]/180.0 * PI;
    //    fprintf(stderr,"dtheta = %lf, a0[type] = %lf, acos(c) = %lf at time step %ld\n", dtheta, a0[type],acos(c), update->ntimestep);
    	tk = ka[type] * dtheta;

      if (eflag) eangle = tk*dtheta;

    	a = -2.0 * tk * s;
    	a11 = a*c / rsq1;
    	a12 = -a / (r1*r2);
    	a22 = a*c / rsq2;

    	f1[0] = a11*delx1 + a12*delx2;
    	f1[1] = a11*dely1 + a12*dely2;
    	f1[2] = a11*delz1 + a12*delz2;
    	f3[0] = a22*delx2 + a12*delx1;
    	f3[1] = a22*dely2 + a12*dely1;
    	f3[2] = a22*delz2 + a12*delz1;

    	// apply force to each of 3 atoms

    	if (newton_bond || i1 < nlocal) {
      	    f[i1][0] += f1[0];
            f[i1][1] += f1[1];
            f[i1][2] += f1[2];
    	}

    	if (newton_bond || i2 < nlocal) {
            f[i2][0] -= f1[0] + f3[0];
            f[i2][1] -= f1[1] + f3[1];
            f[i2][2] -= f1[2] + f3[2];
    	}

       if (newton_bond || i3 < nlocal) {
      	    f[i3][0] += f3[0];
            f[i3][1] += f3[1];
            f[i3][2] += f3[2];
       }

    if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                         delx1,dely1,delz1,delx2,dely2,delz2);
    }
  }
}

/* ---------------------------------------------------------------------- */

void AngleAreaVolume::allocate()
{
  int i; 

  allocated = 1;
  int n = atom->nangletypes;
  init_on = 0;

  ka = (double *) memory->smalloc((n+1)*sizeof(double),"angle:ka");
  a0 = (double *) memory->smalloc((n+1)*sizeof(double),"angle:a0");
  kv = (double *) memory->smalloc((n+1)*sizeof(double),"angle:kv");
  v0 = (double *) memory->smalloc((n+1)*sizeof(double),"angle:v0");
  kl = (double *) memory->smalloc((n+1)*sizeof(double),"angle:kl");
  t_remod = (double *) memory->smalloc((n+1)*sizeof(double),"angle:t_remod");

  setflag = (int *) memory->smalloc((n+1)*sizeof(int),"angle:setflag");
  for (i = 1; i <= n; i++) setflag[i] = 0;

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleAreaVolume::coeff(int which, int narg, char **arg)
{
  int i;

  if (which != 0) error->all("Illegal coeffs for this angle style");
  if (narg != 7) error->all("Incorrect args in angle_coeff command");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nangletypes,ilo,ihi);

  double ka_one = atof(arg[1]);
  double a0_one = atof(arg[2]);
  double kv_one = atof(arg[3]);
  double v0_one = atof(arg[4]);
  double kl_one = atof(arg[5]);
  double t_remod_one = atof(arg[6]);

  int count = 0;
  for (i = ilo; i <= ihi; i++) {
    ka[i] = ka_one;
    a0[i] = a0_one;
    kv[i] = kv_one;
    v0[i] = v0_one;
    kl[i] = kl_one;
    t_remod[i] = t_remod_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all("Incorrect args in angle_coeff command");

}

/* ---------------------------------------------------------------------- */

double AngleAreaVolume::equilibrium_angle(int i)
{
  return -1;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleAreaVolume::write_restart(FILE *fp)
{
  fwrite(&ka[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&a0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kv[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&v0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kl[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&t_remod[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them 
------------------------------------------------------------------------- */

void AngleAreaVolume::read_restart(FILE *fp)
{
  int i;

  allocate();

  if (comm->me == 0){ 
    fread(&ka[1],sizeof(double),atom->nangletypes,fp);
    fread(&a0[1],sizeof(double),atom->nangletypes,fp);
    fread(&kv[1],sizeof(double),atom->nangletypes,fp);
    fread(&v0[1],sizeof(double),atom->nangletypes,fp);
    fread(&kl[1],sizeof(double),atom->nangletypes,fp);
    fread(&t_remod[1],sizeof(double),atom->nangletypes,fp);  
  }
  MPI_Bcast(&ka[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kv[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&v0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kl[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&t_remod[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ---------------------------------------------------------------------- */

double AngleAreaVolume::single(int type, int i1, int i2, int i3)
{
    return -1;
}

/* ---------------------------------------------------------------------- */

double AngleAreaVolume::cell_remodel(int i1, int i2, int i3, double tau)
{
  if (atom->molecular == 0) return 1.0;
  if (molecule == NULL){
    int ifix;
    for (ifix = 0; ifix < modify->nfix; ifix++)
      if (strcmp("molecule",modify->fix[ifix]->style) == 0) break;
  	if (ifix == modify->nfix) error->one("AngleAreaVolume: FixMolecule style is not defined.");
 	 	molecule = dynamic_cast<FixMolecule *> (modify->fix[ifix]);
  }

  double tact = 0.0;
  if (molecule)
    tact = molecule->list_tact[ atom->molecule[i1] ];
  if (tact == 0.0) return 1.0;

  double coeff = exp( -(update->dt*update->ntimestep-tact) / tau );
  coeff = (coeff > 0.1) ? coeff : 0.1;

  return coeff;
}

/* ----------------------------------------------------------------------*/
