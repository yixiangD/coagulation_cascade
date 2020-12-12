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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_stat.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "error.h"
#include "math.h"
#include "memory.h"
#include "common.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStat::FixStat(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
	int i;
	double alpha, len, theta, u;

  if (strcmp(arg[3],"cyl") && strcmp(arg[3],"box") && strcmp(arg[3],"sten"))
  	error->all("Illegal fix stat command");

	if (strcmp(arg[2],"stat/velT") == 0 && strcmp(arg[3],"box") == 0){
    if (narg != 12 && narg != 18) error->all("Illegal fix stat command");
    nx = atoi(arg[4]);
    ny = atoi(arg[5]);
    nz = atoi(arg[6]);
    st_start = atoi(arg[7]);
    nevery = atoi(arg[8]);
    dump_each = atoi(arg[9]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    if (narg == 18){
      xlo = atof(arg[10]);
      xhi = atof(arg[11]);
      ylo = atof(arg[12]);
      yhi = atof(arg[13]);
      zlo = atof(arg[14]);
      zhi = atof(arg[15]);
      sprintf(fname,arg[16]);
      index = atoi(arg[17]);
    }
    else{
      xlo = domain->boxlo[0];
      xhi = domain->boxhi[0];
      ylo = domain->boxlo[1];
      yhi = domain->boxhi[1];
      zlo = domain->boxlo[2];
      zhi = domain->boxhi[2];
      sprintf(fname,arg[10]);
      index = atoi(arg[11]);
    }
		if (index < 1 || index > CTYPES) error->all("Species number out of bound");
	}
	if ((strcmp(arg[2],"stat/dens") == 0 || strcmp(arg[2],"stat/vel") == 0) && strcmp(arg[3],"box") == 0){
  	if (narg != 11 && narg != 17) error->all("Illegal fix stat command");
  	nx = atoi(arg[4]);
  	ny = atoi(arg[5]);
  	nz = atoi(arg[6]);
  	st_start = atoi(arg[7]);
  	nevery = atoi(arg[8]);
  	dump_each = atoi(arg[9]);
  	if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
  	if (narg == 17){
    	xlo = atof(arg[10]);
    	xhi = atof(arg[11]);
    	ylo = atof(arg[12]);
    	yhi = atof(arg[13]);
    	zlo = atof(arg[14]);
    	zhi = atof(arg[15]);
    	sprintf(fname,arg[16]);
  	}
  	else{
    	xlo = domain->boxlo[0];
    	xhi = domain->boxhi[0];
    	ylo = domain->boxlo[1];
    	yhi = domain->boxhi[1];
    	zlo = domain->boxlo[2];
    	zhi = domain->boxhi[2];
    	sprintf(fname,arg[10]);
  	}
	}
  if (strcmp(arg[2],"stat/velT/cyl") == 0 && strcmp(arg[3],"cyl") == 0){
    if (narg != 17 && narg != 23) error->all("Illegal fix stat command");
    sprintf(dir,arg[4]);
    R = atof(arg[5]);
    centx = atof(arg[6]);
    centy = atof(arg[7]);
    centz = atof(arg[8]);
    nx = atoi(arg[9]);
    ny = atoi(arg[10]);
    nz = atoi(arg[11]);
    st_start = atoi(arg[12]);
    nevery = atoi(arg[13]);
    dump_each = atoi(arg[14]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    if (narg == 23){
      xlo = atof(arg[15]);
      xhi = atof(arg[16]);
      ylo = atof(arg[17]);
      yhi = atof(arg[18]);
      zlo = atof(arg[19]);
      zhi = atof(arg[20]);
      sprintf(fname,arg[21]);
			index = atoi(arg[22]);
    }
    else{
      xlo = domain->boxlo[0];
      xhi = domain->boxhi[0];
      ylo = domain->boxlo[1];
      yhi = domain->boxhi[1];
      zlo = domain->boxlo[2];
      zhi = domain->boxhi[2];
      sprintf(fname,arg[15]);
			index = atoi(arg[16]);
    }
		if (index < 1 || index > CTYPES) error->all("Species number out of bound");
  }
	if ((strcmp(arg[2],"stat/dens/cyl") == 0 || strcmp(arg[2],"stat/vel/cyl") == 0) && strcmp(arg[3],"cyl") == 0){
    if (narg != 16 && narg != 22) error->all("Illegal fix stat command");
		sprintf(dir,arg[4]);
		R = atof(arg[5]);
		centx = atof(arg[6]);
		centy = atof(arg[7]);
		centz = atof(arg[8]);
    nx = atoi(arg[9]);
    ny = atoi(arg[10]);
    nz = atoi(arg[11]);
    st_start = atoi(arg[12]);
    nevery = atoi(arg[13]);
    dump_each = atoi(arg[14]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    if (narg == 22){
      xlo = atof(arg[15]);
      xhi = atof(arg[16]);
      ylo = atof(arg[17]);
      yhi = atof(arg[18]);
      zlo = atof(arg[19]);
      zhi = atof(arg[20]);
      sprintf(fname,arg[21]);
    }
    else{
      xlo = domain->boxlo[0];
      xhi = domain->boxhi[0];
      ylo = domain->boxlo[1];
      yhi = domain->boxhi[1];
      zlo = domain->boxlo[2];
      zhi = domain->boxhi[2];
      sprintf(fname,arg[15]);
    }
	}
  if (strcmp(arg[2],"stat/velT/cyl") == 0 && strcmp(arg[3],"sten") == 0){
    if (narg != 25) error->all("Illegal fix stat command");
    sprintf(dir,arg[4]);
    R = atof(arg[5]);
    alpha = atof(arg[6]);
    len = atof(arg[7]);
    centx = atof(arg[8]);
    centy = atof(arg[9]);
    centz = atof(arg[10]);
    nx = atoi(arg[11]);
    ny = atoi(arg[12]);
    nz = atoi(arg[13]);
    st_start = atoi(arg[14]);
    nevery = atoi(arg[15]);
    dump_each = atoi(arg[16]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    xlo = atof(arg[17]);
    xhi = atof(arg[18]);
    ylo = atof(arg[19]);
    yhi = atof(arg[20]);
    zlo = atof(arg[21]);
    zhi = atof(arg[22]);
    sprintf(fname,arg[23]);
		index = atoi(arg[24]);
		if (index < 1 || index > CTYPES) error->all("Species number out of bound");
  }
  if ((strcmp(arg[2],"stat/dens/cyl") == 0 || strcmp(arg[2],"stat/vel/cyl") == 0) && strcmp(arg[3],"sten") == 0){
    if (narg != 24) error->all("Illegal fix stat command");
    sprintf(dir,arg[4]);
    R = atof(arg[5]);
    alpha = atof(arg[6]);
    len = atof(arg[7]);
    centx = atof(arg[8]);
    centy = atof(arg[9]);
    centz = atof(arg[10]);
    nx = atoi(arg[11]);
    ny = atoi(arg[12]);
    nz = atoi(arg[13]);
    st_start = atoi(arg[14]);
    nevery = atoi(arg[15]);
    dump_each = atoi(arg[16]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    xlo = atof(arg[17]);
    xhi = atof(arg[18]);
    ylo = atof(arg[19]);
    yhi = atof(arg[20]);
    zlo = atof(arg[21]);
    zhi = atof(arg[22]);
    sprintf(fname,arg[23]);
  }
  if (strcmp(arg[2],"stat/all") == 0){
    if (narg != 10 && narg != 16) error->all("Illegal fix stat command");
    nx = atoi(arg[4]);
    ny = atoi(arg[5]);
    nz = atoi(arg[6]);
    st_start = atoi(arg[7]);
    nevery = atoi(arg[8]);
    dump_each = atoi(arg[9]);
    if ((nx < 1)||(ny < 1)||(nz < 1)) error->all("Illegal division for statistics");
    if (narg == 17){
      xlo = atof(arg[10]);
      xhi = atof(arg[11]);
      ylo = atof(arg[12]);
      yhi = atof(arg[13]);
      zlo = atof(arg[14]);
      zhi = atof(arg[15]);
    }
    else{
      xlo = domain->boxlo[0];
      xhi = domain->boxhi[0];
      ylo = domain->boxlo[1];
      yhi = domain->boxhi[1];
      zlo = domain->boxlo[2];
      zhi = domain->boxhi[2];
    }
  }

  if ((xlo >= xhi)||(ylo >= yhi)||(zlo >= zhi)) error->all("Illegal coordinates for statistics");
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

  if (strcmp(arg[3],"cyl") == 0 || strcmp(arg[3],"sten") == 0){
	  if (strcmp(dir,"x") == 0){
    	cent2 = centy; cent3 = centz;
    	low = xlo; high = xhi;
    	n1 = nx; n2 = ny; n3 = nz;
  	}else if (strcmp(dir,"y") == 0){
    	cent2 = centx; cent3 = centz;
    	low = ylo; high = yhi;
    	n1 = ny; n2 = nx; n3 = nz;
  	}else if (strcmp(dir,"z") == 0){
    	cent2 = centx; cent3 = centy;
    	low = zlo; high = zhi;
    	n1 = nz; n2 = nx; n3 = ny;
  	}
	}
	if (strcmp(arg[3],"sten") == 0){
		profile = (double *) memory->smalloc(n1*sizeof(double),"fix_stat:profile");
		u = -len;
		for (i = 0; i < n1; i++){
			theta = M_PI*u/len;
			profile[i] = R * (1.0 - alpha*(1.0+cos(theta))/2.0);
			u += 2.0*len/static_cast<double> (n1-1);
		}
	}
	sprintf(stype,arg[3]);
}

/* ---------------------------------------------------------------------- */

FixStat::~FixStat()
{
	if (strcmp(stype,"sten") == 0){
		memory->sfree(profile);
	}
}

/* ---------------------------------------------------------------------- */

int FixStat::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

int FixStat::map_index(double x, double y, double z)
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
    is = static_cast<int> ((x - xlo)*nx/xs);
    js = static_cast<int> ((y - ylo)*ny/ys);
    ks = static_cast<int> ((z - zlo)*nz/zs); 
    ind = 1;
    if (is > nx-1) is--;
    if (js > ny-1) js--;
    if (ks > nz-1) ks--;
  }
  return ind;
}

/* ---------------------------------------------------------------------- */

int FixStat::map_index_cyl(double x, double y, double z)
{
  double theta, theta1, rr;
	double x1, x2, x3;
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
  
	if (strcmp(dir,"x") == 0){
		x1 = x; x2 = y; x3 = z;
	}else if (strcmp(dir,"y") == 0){
    x1 = y; x2 = x; x3 = z;
	}else if (strcmp(dir,"z") == 0){
    x1 = z; x2 = x; x3 = y;
	}
	
  rr = sqrt((x2-cent2)*(x2-cent2) + (x3-cent3)*(x3-cent3));

  if (x1>=low && x1<high){
    theta = acos((x2-cent2)/rr);
    theta1 = asin((x3-cent3)/rr);
    if (theta1 < 0.0)
      theta = 2.0*M_PI - theta;
    is = static_cast<int> ((x1 - low)*n1/(high - low));
    if (is > n1-1) is--;
		if (strcmp(stype,"cyl") == 0)
			js = static_cast<int> (rr*n2/R);
		else if (strcmp(stype,"sten") == 0)
			js = static_cast<int> (rr*n2/profile[is]);
    ks = static_cast<int> (0.5*theta*n3/M_PI);
    ind = 1;
    if (js > n2-1) js--;
    if (ks > n3-1) ks--;
  }
  return ind;
} 

/* ---------------------------------------------------------------------- */
/*
int FixStat::map_index_cyl(double x, double y, double z)
{
  double theta, theta1,rr;
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
   
  rr = sqrt((y-ylo)*(y-ylo) + (x-xlo)*(x-xlo));
  
  if (z>=zlo && z<zhi && rr<xhi){ 
    theta = acos((x-xlo)/rr);
    theta1 = asin((y-ylo)/rr);
    if (theta1 < 0.0)
      theta = 2.0*M_PI - theta; 
    is = static_cast<int> ((z - zlo)*nz/zs);
    js = static_cast<int> (rr*nx/xhi);
    ks = static_cast<int> (0.5*theta*ny/M_PI);
    ind = 1;
    if (is > nx-1) is--;
    if (js > ny-1) js--;
    if (ks > nz-1) ks--;
  }
  return ind;
} 
*/
/* ---------------------------------------------------------------------- */
