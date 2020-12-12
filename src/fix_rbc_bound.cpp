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
#include "stdio.h"
#include "string.h"
#include "fix_rbc_bound.h"
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "domain.h"
#include "update.h"
#include "stdlib.h"
#include "time.h"
#define BUFFACTOR 1.5
#define BUFMIN 1000
#define BUFEXTRA 100
#define EPS 1.0e-32
#define NUMH 500

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
using namespace LAMMPS_NS;
/* ---------------------------------------------------------------------- */

FixBoundRBC::FixBoundRBC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  int igroup;
  char gn[50], gn1[50];
  
  if (narg != 8) error->all("Illegal fix bound/rbc command");

  sprintf(gn,arg[3]);
  sprintf(gn1,arg[4]); 
  ind_bounce = atoi(arg[5]);   // 0 - bounce-back old; 1 - bounce-back new; 2 - bounce-forward 
  d_cut = atof(arg[6]);        
  comm_cut = atof(arg[7]);  
  d_cut_sq = d_cut*d_cut; 
  size_border = 7;
  size_comm = 7;
  max_count = 3;
  igroup = group->find(gn);
  groupbit_inner = group->bitmask[igroup];
  igroup = group->find(gn1);
  groupbit_comm = group->bitmask[igroup];
  if (domain->xperiodic)
    prd[0] = domain->xprd;
  else
    prd[0] = 0;
  if (domain->yperiodic)
    prd[1] = domain->yprd;
  else
    prd[1] = 0;
  if (domain->zperiodic)
    prd[2] = domain->zprd;
  else
    prd[2] = 0;

  ind_ntri = ind_h = tri_mol = ind_norm = NULL;
  ind_tri = tri_ind = NULL;
  ind_x = ind_v = NULL;
  norm = edge1 = edge2 = dif1 = dif2 = dd12 = NULL;
  num_faces = NULL;
  ind_faces = NULL;
  sendnum = recvnum = size_comm_send = size_comm_recv = NULL; 
  maxsendlist = NULL;
  sendlist = NULL;
  buf_send = buf_recv = NULL; 
  slablo = slabhi = NULL;
  
  max_faces = 40;
  nnmax = atom->nlocal;
  ind_faces = memory->create_2d_int_array(nnmax,max_faces,"fix_bound_rbc:ind_faces");  
  num_faces = (int *) memory->smalloc(nnmax*sizeof(int),"fix_bound_rbc:num_faces");
}

/* ---------------------------------------------------------------------- */

FixBoundRBC::~FixBoundRBC()
{
  int i;

  memory->sfree(ind_ntri);
  memory->sfree(ind_h);
  memory->sfree(tri_mol);
  memory->sfree(tri_check);
  memory->sfree(point_check);
  memory->destroy_2d_int_array(ind_tri);
  memory->destroy_2d_int_array(tri_ind);
  memory->destroy_2d_double_array(ind_x);
  memory->destroy_2d_double_array(ind_v); 

  memory->sfree(sendnum);
  memory->sfree(recvnum);
  memory->sfree(size_comm_send);
  memory->sfree(size_comm_recv);
  memory->sfree(slablo);
  memory->sfree(slabhi);    

  memory->sfree(maxsendlist);
  if (sendlist) 
    for (i = 0; i < maxswap; i++) 
      memory->sfree(sendlist[i]);
  memory->sfree(sendlist);

  memory->sfree(buf_send);
  memory->sfree(buf_recv);

  memory->sfree(ind_norm);   
  memory->destroy_2d_double_array(norm);
  memory->destroy_2d_double_array(edge1);
  memory->destroy_2d_double_array(edge2);
  memory->destroy_2d_double_array(dif1);
  memory->destroy_2d_double_array(dif2);
  memory->destroy_2d_double_array(dd12); 
  memory->destroy_2d_int_array(ind_faces);
  memory->sfree(num_faces);
}

/* ---------------------------------------------------------------------- */

int FixBoundRBC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}


/* ---------------------------------------------------------------------- */

void FixBoundRBC::setup(int vflag)
{
  int i,j,k,offset,max_tri;
  int iswap, dim, ineed, nbox;
  int *tag = atom->tag; 
  int *molecule = atom->molecule; 
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int *displs, *rcounts;  

  rcounts = (int *) memory->smalloc(comm->nprocs*sizeof(int),"fix_bound_rbc:rcounts");
  displs = (int *) memory->smalloc(comm->nprocs*sizeof(int),"fix_bound_rbc:displs");   
  if (comm->nprocs > 1){
    MPI_Gather(&nanglelist, 1, MPI_INT, rcounts, 1, MPI_INT, 0, world);
    if (comm->me == 0){
      ntri = 0;  
      offset = 0;
      for (i = 0; i < comm->nprocs; i++) {
        ntri += rcounts[i];  
        displs[i] = offset;
        rcounts[i] *= 1;
        offset += rcounts[i];
      }
    }
    MPI_Bcast(&ntri,1,MPI_INT,0,world); 
  } else{
    ntri = nanglelist; 
  }
  
  //int bff[ntri*4], bff1[ntri*4];
  int *bff, *bff1, *bff2;
  bff = new int[4*ntri];
  bff1 = new int[ntri];
  bff2 = new int[ntri];
  /*
  if (comm->nprocs > 1){
    for (i = 0; i < nanglelist; i++) {
      bff1[4*i] = tag[anglelist[i][0]];
      bff1[4*i+1] = tag[anglelist[i][1]];
      bff1[4*i+2] = tag[anglelist[i][2]];
      bff1[4*i+3] = molecule[anglelist[i][0]];
    }
    printf("fixrbcbound:step:ffff, cpu is %d\n", comm->me);
    MPI_Barrier(world);
    MPI_Gatherv(bff1, 4*nanglelist, MPI_INT, bff, rcounts, displs, MPI_INT, 0, world);
    MPI_Barrier(world);
    printf("fixrbcbound:step:gggg, cpu is %d, 4*nanglelist is %d\n", comm->me, 4*nanglelist);
    MPI_Barrier(world);
    MPI_Bcast(bff,4*ntri,MPI_INT,0,world);
    printf("fixrbcbound:step:hhhh, cpu is %d\n", comm->me);
    MPI_Barrier(world);
  }*/ 
  if (comm->nprocs > 1){
    for (i = 0; i < nanglelist; i++) 
      bff1[i] = tag[anglelist[i][0]];
    
    MPI_Gatherv(bff1, nanglelist, MPI_INT, bff2, rcounts, displs, MPI_INT, 0, world);
    for(i = 0; i < ntri; ++i)
      bff[4*i] = bff2[i];

    for (i = 0; i < nanglelist; i++) 
      bff1[i] = tag[anglelist[i][1]];
    MPI_Gatherv(bff1, nanglelist, MPI_INT, bff2, rcounts, displs, MPI_INT, 0, world);
    for(i = 0; i < ntri; ++i)
      bff[4*i+1] = bff2[i];

    for (i = 0; i < nanglelist; i++) 
      bff1[i] = tag[anglelist[i][2]];
    MPI_Gatherv(bff1, nanglelist, MPI_INT, bff2, rcounts, displs, MPI_INT, 0, world);
    for(i = 0; i < ntri; ++i)
      bff[4*i+2] = bff2[i];

    for (i = 0; i < nanglelist; i++) 
      bff1[i] = molecule[anglelist[i][0]];
    MPI_Gatherv(bff1, nanglelist, MPI_INT, bff2, rcounts, displs, MPI_INT, 0, world);
    for(i = 0; i < ntri; ++i)
      bff[4*i+3] = bff2[i];

    MPI_Bcast(bff,4*ntri,MPI_INT,0,world);
  } 
  else{
    for (i = 0; i < ntri; i++) {
      bff[4*i] = tag[anglelist[i][0]];
      bff[4*i+1] = tag[anglelist[i][1]];
      bff[4*i+2] = tag[anglelist[i][2]];
      bff[4*i+3] = molecule[anglelist[i][0]];
    } 
  }
  
  tri_check = (int *) memory->smalloc(ntri*sizeof(int),"fix_bound_rbc:tri_check");
  tri_ind = memory->create_2d_int_array(ntri,3,"fix_bound_rbc:tri_ind");
  tri_mol = (int *) memory->smalloc(ntri*sizeof(int),"fix_bound_rbc:tri_mol"); 
  min_ind = 999999999;
  max_ind = -1;
  for (i = 0; i < ntri; i++){
    tri_mol[i] = bff[4*i+3]; 
    for (j = 0; j < 3; j++){
      tri_ind[i][j] = bff[4*i+j];
      if (tri_ind[i][j] < min_ind) min_ind = tri_ind[i][j];
      if (tri_ind[i][j] > max_ind) max_ind = tri_ind[i][j];   
    } 
  }
  ind_tot = max_ind - min_ind + 1;

  point_check = (int *) memory->smalloc(ind_tot*sizeof(int),"fix_bound_rbc:point_check");
  ind_ntri = (int *) memory->smalloc(ind_tot*sizeof(int),"fix_bound_rbc:ind_ntri");
  ind_h = (int *) memory->smalloc(ind_tot*sizeof(int),"fix_bound_rbc:ind_h");
  ind_x = memory->create_2d_double_array(ind_tot,3,"fix_bound_rbc:ind_x");
  ind_v = memory->create_2d_double_array(ind_tot,3,"fix_bound_rbc:ind_v"); 
  
  for (i = 0; i < ind_tot; i++)
    ind_ntri[i] = 0;

  for (i = 0; i < ntri; i++)
    for (j = 0; j < 3; j++){
      tri_ind[i][j] -= min_ind;
      ind_ntri[tri_ind[i][j]]++;   
    }
  max_tri = -1;
  for (i = 0; i < ind_tot; i++){
    if (ind_ntri[i] > max_tri) max_tri = ind_ntri[i];
    ind_ntri[i] = 0;
  }

  ind_tri = memory->create_2d_int_array(ind_tot,max_tri,"fix_bound_rbc:ind_tri"); 
  for (i = 0; i < ntri; i++)
    for (j = 0; j < 3; j++){
      k = tri_ind[i][j];
      ind_tri[k][ind_ntri[k]] = i;
      ind_ntri[k]++;   
    }  

  /*char f_name[FILENAME_MAX]; 
  FILE* out_stat;
  sprintf(f_name,"output_%d.plt",comm->me);
  out_stat=fopen(f_name,"w");
  fprintf(out_stat,"min_ind=%d, max_ind=%d, ind_tot=%d  \n", min_ind, max_ind, ind_tot);
  fprintf(out_stat,"ntri=%d, max_tri=%d  \n", ntri, max_tri); 
  fprintf(out_stat,"triangles  \n"); 
  for (i = 0; i < ntri; i++)
    fprintf(out_stat,"%d    %ld %ld %ld mol=%d \n", i, tri_ind[i][0], tri_ind[i][1], tri_ind[i][2],tri_mol[i]);
  fprintf(out_stat,"indexes  \n"); 
  for (i = 0; i < ind_tot; i++){
    fprintf(out_stat,"%d   %ld    ", i, ind_ntri[i]);
    for (j = 0; j < ind_ntri[i]; j++)  
      fprintf(out_stat,"%ld ", ind_tri[i][j]);
    fprintf(out_stat," \n"); 
  }
  fclose(out_stat);*/  

   
  norm = memory->create_2d_double_array(ntri,4,"fix_bound_rbc:norm");
  edge1 = memory->create_2d_double_array(ntri,3,"fix_bound_rbc:edge1");
  edge2 = memory->create_2d_double_array(ntri,3,"fix_bound_rbc:edge2"); 
  dif1 = memory->create_2d_double_array(ntri,3,"fix_bound_rbc:dif1");
  dif2 = memory->create_2d_double_array(ntri,3,"fix_bound_rbc:dif2");
  dd12 = memory->create_2d_double_array(ntri,3,"fix_bound_rbc:dd12"); 
  ind_norm = (int *) memory->smalloc(ntri*sizeof(int),"fix_bound_rbc:ind_norm"); 
 
  maxsend = BUFMIN;
  buf_send = (double *) 
    memory->smalloc((maxsend+BUFEXTRA)*sizeof(double),"fix_bound_rbc:buf_send");
  maxrecv = BUFMIN;
  buf_recv = (double *) 
    memory->smalloc(maxrecv*sizeof(double),"fix_bound_rbc:buf_recv");
 
  maxswap = comm->maxswap;
  sendlist = (int **) memory->smalloc(maxswap*sizeof(int *),"sendlist");
  maxsendlist = (int *) memory->smalloc(maxswap*sizeof(int),"maxsendlist");
  for (int i = 0; i < maxswap; i++) {
    maxsendlist[i] = BUFMIN;
    sendlist[i] = (int *) memory->smalloc(BUFMIN*sizeof(int),"sendlist[i]");
  } 

  sendnum = (int *) memory->smalloc(maxswap*sizeof(int),"fix_bound_rbc:sendnum");
  recvnum = (int *) memory->smalloc(maxswap*sizeof(int),"fix_bound_rbc:recvnum");
  size_comm_send = (int *) memory->smalloc(maxswap*sizeof(int),"fix_bound_rbc:size");
  size_comm_recv = (int *) memory->smalloc(maxswap*sizeof(int),"fix_bound_rbc:size");
  slablo = (double *) memory->smalloc(maxswap*sizeof(double),"fix_bound_rbc:slablo");
  slabhi = (double *) memory->smalloc(maxswap*sizeof(double),"fix_bound_rbc:slabhi");     

  iswap = 0;
  for (dim = 0; dim < 3; dim++) 
    for (ineed = 0; ineed < 2*comm->need[dim]; ineed++) {
      if (ineed % 2 == 0) {
	nbox = comm->myloc[dim] + ineed/2;
	slablo[iswap] = domain->boxlo[dim] + 
	  nbox * domain->prd[dim] / comm->procgrid[dim];
	slabhi[iswap] = domain->sublo[dim] + comm_cut;
	slabhi[iswap] = MIN(slabhi[iswap],domain->boxlo[dim] + 
			    (nbox+1) * domain->prd[dim] / comm->procgrid[dim]);
	if (comm->myloc[dim] == 0 && domain->periodicity[dim] == 0) slabhi[iswap] = slablo[iswap];
      } else {
	nbox = comm->myloc[dim] - ineed/2;
	slablo[iswap] = domain->subhi[dim] - comm_cut;
	slablo[iswap] = MAX(slablo[iswap],domain->boxlo[dim] +
			    nbox * domain->prd[dim] / comm->procgrid[dim]);
	slabhi[iswap] = domain->boxlo[dim] + 
	  (nbox+1) * domain->prd[dim] / comm->procgrid[dim];
	if (comm->myloc[dim] == comm->procgrid[dim]-1 && domain->periodicity[dim] == 0) slabhi[iswap] = slablo[iswap];
      }
      slablo[iswap] -= neighbor->skin;
      slabhi[iswap] += neighbor->skin;
      iswap++;
    }

  memory->sfree(rcounts);
  memory->sfree(displs);


  delete[] bff;
  delete[] bff1;
  delete[] bff2;
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::initial_integrate(int vflag)
{
  int i,j,k,m,kk,ind,cond,bounce,ccount,cplane;
  int i_x[NUMH],shift[NUMH][3];
  double tt,dtt,u_inv,u0,u1,uu,uun; 
  double dl[3],t_x[NUMH],x_x[3];
  double vv[3],xd[3],xc[3],cf[4],a1[3],a2[3],xp[3],nnp[4];
  double dtv = update->dt; 

  double **xn = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  //double **i_v = atom->i_v;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *tag = atom->tag; 

  time_t start,end;
  double tdiff;

  if (neighbor->ago == 0) { 
    local_faces(0);
    if (comm->nprocs > 1) borders();
    
    time (&start);
    face_decide();
    time (&end);
    tdiff = difftime (end,start);
    //printf("time between facedecide is %f, cpu is %d\n", tdiff, comm->me);

  } else{
    local_faces(1);
    
    if (comm->nprocs > 1) communicate();

  }  
 
  for(i = 0; i < ntri; i++)
    ind_norm[i] = 0;

  time(&start);
  for(i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      cond = 1;
      ccount = 0;
      dtt = dtv;
      cplane = -1;
      while (cond && ccount < max_count){
        for(j = 0; j < 3; j++){
          //dl[j] = i_v[i][j]*dtt;
          dl[j] = v[i][j]*dtt;
          xp[j] = xn[i][j] - dl[j];
	}
        ind = 0;
        for(k = 0; k < num_faces[i]; k++){
          tt = -1.0;
          m = ind_faces[i][k];
          if (m != cplane){
            if (ind_norm[m] == 0) calc_norm(m);
            for(j = 0; j < 3; j++){
              x_x[j] = ind_x[tri_ind[m][0]][j] - xp[j]; 
              if (fabs(x_x[j]) > 0.5*prd[j]){
                if (x_x[j] > 0)
                  shift[ind][j] = 1;
                else
                  shift[ind][j] = -1;
	      } else
                shift[ind][j] = 0; 
              a1[j] = edge1[m][j] - dif1[m][j]*dtt;
              a2[j] = edge2[m][j] - dif2[m][j]*dtt; 
	    }
            nnp[0] = a1[1]*a2[2] - a1[2]*a2[1];
            nnp[1] = a1[2]*a2[0] - a1[0]*a2[2];
            nnp[2] = a1[0]*a2[1] - a1[1]*a2[0];
            nnp[3] = - nnp[0]*(ind_x[tri_ind[m][0]][0] - ind_v[tri_ind[m][0]][0]*dtt) - nnp[1]*(ind_x[tri_ind[m][0]][1] - ind_v[tri_ind[m][0]][1]*dtt) - nnp[2]*(ind_x[tri_ind[m][0]][2] - ind_v[tri_ind[m][0]][2]*dtt);    
            u0 = nnp[0]*(xp[0]+prd[0]*shift[ind][0]) + nnp[1]*(xp[1]+prd[1]*shift[ind][1]) 
              + nnp[2]*(xp[2]+prd[2]*shift[ind][2]) + nnp[3];
          if ((u0 <= 0.0 && mask[i]&groupbit_inner) || (u0 >= 0.0 && !(mask[i]&groupbit_inner))){
            u1 = norm[m][0]*(xn[i][0]+prd[0]*shift[ind][0]) + norm[m][1]*(xn[i][1]+prd[1]*shift[ind][1]) 
              + norm[m][2]*(xn[i][2]+prd[2]*shift[ind][2]) + norm[m][3];
            if (u0*u1 <= 0.0){ 
              for(j = 0; j < 3; j++){
                vv[j] = dl[j] - dtt*ind_v[tri_ind[m][0]][j];
                xd[j] = xp[j]+prd[j]*shift[ind][j] - ind_x[tri_ind[m][0]][j] + dtt*ind_v[tri_ind[m][0]][j];
	      }
              xc[0] = dtt*(a1[1]*dif2[m][2] - a1[2]*dif2[m][1] + dif1[m][1]*a2[2] - dif1[m][2]*a2[1]);
              xc[1] = dtt*(a1[2]*dif2[m][0] - a1[0]*dif2[m][2] + dif1[m][2]*a2[0] - dif1[m][0]*a2[2]);
              xc[2] = dtt*(a1[0]*dif2[m][1] - a1[1]*dif2[m][0] + dif1[m][0]*a2[1] - dif1[m][1]*a2[0]);
              cf[0] = dtt*dtt*(dd12[m][0]*vv[0] + dd12[m][1]*vv[1] + dd12[m][2]*vv[2]);
              cf[1] = dtt*dtt*(dd12[m][0]*xd[0] + dd12[m][1]*xd[1] + dd12[m][2]*xd[2]) + xc[0]*vv[0] + xc[1]*vv[1] + xc[2]*vv[2]; 
              cf[2] = xc[0]*xd[0] + xc[1]*xd[1] + xc[2]*xd[2] + nnp[0]*vv[0] + nnp[1]*vv[1] + nnp[2]*vv[2];
              cf[3] = nnp[0]*xd[0] + nnp[1]*xd[1] + nnp[2]*xd[2];  
              tt = find_root(cf);
              //tt =  solve_quadratic(cf);
	    }
          }
	    if (tt>=0.0 && tt<=1.0){  
              t_x[ind] = tt;
	      i_x[ind] = m;
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
            m = i_x[kk];
            t_x[kk] = t_x[ind-j];
            i_x[kk] = i_x[ind-j];
            t_x[ind-j] = tt;
            i_x[ind-j] = m; 
            for(k = 0; k < 3; k++){
              m = shift[kk][k];
              shift[kk][k] = shift[ind-j][k];
              shift[ind-j][k] = m;
	    }  
          }
	  while(ind){
            bounce = 1;
            m = ind - 1;
            tt = t_x[m];
	    kk = i_x[m];
            for (k = 0; k < 3; k++){
	      x_x[k] = xp[k] + prd[k]*shift[m][k] + tt*dl[k];
	      xc[k] = x_x[k] - ind_x[tri_ind[kk][0]][k] + dtt*(1.0-tt)*ind_v[tri_ind[kk][0]][k];
              a1[k] = edge1[kk][k] - dif1[kk][k]*dtt*(1.0-tt);
              a2[k] = edge2[kk][k] - dif2[kk][k]*dtt*(1.0-tt);                	              
	    }
	    xd[0] = a1[0]*xc[0] + a1[1]*xc[1] + a1[2]*xc[2];
            xd[1] = a2[0]*xc[0] + a2[1]*xc[1] + a2[2]*xc[2];
            vv[0] = a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2];
            vv[1] = a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2];
            vv[2] = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]; 
            u_inv = vv[0]*vv[1] - vv[2]*vv[2];
            u0 = (xd[0]*vv[1] - xd[1]*vv[2])/u_inv;
            u1 = (xd[1]*vv[0] - xd[0]*vv[2])/u_inv;  
	    if (u0 <= 0.0 || u1 <= 0.0 || u0+u1 >= 1.0)
              bounce = 0;

	    if (bounce){
              cplane = kk; 
              ccount++;
              for (k = 0; k < 3; k++)
                vv[k] = (1.0-u0-u1)*ind_v[tri_ind[kk][0]][k] + u0*ind_v[tri_ind[kk][1]][k] + u1*ind_v[tri_ind[kk][2]][k];  
              if (ind_bounce == 2){
                nnp[0] = a1[1]*a2[2] - a1[2]*a2[1];
                nnp[1] = a1[2]*a2[0] - a1[0]*a2[2];
                nnp[2] = a1[0]*a2[1] - a1[1]*a2[0];
                uun = sqrt(nnp[0]*nnp[0] + nnp[1]*nnp[1] + nnp[2]*nnp[2]);
                //uu = 2.0 * ((i_v[i][0]-vv[0])*nnp[0] + (i_v[i][1]-vv[1])*nnp[1] + (i_v[i][2]-vv[2])*nnp[2])/uun;
                uu = 2.0 * ((v[i][0]-vv[0])*nnp[0] + (v[i][1]-vv[1])*nnp[1] + (v[i][2]-vv[2])*nnp[2])/uun;
                for (k = 0; k < 3; k++){
                  xp[k] = x_x[k] - prd[k]*shift[m][k];        
                  //xn[i][k] = xp[k] + (1.0-tt)*dtt*(i_v[i][k] - uu*nnp[k]/uun);      // verlet 
                  xn[i][k] = xp[k] + (1.0-tt)*dtt*(v[i][k] - uu*nnp[k]/uun);      // verlet 
	        }
                for (k = 0; k < 3; k++){
                  //i_v[i][k] = -i_v[i][k] + 2.0*vv[k];     
                  v[i][k] = -v[i][k] + 2.0*vv[k];     
                  f[i][k] = -f[i][k];  
	        }
	      }  
              else{
                for (k = 0; k < 3; k++){
                  //i_v[i][k] = -i_v[i][k] + 2.0*vv[k];     
                  v[i][k] = -v[i][k] + 2.0*vv[k];     
	        }
                for (k = 0; k < 3; k++){
                  xp[k] = x_x[k] - prd[k]*shift[m][k];      
                  //xn[i][k] = xp[k] + (1.0-tt)*dtt*i_v[i][k];      // verlet 
                  xn[i][k] = xp[k] + (1.0-tt)*dtt*v[i][k];      // verlet 
	        }
                if (ind_bounce == 1)
                  for (k = 0; k < 3; k++)
                    f[i][k] = -f[i][k];
	      }
              if (mask[i] & groupbit_comm){
                m = tag[i] - min_ind;
                for (k = 0; k < 3; k++){
                  ind_x[m][k] = xn[i][k];
                  //ind_v[m][k] = i_v[i][k];         
                  ind_v[m][k] = v[i][k];         
	        }
                for(k = 0; k < ind_ntri[m]; k++)
                  ind_norm[ind_tri[m][k]] = 0;     
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
  time (&end);
  tdiff = difftime (end,start);
  //if(!comm->me)
  //  printf("time between reflection is %f\n", tdiff);
  
}

/* ---------------------------------------------------------------------- */
void FixBoundRBC::final_integrate()
{
  if (atom->nlocal > nnmax){
    nnmax = atom->nlocal;
    ind_faces = memory->grow_2d_int_array(ind_faces,nnmax,max_faces,"fix_bound_rbc:ind_faces");
    num_faces = (int *) memory->srealloc(num_faces,nnmax*sizeof(int),"fix_bound_rbc:num_faces");    
  }
}


/* ---------------------------------------------------------------------- */

void FixBoundRBC::borders()
{
  int i,m,iswap,dim,ineed,nsend,nrecv;
  double lo,hi;
  double *buf;
  MPI_Request request;
  MPI_Status status; 

  iswap = 0;
  
  for (dim = 0; dim < 3; dim++) {
    for (ineed = 0; ineed < 2*comm->need[dim]; ineed++) {

      // find all gost faces within slab boundaries lo/hi
      // store face indices in list for use in future timesteps

      lo = slablo[iswap];
      hi = slabhi[iswap];
      m = nsend = 0;
      
      for (i = 0; i < ind_tot; i++) {
	if (ind_h[i] && ind_x[i][dim] >= lo && ind_x[i][dim] < hi) {
	  if (m > maxsend) grow_send(m);
	  m += pack_border(i,&buf_send[m],comm->pbc_flag[iswap], comm->pbc[iswap]);
	  if (nsend == maxsendlist[iswap]) grow_list(iswap,nsend);
	  sendlist[iswap][nsend++] = i;
	}
      }
      // swap faces with other proc
      // put incoming ghosts at end of my faces arrays
      // if swapping with self, simply do nothing

      if (comm->sendproc[iswap] != comm->me) {
	MPI_Send(&nsend,1,MPI_INT,comm->sendproc[iswap],0,world);
	MPI_Recv(&nrecv,1,MPI_INT,comm->recvproc[iswap],0,world,&status);
	if (nrecv*size_border > maxrecv) 
	  grow_recv(nrecv*size_border);
	MPI_Irecv(buf_recv,nrecv*size_border,MPI_DOUBLE,
		  comm->recvproc[iswap],0,world,&request);
	MPI_Send(buf_send,nsend*size_border,MPI_DOUBLE,
		 comm->sendproc[iswap],0,world);
	MPI_Wait(&request,&status);
	buf = buf_recv; 
      } else {
        nsend = 0;
	nrecv = 0;
      }
      
      // unpack buffer

      m = 0;
      for (i = 0; i < nrecv; i++) m += unpack_border(&buf[m]);
     
      // set all pointers & counters

      sendnum[iswap] = nsend;
      recvnum[iswap] = nrecv;
      size_comm_send[iswap] = nsend * size_comm;
      size_comm_recv[iswap] = nrecv * size_comm;

      iswap++;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::communicate()
{
  int iswap;
  double *buf;
  MPI_Request request;
  MPI_Status status;
  int dummy;
  for (iswap = 0; iswap < comm->nswap; iswap++) {

    // pack buffer

    dummy = pack_comm(sendnum[iswap],sendlist[iswap],buf_send, comm->pbc_flag[iswap], comm->pbc[iswap]);

    // exchange with another proc
    // if self, set recv buffer to send buffer

    if (comm->sendproc[iswap] != comm->me) {
      MPI_Irecv(buf_recv,size_comm_recv[iswap],MPI_DOUBLE,
		comm->recvproc[iswap],0,world,&request);
      MPI_Send(buf_send,size_comm_send[iswap],MPI_DOUBLE,
	       comm->sendproc[iswap],0,world);
      MPI_Wait(&request,&status);
      buf = buf_recv;
    } else buf = buf_send;

    // unpack buffer

    unpack_comm(recvnum[iswap],0,buf);
  }
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::local_faces(int n)
{
  int i,j,k;
  double **x = atom->x;
  //double **i_v = atom->i_v;
  double **v = atom->v;
  int *tag = atom->tag;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  if (n){
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit_comm) {
        k = tag[i] - min_ind;
        for (j = 0; j < 3; j++){
          ind_x[k][j] = x[i][j];
          //ind_v[k][j] = i_v[i][j];         
          ind_v[k][j] = v[i][j];         
	}            
      }
  } else{
    for (i = 0; i < ind_tot; i++)
      ind_h[i] = 0;
 
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit_comm) {
        k = tag[i] - min_ind;
        ind_h[k] = 1;
        for (j = 0; j < 3; j++){
          ind_x[k][j] = x[i][j];
          //ind_v[k][j] = i_v[i][j];         
          ind_v[k][j] = v[i][j];         
	}            
      } 
  }
}


/* ---------------------------------------------------------------------- */

void FixBoundRBC::grow_send(int n)
{
  maxsend = static_cast<int> (BUFFACTOR * n);
  buf_send = (double *) 
    memory->srealloc(buf_send,(maxsend+BUFEXTRA)*sizeof(double),
		     "fix_bound_rbc:buf_send");
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::grow_list(int iswap, int n)
{
  maxsendlist[iswap] = static_cast<int> (BUFFACTOR * n);
  sendlist[iswap] = (int *) 
    memory->srealloc(sendlist[iswap],maxsendlist[iswap]*sizeof(int),
		     "fix_bound_rbc:sendlist[iswap]");
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::grow_recv(int n)
{
  maxrecv = static_cast<int> (BUFFACTOR * n);
  memory->sfree(buf_recv);
  buf_recv = (double *) memory->smalloc(maxrecv*sizeof(double),
					"fix_bound_rbc:buf_recv");
}

/* ---------------------------------------------------------------------- */

int FixBoundRBC::pack_border(int i, double *buf, int pbc_flags, int *pbc)
{
  int m = 0;

  buf[m++] = static_cast<double> (i);
  if (pbc_flags == 0) {
    buf[m++] = ind_x[i][0];
    buf[m++] = ind_x[i][1];
    buf[m++] = ind_x[i][2];
  } else {
    buf[m++] = ind_x[i][0] + pbc[0]*prd[0];
    buf[m++] = ind_x[i][1] + pbc[1]*prd[1];
    buf[m++] = ind_x[i][2] + pbc[2]*prd[2];
  }
  buf[m++] = ind_v[i][0];
  buf[m++] = ind_v[i][1];
  buf[m++] = ind_v[i][2]; 

  return m;
}

/* ---------------------------------------------------------------------- */

int FixBoundRBC::unpack_border(double *buf)
{
  int k;
  int m = 0;

  k = static_cast<int> (buf[m++]);
  ind_h[k] = 1;
  ind_x[k][0] = buf[m++];
  ind_x[k][1] = buf[m++];
  ind_x[k][2] = buf[m++];
  ind_v[k][0] = buf[m++];
  ind_v[k][1] = buf[m++];
  ind_v[k][2] = buf[m++];
  
  return m;
}

/* ---------------------------------------------------------------------- */

int FixBoundRBC::pack_comm(int n, int *list, double *buf, int pbc_flags, int *pbc)
{
  int i,j,m;

  m = 0;
  if (pbc_flags == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<double> (j);
      buf[m++] = ind_x[j][0];
      buf[m++] = ind_x[j][1];
      buf[m++] = ind_x[j][2];
      buf[m++] = ind_v[j][0];
      buf[m++] = ind_v[j][1];
      buf[m++] = ind_v[j][2];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = static_cast<double> (j);
      buf[m++] = ind_x[j][0] + pbc[0]*prd[0];
      buf[m++] = ind_x[j][1] + pbc[1]*prd[1];
      buf[m++] = ind_x[j][2] + pbc[2]*prd[2];       
      buf[m++] = ind_v[j][0];
      buf[m++] = ind_v[j][1];
      buf[m++] = ind_v[j][2];
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixBoundRBC::unpack_comm(int n, int dummy, double *buf)
{
  int i,m,k;

  m = 0;
  for (i = 0; i < n; i++) {
    k = static_cast<int> (buf[m++]);
    ind_x[k][0] = buf[m++];
    ind_x[k][1] = buf[m++];
    ind_x[k][2] = buf[m++];    
    ind_v[k][0] = buf[m++];
    ind_v[k][1] = buf[m++];
    ind_v[k][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
/*
void FixBoundRBC::face_decide()
{
  int i,j,k,l,m,n;
  double dd[3],rr;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int **ind_tmp;
    

  for (i=0; i<atom->nlocal; i++)
    if (mask[i] & groupbit){
      for (j=0; j<ntri; j++)
        ind_norm[j] = 0; 
      num_faces[i] = 0;
      for (j = 0; j < ind_tot; j++) 
        if (ind_h[j] && molecule[i] != tri_mol[ind_tri[j][0]]){
          for (k = 0; k < 3; k++)
            dd[k] = x[i][k] - ind_x[j][k];
          domain->minimum_image(dd);
          rr = dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2];
	  if (rr < d_cut_sq)
            for (n = 0; n < ind_ntri[j]; n++){
              k = ind_tri[j][n];
              if (ind_norm[k] == 0 && ind_h[tri_ind[k][0]] &&  ind_h[tri_ind[k][1]] && ind_h[tri_ind[k][2]]) {
                ind_norm[k] = 1;
                ind_faces[i][num_faces[i]] = k;
                num_faces[i]++;
	      }
              else 
	      if (ind_norm[k] == 0) printf("It might make sense to increase comm_cut, because some triangles are omitted!!! \n");
              if (num_faces[i] == max_faces) {
                m = max_faces;
                ind_tmp = memory->create_2d_int_array(nnmax,m,"fix_bound_rbc:ind_tmp");
                for (k=0; k<nnmax; k++)
                  for (l=0; l<m; l++)
                    ind_tmp[k][l] = ind_faces[k][l];   
                max_faces += 5;
                ind_faces = memory->grow_2d_int_array(ind_faces,nnmax,max_faces,"fix_bound_rbc:ind_faces");
                for (k=0; k<nnmax; k++)
                  for (l=0; l<m; l++)
                    ind_faces[k][l] = ind_tmp[k][l];
                memory->destroy_2d_int_array(ind_tmp); 
	      }     
	    }
	}
    } 
}
*/



void FixBoundRBC::face_decide()
{
  int i,j,k,l,m,n,ind;
  double dd[3],rr;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int **ind_tmp;   

  time_t t1a, t1b, t2a, t2b, t3a, t3b;
  double tdiff1, tdiff2, tdiff3;
  time(&t1a);
  
  n_tri_check = 0;
  for (i=0; i<ntri; i++)
    if (ind_h[tri_ind[i][0]] &&  ind_h[tri_ind[i][1]] && ind_h[tri_ind[i][2]]){
      tri_check[n_tri_check] = i;
      n_tri_check++;
    }    
  time(&t1b);
  tdiff1 = difftime(t1b, t1a);
  tdiff2 = 0;
  tdiff3 = 0;
   
  time(&t2a);
  for (i=0; i<atom->nlocal; i++)
    if (mask[i] & groupbit){ 
      num_faces[i] = 0;
      for (j = 0; j < n_tri_check; j++){
        n = tri_check[j]; 
        if (molecule[i] != tri_mol[n]){
          ind = 0;
          for (l = 0; l < 3; l++){  
            for (k = 0; k < 3; k++)
              dd[k] = x[i][k] - ind_x[tri_ind[n][l]][k];
            domain->minimum_image(dd);
            rr = dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2];
            if (rr < d_cut_sq){
              ind = 1;
              break;      
            }
          }  
	  if (ind){
            ind_faces[i][num_faces[i]] = n;
            num_faces[i]++;
	    if (num_faces[i] == max_faces) {
              m = max_faces;
              ind_tmp = memory->create_2d_int_array(nnmax,m,"fix_bound_rbc:ind_tmp");
              for (k=0; k<nnmax; k++)
                for (l=0; l<m; l++)
                  ind_tmp[k][l] = ind_faces[k][l];   
              max_faces += 5;
              ind_faces = memory->grow_2d_int_array(ind_faces,nnmax,max_faces,"fix_bound_rbc:ind_faces");
              for (k=0; k<nnmax; k++)
                for (l=0; l<m; l++)
                  ind_faces[k][l] = ind_tmp[k][l];
              memory->destroy_2d_int_array(ind_tmp); 
	    }     
	  }
        }
      }
    }
  time(&t2b);
  tdiff2 = difftime(t2b,t2a);
  if(tdiff2 > 10.0)
    printf("bb:step is %d time between facedecide three stage is %f %f, total is %f, cpu is %d\n", update->ntimestep, tdiff1, tdiff2, tdiff1 + tdiff2,  comm->me);
}
/*
void FixBoundRBC::face_decide()
{
  int i,j,k,l,m,n,jj;
  double dd[3],rr;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int **ind_tmp;

  time_t t1a, t1b, t2a, t2b, t3a, t3b;
  double tdiff1, tdiff2, tdiff3;
  time(&t1a);
  n_point_check = 0;
  for (i=0; i<ind_tot; i++)
    if (ind_h[i]){
      point_check[n_point_check] = i;
      n_point_check++;
    }
  time(&t1b);
  tdiff1 = difftime(t1b, t1a);
  tdiff2 = 0;
  tdiff3 = 0;
  for (i=0; i<atom->nlocal; i++)
    if (mask[i] & groupbit){
      time(&t2a);
      for (j=0; j<ntri; j++)
        ind_norm[j] = 0; 
      time(&t2b);
      tdiff2 += difftime(t2b, t2a);
      time(&t3a);
      num_faces[i] = 0;
      for (jj = 0; jj < n_point_check; jj++){
        j = point_check[jj];
        if (molecule[i] != tri_mol[ind_tri[j][0]]){
          for (k = 0; k < 3; k++)
            dd[k] = x[i][k] - ind_x[j][k];
          domain->minimum_image(dd);
          rr = dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2];
	  if (rr < d_cut_sq)
            for (n = 0; n < ind_ntri[j]; n++){
              k = ind_tri[j][n];
              if (ind_norm[k] == 0 && ind_h[tri_ind[k][0]] &&  ind_h[tri_ind[k][1]] && ind_h[tri_ind[k][2]]) {
                ind_norm[k] = 1;
                ind_faces[i][num_faces[i]] = k;
                num_faces[i]++;
	      }
              else 
	      if (ind_norm[k] == 0) printf("It might make sense to increase comm_cut, because some triangles are omitted!!! \n");
              if (num_faces[i] == max_faces) {
                m = max_faces;
                ind_tmp = memory->create_2d_int_array(nnmax,m,"fix_bound_rbc:ind_tmp");
                for (k=0; k<nnmax; k++)
                  for (l=0; l<m; l++)
                    ind_tmp[k][l] = ind_faces[k][l];   
                max_faces += 5;
                ind_faces = memory->grow_2d_int_array(ind_faces,nnmax,max_faces,"fix_bound_rbc:ind_faces");
                for (k=0; k<nnmax; k++)
                  for (l=0; l<m; l++)
                    ind_faces[k][l] = ind_tmp[k][l];
                memory->destroy_2d_int_array(ind_tmp); 
	      }     
	    }
        }
      }
      time(&t3b);
      tdiff3 += difftime(t3b, t3a);
    }
    
  printf("time between facedecide three stage is %f %f %f, total is %f, cpu is %d\n", tdiff1, tdiff2, tdiff3, tdiff1 + tdiff2 + tdiff3,  comm->me);

}
*/

/* ---------------------------------------------------------------------- */

void FixBoundRBC::calc_norm(int m)
{
  int k;

  for (k = 0; k < 3; k++){
    edge1[m][k] = ind_x[tri_ind[m][1]][k] - ind_x[tri_ind[m][0]][k];
    edge2[m][k] = ind_x[tri_ind[m][2]][k] - ind_x[tri_ind[m][0]][k];
    dif1[m][k] = ind_v[tri_ind[m][1]][k] - ind_v[tri_ind[m][0]][k];
    dif2[m][k] = ind_v[tri_ind[m][2]][k] - ind_v[tri_ind[m][0]][k]; 
  }
  domain->minimum_image(edge1[m]);
  domain->minimum_image(edge2[m]);

  norm[m][0] = edge1[m][1]*edge2[m][2] - edge1[m][2]*edge2[m][1];
  norm[m][1] = edge1[m][2]*edge2[m][0] - edge1[m][0]*edge2[m][2];
  norm[m][2] = edge1[m][0]*edge2[m][1] - edge1[m][1]*edge2[m][0];
  norm[m][3] = - norm[m][0]*ind_x[tri_ind[m][0]][0] - norm[m][1]*ind_x[tri_ind[m][0]][1] - norm[m][2]*ind_x[tri_ind[m][0]][2];

  dd12[m][0] = dif1[m][1]*dif2[m][2] - dif1[m][2]*dif2[m][1];
  dd12[m][1] = dif1[m][2]*dif2[m][0] - dif1[m][0]*dif2[m][2];
  dd12[m][2] = dif1[m][0]*dif2[m][1] - dif1[m][1]*dif2[m][0];    

  ind_norm[m] = 1;
}

/* ------------------------------------------------------------------------- */

double FixBoundRBC::find_root(double *cc)
{
  double dot,tt,t1,t2;
  double cp[3];

  cp[0] = 3.0*cc[0];
  cp[1] = 2.0*cc[1];
  cp[2] = cc[2];  

  dot = cp[1]*cp[1] - 4.0*cp[0]*cp[2];
  if (dot < 0){
    tt = newton(0.0,cc,cp);
  }
  else {
    t1 = 0.5*(-cp[1]+sqrt(dot))/cp[0];
    t2 = 0.5*(-cp[1]-sqrt(dot))/cp[0];
    if (t2 < t1) {
      tt = t2;
      t2 = t1;
      t1 = tt;
    }
    if (t2<0.0 || t1>1.0 || (t1<0.0 && t2>1.0)){
      tt = newton(0.0,cc,cp); 
    }
    else{
      if (t1>=0.0 && t2<=1.0){
        tt = newton(0.0,cc,cp);
        if (tt<0.0 && tt>1.0)
          tt = newton(1.0,cc,cp);
        if (tt<0.0 && tt>1.0)
          tt = newton(0.5*(t1+t2),cc,cp);  
      }else{
        tt = newton(0.0,cc,cp);
        if (tt<0.0 && tt>1.0)
          tt = newton(1.0,cc,cp); 
      }  
    }
  }

  return tt;
}

/* ------------------------------------------------------------------------- */

double FixBoundRBC::newton(double tt, double *cc, double *cp)
{
  int nn;
  double tt1,val,valp;
  int it_max = 50;

  nn = 0;
  val = cc[0]*tt*tt*tt + cc[1]*tt*tt + cc[2]*tt + cc[3];
  while (fabs(val)>EPS && nn<it_max){
    valp = cp[0]*tt*tt + cp[1]*tt + cp[2];
    tt1 = tt - val/valp;
    nn++;
    tt = tt1;
    val = cc[0]*tt*tt*tt + cc[1]*tt*tt + cc[2]*tt + cc[3];
  }

  return tt;
}

/* ------------------------------------------------------------------------- */

double FixBoundRBC::solve_quadratic(double *cc)
{
  double tt,dot,t1,t2,uv;

  tt = -1;
  dot = cc[2]*cc[2] - 4.0*cc[1]*cc[3];
  if (dot >= 0){
    t1 = 0.5*(-cc[2]+sqrt(dot))/cc[1];
    t2 = 0.5*(-cc[2]-sqrt(dot))/cc[1];
    if (t2 < t1) {
      uv = t2;
      t2 = t1;
      t1 = uv;
    }
    if (t2>=0.0 && t2<=1.0)
      tt = t2;
    if (t1>=0.0 && t1<=1.0)
      tt = t1;
  }
  
  return tt;
}

/* ------------------------------------------------------------------------- */
