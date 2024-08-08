/*
* PACK - A Monte Carlo structure searching program
* Copyright (C) 2008 Eric H. Majzoub and Vidvuds Ozolins
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "global_defs.h"
#include "packlib.h"

extern int debug;
extern int numat;

extern double LJrer[MAX_TABLE_SQR];
extern double *SShrd,*SSmin;
extern double *SShrd_2,*SShrd_6,*SShrd12;
extern struct obj *object;
extern struct atom *s;


/****************************************/
/*                                      */
/*         Set Hard Distances           */
/*                                      */
/****************************************/
void get_hrd_dist(void)
{
  int i,j,k,l;
  int at1=0,at2=0,idx=0;

  double dist,t;

  /* kittel (R/sigma) = 2^(1/6) = 1.12, recip= 0.8908...) */
  static double SGMA = 0.8908987184;


  if ( debug>5 ){
    printf("* %5s%5s%5s%5s%5s%7s%7s%5s%5s%10s%10s\n",
	   "idx","i","j","k","l","Rij","Rkl","at1","at2","SShrd","SShrd_2");
  }
  /* set the hard distances for SS repulsion */
  for(i=0, at1=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      at1++;
      for(k=0, at2=0; k<MAX_OBJS && object[k].used==1; k++){
	for(l=0; l<MAX_VERTS && object[k].vused[l]==1; l++){
	  at2++;

	  idx = (at1-1) * numat + (at2-1);

	  dist = LJrer[ object[i].Z[j] * MAX_TABLE + object[k].Z[l] ];
	  if ( !tol( dist, RER_DEFAULT, 1e-6 ) ){
	    t=SShrd[ idx ] = dist / SGMA;
	  } else {
	    t=SShrd[ idx ] = object[i].R[j]+object[k].R[l];
	  }
	  SShrd_2[ idx ] = pow(t,2);
	  SShrd_6[ idx ] = pow(t*SGMA,6);
	  SShrd12[ idx ] = pow(t*SGMA,12);

	  if ( debug>5 ){
	    printf("* %5d%5d%5d%5d%5d%7.3f%7.3f%5d%5d%10.5f%10.5f\n",
		   idx,i,j,k,l,object[i].R[j],object[k].R[l],at1-1,at2-1,SShrd[ idx ],SShrd_2[ idx ]);
	  }

	}
      }
    }
  }

  
  return;
}

/****************************************/
/*                                      */
/*      Get Minimum Distances           */
/*                                      */
/****************************************/
void get_min_distances(void)
{

  int i,j,k,l;
  int at1=0,at2=0,idx=0;
  extern int debug;

  double d_ij=0,d_hd=0;

  struct matrix L;

  extern struct cellprm cell;

  get_hrd_dist( );

  L = L_e3(&cell);


  if ( debug > 2 ){
    printf("*\n");
    printf("* Distances (PBC only):\n");
    printf("* %8s%8s%15s%15s%15s\n","at1","at2","min","hard","diff");
  }

  /* find minimum distances */

  for(i=0, at1=0; i<MAX_OBJS && object[i].used==1 ; i++){
    for(j=1; j<MAX_VERTS && object[i].vused[j]==1 ; j++){
      at1++;

      for(k=0, at2=0; k<MAX_OBJS && object[k].used==1 ; k++){
	for(l=1; l<MAX_VERTS && object[k].vused[l]==1 ; l++){
	  at2++;

	  idx = (at1-1) * numat + (at2-1);

	  d_ij = dist_pbc(&L,&object[i].v[j],&object[k].v[l]);
	  SSmin[ idx ] = d_ij;
	  d_hd = SShrd[ idx ];
	  if ( debug > 2 ){
	    printf("* %8d%8d%15.6e%15.6e%15.6e\n",
		   at1-1,at2-1,d_ij,d_hd,d_hd-d_ij);
	  }
	}
      }
    }
  }

  if ( debug > 2 ){
    printf("*\n");
  }

  return; 
}
