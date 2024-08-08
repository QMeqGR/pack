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
extern int efunc;
extern int numat;
extern double rep_epsilon;
extern double *etable;
extern double LJeps[MAX_TABLE_SQR];
extern double LJrer[MAX_TABLE_SQR];
extern double *SShrd,*SShrd_2,*SShrd_6,*SShrd12;
extern struct cellprm cell;
extern struct obj *object;
extern struct atom *s;

/************************************************************/
/*                                                          */
/*                 Soft Sphere and LJ                       */
/*             Repulsive Energy Functional                  */
/*                                                          */
/************************************************************/
double Energy_rep(double orth)
{

  int i,j,k,l;
  int idx,n1,n2,o1,o2;
  int Z1,Z2;
  int same_atom=0;

  static int called_flag=0;

  /* kittel (R/sigma) = 2^(1/6) = 1.12, recip= 0.8908...) */
  static double SGMA = 0.8908987184;

  static double epsilon;
  static double min_val,OVLP;
  double energy=0;
  static double temp=0,distsq=0;
  static double D6,D12;
  static double LJ_D6,LJ_D12;
  static double ljeps,ljrer,ljrersq=0,ljtemp;

  struct vector v1,v2;

  struct matrix L;


  if ( !called_flag ) {

    if ( efunc == 0 ) epsilon = rep_epsilon;
    if ( efunc == 1 ) epsilon = 4*rep_epsilon;

    /* set LJrer array by hand if not set in the input file */
    /* This assumes energ init was called */
    for(i=0;i<numat;i++){
      for(j=0;j<numat;j++){
	if ( tol( RER_DEFAULT , LJrer[ s[i].Z * MAX_TABLE + s[j].Z ] , 1e-5 ) )
	  LJrer[ s[i].Z * MAX_TABLE + s[j].Z ] = SGMA*( s[i].R + s[j].R );
      }
    }

    get_hrd_dist( );

    /* If rep_epsilon was set to reduce the LJ energies overall, then reset them here */
      for(i=0;i<MAX_TABLE;i++){
	for(j=0;j<MAX_TABLE;j++){
	  /* !!!!! note rep_epsilon should be default 1.0 for true LJ !!!!! */
	  LJeps[i*MAX_TABLE+j] = LJeps[i*MAX_TABLE+j] * rep_epsilon;
	}
      }

    /* Print lennard-jones interaction energy table */
    if ( debug > 4 ) {
      printf("* Lennard-Jones eps table data (if different from %f)\n",EPS_DEFAULT);
      for(i=0;i<MAX_TABLE;i++){
	for(j=0;j<MAX_TABLE;j++){
	  ljeps = LJeps[i*MAX_TABLE+j];
	  if ( !tol( EPS_DEFAULT , ljeps , 1e-6) ) printf("* en_rep eps: [%3d] [%3d] : %5.3f\n",i,j,ljeps);
	}
      }
      printf("* Lennard-Jones rer table data (if different from %f)\n",RER_DEFAULT);
      for(i=0;i<MAX_TABLE;i++){
	for(j=0;j<=i;j++){
	  ljrer = LJrer[i*MAX_TABLE+j];
	  if ( !tol( RER_DEFAULT , ljrer , 1e-6) ) printf("* en_rep rer: [%3d] [%3d] : %5.3f (* SGMA in table)\n",i,j,ljrer/SGMA);
	}
      }
    }

    called_flag = 1;
  }

  /*************************************************************/
  /*************************************************************/

  /* calculate the L matrix for conversion to cartesian coords */
  L = L_e3(&cell);

  if ( debug > 5 ) {
    printf("* SS energy run:\n*\n");
    printf("*%18s%9s%8s%8s%14s%8s%8s%4s%4s%4s%4s%4s%4s%10s%10s\n",
	   "d_ij","x1","y1","z1","x2","y2","z2","i","k","j","l","at1","at2","OVLP","U_ij");
  }

  /*------------------------------------------------------------------------------------
    OBJ-OBJ interaction 
    ------------------------------------------------------------------------------------*/

  for (i=0, n1=0; i<MAX_OBJS && object[i].used==1; i++) {
    for (j=0; j<MAX_VERTS && object[i].vused[j]==1; j++) {
      n1++;
      o1 = n1-1;
      
      for (k=0, n2=0; k<MAX_OBJS && object[k].used==1; k++) {
	for (l=0; l<MAX_VERTS && object[k].vused[l]==1; l++) {
	  n2++;
	  o2 = n2-1;

	  idx = o1 * numat + o2;
	  temp = 0;

	  if ( o2 > o1 ) continue;
	  if ( i==k && j==l ) same_atom=1;
	  else same_atom=0;

	  if ( s[o1].moved==0 && s[o2].moved==0 ) continue;

	  if ( efunc == 1 ) { /* efunc=1 is L-J */
	    Z1 = s[o1].Z;
	    Z2 = s[o2].Z;
	    ljeps = LJeps[ Z1 * MAX_TABLE + Z2 ];
	    ljrer = LJrer[ Z1 * MAX_TABLE + Z2 ];
	    epsilon = ljeps;
	    ljrersq = ljrer*ljrer;
	  }

	  v1 = makevec( s[o1].x, s[o1].y, s[o1].z );
	  v2 = makevec( s[o2].x, s[o2].y, s[o2].z );

	  distsq = distsq_ssrep(&L,&v1,&v2,&cell,orth,same_atom);
	  
	  if ( efunc == 0 ) {
	    D6  = 1.0 / (distsq*distsq*distsq);
	    D12 = D6*D6;
	  } else if ( efunc == 1 ) {
	    ljtemp = ljrersq / distsq;
	    LJ_D6  = ljtemp*ljtemp*ljtemp;
	    LJ_D12 = LJ_D6*LJ_D6; 
	  }

	  if ( efunc == 0 ){

	    min_val = SShrd_2[ idx ];
	    if ( (distsq+TINY)<min_val ) temp = epsilon * ( SShrd12[ idx ] * D12 );
	    
	  } else if ( efunc == 1 ) {
	    
	    temp = epsilon * ( LJ_D12 - LJ_D6 );
	    
	  }
	  
	  etable[ idx ] = temp;
	  
	  if ( (debug>5 && temp>0) || debug>6 ) {
	    printf("* ob-ob : (%7.4g) ", sqrt(distsq) );
	    printf("%8.3f%8.3f%8.3f%14.3f%8.3f%8.3f",
		   v1.x,v1.y,v1.z,
		   v2.x,v2.y,v2.z);
	    OVLP = sqrt(distsq) - sqrt(min_val);
	    if ( OVLP < 0.0 ) printf("%4d%4d%4d%4d%4d%4d%10.2e%10.2e\n",
				     i,k,j,l,o1,o2,fabs(OVLP),temp);
	    else printf("%4d%4d%4d%4d%4d%4d%10s%10.2e\n",
			i,k,j,l,o1,o2,"---",temp);

	  }

	}
      }
    }
  }


  for(i=0;i<numat;i++){
    for(j=0;j<=i;j++){
      energy += etable[ numat * i + j ];
    }
  }

  /* debug printing of etable */
  if ( debug>6 ){
    printf("etable: e= %f\n",energy);
    for(i=0;i<numat;i++){
      printf("%3d: ",i);
      for(j=0;j<=i;j++){
	printf("%9.5f ",etable[ numat * i + j ]);
	if ( j==i ) printf("\n");
      }
    }
    printf("\n");
  }
  
  return(energy);
}

