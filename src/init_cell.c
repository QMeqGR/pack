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

extern int restart,forces;
extern int debug;
extern int numat;
extern int non_per_ewald;

extern double cvd_max;
extern double wl_cell_init_dx;

extern struct cellprm cell;
extern struct obj *object;

/*****************************************/
/*                                       */
/*    Initialize the cell structure      */
/*                                       */
/*****************************************/
void cell_init(void)
{

  int i,j;
  int viol_flag=0,place_object_trial_number=0;
  int pathological_flag=0;
  int make_obj=0;
  static int been_called=0;

  double temp=0,phi=0,the=0,psi=0;
  double max_cvd = cvd_max;
  static double orig_a=0,orig_b=0,orig_c=0;
  static double orig_alph=0,orig_beta=0,orig_gamm=0;
  static double max_R=0;

  struct vector tc;
  struct matrix L,Linv;

  if ( ! been_called ){

    orig_a = cell.a;
    orig_b = cell.b;
    orig_c = cell.c;
    orig_alph = cell.alph;
    orig_beta = cell.beta;
    orig_gamm = cell.gamm;

    for (i=0;i<MAX_OBJS && object[i].used==1 ;i++){
      for (j=0;j<MAX_VERTS && object[i].vused[j]==1 ;j++){
	if ( object[i].vfxd[j]==1 ) continue; /* object could be a surface */
	max_R   = fmax( object[i].R[j] , max_R );
      }
    }

    /* set constraint flag if there are constraints */
    for (i=0;i<MAX_OBJS && object[i].used==1 ;i++){
      object[i].vcnstr[0].constrained=0;
      if ( object[i].vcnstr[0].xmin > CONSTR_XMIN ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].xmax < CONSTR_XMAX ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].ymin > CONSTR_YMIN ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].ymax < CONSTR_YMAX ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].zmin > CONSTR_ZMIN ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].zmax < CONSTR_ZMAX ) object[i].vcnstr[0].constrained=1;
      if ( object[i].vcnstr[0].constrained==1 &&
	   debug > 1 ) {
	printf("* Found constrained object: i=%d\n",i);
	printf("* xmin = %10.3f   xmax = %10.3f\n",object[i].vcnstr[0].xmin,object[i].vcnstr[0].xmax);
	printf("* ymin = %10.3f   ymax = %10.3f\n",object[i].vcnstr[0].ymin,object[i].vcnstr[0].ymax);
	printf("* zmin = %10.3f   zmax = %10.3f\n",object[i].vcnstr[0].zmin,object[i].vcnstr[0].zmax);
      }
    }

    if ( debug > 2 ){

      /* print out constrained object table */
      printf("* Constrained objects:\n");
      for (i=0;i<MAX_OBJS && object[i].used==1 ;i++){
	printf("* object %3d constrained= %d\n",i,object[i].vcnstr[0].constrained);
      }

      printf("* cell init: orig_ a b c = %15.10f%15.10f%15.10f\n",orig_a,orig_b,orig_c);
      printf("* cell init: orig_ angle = %15.10f%15.10f%15.10f\n",orig_alph,orig_beta,orig_gamm);
      printf("* cell init: max_R = %15.10f   max_cvd = %15.10f\n",max_R,max_cvd);
    }
    
    been_called = 1;
  }
  


  /* initialize the cell structure */
  cell.a = orig_a;
  cell.b = orig_b;
  cell.c = orig_c;
  cell.alph = orig_alph;
  cell.beta = orig_beta;
  cell.gamm = orig_gamm;
  pathological_flag = calc_cell_basis(&cell);
  if ( pathological_flag ){
    printf("* ERROR! cell basis problems! pathological_flag= %d\n",
	   pathological_flag);
  }
  L = L_e3(&cell);
  

  /*       Initialize Random Locations        */
  /*
    First, we randomly pick the locations of the centers
    of the objects to be placed in the unit cell.  Note that
    the centers are already going to be in fractional coordinates.
    The functions make_object, will put the vertices
    in fractional coordinates for the unit cell.
  */
  if ( !restart ) {

    /* initialize the object locations */
    for( i=0; i<MAX_OBJS && object[i].used==1 ;i++ ){    
      place_object_trial_number=0;      

      /* fixed or constrained?, just make it where it is */
      if ( object[i].vfxd[0]==1 || 
	   object[i].vcnstr[0].constrained==1 ){ 
	if ( debug >1 ) printf("* Found fixed or constrained object: i=%d Placing with input coords.\n",i);
	tc = object[i].v[0]; /* need to put this in frac coords */
	L = L_e3( &cell );
	Linv = Linverse( &L );
	tc = UTMvmult( &Linv , &tc );
	phi=0,the=0,psi=0;
	make_obj = make_object(i,tc,phi,the,psi);
	if ( make_obj == 0 ){
	  printf("* Fatal error making object\n");
	  exit_pack(0);
	}
	continue;
      }

      if ( i==0 ) tc = get_rand_pos( non_per_ewald );

      /* make sure the objects start out at least 2*cvd apart */
      if ( i>0 ) {
	
	do {
	  place_object_trial_number++;
	  tc = get_rand_pos( non_per_ewald );
	  
	  for(j=0,viol_flag=0; j<i; j++) {
	
	    if ( non_per_ewald ) temp = dist_NOpbc(&L,&object[j].v[0],&tc);
	    else temp = dist_pbc(&L,&object[j].v[0],&tc);
	    
	    if ( debug>1 ) {
	      printf("* cell init: obj i=%3d  j=%3d  [%8.4f%8.4f%8.4f]\n",
		     i,j,tc.x,tc.y,tc.z);
	      printf("* cell init: obj                 temp=%8.4f   2*max_cvd=%8.4f   trial= %d\n",
		     temp,2*max_cvd,place_object_trial_number);
	    }
	    if ( place_object_trial_number > CELL_INIT_MAX_TRY ){
	      printf("* ERROR: initializing object positions i= %d  j= %d\n",i,j);
	      printf("* ERROR: can not init cell. Try increasing cell size or wl_cell_init_dx?\n");
	      printf("* place_object_trial_number = %d\n",place_object_trial_number);
	      exit(0);
	    }
	    if ( temp<(2*max_cvd) ) { viol_flag=1; break; }

	  } /* j loop */

	} while ( viol_flag==1 );
	
      }
      
      phi = 2*PI * drand48();
      the =   PI * drand48();
      psi = 2*PI * drand48();
      if ( debug>1 ) printf("* cell init: obj i=%3d  %8.4f%8.4f%8.4f\n",i,tc.x,tc.y,tc.z);
      make_obj = make_object(i,tc,phi,the,psi);
      if ( make_obj == 0 ){
	printf("Fatal error making object\n");
	exit_pack(0);
      }
      
    }
    
  } /* matches !restart */
  
  return;
}


struct vector get_rand_pos(int cluster)
{
  struct vector randpos;

  double ddx,ddy,ddz,DD,shift;

  if ( cluster==1 ) DD = wl_cell_init_dx;
  else DD = 1;

  ddx = DD*drand48();
  ddy = DD*drand48();
  ddz = DD*drand48();

  if ( cluster==1 ){
    shift=0.5-DD/2;
    randpos = vadd( makevec( shift , shift , shift ), makevec( ddx ,ddy , ddz ) );
  }
  else randpos = makevec( ddx , ddy , ddz );

  return randpos;
}
