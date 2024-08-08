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

extern double dplane;

extern struct obj *object;
extern struct cellprm cell;


/**************************************/
/*                                    */
/* translate accept                   */
/*                                    */
/**************************************/
int trans_accept(int obnum, struct vector v, int non_per_ewald)
{

  static int been_here=0;
  static int cnstr_viol_count=0;
  static double wall_bot=0;
  static double wall_top=0;
  static double ACCEPT_CUSHION=ACCEPT_NPE_DPLANE_WALL_CUSHION;

  if ( been_here==0 ){

    wall_bot = ( ACCEPT_CUSHION * dplane ) / cell.a;
    wall_top = 1.0 - wall_bot;

    if ( debug > 2 ){
      printf("* trans_accept:  CUSHION wall_bot= %.4f      wall_top= %.4f\n",wall_bot,wall_top);
    }

    if ( wall_bot < 0.0 || wall_bot > 1.0 ||
	 wall_top < 0.0 || wall_top > 1.0 ||
	 wall_top <= wall_bot ){
      printf("* ERROR: wall_bot= %f  wall_top= %f\n",wall_bot,wall_top);
      exit_pack(0);
    }

    been_here=1;
  }

  /* check for constraint violations and reject the move if any are violated */
  if ( obnum != -1 && object[obnum].vcnstr[0].constrained==1 ){
    if ( check_constr_violation( obnum , v ) ){
      cnstr_viol_count++;
      if ( debug > 3 ) printf("* Constraint violation ob %2d. Returning 0 to trans_accept. Violation count= %d !!! Swapping with unconstrained objects will hang here!!!\n",obnum,cnstr_viol_count);
      if ( cnstr_viol_count > 1e8 ){
	printf("* WARN: constraint violation count is large. Are you swapping with unconstrained objects?\n");
	cnstr_viol_count=0;
      }
      return(0);
    }
  }

  if ( non_per_ewald ) {
    if ( v.x < wall_bot ) return(0);
    if ( v.x > wall_top ) return(0);
    if ( v.y < wall_bot ) return(0);
    if ( v.y > wall_top ) return(0);
    if ( v.z < wall_bot ) return(0);
    if ( v.z > wall_top ) return(0);
  }
  
  return(1);

}

/**************************************/
/*                                    */
/* translate all objects              */
/*                                    */
/**************************************/
void trans_all(struct vector a)
{
  int i;
  struct matrix L;

  a = vsmult( -1.0 , a );

  L = L_e3(&cell);
  for(i=0;i<MAX_OBJS && object[i].used==1;i++) trans_object( i , &L , a );

  return;
}

/**************************************/
/*                                    */
/* translate an object                */
/*                                    */
/**************************************/
void trans_object(int m,
		  struct matrix *L,
		  struct vector a)
{
  int debug=0;
  int i;
  int nv;
  struct vector ctrans;
  struct vector center_old={0},center_new={0};

  nv = object[m].num_vert;

  for(i=0;i<MAX_VERTS && object[m].vused[i]==1; i++){
    object[m].vmoved[i]=1;
  }

  center_old = object[m].v[0];
  center_new = object[m].v[0] = vadd( center_old, a );

  /* calculate the cartesian trans vector for the vertices */
  ctrans = Ctrans( L, L, &center_old, &center_new );

  if ( debug ){
    for (i=1; i<nv; i++) {
      printf("trans_ob: before [vertex %d]: %15.10f%15.10f%15.10f\n",
	     i,
	     object[m].v[i].x,
	     object[m].v[i].y,
	     object[m].v[i].z);
    }
    printf("trans_objec: center old: %15.10f%15.10f%15.10f\n",
	   center_old.x,
	   center_old.y,
	   center_old.z);
    printf("trans_objec: center new: %15.10f%15.10f%15.10f\n",
	   center_new.x,
	   center_new.y,
	   center_new.z);
    printf("trans_objec: ctrans vec: %15.10f%15.10f%15.10f (cartesian)\n",
	   ctrans.x,
	   ctrans.y,
	   ctrans.z);
  }

  
  for (i=1; i<nv; i++) object[m].v[i] = vertex_new( L, L, &ctrans, &object[m].v[i] );

  if ( debug ){
    for (i=1; i<nv; i++) {
      printf("trans_anion: after [vertex %d]: %15.10f%15.10f%15.10f\n",i,
	     object[m].v[i].x,
	     object[m].v[i].y,
	     object[m].v[i].z);
    }
  }

  return;
}

void mark_all_objects_moved(void)
{

  int i,j;

  /* set 'moved' flags to one */
  for(i=0;i<MAX_OBJS && object[i].used==1;i++){
    for(j=0;j<MAX_VERTS && object[i].vused[j]==1;j++){
      object[i].vmoved[j]=1;
    }
  }
  
  return;
  
}
