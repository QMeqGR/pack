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
extern int efunc,num_ob;
extern double rep_epsilon;
extern double dplane;

extern struct cellprm cell;
extern struct atom *p;
extern struct atom *s;
extern struct obj *object;


void check_constraints(void)
{
  int chk_cnstr_err= -1;

  chk_cnstr_err=check_all_object_constraints();
  if ( chk_cnstr_err > -1 ){
    printf("* Error: Input constraints are violated. Check input file coords for object %d. Exiting.\n",
	   chk_cnstr_err);
    exit_pack(0);
  }

}

int check_all_object_constraints(void)
{
  int i;
  int cnstr_err_flag= -1;
  struct vector v;

  for (i=0; i<num_ob; i++){
    if ( object[i].vcnstr[0].constrained==1 ){
      printf("in check all object constraints\n");
      v = object[i].v[0];

      printf("ob %2d %10.3f%10.3f%10.3f    %10.3f%10.3f%10.3f\n",
	   i,
	   v.x,v.y,v.z,
	   object[i].v[0].x,object[i].v[0].y,object[i].v[0].z);
      if ( check_constr_violation( i , v ) ) return(i);
    }
  }
  return(-1);
}


/* This function checks if there is a constraint violation
   within some buffer zone.  See also the function trans_accept */

int check_constr_violation(int obnum, struct vector v)
{

  int debug=0;

  if ( debug )
    printf("ob %2d %10.3f%10.3f%10.3f    %10.3f%10.3f%10.3f\n",
	   obnum,
	   v.x,v.y,v.z,
	   object[obnum].v[0].x,object[obnum].v[0].y,object[obnum].v[0].z);
  
  /* Currently not checking every vertex. Checking only the vertex center for speed. */
  if ( v.z < object[obnum].vcnstr[0].zmin ) return(1);
  if ( v.z > object[obnum].vcnstr[0].zmax ) return(1);
  
  if ( v.y < object[obnum].vcnstr[0].ymin ) return(1);
  if ( v.y > object[obnum].vcnstr[0].ymax ) return(1);
  
  if ( v.x < object[obnum].vcnstr[0].xmin ) return(1);
  if ( v.x > object[obnum].vcnstr[0].xmax ) return(1);

  return(0);
}

