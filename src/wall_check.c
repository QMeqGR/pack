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


/* This function checks if there is a wall violation
   within some buffer zone.  See also the function trans_accept */

int check_wall_violation(void)
{
  int i;
  int been_here=0;

  static double np_wall_check_hi=0;
  static double np_wall_check_lo=0;

  if ( been_here == 0 ){

    np_wall_check_lo = ( VIOLAT_NPE_DPLANE_WALL_CUSHION * dplane ) / cell.a;
    np_wall_check_hi = 1.0 - np_wall_check_lo;


    if ( debug > 2 ) printf("* check_wall_violation: np_wall_check_lo= %.4f      np_wall_check_hi= %.4f\n",
			    np_wall_check_lo,np_wall_check_hi);
    
    if ( np_wall_check_hi > 1.0 || np_wall_check_hi < 0.0 || 
	 np_wall_check_lo > 1.0 || np_wall_check_lo < 0.0 || 
	 np_wall_check_lo > np_wall_check_hi ){
      printf("* ERROR: wall_check init: np_wall_check_hi= %f   np_wall_check_lo= %f\n",
	     np_wall_check_hi, np_wall_check_lo);
      exit_pack(0);
    }
    
    been_here=1;
  }

  for(i=0;i<num_ob;i++) {

    if ( object[i].v[0].x > np_wall_check_hi ) { wall_check_error(10,i,object[i].v[0].x); return(10); }
    if ( object[i].v[0].x < np_wall_check_lo ) { wall_check_error(10,i,object[i].v[0].x); return(10); }

    if ( object[i].v[0].y > np_wall_check_hi ) { wall_check_error(10,i,object[i].v[0].y); return(10); }
    if ( object[i].v[0].y < np_wall_check_lo ) { wall_check_error(10,i,object[i].v[0].y); return(10); }

    if ( object[i].v[0].z > np_wall_check_hi ) { wall_check_error(10,i,object[i].v[0].z); return(10); }
    if ( object[i].v[0].z < np_wall_check_lo ) { wall_check_error(10,i,object[i].v[0].z); return(10); }

  }

  return(0);
}

void wall_check_error(int i,int j, double d)
{

  if ( i==1  ) printf("* ERROR! wall_check_error: cation number %d coordinate %f\n",j,d);
  if ( i==10 ) printf("* ERROR! wall_check_error: anion  number %d coordinate %f\n",j,d);

}
