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
extern double rep_epsilon;

extern struct cellprm cell;
extern struct atom *p;
extern struct atom *s;
extern struct obj *object;


struct vector geometric_center(void)
{
  int i,j;
  struct vector gc={0};

  for(i=0;i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      gc = vsmult( 0.5 , vadd( gc, object[i].v[j] ) );
    }
  }

  return gc;
}


struct vector center_of_mass(void)
{
  int i,j,nv;
  static int been_called=0;
  static double total_mass=0;
  struct vector cm={0};

  if ( been_called == 0 ) {
    
    for(i=0;i<MAX_OBJS && object[i].used==1; i++){
      for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
	total_mass += molwt( object[i].Z[j] );
      }
    }
    
    if ( debug > 4 ) printf("Total mass of all atoms: %10.4f [g/mol]\n",total_mass);
    been_called = 1;
  }


  for(i=0;i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      cm = vadd( cm , vsmult( molwt(object[i].Z[j]) , object[i].v[j] ) );
    }
  }

  cm = vsmult( (double)1/total_mass , cm );

  return cm;
}
