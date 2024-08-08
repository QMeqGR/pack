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

extern struct obj *object;

/**************************************/
/*                                    */
/* Check that atoms are within bounds */
/*                                    */
/**************************************/

void check_objects(void)
{
  int i,j;
  double x=1e10,y=1e10,z=1e10;

  for (i=0; i<MAX_OBJS && object[i].used==1; i++){
    for (j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      x = (*(object+i)).v[j].x; y = (*(object+i)).v[j].y; z = (*(object+i)).v[j].z;
      if ( j==0 ){ /* only checks the center?? */
	if ( x > 1.0 || x < 0.0 ) printf("* WARN: obj [%d] v[%d]: out of bounds X %15.9f\n",i,j,x);
	if ( y > 1.0 || y < 0.0 ) printf("* WARN: obj [%d] v[%d]: out of bounds Y %15.9f\n",i,j,y);
	if ( z > 1.0 || z < 0.0 ) printf("* WARN: obj [%d] v[%d]: out of bounds Z %15.9f\n",i,j,z);
      }
    }
  }
  return;
}
