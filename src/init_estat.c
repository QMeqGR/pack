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

extern int debug,numat;

extern struct atom *p;
extern struct obj *object;
extern struct cellprm cell;

/* 
 * This function sets up the p structure which contains all the atoms and
 * their charges. Eion needs the atoms in CARTESIAN coordinates!!!!!!
 */
void estat_init(void)
{

  int i,j,idx=0,iidx=0;
  static int called_flag=0;
  struct vector T;
  struct matrix L;

  if ( !called_flag ) {

    if ( debug > 2 ) printf("* estat_init: ASSIGN CHARGES, Z, and TYPE in p struct\n");

    /* set object charges and Z */
    for( i=0, idx=0; i<MAX_OBJS && object[i].used==1; i++ ){
      for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
	object[i].p_num[j] = idx;
	p[idx].chrg = object[i].vchrg[j];
	p[idx].Z    = object[i].Z[j];
	p[idx].R    = object[i].R[j];
	p[idx].typ  = object[i].typ;
	p[idx].moved  = 1;
	if ( debug > 2 ) printf("* idx= %4d  obj %2d vrt %2d chrg= %+8.3f  Z= %3d  typ= %3d\n",
				idx, i, j, p[idx].chrg, p[idx].Z , p[idx].typ);
	idx++;
      }
    }

    called_flag=1;
  }


  /*********************************************/
  /* fill all the atom positions with new info */
  /*********************************************/

  L = L_e3(&cell);

  for( i=0, idx=0; i<MAX_OBJS && object[i].used==1; i++ ){
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      idx++;
      iidx = idx - 1;

      T = rezone( (*(object+i)).v[j] );
      T = cart(&L,&T);
      p[iidx].x = T.x;
      p[iidx].y = T.y;
      p[iidx].z = T.z;
      p[iidx].moved = object[i].vmoved[j];

    }
  }

  if ( debug>4 ) {
    printf("* estat_init: atom info:\n");
    printf("* [ n]%10s%10s%10s%8s%8s%8s%4s%6s%4s\n",
	   "x","y","z","w","R","chrg","Z","typ","mv");
    for(i=0; i<numat; i++){
      printf("* [%2d]%10.5f%10.5f%10.5f%8.3f%8.3f%8.3f%4d%6d%4d\n",
	     i,p[i].x,p[i].y,p[i].z,p[i].w,p[i].R,p[i].chrg,p[i].Z,p[i].typ,p[i].moved);
    }
  }

  return;
}
