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

extern struct obj *object;
extern struct atom *s;

/*
 * The function energ_init will initialize the s atom structure with
 * all the cations and anions.  Call this function before calling the
 * energy function Energy_rep.  This takes the atoms in FRAC coords
 * and fills the structure for Energy_rep with atoms in FRAC coords.
 * The Energy_rep function will get the cartesian positions!!!!!!!
 *
 */

void energ_init(void)
{

  int i,j,idx=0,iidx=0;
  static int called_flag=0;
  struct vector T;

  /* set all radii ONCE, the first time the function is called */
  if ( !called_flag ) {

    if ( debug > 3 ) printf("* energ_init: SETTING RADII and TYPE\n");

    /* set radii and Z */
    for(i=0, idx=0; i<MAX_OBJS && object[i].used==1; i++){
      for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
	(*(s+idx)).R = object[i].R[j];
	(*(s+idx)).Z = object[i].Z[j];
	(*(s+idx)).chrg = object[i].vchrg[j];
	(*(s+idx)).typ = object[i].typ;
	(*(s+idx)).moved = 1;
	if ( debug > 2 ) printf("* idx= %4d  obj %2d vrt %2d R= %+6.3f  Z= %3d  typ= %3d\n",
				idx, i, j, s[idx].R, s[idx].Z, s[idx].typ );
	idx++;
      }
    }
    
    called_flag=1;
  }
  
  /*****************************************************/
  /*****************************************************/

  /* fill all the necessary positions with new info */

  /* set radii and Z */
  for(i=0, idx=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      idx++;
      iidx = idx-1;

      T = rezone( (*(object+i)).v[j] );
      
      (*(s+iidx)).x = T.x;
      (*(s+iidx)).y = T.y;
      (*(s+iidx)).z = T.z;
      (*(s+iidx)).moved = object[i].vmoved[j];
      if ( debug >3 ) printf("* energ_init: s[%3d] object [%3d %3d] %6.3f %6.3f %6.3f   moved= %d\n",
			     iidx,i,j,(*(s+iidx)).x,(*(s+iidx)).y,(*(s+iidx)).z, (*(s+iidx)).moved );

    }
  }
  
  if ( debug > 4 ) {
    printf("* energ_init: atom info:\n");
    printf("* [ n]%10s%10s%10s%8s%8s%8s%4s%6s%4s\n",
	   "x","y","z","w","R","chrg","Z","typ","mv");
    for(i=0; i<numat; i++){
      printf("* [%2d]%10.5f%10.5f%10.5f%8.3f%8.3f%8.3f%4d%6d%4d\n",
	     i,s[i].x,s[i].y,s[i].z,s[i].w,s[i].R,s[i].chrg,s[i].Z,s[i].typ,s[i].moved);
    }
  }

  return;
}
