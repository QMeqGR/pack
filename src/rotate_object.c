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

extern double c2vtol;

extern struct obj *object;

/***************************************/
/*                                     */
/*          Rotate an object           */
/*                                     */
/***************************************/
void rot_object(int n,
		struct matrix *L,
		struct cellprm *cell,
		struct matrix R)
{
  
  /* plan:
     0. get cartesian coords for each atom
     1. translate the object to the origin
     2. perform rotation
     3. check for distance violations
     4. translate object back to location it came from
     5. put back in lattice vector coordinates
  */

  int i,nv;
  int debug=0;

  struct vector c={0};
  struct matrix Linv;

  nv = object[n].num_vert;
  /* mark all vertices as moved */
  for(i=0;i<MAX_VERTS && object[n].vused[i]==1;i++) object[n].vmoved[i] = 1;

  if ( nv==1 ) return; /* no need to rotate a single atom */
  c = cart( L, &object[n].v[0] );

  Linv = Linverse(L);

  if ( debug ) {
    printf("rot_objec: cell:  a   b   c %15.10f%15.10f%15.10f\n",
	   (*cell).a,
	   (*cell).b,
	   (*cell).c);
    printf("rot_objec: cell:  angles    %15.10f%15.10f%15.10f\n",
	   (*cell).alph,
	   (*cell).beta,
	   (*cell).gamm);
    printf("rot_objec: center cartesian %15.10f%15.10f%15.10f\n",c.x,c.y,c.z);
    printf("rot_objec: L1_ %15.10f%15.10f%15.10f\n",(*L).r11,(*L).r12,(*L).r13);
    printf("rot_objec: L2_ %15.10f%15.10f%15.10f\n",(*L).r21,(*L).r22,(*L).r23);
    printf("rot_objec: L3_ %15.10f%15.10f%15.10f\n",(*L).r31,(*L).r32,(*L).r33);
  }

  /*******************************************************/
  /* 0. get cartesian coords for each atom               */

  for (i=0; i<nv; i++) object[n].v[i] = cart( L, &object[n].v[i] );

  /* debug printing */
  if (debug) for (i=0; i<nv; i++) printf("rot_objec: pre trans ob[%d] v[%d]: atom C %15.9lf%15.9lf%15.9lf\n",
					 n,i,object[n].v[i].x,object[n].v[i].y,object[n].v[i].z);

  /*******************************************************/
  /* 1. translate the object to the origin */
  for (i=0; i<nv; i++)  object[n].v[i] = vsub( object[n].v[i], c );

  /* debug printing */
  if (debug) for (i=0; i<nv; i++) printf("rot_objec: trans ob[%d] v[%d]: atom C %15.9lf%15.9lf%15.9lf\n",
					 n,i,object[n].v[i].x,object[n].v[i].y,object[n].v[i].z);

  /*******************************************************/
  /* 2. Rotate the vertex points about the origin */
  for (i=1; i<nv; i++) object[n].v[i] = multiply(&R,&object[n].v[i]);

  if (debug) for (i=0; i<nv; i++) printf("rot_objec: rotated object[%d] v[%d]: atom C %15.9lf%15.9lf%15.9lf\n",
					 n,i,object[n].v[i].x,object[n].v[i].y,object[n].v[i].z);

  /*******************************************************/
  /* 3. Check for distance violations                    */
  for (i=1; i<nv; i++) {
    if ( !tol( vmagsq( object[n].v[i] ) , object[n].cvdsq[i], c2vtol ) ) {
      printf("# rot_objec: error: objec [%d] vertex [%d]\n",n,i);
      printf("# rot_objec: step 3; objec center at origin\n");
      printf("# rot_objec: expecting d= %15.10f, found %15.10f   c2vtol= %f\n", object[n].cvd[i] , vmag( object[n].v[i] ) , c2vtol );
      printf("# rot_objec: vertex [%d]: %15.10f%15.10f%15.10f   (scaled coords)\n",
	     i,
	     object[n].v[i].x,
	     object[n].v[i].y,
	     object[n].v[i].z);
      printf("# rot_objec: cell: %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
	     (*cell).a,(*cell).b,(*cell).c,R2D*(*cell).alph,R2D*(*cell).beta,R2D*(*cell).gamm);
      printf("# rot_objec: exiting.\n");
      exit(0);
    }
  }

  /*******************************************************/
  /* 4. translate back to where it came from */
  for (i=0; i<nv; i++) object[n].v[i] = vadd( object[n].v[i], c );

  if (debug) for (i=0; i<nv; i++) printf("rot_objec: trans object[%d] v[%d]: atom C %15.9lf%15.9lf%15.9lf\n",
					 n,i,object[n].v[i].x,object[n].v[i].y,object[n].v[i].z);

  /*******************************************************/
  /* 5. change representation to lattice coordinates     */
  for (i=0; i<nv; i++) object[n].v[i] = UTMvmult( &Linv, &object[n].v[i] );

  if (debug) for (i=0; i<nv; i++) printf("rot_objec: frac object[%d] v[%d]: atom C %15.9lf%15.9lf%15.9lf\n",
					 n,i,object[n].v[i].x,object[n].v[i].y,object[n].v[i].z);

  return;
}

