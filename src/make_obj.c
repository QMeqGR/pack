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


extern struct cellprm cell;
extern struct obj *object,*object_orig;

/************************
 * make_object
 ************************/
int make_object(int N,
		struct vector c,
		double phi,
		double the,
		double psi)
{
  int i,j;
  int nv=0;
  int debug=0;    /* only for very low level debugging */

  struct vector w[MAX_VERTS]={0};
  struct vector centroid={0},v={0},vnew={0};
  struct matrix R,L,Linv;

  /* input debug */
  /* printf(" c = %lf %lf %lf\n",c.x,c.y,c.z); */

  if ( debug ) printf("f:make_object: MAX_VERTS = %d\n",MAX_VERTS);


  /*----------------------------------------
    General Object
    ---------------------------------------*/
  nv = object[N].num_vert;
  for( j=0; j<MAX_VERTS && object_orig[N].vused[j]==1; j++ ){
    w[j] = object_orig[N].v[j];
  }
  centroid = w[0];
  
  /* debug printing */
  if (debug) {
    for (i=0; i<nv; i++) printf("f:make_object: pre  tran atom [%2d] %15.9lf%15.9lf%15.9lf cart\n",
				i,w[i].x,w[i].y,w[i].z);
    printf("f:make_object: ---\n");
  }

  /* translate the centroid (and whole object) to the origin */
  for (i=0; i<nv; i++) {
    w[i].x -= centroid.x;
    w[i].y -= centroid.y;
    w[i].z -= centroid.z;
  }

  /* debug printing */
  if (debug) {
    for (i=0; i<nv; i++) printf("f:make_object: pre  rot  atom [%2d] %15.9lf%15.9lf%15.9lf cart\n",
				i,w[i].x,w[i].y,w[i].z);
    printf("f:make_object: ---\n");
  }
  
  /* orientation decided by random 3 numbers for phi,theta,psi */
  R = R_ptp(phi,the,psi);
  
  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=0; i<nv; i++) w[i] = multiply(&R,&w[i]);

  /* debug printing */
  if (debug) {
    for (i=0; i<nv; i++) printf("f:make_object: post  rot  atom [%2d] %15.9lf%15.9lf%15.9lf cart\n",
				i,w[i].x,w[i].y,w[i].z);
    printf("f:make_object: ---\n");
  }

  /* put the coordinates in fractional coords */
  for (i=1; i<nv; i++) {
    v = makevec(w[i].x,w[i].y,w[i].z);
    L = L_e3( &cell );
    Linv = Linverse( &L );
    vnew = UTMvmult( &Linv , &v );
    w[i].x = vnew.x;
    w[i].y = vnew.y;
    w[i].z = vnew.z;
  }
  
  /* Translate the object */
  for (i=0; i<nv; i++) w[i] = vadd(w[i],c);


  /* debug printing */
  if (debug) {
    for (i=0; i<nv; i++) printf("f:make_object: post tran atom [%2d] %15.9lf%15.9lf%15.9lf frac\n",
				i,w[i].x,w[i].y,w[i].z);
    printf("f:make_object: ---\n");
  }

  /* fill the object vertices with coordinates */
  for (i=0; i<nv;  i++) object[N].v[i] = w[i];

  return(1);

}

