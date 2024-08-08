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


/*
 * make_tetr
 * 
 * Given the center vector, two vectors defining another frame in E^3, and the
 * distance to the verices from the center, will return a 'tetrahedron'
 * rotated with the bottom face of the 'standard' tetrahedron pointing
 * toward z-hat direction of the new frame.
 *
 * Some of this code is taken from my awk program tetrahedron.awk
 *
 * 1. take a tetrahedron (standard defined in the code)
 * 2. translate the center to the origin
 * 3. scale the center to vertex distances (given d)
 * 3.5 scale the distances of the vertexes to give fractional coords
 * 4. rotate it approprately (given n1,n2)
 * 5. translate it out to the position (given c)
 * 6. return the tetrahedron
 *
 */
struct tetrahedron make_tetr(struct vector c,
			     double phi,
			     double the,
			     double psi,
			     double d,
			     struct cellprm cell)
{
  int i;
  int num_vert=5;
  int debug=0;    /* only for very low level debugging */

  struct vector w[5];
  struct tetrahedron tetr;
  
  struct vector centroid;
  
  struct matrix R;

  /* input debug */
  /* printf(" c = %lf %lf %lf\n",c.x,c.y,c.z); */
  
  /* initialize the default tetrahedron, from diffraction book 1, 1997 */
  w[1].x = 0;
  w[1].y = 1/SRT3;
  w[1].z = 0;
  
  w[2].x = -0.5;
  w[2].y = -1/(2*SRT3);
  w[2].z = 0;
  
  w[3].x = 0.5;
  w[3].y = -1/(2*SRT3);
  w[3].z = 0;
  
  w[4].x = 0;
  w[4].y = 0;
  w[4].z = SRT2/SRT3;
  
  /* the centroid */
  centroid.x = w[0].x = 0;
  centroid.y = w[0].y = 0;
  centroid.z = w[0].z = 1/(2*SRT6);
  
  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: pre tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* translate the centroid (and whole tetr) to the origin */
  for (i=0; i<num_vert; i++) {
    w[i].x -= centroid.x;
    w[i].y -= centroid.y;
    w[i].z -= centroid.z;
  }
  /* scale the tetrahedral vertex lengths (from centroid) */
  for (i=0; i<num_vert; i++) {
    w[i].x *= (d * 2*SRT2/SRT3);
    w[i].y *= (d * 2*SRT2/SRT3);
    w[i].z *= (d * 2*SRT2/SRT3);
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: pre rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }
  
  /*
   * rotate the tetrahedron so the 'bottom' face normal is pointing along
   * the direction of z-hat in the new frame.  This will be the default
   * since the face normal for the 'bottom' face points along -zhat.
  */


  /* orientation decided by random 3 numbers for phi,theta,psi */
  R = R_ptp(phi,the,psi);
  
  /* Rotate the vertex points on the centered polygon into the nhat frame */
  for (i=0; i<num_vert; i++) w[i] = multiply(&R,&w[i]);

  /* scale the tetrahedral vertex lengths (to frac coords) */
  for (i=1; i<num_vert; i++) {
    w[i].x /= cell.a;
    w[i].y /= cell.b;
    w[i].z /= cell.c;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: post rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* Translate the polygon */
  for (i=0; i<num_vert; i++) w[i] = vadd(w[i],c);

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_tetr: post tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* make the tetrahedron structure */
  for (i=0; i<num_vert;  i++) tetr.v[i] = w[i];

  return(tetr);

}

