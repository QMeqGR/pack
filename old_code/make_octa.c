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
 * make_octa
 * 
 * Given the center vector, two vectors defining another frame in E^3, and the
 * distance to the verices from the center, will return an 'octahedron'.
 *
 *
 * 1. take an octahedron (standard defined in the code)
 * 2. translate the center to the origin
 * 3. scale the center to vertex distances (given d)
 * 3.5 scale the distances of the vertexes to give fractional coords
 * 4. rotate it approprately (given n1,n2)
 * 5. translate it out to the position (given c)
 * 6. return the octaahedron
 *
 */
struct octahedron make_octa(struct vector c,
			    double phi,
			    double the,
			    double psi,
			    double d,
			    struct cellprm cell)
{
  int i;
  int num_vert=7;
  int debug=0;    /* only for very low level debugging */

  struct vector w[7];
  struct octahedron octa;
  
  struct vector centroid;

  struct matrix R;

  /* input debug */
  /* printf(" c = %lf %lf %lf\n",c.x,c.y,c.z); */
  
  /* initialize the default octahedron */
  /* the centroid */
  centroid.x = w[0].x = 0;
  centroid.y = w[0].y = 0;
  centroid.z = w[0].z = 0;

  w[1].x =  1;
  w[1].y =  0;
  w[1].z =  0;
  
  w[2].x = -1;
  w[2].y =  0;
  w[2].z =  0;
  
  w[3].x =  0;
  w[3].y =  1;
  w[3].z =  0;
  
  w[4].x =  0;
  w[4].y = -1;
  w[4].z =  0;

  w[5].x =  0;
  w[5].y =  0;
  w[5].z =  1;
  
  w[6].x =  0;
  w[6].y =  0;
  w[6].z = -1;
  
  
  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: pre tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* translate the centroid (and whole polygon) to the origin */
  for (i=0; i<num_vert; i++) {
    w[i].x -= centroid.x;
    w[i].y -= centroid.y;
    w[i].z -= centroid.z;
  }
  /* scale the polygon vertex lengths (from centroid) */
  for (i=0; i<num_vert; i++) {
    w[i].x *= d;
    w[i].y *= d;
    w[i].z *= d;
  }

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: pre rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }
  
  /*
   * rotate the octahedron 
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
    printf("f:make_octa: post rot atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* Translate the polygon  */
  for (i=0; i<num_vert; i++) w[i] = vadd(w[i],c);

  /* degub printing */
  if (debug) for (i=0; i<num_vert; i++) {
    printf("f:make_octa: post tran atom C %15.9lf%15.9lf%15.9lf\n",w[i].x,w[i].y,w[i].z);
  }

  /* make the octahedron structure */
  for (i=0; i<num_vert;  i++) octa.v[i] = w[i];

  return(octa);

}
