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

#define ORTH_MINIMUM 0.7

/**************************************/
/*                                    */
/* distance between two points        */
/*                                    */
/**************************************/
double dist(struct vector *a,struct vector *b)
   {
   double c;
   double ax,ay,az,bx,by,bz;

   ax = (*a).x; ay = (*a).y; az = (*a).z;
   bx = (*b).x; by = (*b).y; bz = (*b).z;

   c = sqrt((ax-bx)*(ax-bx)+
	    (ay-by)*(ay-by)+
	    (az-bz)*(az-bz));
   return(c);
   }

/**************************************/
/*                                    */
/* distance between two points        */
/* NO Periodic Boundary Conditions    */
/*                                    */
/**************************************/
double dist_NOpbc(struct matrix *L,struct vector *v,struct vector *w)
{
  return( vmag( vsub( cart(L,v), cart(L,w) ) ) );
}

/**************************************/
/*                                    */
/* distance between two points        */
/* Periodic Boundary Conditions       */
/*                                    */
/**************************************/
double dist_pbc(struct matrix *L,struct vector *v,struct vector *w)
{
  double d;

  d = sqrt( distsq_pbc(L,v,w) );

  return(d);
}

/* This function was specificially written for Energy_rep() */
double distsq_pbc(struct matrix *L,struct vector *v,struct vector *w)
{
  static struct vector t1;

  t1 = vsub( rezone( *v ) , rezone( *w ) );
  t1 = rezone( t1 );

  if ( t1.x > 0.5 ) t1.x -= 1.0;
  if ( t1.y > 0.5 ) t1.y -= 1.0;
  if ( t1.z > 0.5 ) t1.z -= 1.0;

  return ( vmagsq( cart(L,&t1) ) );
}


/*
 * This function was written for Energy_rep
 * to take care of self-interaction of
 * periodic images.
 *
 *
*/
double distsq_ssrep(struct matrix *L,
		    struct vector *v,
		    struct vector *w,
		    struct cellprm *cell,
		    double orth,
		    int same_atom)
{

  int a1,a2,a3;
  int n_a1=1,n_a2=1,n_a3=1;
  extern int debug;
  extern int non_per_ewald;

  double d=0,dmin=1e12;

  struct vector ta1,ta2,Rcel,T;
  struct vector A,B,C;

  A = (*cell).bas.A;
  B = (*cell).bas.B;
  C = (*cell).bas.C;

  if ( orth < ORTH_MINIMUM ) { n_a1=2; n_a2=2; n_a3=2; }
  if ( non_per_ewald == 1  ) { n_a1=0; n_a2=0; n_a3=0; }

  for (a1=-n_a1; a1<(n_a1+1); a1++) {
    ta1 = vsmult(a1,A);
    for (a2=-n_a2; a2<(n_a2+1); a2++) {
      ta2 = vsmult(a2,B);
      for (a3=-n_a3; a3<(n_a3+1); a3++) {
	
	if ( a1==0 && a2==0 && a3==0 && same_atom ) continue;
	
	Rcel = vadd( ta1 , vadd( ta2 , vsmult(a3,C) ) );

	T = vadd( Rcel , cart(L,w) );

	d = vmagsq( vsub( T, cart(L,v) ) );

	if ( debug > 7 ) {
	  printf("*\n");
	  printf("* distsq_ssrep: a1 a2 a3 = %+2d %+2d %+2d\n",a1,a2,a3);
	  if ( debug > 8 ){
	    printf("* distsq_ssrep: v: %10.4f%10.4f%10.4f (frac)\n",(*v).x,(*v).y,(*v).z);
	    printf("* distsq_ssrep: w: %10.4f%10.4f%10.4f (frac)\n",(*w).x,(*w).y,(*w).z);
	  }
	  printf("* distsq_ssrep: v: %10.4f%10.4f%10.4f\n",
		 cart(L,v).x,cart(L,v).y,cart(L,v).z);
	  printf("* distsq_ssrep: w: %10.4f%10.4f%10.4f\n",
		 cart(L,w).x,cart(L,w).y,cart(L,w).z);
	  printf("* distsq_ssrep: T: %10.4f%10.4f%10.4f\n",T.x,T.y,T.z);
	  if ( debug > 8 ){
	    printf("* distsq_ssrep: cell. a b c = %10.4f%10.4f%10.4f\n",(*cell).a,(*cell).b,(*cell).c);
	    printf("* distsq_ssrep: A  .  x y z = %10.4f%10.4f%10.4f\n",A.x,A.y,A.z);
	    printf("* distsq_ssrep: ta1.  x y z = %10.4f%10.4f%10.4f\n",ta1.x,ta1.y,ta1.z);
	    printf("* distsq_ssrep: B  .  x y z = %10.4f%10.4f%10.4f\n",B.x,B.y,B.z);
	    printf("* distsq_ssrep: ta2.  x y z = %10.4f%10.4f%10.4f\n",ta2.x,ta2.y,ta2.z);
	    printf("* distsq_ssrep: C  .  x y z = %10.4f%10.4f%10.4f\n",C.x,C.y,C.z);
	    printf("* distsq_ssrep: Rcel  x y z = %10.4f%10.4f%10.4f\n",Rcel.x,Rcel.y,Rcel.z);
	  }
	  printf("* distsq_ssrep: dsq = %15.10e\n",d);
	  printf("* distsq_ssrep: d   = %15.10e ( T-v dist )\n*\n",sqrt(d) );
	}

	if ( d < dmin ) dmin=d;
      }
    }
  }

  return(dmin);
}
