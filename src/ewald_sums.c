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
//#include <mathimf.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "global_defs.h"
#include "packlib.h"

int ngs = 0;
int nrs = 0;
double *vgx, *vgy, *vgz, *vgexp;
double *vrx, *vry, *vrz, *vrsq;

extern struct plane surf_plane;


/********************************************************/
/********************************************************/
/*             Ewald summation routines                 */
/********************************************************/
/********************************************************/
/* These routines were adapted from V.Ozolins fortran code */
/* They were not f2c'd. */

/*
 * Recip: return reciprocal space vectors
 */
struct basis recip(struct basis bas)
{
  struct basis rec;
  double vol,tpi;

  /*  tpi = 8.0*atan(1.0); */
  tpi = 1.0;
  vol = vtriple( bas.A, bas.B, bas.C );

  /*
  printf("vol=%f\n",vol);
  printf("%f %f %f\n",bas.A.x,bas.A.y,bas.A.z);
  printf("%f %f %f\n",bas.B.x,bas.B.y,bas.B.z);
  printf("%f %f %f\n",bas.C.x,bas.C.y,bas.C.z);
  */

  rec.A = vsmult( tpi/vol, crossprod( bas.B, bas.C ) );
  rec.B = vsmult( tpi/vol, crossprod( bas.C, bas.A ) );
  rec.C = vsmult( tpi/vol, crossprod( bas.A, bas.B ) );

  /*
  printf("%f %f %f\n",rec.A.x,rec.A.y,rec.A.z);
  printf("%f %f %f\n",rec.B.x,rec.B.y,rec.B.z);
  printf("%f %f %f\n",rec.C.x,rec.C.y,rec.C.z);
  */

  return(rec);
}

double find_Ewald_eta(struct basis bas,double errlim,int nat,struct atom *at)
{
  int i,j;
  int iter;
  int MAXIT=100;
  int debug=0;

  static double eta=0;
  double pi,tpi,tpi2;
  double taumax,rsq;
  double FACT;
  double x,y,rmax,ekmax,etaprev;

  struct vector dtau;
  struct basis rec;
  struct intvector K,R;

  if ( errlim > 1.0 ) {
    printf("errlim in find_Ewald_eta is too large!\n");
    exit(0);
  }

  pi   = 4*atan(1.0);
  tpi  = 2*pi;
  tpi2 = tpi*tpi;

  rec = recip(bas);
  taumax = 0.0;
  for(i=0; i<nat; i++){
    for(j=0; j<=i; j++){
      dtau = makevec( at[i].x - at[j].x,
		      at[i].y - at[j].y,
		      at[i].z - at[j].z );
      rsq = vdotprod( dtau , dtau );
      taumax = fmax( taumax , rsq );
    }
  }

  if ( fabs(eta) < 1e-5 ) eta = 1.0;
  FACT = 2.0;
  etaprev = eta;

  for(iter=0; iter<MAXIT; iter++) {
    x = 2.0 * eta * sqrt( -log( errlim ) );
    ekmax = x*x / 2.0;
    K = gbox( ekmax, 1.0, recip(bas) );

    x = sqrt( -log(errlim) ) / eta;
    y = sqrt(taumax);
    rmax = (x+y)*(x+y);
    ekmax = tpi2 * rmax / 2.0;
    R = gbox( ekmax, 1.0, bas );

    if ( K.n1*K.n2*K.n3 > R.n1*R.n2*R.n3 ){
      
      if ( iter == MAXIT ) {
	printf("WARNING: possibly suboptimal value of eta in Ewald sums!\n");
	printf("Real space sums over %d %d %d\n",R.n1,R.n2,R.n3);
	printf("K space sums over %d %d %d\n",K.n1,K.n2,K.n3);
      } else {

	if ( eta <= etaprev ) {
	  x = eta / FACT;
	} else {
	  x = eta - (eta-etaprev) / 3.0;
	}
	etaprev = eta;
	eta = x;
      }
	
    } else if ( 2*K.n1*K.n2*K.n3 < R.n1*R.n2*R.n3 ) {

	if ( iter == MAXIT ) {
	  printf("WARNING: possibly suboptimal value of eta in Ewald sums!\n");
	  printf("Real space sums over %d %d %d\n",R.n1,R.n2,R.n3);
	  printf("K space sums over %d %d %d\n",K.n1,K.n2,K.n3);
	} else {

	  if ( eta >= etaprev ) {
	    x = eta * FACT;
	  } else {
	    x = eta + (etaprev-eta) / 3.0;
	  }
	  etaprev = eta;
	  eta = x;
	}
    } else {
      if ( debug ) printf("iter=%d Found an optimal value of eta = %f\n",iter,eta);
    }

  } /* end for loop */
	  
  return(eta);
}


double Eion(struct cellprm cell, struct atom *at,int nat,double errlim,double eta)
{
  int i, i1,i2,i3;
  int iat1,iat2;
  int m,npairs;

  double pi,tpi,fpi,tpi2,pilog, etasq;
  double rmax,ekmax, omega, logerrlim;
  double x, gsq, rsq, taumax;
  double xsum, xchrg2;
  double *zprod;
  double ecc=0;

  extern int ngs, nrs;
  extern double *vgx, *vgy, *vgz, *vgexp;
  extern double *vrx, *vry, *vrz, *vrsq;

  struct intvector K,R;
  struct basis qb,rb;
  struct vector dtau;
  struct vector *difftau;

  pi   = 4*atan(1.0);
  pilog = log(pi);
  tpi  = 2*pi;
  fpi  = 4*pi;
  tpi2 = tpi*tpi;

  rb = cell.bas;
  omega = vtriple( rb.A, rb.B, rb.C );
  qb = recip( rb );
  etasq = eta*eta;
  logerrlim = log(errlim);

  // Calculate the sum of charges squared
  xchrg2 = 0.0;
  for(iat1=0; iat1<nat; iat1++)
    xchrg2 += at[iat1].chrg * at[iat1].chrg;

  /* Calculate temporary arrays for atomic position differences
     and charge products */

  npairs = nat*(nat-1)/2;
  difftau = (struct vector *)malloc( npairs * sizeof(struct vector) );
  zprod = (double *)malloc( npairs * sizeof(double) );

  m = 0;
  taumax = 0.0;
  for(iat1=0; iat1<nat; iat1++){
    for(iat2=0; iat2<iat1; iat2++){
      difftau[m].x = at[iat1].x-at[iat2].x;
      difftau[m].y = at[iat1].y-at[iat2].y;
      difftau[m].z = at[iat1].z-at[iat2].z;
      zprod[m] = 2 * at[iat1].chrg * at[iat2].chrg;
      rsq = difftau[m].x*difftau[m].x + 
	    difftau[m].y*difftau[m].y + 
	    difftau[m].z*difftau[m].z;
      if( rsq > taumax)   taumax = rsq;
      m++;
    }
  }
  taumax = sqrt( taumax );

  x = 2.0 * eta * sqrt( -logerrlim );
  ekmax = x * x / 2.0;
  K = gbox(ekmax, 1.0, qb);

//  printf("K: %d %d %d\n",K.n1,K.n2,K.n3);


  // First, set up the list of reciprocal space vectors that
  // we must sum over. This list is stored globally and accessed by
  // the deltaEion subroutine.

  if( vgx ) free (vgx);
  if( vgy ) free (vgy);
  if( vgz ) free (vgz);
  if( vgexp ) free (vgexp);

  vgx = (double *)malloc(8*K.n1*K.n2*K.n3 * sizeof(double));
  vgy = (double *)malloc(8*K.n1*K.n2*K.n3 * sizeof(double));
  vgz = (double *)malloc(8*K.n1*K.n2*K.n3 * sizeof(double));
  vgexp = (double *)malloc(8*K.n1*K.n2*K.n3 * sizeof(double));

  ngs = 0;
  ecc = 0.0;
  for(i1=-K.n1; i1<(K.n1+1); i1++){
    for(i2=-K.n2; i2<(K.n2+1); i2++){
      for(i3=-K.n3; i3<(K.n3+1); i3++){

	if ( i1==0 && i2==0 && i3==0 ) continue;

	vgx[ngs] = tpi * ( i1*qb.A.x + i2*qb.B.x + i3*qb.C.x );
	vgy[ngs] = tpi * ( i1*qb.A.y + i2*qb.B.y + i3*qb.C.y );
	vgz[ngs] = tpi * ( i1*qb.A.z + i2*qb.B.z + i3*qb.C.z );
	gsq = vgx[ngs]*vgx[ngs] + vgy[ngs]*vgy[ngs] + vgz[ngs]*vgz[ngs]; 

	x = exp(-gsq/(4.0*etasq));
	if ( x >= errlim ) {
	  /*  vgexp now stores the exponential factors:
	   (fpi/omega) * exp[-G^2/(4.0*eta^2)] / G^2;	  */
	  vgexp[ngs] = fpi * x / (gsq*omega);
	  ecc += xchrg2 * vgexp[ngs];
	  ngs++;
	}
      }
    }
  }
//  printf(" Number of G vectors %i\n",ngs);

  //Calculate size of real space box
  x = sqrt( -logerrlim ) / eta;
  rmax = (x+taumax) * (x+taumax);
  ekmax = tpi2 * rmax / 2.0;
  R = gbox(ekmax, 1.0, rb);

//  printf("YY R = %d %d %d\n",R.n1,R.n2,R.n3);

  // Set up real-space vectors in the Ewald sum.
  if( vrx ) free (vrx);
  if( vry ) free (vry);
  if( vrz ) free (vrz);

  vrx = (double *)malloc(8*R.n1*R.n2*R.n3 * sizeof(double));
  vry = (double *)malloc(8*R.n1*R.n2*R.n3 * sizeof(double));
  vrz = (double *)malloc(8*R.n1*R.n2*R.n3 * sizeof(double));

  nrs=0;
  for(i1=-R.n1; i1<(R.n1+1); i1++){
    for(i2=-R.n2; i2<(R.n2+1); i2++){
      for(i3=-R.n3; i3<(R.n3+1); i3++){

	vrx[nrs] = i1*rb.A.x + i2*rb.B.x + i3*rb.C.x;
	vry[nrs] = i1*rb.A.y + i2*rb.B.y + i3*rb.C.y;
	vrz[nrs] = i1*rb.A.z + i2*rb.B.z + i3*rb.C.z;

	x = sqrt( vrx[nrs]*vrx[nrs] + vry[nrs]*vry[nrs] + vrz[nrs]*vrz[nrs]);
	rsq = x;

	if ( x > taumax ) {
	  x -= taumax;
	  if( erfc( eta * x ) / x > errlim ) {
	    if ( i1 !=0 || i2 !=0 || i3 != 0 )
	      ecc += xchrg2 * erfc( eta * rsq ) / rsq;
	    nrs++;
	  }
	}
	else {
	  if ( i1 !=0 || i2 !=0 || i3 != 0 )
	    ecc += xchrg2 * erfc( eta * rsq ) / rsq;
	  nrs++;
	}
      }
    }
  }

//  printf(" Number of R vectors %i\n",nrs);

  /* These are the main loops in the subroutine.
     The first loop is the reciprocal space sum,
     while the second loop is the real-space sum.
     Both must be vectorized to the maximum. */

  for(m=0; m<npairs; m++){
    for (i=0; i<ngs; i++) {
      ecc += zprod[m] * vgexp[i] * 
	cos( vgx[i]*difftau[m].x +
	     vgy[i]*difftau[m].y +
	     vgz[i]*difftau[m].z );
    }
    for (i=0; i<nrs; i++) {
      dtau.x = vrx[i] + difftau[m].x;
      dtau.y = vry[i] + difftau[m].y;
      dtau.z = vrz[i] + difftau[m].z;
      rsq = sqrt(dtau.x*dtau.x + dtau.y*dtau.y + dtau.z*dtau.z);
      ecc += zprod[m] * erfc( eta * rsq ) / rsq;
    }
  }

  //  printf("YY ecc at 3: %g\n",ecc);

  xsum = 0.0;
  for(m=0; m<npairs; m++){
    xsum += zprod[m];
  }
  ecc -= pi * (xsum + xchrg2) / ( etasq * omega );
  ecc -= 2.0 * eta * xchrg2 / sqrt(pi);

  if ( difftau )  free (difftau);
  if ( zprod )    free (zprod);

  //  printf("YY ecc at 4: %g\n",ecc);
  return(ecc);
}

double deltaEion(int nmoved, int *ixmoved, struct cellprm cell, 
		 struct atom *at_old, struct atom *at_new, int nat,
		 double errlim, double eta)
{
  int i;
  int iat1,iat2,im;
  int n,m,npairs;
  int tbreak;

  extern int ngs, nrs;
  extern double *vgx, *vgy, *vgz, *vgexp;
  extern double *vrx, *vry, *vrz, *vrsq;

  double ecc=0.0;
  double xsum_old,xsum_new,rsq;
  double *zprod;

  struct vector dtau;
  struct vector *difftau_new;
  struct vector *difftau_old;

  /* Calculate temporary arrays for atomic position differences
     and charge products */

  npairs = (nat-1)*nmoved - nmoved*(nmoved-1)/2;
  difftau_new = (struct vector *)malloc( npairs * sizeof(struct vector) );
  difftau_old = (struct vector *)malloc( npairs * sizeof(struct vector) );
  zprod = (double *)malloc( npairs * sizeof(double) );

  m = 0;
  for(im=0; im<nmoved; im++){
    iat1=ixmoved[im];
    for(iat2=0; iat2<nat; iat2++){
      tbreak=0;
      for(i=0; i<=im; i++) {
	if( iat2 == ixmoved[i] ) {
	  tbreak=1;
	  continue;
	}
      }
      if( ! tbreak ) {
	difftau_old[m].x = at_old[iat1].x-at_old[iat2].x;
	difftau_old[m].y = at_old[iat1].y-at_old[iat2].y;
	difftau_old[m].z = at_old[iat1].z-at_old[iat2].z;

	difftau_new[m].x = at_new[iat1].x-at_new[iat2].x;
	difftau_new[m].y = at_new[iat1].y-at_new[iat2].y;
	difftau_new[m].z = at_new[iat1].z-at_new[iat2].z;

	zprod[m] = 2 * at_old[iat1].chrg * at_old[iat2].chrg;
	m++;
      }
    }
  }

  if (m != npairs ) {
    printf("Wrong number of changed pairs! m=%d, npairs=%d\n",m, npairs);
    exit;
  }

  // Calculate change to the reciprocal space sum.
  for(m=0; m<npairs; m++){

    xsum_old = 0.0;
    xsum_new = 0.0;
    for (i=0; i<ngs; i++) {
      xsum_old += vgexp[i] * cos( vgx[i]*difftau_old[m].x +
				  vgy[i]*difftau_old[m].y +
				  vgz[i]*difftau_old[m].z );
      xsum_new += vgexp[i] * cos( vgx[i]*difftau_new[m].x +
				  vgy[i]*difftau_new[m].y +
				  vgz[i]*difftau_new[m].z );
    }
    ecc += zprod[m] * (xsum_new - xsum_old);
  }

  //  printf("deltaE ecc at 1: %g\n",ecc);

  // Calculate change to the real space sum.
  for(m=0; m<npairs; m++){

    xsum_old = 0.0;
    xsum_new = 0.0;
    for (i=0; i<nrs; i++) {
      dtau.x = vrx[i] + difftau_old[m].x;
      dtau.y = vry[i] + difftau_old[m].y;
      dtau.z = vrz[i] + difftau_old[m].z;
      rsq = sqrt(dtau.x*dtau.x + dtau.y*dtau.y + dtau.z*dtau.z);
      xsum_old += erfc( eta * rsq ) / rsq;

      dtau.x = vrx[i] + difftau_new[m].x;
      dtau.y = vry[i] + difftau_new[m].y;
      dtau.z = vrz[i] + difftau_new[m].z;
      rsq = sqrt(dtau.x*dtau.x + dtau.y*dtau.y + dtau.z*dtau.z);
      xsum_new += erfc( eta * rsq ) / rsq;
    }
    ecc += zprod[m] * (xsum_new - xsum_old);
  }

  if ( difftau_new )  free (difftau_new);
  if ( difftau_old )  free (difftau_old);
  if ( zprod )    free (zprod);
  
  return(ecc);
}

void freeEionArrays() {

  extern int ngs, nrs;
  extern double *vgx, *vgy, *vgz, *vgexp;
  extern double *vrx, *vry, *vrz, *vrsq;

  if( vrx ) free (vrx);
  if( vry ) free (vry);
  if( vrz ) free (vrz);
  if( vgx ) free (vgx);
  if( vgy ) free (vgy);
  if( vgz ) free (vgz);
  if( vgexp ) free (vgexp);
  ngs = 0;
  nrs = 0;  
  
  return;
}

/* c============================================================ */
/* c     Find the max. dimensions of the reciprocal-space */
/* c     box which contains all plane waves with Ekin < Emax. */
/* c */
/* c  DATE: Wed Oct 16 14:54:34 MDT 1996 */
/* c============================================================ */
struct intvector gbox(double emax, double alat, struct basis bas)
{
  double pi,tpi,tpiba,Gmax;
  double q11,q12,q13,q22,q23,q33,rvol,range1,range2,range3;
  struct intvector N;
  
  pi    = 4*atan(1.0);
  tpi   = 2*pi;
  tpiba = tpi / alat;
  Gmax  = sqrt( 2 * emax ) / tpiba;

  /* can generalize this to non-orthorhombic later */
  q11 = vdotprod(bas.A,bas.A);
  q12 = vdotprod(bas.A,bas.B);
  q13 = vdotprod(bas.A,bas.C);
  q22 = vdotprod(bas.B,bas.B);
  q23 = vdotprod(bas.B,bas.C);
  q33 = vdotprod(bas.C,bas.C);

  rvol = vtriple( bas.A, bas.B, bas.C );
  
  range1 = rvol / sqrt(q22*q33 - q23*q23);
  range2 = rvol / sqrt(q11*q33 - q13*q13);
  range3 = rvol / sqrt(q11*q22 - q12*q12);
    
  N.n1 = (int)( Gmax / range1 + 1.0 );
  N.n2 = (int)( Gmax / range2 + 1.0 );
  N.n3 = (int)( Gmax / range3 + 1.0 );

  /*
  printf("q11=%f\n",q11);
  printf("q12=%f\n",q12);
  printf("q13=%f\n",q13);
  printf("q22=%f\n",q22);
  printf("q23=%f\n",q23);
  printf("q33=%f\n",q33);
  printf("gmax=%f\n",Gmax);
  printf("N: %d %d %d\n",N.n1,N.n2,N.n3);
  */

  return(N);
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Non-Peroidic Code %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
double Eion_nopbc(struct cellprm cell,
		  struct atom *at,
		  int nat,
		  double errlim,
		  double eta)
{
  int iat1,iat2;

  double ecc=0;
  double logerrlim;
  double rsq;

  extern double *ecctable;

  struct vector dtau;

  /* our cell is orthorhombic, so this is easy 
  rblen = cell.a;
  rblen = fmax(rblen, cell.b);
  rblen = fmax(rblen, cell.c);

  rmax = nat * rblen;
  ekmax = tpi2 * rmax / 2.0;
  */

  /*
  printf("YY cell:\n");
  printf("YY %15.10f%15.10f%15.10f\n",cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
  printf("YY %15.10f%15.10f%15.10f\n",cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
  printf("YY %15.10f%15.10f%15.10f\n",cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);

  printf("YY atoms:  nat=%d\n",nat);
  for(i1=0;i1<nat;i1++){
    printf("YY %15.12f%15.12f%15.12f%15.12f\n",at[i1].x,at[i1].y,at[i1].z,at[i1].chrg);
  }
  */


  ecc = 0.0;

  logerrlim = log(errlim);

  for(iat1=0; iat1<nat; iat1++){
    for(iat2=0; iat2<=iat1; iat2++){
      dtau = makevec( at[iat1].x-at[iat2].x,
		      at[iat1].y-at[iat2].y,
		      at[iat1].z-at[iat2].z );
      rsq = sqrt( vdotprod( dtau, dtau ) );
      if ( iat1 == iat2 ) continue;
      if ( at[iat1].moved==0 && at[iat2].moved==0 ) continue;

      
      /* Uncomment this section if the value of ecc is getting too large
	 and you wish to debug the atom locations           */
      /*
	if ( rsq < 1e-10 ) {
	printf("* WARN: Eion: atomic positions coincide!\n");
	printf("* i1 i2 i3: %d %d %d\n",i1,i2,i3);
	printf("* at1 at2: %d %d\n",iat1,iat2);
	printf("* ch1 ch2: %f %f\n",at[iat1].chrg,at[iat2].chrg);
	printf("* at1:   %+10.6e %+10.6e %+10.6e\n",at[iat1].x,at[iat1].y,at[iat1].z);
	printf("* at2:   %+10.6e %+10.6e %+10.6e\n",at[iat2].x,at[iat2].y,at[iat2].z);
	printf("* a1-a2: %+10.6e %+10.6e %+10.6e\n",at[iat1].x-at[iat2].x,at[iat1].y-at[iat2].y,at[iat1].z-at[iat2].z);
	printf("* vr:    %+10.6e %+10.6e %+10.6e\n",vr.x,vr.y,vr.z);
	printf("* dtau:  %+10.6e %+10.6e %+10.6e\n",dtau.x,dtau.y,dtau.z);
	printf("* eta = %20.18e\n",eta);
	printf("* rsq = %20.18e\n",rsq);
	printf("* erfc(rsq) = %e\n",erfc(rsq));
	printf("* x = %20.18e\n",x);
	printf("* \n");
	}
      */
      
      ecctable[ nat * iat1 + iat2 ] = at[iat1].chrg * at[iat2].chrg / rsq;
    }
  }
  
  /*  printf("YY ecc at 2: %g\n",ecc); */

  for(iat1=0; iat1<nat; iat1++){
    for(iat2=0; iat2<=iat1; iat2++){
      ecc += ecctable[ nat * iat1 + iat2 ];
    }
  }
  
  return(ecc);
}


