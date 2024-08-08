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



#include <math.h>
#define NRANSI
#define float double
#include <stdio.h>
#include <stdlib.h>
#include "packlib.h"
#include "global_defs.h"
#define PI (3.14159265358979323846264338327950288419716939937510)
#define R2D ( 180.0 / PI )
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}


double aspect(void);
void amoeba_ehm(double **p, double y[], int ndim, double ftol,
		struct Energy (*funk)(void), int *nfunk)
{
  double amotry_ehm(double **p, double y[], double psum[], int ndim,
		   struct Energy (*funk)(void), int ihi, double fac);
  int i,ihi=0,ilo=0,inhi=0,j,mpts=ndim+1;
  int k;
  int n_elo=0;
  int chunks=0;
  int dchunks=2000;
  int loop_count=0;
  extern int ecol;
  extern int debug;
  extern int simp_nmax;
  extern double rtol;
  double sum,swap,ysave,ytry;
  double elo=1e10;
  extern double *psum;
  struct Energy E;
  
  *nfunk=0;

  /* GET PSUM */
  for (j=0;j<ndim;j++) {		
    for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];
    psum[j]=sum;
  }

  /* start the main simplex loop */
  for (;*nfunk<simp_nmax;) {

    ilo=0;
    ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);

    for (i=0;i<mpts;i++) {
      if (y[i] <= y[ilo]) ilo=i;
      if (y[i] > y[ihi]) {
	inhi=ihi;
	ihi=i;
      } else if (y[i] > y[inhi] && i != ihi) inhi=i;
    }


    rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

    /* print some info to the screen */
    if ( *nfunk > 10 && *nfunk > chunks-1 ) {
      simplex_trans(p[ilo]);
      E = (*funk)();
      if ( E.ecc < -1e8 ){
	printf("* WARN: amoeba: ecc getting too large: %10.6e   Halting Simplex.\n",E.ecc);
	simp_nmax = 0;
      }
      if ( ecol && debug > -1 ) {
	printf("*X------------------------------------------------------------\n");
	printf("*X %20s%12d%25s\n","nfunk  =",*nfunk,"en function calls");
	printf("*X %20s%12.3e%25s\n","rtol  =",rtol,"simplex variance");
	printf("*X %20s%12.5f%25s\n","pf  =",E.pf,"packing fraction");
	printf("*X %20s%12.5f%25s\n","vol  =",E.vol,"cell volume");
	printf("*X %20s%12.5f%25s\n","ortho  =",E.orth,"cv/ortho");
	printf("*X %20s%12.3f%25s\n","lat a  =",E.a,"lattice vector");
	printf("*X %20s%12.3f%25s\n","lat b  =",E.b,"lattice vector");
	printf("*X %20s%12.3f%25s\n","lat c  =",E.c,"lattice vector");
	printf("*X %20s%12.3f%25s\n","alph  =",E.alph*R2D,"angle 2--3");
	printf("*X %20s%12.3f%25s\n","beta  =",E.beta*R2D,"angle 1--3");
	printf("*X %20s%12.3f%25s\n","gamm  =",E.gamm*R2D,"angle 1--2");
	printf("*X %20s%12.5f%25s\n","ar  =",aspect(),"aspect ratio");
	printf("*X %20s%12.3e%25s\n","ecc  =",E.ecc,"electrostatic");
	printf("*X %20s%12.3e%25s\n","eng  =",E.eng,"repulsive");
	printf("*X %20s%12.3e%25s\n","bnd  =",E.bnd,"boundary");
	printf("*X %20s%12.3e%25s\n","tot  =",E.tot,"total energy");
      } else if ( debug > -1 ){
	printf("%2s%10d%15.4e%10.5f%10.5f%10.5f%18.7e%18.7e%18.7e\n","*X",
	       *nfunk,rtol,E.pf,aspect(),E.orth,E.ecc,E.eng,E.tot);
      }
      fflush(stdout);
      chunks += dchunks;
    }


    if (rtol < ftol) {
      SWAP(y[0],y[ilo])
	for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i])
	 
      simplex_trans(p[0]); /* move to the lowest energy state */
      break;
    }

    *nfunk += 2;

    if ( debug > 4 ) {
      printf("amoeba: loop count: %d\n",loop_count++);
      printf("psum: ");
      for(k=0;k<ndim;k++) printf(" %f",psum[k]);
      printf("\n");
      fflush(stdout);
    }

    ytry=amotry_ehm(p,y,psum,ndim,funk,ihi,-1.0);
    if ( debug > 5 ) {
      E = (*funk)();
      printf("%2s%15s%15.8f%15.8f%18.7e%18.7e%18.7e\n","* ","A",rtol,E.pf,E.ecc,E.eng,E.tot);
      fflush(stdout);
    }
    if (ytry <= y[ilo]) {
      ytry=amotry_ehm(p,y,psum,ndim,funk,ihi,2.0);
      if ( debug > 5 ) {
	E = (*funk)();
	printf("%2s%15s%15.8f%15.8f%18.7e%18.7e%18.7e\n","* ","B",rtol,E.pf,E.ecc,E.eng,E.tot);
	fflush(stdout);
      }
    }
    else if (ytry >= y[inhi]) {

      ysave=y[ihi];
      ytry=amotry_ehm(p,y,psum,ndim,funk,ihi,0.5);
      if ( debug > 5 ) {
	E = (*funk)();
	printf("%2s%15s%15.8f%15.8f%18.7e%18.7e%18.7e\n","* ","C",rtol,E.pf,E.ecc,E.eng,E.tot);
	fflush(stdout);
      }
      if (ytry >= ysave) {

	for (i=0;i<mpts;i++) {

	  if (i != ilo) {

	    for (j=0;j<ndim;j++)
	      p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);

	    simplex_trans(psum);
	    E = (*funk)();  /* arg was psum */
	    y[i] = E.tot;

	  }
	}
	*nfunk += ndim;

	for (j=0;j<ndim;j++) {		
	  for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];
	  psum[j]=sum;
	}

      }
    }
    else --(*nfunk);
  } /* end the main loop */

  /* if we exceed the max iterations, just return the lowest point found */
  for (i=0;i<ndim+1;i++) {
    simplex_trans(p[i]);
    E = (*funk)();
    if ( E.tot < elo ) { elo=E.tot; n_elo=i; }
  }
  simplex_trans(p[n_elo]);

}
#undef SWAP
#undef NRANSI
