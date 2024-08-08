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


/*********************/
/* MATRIX OPERATIONS */
/*********************/
double det3x3(struct matrix R)
{
  return( R.r11*(R.r22*R.r33-R.r23*R.r32)-
	  R.r12*(R.r21*R.r33-R.r23*R.r31)+
	  R.r13*(R.r21*R.r32-R.r31*R.r22)  );
}

/* create a rotation matrix for rotaions away from the standard */
/* cartesian coordinate system.  These follow Byron and Fuller's convention */
/* ranges: theta:0..Pi   phi:0..2Pi  psi:0..2Pi */
struct matrix R_ptp(double phi,double the,double psi)
{
  double sphi,cphi,sthe,cthe,spsi,cpsi;
  struct matrix R;

  sphi=sin(phi); cphi=cos(phi);
  sthe=sin(the); cthe=cos(the);
  spsi=sin(psi); cpsi=cos(psi);

  R.r11 =  cphi * cthe * cpsi - sphi * spsi;
  R.r12 =  sphi * cthe * cpsi + cphi * spsi;
  R.r13 = -sthe * cpsi;

  R.r21 = -cphi * cthe * spsi - sphi * cpsi;
  R.r22 = -sphi * cthe * spsi + cphi * cpsi;
  R.r23 =  sthe * spsi;

  R.r31 = cphi * sthe;
  R.r32 = sphi * sthe;
  R.r33 = cthe;

  return R;
}
struct matrix R_ptp_inv(struct matrix Rf)
{
  double r[9];
  struct matrix R;

  r[0] = Rf.r11;
  r[1] = Rf.r12;
  r[2] = Rf.r13;

  r[3] = Rf.r21;
  r[4] = Rf.r22;
  r[5] = Rf.r23;

  r[6] = Rf.r31;
  r[7] = Rf.r32;
  r[8] = Rf.r33;

  invt_matrx(3,r);

  R.r11 = r[0];
  R.r12 = r[1];
  R.r13 = r[2];

  R.r21 = r[3];
  R.r22 = r[4];
  R.r23 = r[5];

  R.r31 = r[6];
  R.r32 = r[7];
  R.r33 = r[8];

  return R;
}


/**************************************/
/*                                    */
/*    Matrix on Vector  multiply      */
/*                                    */
/**************************************/
struct vector multiply(struct matrix *R,struct vector *v)
{
  double vx,vy,vz;
  double R11,R12,R13,R21,R22,R23,R31,R32,R33;
  struct vector temp;
  
  vx = (*v).x; vy = (*v).y; vz = (*v).z;
  
  R11 = (*R).r11; R12 = (*R).r12; R13 = (*R).r13;
  R21 = (*R).r21; R22 = (*R).r22; R23 = (*R).r23;
  R31 = (*R).r31; R32 = (*R).r32; R33 = (*R).r33;
  
  temp.x = R11*vx + R12*vy + R13*vz;
  temp.y = R21*vx + R22*vy + R23*vz;
  temp.z = R31*vx + R32*vy + R33*vz;
  
  return(temp);
}

/**************************************/
/*                                    */
/*  NxN Matrix on Matrix  multiply    */
/*                                    */
/**************************************/
/*
 * This function takes:
 *
 * int N, size of NxN matrices
 * double *mat1, pointer to mat1
 * double *mat2, pointer to mat2
 * double *mat3, pointer to mat3
 */
void NxNmult(int N,double *mat1,double *mat2,double *mat3)
{
  int k,l,r;
  double sum;

  printf("Test this function before use:\n");
  for(k=0;k<N;k++){      /* row of mat1 */
    for(l=0;l<N;l++){    /* col of mat2 */
      for(sum=0,r=0;r<N;r++){
	sum += mat1[N*k+r] * mat2[N*r+l];
      }
      mat3[N*k+l]=sum;
    }
  }

  return;
}
/**************************************/
/*                                    */
/*  NxM Matrix on N Vector multiply   */
/*                                    */
/**************************************/
void NxMmat_Nvec_mult(int N,int M,double *mat,double *v)
{
  int i,j;
  double sum;
  double *tmp=NULL;

  tmp = (double *)malloc(N*sizeof(double));
  if ( !tmp ){
    printf("Error mallocing for tmp.\n"); exit(0);
  }

  for(i=0;i<N;i++){
    for(sum=0,j=0;j<M;j++){
      sum += mat[M*i+j] * v[j];
    }
    tmp[i]=sum;
  }

  for(i=0;i<N;i++) v[i]=tmp[i];

  free(tmp);
  return;
}

/**************************************/
/*                                    */
/*  Create matrix L lattrep in E3     */
/*                                    */
/**************************************/
struct matrix L_e3(struct cellprm *cell)
{
  /* make sure that calc_cell_basis() has been called before
     calling this function!! */

  struct matrix L;

  calc_cell_basis(cell);

  /*  NOT WHAT WE NEED
  L.r11 = (*cell).bas.A.x;
  L.r12 = (*cell).bas.A.y;
  L.r13 = (*cell).bas.A.z;
  
  L.r21 = (*cell).bas.B.x;
  L.r22 = (*cell).bas.B.y;
  L.r23 = (*cell).bas.B.z;

  L.r31 = (*cell).bas.C.x;
  L.r32 = (*cell).bas.C.y;
  L.r33 = (*cell).bas.C.z;
  */

  /* we really need the transpose of the matrix above */

  L.r11 = (*cell).bas.A.x;
  L.r12 = (*cell).bas.B.x;
  L.r13 = (*cell).bas.C.x;
  
  L.r21 = (*cell).bas.A.y;
  L.r22 = (*cell).bas.B.y;
  L.r23 = (*cell).bas.C.y;

  L.r31 = (*cell).bas.A.z;
  L.r32 = (*cell).bas.B.z;
  L.r33 = (*cell).bas.C.z;

  /*
  printf("L_e3: cell   = %15.10f %15.10f %15.10f\n",(*cell).a,(*cell).b,(*cell).c);
  printf("L_e3: angles = %15.10f %15.10f %15.10f\n",
	 (*cell).alph,
	 (*cell).beta,
	 (*cell).gamm);
  printf("L_e3: L1_ %15.10f%15.10f%15.10f\n",L.r11,L.r12,L.r13);
  printf("L_e3: L2_ %15.10f%15.10f%15.10f\n",L.r21,L.r22,L.r23);
  printf("L_e3: L3_ %15.10f%15.10f%15.10f\n",L.r31,L.r32,L.r33);
  */
  
  return (L);
}

struct matrix Linverse(struct matrix *L)
{
  struct matrix Li;

  Li.r11 = 1.0 / (*L).r11;
  Li.r12 = - ( (*L).r12 / ((*L).r22*(*L).r11) );
  Li.r13 = ((*L).r12*(*L).r23-(*L).r22*(*L).r13) / ((*L).r11*(*L).r22*(*L).r33);

  Li.r21 = 0.0;
  Li.r22 = 1.0 / (*L).r22;
  Li.r23 = - ( (*L).r23 / ((*L).r33*(*L).r22) );

  Li.r31 = 0.0;
  Li.r32 = 0.0;
  Li.r33 = 1.0 / (*L).r33;
  
  return (Li);
}

/* multiplication of a vector by an upper triangular matrix */
struct vector UTMvmult(struct matrix *L, struct vector *v)
{
  struct vector vp;
  vp.x = (*L).r11 * (*v).x + (*L).r12 * (*v).y + (*L).r13 * (*v).z;
  vp.y = (*L).r22 * (*v).y + (*L).r23 * (*v).z;
  vp.z = (*L).r33 * (*v).z;
  return (vp);
}

struct vector cart(struct matrix *L, struct vector *v)
{
  return( UTMvmult(L,v) );
}

/* Ctrans serves two purposes:
   1. if L_old = L_new, and the center moves, it gives the translation
      vector for anion translations in cartesian coordinates. Use this
      for anion translations.
   2. if the center is fixed (fract coords), and L changes, then it
      gives the cartesian 'effective' translation vector for the vertices
      during a lat parm change. Use this for rescaling after lat parm changes.
*/
struct vector Ctrans(struct matrix *L_old,
		     struct matrix *L_new,
		     struct vector *center_old,
		     struct vector *center_new)
{
  static struct vector t;
  static struct vector t1,t2;
  t1 = cart( L_new, center_new );
  t2 = cart( L_old, center_old );
  /* printf("Ctrans: t1: %10.4f%10.4f%10.4f    t2: %10.4f%10.4f%10.4f\n",t1.x,t1.y,t1.z,t2.x,t2.y,t2.z); */
  t = vsub( t1 , t2 );
  return (t);
}

/* vertex_new serves two purposes:
   1. if L_old = L_new, and ctrans is the center translation,
      it gives the new vertex location.
      Use for anion translations.
   2. if L_old != L_new, and ctrans is the effective translation of
      the center of the anion, then we get the new vertex location
      in terms of the new basis vectors.
      Use this for rescaling after lat parm changes.
*/
struct vector vertex_new(struct matrix *L_old,
			 struct matrix *L_new,
			 struct vector *ctrans,
			 struct vector *v_old)
{
  struct vector vcart_new;
  struct vector v_new;
  struct matrix Linv_new;

  vcart_new = vadd( cart( L_old, v_old ), *ctrans );
  Linv_new = Linverse( L_new );
  v_new = UTMvmult( &Linv_new , &vcart_new );
  return (v_new);
}

/**************************************/
/*                                    */
/*  Create Rotation Matrix Forward    */
/*                                    */
/**************************************/
struct matrix matrix_for(struct vector n1,struct vector n2)
{

  /* Note: This requires only two orthogonal vectors, it makes the third
     from the first two by the cross product. */
  
  struct vector e1hat,e2hat,e3hat;
  struct vector n1hat,n2hat,n3hat;
  struct matrix R;

  e1hat.x = 1;   e2hat.x = 0;   e3hat.x = 0; 
  e1hat.y = 0;   e2hat.y = 1;   e3hat.y = 0;
  e1hat.z = 0;   e2hat.z = 0;   e3hat.z = 1;

  /* orientation decided by the input vectors definition of a new frame. */
  /* use grahm-schmidt */
  n1hat = normvec(n1);
  n2hat = normvec( vsub( n2, vsmult( vdotprod(n1,n2), n1 ) ) );
  n3hat = normvec( crossprod(n1hat,n2hat) );  
  
  /* define rotation matrix elements */
  R.r11 = vdotprod(e1hat,n1hat);
  R.r12 = vdotprod(e1hat,n2hat);
  R.r13 = vdotprod(e1hat,n3hat);
  
  R.r21 = vdotprod(e2hat,n1hat);
  R.r22 = vdotprod(e2hat,n2hat);
  R.r23 = vdotprod(e2hat,n3hat);
  
  R.r31 = vdotprod(e3hat,n1hat);
  R.r32 = vdotprod(e3hat,n2hat);
  R.r33 = vdotprod(e3hat,n3hat);
  
  return(R);
}


/*
          This function takes an int, and the pointer to the matrix to
          be inverted.  It returns void.

          This program reduces the given matrix to identity, while at
          the same time performing the same operation to an identity.
          When the given reaches identity, the identity reaches givens
          inverse.
   
          This program only uses partial pivoting (of the rows).
          Full pivoting may be added later along with other
          modifications.
          Remember to free id after you use it in whatever program
          calls this function.
 */
void invt_matrx(int n, double *matrx)
{
    int i, j, k, l, m, p, r, s;
    double *tempm=NULL, *tempi=NULL,*id=NULL;

    tempm = (double *)malloc( n * sizeof( double ) );
    if ( !tempm ) { printf("error: invt_matrx tempm malloc\n"); exit(0); }

    tempi = (double *)malloc( n * sizeof( double ) );
    if ( !tempi ) { printf("error: invt_matrx tempi malloc\n"); exit(0); }

    id = (double *)malloc( n * n * sizeof( double ) );
    if ( !id ) { printf("error: invt_matrx id malloc\n"); exit(0); }

    /* create nxn identity */
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    if (i == j) {
		id[n * i + j] = 1;
	    } else {
		id[n * i + j] = 0;
	    }
	}
    }

    for (l = 0; l < n; l++) {
	k = l;

	/* in case diagonal element is 0 */
	if (matrx[n * k + l] == 0) {
	    for (r = 0; r < n; r++) {
		if (matrx[n * r + l] != 0) {
		    break;
		}
	    }
	    for (s = 0; s < n; s++) {
		matrx[n * k + s] += matrx[n * r + s];
		id[n * k + s] += id[n * r + s];
	    }
	}
	/* get a 1 in row l, column l */
	for (m = 0; m < n; m++) {
	    tempm[m] = matrx[n * k + m] / matrx[n * k + l];
	    tempi[m] = id[n * k + m] / matrx[n * k + l];
	}
	for (m = 0; m < n; m++) {
	    matrx[n * k + m] = tempm[m];
	    id[n * k + m] = tempi[m];
	}

	/* eliminate other elements in lth column */
	for (i = 0; i < n; i++) {
	    if (i == k && k == n - 1) {
		break;
	    }
	    if (i == k && k != n - 1) {
		i++;
	    }
	    /* mult kth row by ith column elem */
	    for (p = 0; p < n; p++) {
		tempm[p] = matrx[n * k + p] * matrx[n * i + l];
		tempi[p] = id[n * k + p] * matrx[n * i + l];
	    }

	    /* perform the row op */
	    for (p = 0; p < n; p++) {
		matrx[n * i + p] -= tempm[p];
		id[n * i + p] -= tempi[p];
	    }
	}
    }

    /* copy the new matrix back to the old positions and free memory */
    for (i=0; i<n*n; i++) matrx[i] = id[i];

    free(tempi);
    free(tempm);
    free(id);
    return;
}
