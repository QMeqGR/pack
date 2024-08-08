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


/*
 * Packlib library
 * Eric Majzoub
 * Sandia National Laboratories
 * Livermore, CA
 * 10 October 2001
 *

 May 8 2011
 - Added RandomChoice function

 14 feb 2007
 -fixed NxM on vec mult bug; sum was initialized as int, needed double.

 22 jan 2006
 -begin coding for non-ortho lattices

 21 jan 2006
 -added dimer type functions for pack.

 */

#include <math.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <signal.h>
#include "global_defs.h"
#include "packlib.h"

extern int process_num, process_rank;

/*
 * Tol: tolerance function.
 * Returns: 1 -- if inside tolerance
 *          0 -- outside tolerance
 */
int tol(double a,double b,double t)
{
  if ( fabs(a-b) < t ) return (1);
  else return(0);
}

/*
 * Random sign function.
 * Returns +1 or -1 randomly
 *
 */
int rsign(void)
{
  if ( drand48() < 0.5 ) { return(1); }
  else return (-1);
}

/* dmod: double modulus function */
double dmod( double a, double b)
{

  if ( a<0.0 ) a = -a;
  if ( b<0.0 ) b = -b;

  while( (a-b)>b ){
    a -= b;
  }

  if ( (a-b) < 0.0 ) return(a);
  else if ( (a-b) > 0.0 ) return( (a-b) );
  else return(-1);

}

/************************************/
/*                                  */
/*        Complex Number            */
/*         Functions                */
/*                                  */
/************************************/

/*
 * cmag: magnitude
 *
 */
double cmag(struct cmplx a)
{
  return ( a.real*a.real + a.imag*a.imag );
}

/*
 * cmult: complex multiply
 *
 */
struct cmplx cmult(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real*b.real - a.imag*b.imag;
  c.imag = a.real*b.imag + a.imag*b.real;
  return(c);
}

/*
 * cadd: complex add
 *
 */
struct cmplx cadd(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real + b.real;
  c.imag = a.imag + b.imag;
  return(c);
}

/*
 * cadd: complex sub
 *
 */
struct cmplx csub(struct cmplx a,struct cmplx b)
{
  struct cmplx c;
  c.real = a.real - b.real;
  c.imag = a.imag - b.imag;
  return(c);
}

/*****************************/
/*                           */
/*     Vector                */
/*     Operations            */
/*                           */
/*****************************/


/*
 * makevec: make a vector out of doubles
 */
struct vector makevec(double a, double b, double c)
{
  struct vector v;

  v.x = a;
  v.y = b;
  v.z = c;

  return(v);
}

/*
 * vmag: vector magnitude
 */
double vmag(struct vector a)
{
  return ( sqrt(a.x*a.x + a.y*a.y + a.z*a.z) );
}
/*
 * vmagsq: square of the vector magnitude
 */
double vmagsq(struct vector a)
{
  return ( a.x*a.x + a.y*a.y + a.z*a.z );
}


/**************************************/
/*                                    */
/*            vsub                    */
/*                                    */
/**************************************/
struct vector vsub(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  return(c);
  }

/**************************************/
/*                                    */
/*        vadd                        */
/*                                    */
/**************************************/
struct vector vadd(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  return(c);
  }

/*
 * vector scalar multiply
 * vsmult
 */
struct vector vsmult(double a, struct vector b)
{
  struct vector c;
  
  c.x = a * b.x;
  c.y = a * b.y;
  c.z = a * b.z;
  return(c);
}

/*
 * vdotprod
 */
double vdotprod(struct vector a, struct vector b)
{
  return( a.x * b.x + a.y * b.y + a.z * b.z );
}

/*
 * vcross
 *
 */
struct vector crossprod(struct vector a, struct vector b)
  {
  struct vector c;

  c.x = a.y*b.z - a.z*b.y;
  c.y = a.z*b.x - a.x*b.z;
  c.z = a.x*b.y - a.y*b.x;
  return(c);
  }

/*
 * normvec
 *
 */
struct vector normvec(struct vector a)
  {
  double d;
  struct vector c;

  d = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  c.x = a.x/d;
  c.y = a.y/d;
  c.z = a.z/d;

  return(c);
  }

/*
 * vtriple
 */
double vtriple(struct vector a, struct vector b, struct vector c)
{
  double v;
  v = fabs( vdotprod( c, crossprod( a, b ) ) );
  return(v);
}


/*
 * is cubic cell?
 */
int is_cubic(void)
{

  extern struct cellprm cell;
  
  if ( tol( cell.bas.A.x, cell.bas.B.y, 1e-8 ) &&
       tol( cell.bas.A.x, cell.bas.C.z, 1e-8 ) &&
       tol( cell.bas.A.y, 0.0, 1e-8 ) &&
       tol( cell.bas.A.z, 0.0, 1e-8 ) &&
       tol( cell.bas.B.x, 0.0, 1e-8 ) &&
       tol( cell.bas.B.z, 0.0, 1e-8 ) &&
       tol( cell.bas.C.x, 0.0, 1e-8 ) &&
       tol( cell.bas.C.y, 0.0, 1e-8 ) ) return (1);
  else return(0);
	      
}


/*
 * molwt
 *
 */
double molwt(int Z)
{

struct table {
  char *element;
  double weight;
  char *latexsymbol;
  int Z;
  double radius;
  double R,G,B;
} PeriodicTable[] = {
  {"Du",0.0000,"Dummy",0,0,0,0,0},
  {"H" ,1.0080,"H",1,0.30,0.40,0.36,0.36},
  {"He",4.0026,"He",2,0,0,0,0},
  {"Li",6.941,"Li",3,0,0,0,0},
  {"Be",9.0122,"Be",4,0,0,0,0},
  {"B",10.81,"B",5,0,0,0,0},
  {"C",12.01,"C",6,0.40,0.5,0.5,0},
  {"N",14.007,"N",7,0,0,0,0},
  {"O",15.9994,"O",8,0,0,0,0},
  {"F",18.998,"F",9,0.45,0.1,0.2,0.7},
  {"Ne",20.18,"Ne",10,0,0,0,0},
  {"Na",22.9898,"Na",11,0.60,0.90,0.90,0.00},
  {"Mg",24.305,"Mg",12,0,0,0,0},
  {"Al",26.98154,"Al",13,0.45,0.10,0.90,0.00},
  {"Si",28.0885,"Si",14,0,0,0,0},
  {"P",30.974,"P",15,0,0,0,0},
  {"S",32.064,"S",16,0,0,0,0},
  {"Cl",35.453,"Cl",17,0.5,0.2,0.6,0.2},
  {"Ar",39.948,"Ar",18,0,0,0,0},
  {"K",39.09,"K",19,0,0,0,0},
  {"Ca",40.08,"Ca",20,0,0,0,0},
  {"Sc",44.9559,"Sc",21,0,0,0,0},
  {"Ti",47.90,"Ti",22,0,0,0,0},
  {"V",50.942,"V",23,0,0,0,0},
  {"Cr",51.996,"Cr",24,0,0,0,0},
  {"Mn",54.9380,"Mn",25,0,0,0,0},
  {"Fe",55.847,"Fe",26,0,0,0,0},
  {"Co",58.9332,"Co",27,0,0,0,0},
  {"Ni",58.71,"Ni",28,0,0,0,0},
  {"Cu",63.546,"Cu",29,0,0,0,0},
  {"Zn",65.38,"Zn",30,0,0,0,0},
  {"Ga",69.72,"Ga",31,0,0,0,0},
  {"Ge",72.59,"Ge",32,0,0,0,0},
  {"As",74.922,"As",33,0,0,0,0},
  {"Se",78.96,"Se",34,0,0,0,0},
  {"Br",79.91,"Br",35,0,0,0,0},
  {"Kr",83.80,"Kr",36,0,0,0,0},
  {"Rb",85.47,"Rb",37,0,0,0,0},
  {"Sr",87.62,"Sr",38,0,0,0,0},
  {"Y",88.906,"Y",39,0,0,0,0},
  {"Zr",91.22,"Zr",40,0,0,0,0},
  {"Nb",92.9064,"Nb",41,0,0,0,0},
  {"Mo",95.94, "Mo",42,0,0,0,0},
  {"Tc",98.91, "Tc",43,0,0,0,0},
  {"Ru",101.07,"Ru",44,0,0,0,0},
  {"Rh",102.90,"Rh",45,0,0,0,0},
  {"Pd",106.4, "Pd",46,0,0,0,0},
  {"Ag",107.87,"Ag",47,0,0,0,0},
  {"Cd",112.40,"Cd",48,0,0,0,0},
  {"In",114.82,"In",49,0,0,0,0},
  {"Sn",118.69,"Sn",50,0,0,0,0},
  {"Sb",121.75,"Sb",51,0,0,0,0},
  {"Te",127.60,"Te",52,0,0,0,0},
  {"I", 126.90,"I", 53,0,0,0,0},
  {"Xe",131.30,"Xe",54,0,0,0,0},
  {"Cs",132.91,"Cs",55,0,0,0,0},
  {"Ba",137.34,"Ba",56,0,0,0,0},
  {"La",138.91,"La",57,0,0,0,0},
  {"Ce",140.12,"Ce",58,0,0,0,0},
  {"Pr",140.91,"Pr",59,0,0,0,0},
  {"Nd",144.24,"Nd",60,0,0,0,0},
  {"Pm",145.00,"Pm",61,0,0,0,0},
  {"Sm",150.35,"Sm",62,0,0,0,0},
  {"Eu",151.96,"Eu",63,0,0,0,0},
  {"Gd",157.25,"Gd",64,0,0,0,0},
  {"Tb",158.92,"Tb",65,0,0,0,0},
  {"Dy",162.50,"Dy",66,0,0,0,0},
  {"Ho",164.93,"Ho",67,0,0,0,0},
  {"Er",167.26,"Er",68,0,0,0,0},
  {"Tm",168.93,"Tm",69,0,0,0,0},
  {"Yb",173.04,"Yb",70,0,0,0,0},
  {"Lu",174.97,"Lu",71,0,0,0,0},
  {"Hf",178.49,"Hf",72,0,0,0,0},
  {"Ta",180.95,"Ta",73,0,0,0,0},
  {"W", 183.85,"W", 74,0,0,0,0},
  {"Re",186.2, "Re",75,0,0,0,0},
  {"Os",190.20,"Os",76,0,0,0,0},
  {"Ir",192.22,"Ir",77,0,0,0,0},
  {"Pt",195.09,"Pt",78,0,0,0,0},
  {"Au",196.97,"Au",79,0,0,0,0},
  {"Hg",200.59,"Hg",80,0,0,0,0},
  {"Tl",204.37,"Tl",81,0,0,0,0},
  {"Pb",207.2, "Pb",82,0,0,0,0},
  {"Bi",208.98,"Bi",83,0,0,0,0},
  {"Po",210.00,"Po",84,0,0,0,0},
  {"At",210.00,"At",85,0,0,0,0},
  {"Rn",222.00,"Rn",86,0,0,0,0},
  {"Fr",223.00,"Fr",87,0,0,0,0},
  {"Ra",226.02,"Ra",88,0,0,0,0},
  {"Ac",227.00,"Ac",89,0,0,0,0},
  {"Th",232.04,"Th",90,0,0,0,0},
  {"Pa",231.00,"Pa",91,0,0,0,0},
  {"U", 238.03,"U", 92,0,0,0,0},
  {"Np",237.0, "Np",93,0,0,0,0},
  {"Pu",244.0, "Pu",94,0,0,0,0},
  {"Am",243.0, "Am",95,0,0,0,0}
};


  return PeriodicTable[Z].weight;
}


void kr_reverse(char s[])
{
  int c, i, j;

  for(i=0, j=strlen(s)-1; i<j; i++, j--){
    c=s[i];
    s[i]=s[j];
    s[j]=c;
  }
  return;
}

void kr_itoa(int n, char s[])
{
  int i, sign;
  
  if ((sign = n) < 0)
    n = -n;
  i=0;
  do {
    s[i++] = n % 10 + '0';
  } while ((n /= 10) > 0);
  if ( sign < 0 )
    s[i++] = '-';
  s[i] = '\0';
  kr_reverse(s);

  return;
}

int RandomChoice(double * prob, int n){
	int i;
	double p=drand48();
	for (i=0; i<n; i++) {
		if ((p<=prob[n+i])&&(prob[i]>0)) {
			return(i);
		}
	}
	if (i>=n) {
		printf("error finding a random choice out of prob:");
		for (i=0; i<n; i++) printf(" %8f",prob[i]);
		printf("\nexiting\n");
		exit(0);
	}
}
