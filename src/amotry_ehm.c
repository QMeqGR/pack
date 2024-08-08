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


#define NRANSI
#define float double
#include <stdio.h>
#include <stdlib.h>
#include "packlib.h"

double amotry_ehm(double **p, double y[], double psum[], int ndim,
		  struct Energy (*funk)(), int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry;
  extern double *ptry;
  struct Energy E;

  fac1=(1.0-fac)/ndim;
  fac2=fac1-fac;
  for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

  simplex_trans(ptry);

  E = (*funk)(); /* arg is ptry */
  ytry = E.tot;
  if (ytry < y[ihi]) {
    y[ihi]=ytry;
    for (j=0;j<ndim;j++) {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j]=ptry[j];
    }
  }
  return ytry;
}
#undef NRANSI



