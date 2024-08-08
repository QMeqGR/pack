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


/**************************************/
/*                                    */
/*      Cell functions                */
/*                                    */
/**************************************/
int calc_cell_basis(struct cellprm *cell){

  static double temp;

  /*
  printf("calc_cell_basis: a b c : %15.10f%15.10f%15.10f\n",
	 (*cell).a,
	 (*cell).b,
	 (*cell).c);
  printf("calc_cell_basis: angles: %15.10f%15.10f%15.10f\n",
	 (*cell).alph,
	 (*cell).beta,
	 (*cell).gamm);
  */

  (*cell).bas.A.x = (*cell).a;
  (*cell).bas.A.y = 0;
  (*cell).bas.A.z = 0;
  
  (*cell).bas.B.x = (*cell).b * cos((*cell).gamm);
  (*cell).bas.B.y = (*cell).b * sin((*cell).gamm);
  (*cell).bas.B.z = 0;

  (*cell).bas.C.x = (*cell).c * cos((*cell).beta);
  (*cell).bas.C.y = (*cell).c / sin((*cell).gamm) * ( cos((*cell).alph) - cos((*cell).gamm)*cos((*cell).beta)  );
  temp = ( (*cell).c*(*cell).c - (*cell).bas.C.x*(*cell).bas.C.x - (*cell).bas.C.y*(*cell).bas.C.y );
  (*cell).bas.C.z = sqrt( temp );

  if ( temp < 0.0 ) return 1;
  else return 0;

}
/* cell volume */
double cell_volume(struct cellprm *cell){
  calc_cell_basis(cell);
  return( vtriple((*cell).bas.A, (*cell).bas.B, (*cell).bas.C) );
}

/*************************************/
/* calculate and return aspect ratio */
/*************************************/
double aspect(void)
{
  extern struct cellprm cell;
  return fmax(cell.a,fmax(cell.b,cell.c))/fmin(cell.a,fmin(cell.b,cell.c));
}


/*********** surface area of cell *************/
double surface_area(void)
{
  extern struct cellprm cell;
  return(2*(vmag(crossprod(cell.bas.A,cell.bas.B))+
	    vmag(crossprod(cell.bas.A,cell.bas.C))+
	    vmag(crossprod(cell.bas.B,cell.bas.C))));
}
