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

extern int debug,numat;


/************************************************************/
/*                                                          */
/*                  Pressure Calculation                    */
/*                                                          */
/************************************************************/
struct vector Pressure(void)
{
  int i;

  double da,db,dc,vol;
  double Ei,Ef,dE,dV=PRESS_DVOL;
  double Px=0,Py=0,Pz=0;

  struct vector cell_chng,ang_chng;
  struct vector P;
  struct Energy E;

  extern struct cellprm cell;

  E = Total_Energy();
  Ei = E.ecc_eng;

  ang_chng = makevec(0.0,0.0,0.0);
  vol = cell_volume(&cell);

  for(i=0; i<3; i++){

    if ( debug > 4 ) printf("* (dE/dV): changing coord %d\n",i);
    
    da=0.0; db=0.0; dc=0.0;

    if ( i==0 ) da = dV * cell.a / vol;
    if ( i==1 ) db = dV * cell.b / vol;
    if ( i==2 ) dc = dV * cell.c / vol;

    if ( debug > 5 ) {
      printf("da=%20.19f\n",da);
      printf("db=%20.19f\n",db);
      printf("dc=%20.19f\n",dc);
    }

    cell_chng = makevec(da,db,dc);
    rescale_objects( &cell, cell_chng, ang_chng);

    E = Total_Energy();
    Ef = E.ecc_eng;
    dE = Ef-Ei;

    if ( i==0 ) { Px = dE/dV; }
    if ( i==1 ) { Py = dE/dV; }
    if ( i==2 ) { Pz = dE/dV; }

    /* return things to normal */
    cell_chng = makevec(-da,-db,-dc);
    rescale_objects( &cell, cell_chng, ang_chng);
  }

  P = makevec(Px,Py,Pz);
  return P;
}

