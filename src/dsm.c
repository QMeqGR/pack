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

extern int debug;
extern struct cellprm cell;

/******************************************************************/
/*                                                                */
/*                                                                */
/*                Distance Scaling Method Code                    */
/*                                                                */
/*                                                                */
/******************************************************************/
/*
  -Scale factor will be given as a percent to expand
  -Integer argument +-1 will indicate to expand or contract

*/
void scale_cell(
		double sf,
		int updown,
		int scale_cell)
{
  double da=0,db=0,dc=0;
  double dvol=0;
  struct vector cell_chng,ang_chng;
  static struct vector cell_save;

  if ( scale_cell == 0 ) return;

  /* angles won't be changed here, just scale the lattice vectors */
  ang_chng = makevec(0,0,0);

  if (updown == 1){
    cell_save = makevec(cell.a,cell.b,cell.c);
    da = (sf - 1.0) * cell.a;
    db = (sf - 1.0) * cell.b;
    dc = (sf - 1.0) * cell.c;
  }
  else if (updown == -1){
    da = cell_save.x - cell.a;
    db = cell_save.y - cell.b;
    dc = cell_save.z - cell.c;
  }
  else{
    printf("* FUNCTION: scale_cell\n");
    printf("* !!!!!!!!!! ERROR\n");
    printf("* integer arg must be +1 or -1!\n");
  }

  cell_chng = makevec(da,db,dc);

  if ( debug > 3 ) {
    printf("* FUNCTION: scale_cell\n");
    printf("*  updown  = %d      (1=up , -1=down)\n",updown);
    printf("*  scale_factor = %f\n",sf);
    printf("*  a  b  c = %12.6f%12.6f%12.6f\n",cell.a,cell.b,cell.c);
    printf("* da db dc = %12.6f%12.6f%12.6f\n",da,db,dc);
    printf("*  angles  = %12.6f%12.6f%12.6f\n",cell.alph,cell.beta,cell.gamm);
    printf("* dvol = %12.6f    cell_vol = %12.6f     ratio = %12.6f\n",
	   dvol,cell_volume(&cell),dvol/cell_volume(&cell));
  }
  
  /* rescale the object locations */
  /* the rescaling needs the pre-changed cell parameters and the change params */
  rescale_objects( &cell, cell_chng, ang_chng);
  
  return;
}
