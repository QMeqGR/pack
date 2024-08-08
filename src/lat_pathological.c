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
/*    LAT_PATHOLOGICAL                */
/*                                    */
/*    Check for DPLANE violation      */
/*    Check for sqrt( neg number )    */
/*                                    */
/**************************************/
int lat_pathological(struct cellprm *cell,
		     struct vector lat_chng,
		     struct vector ang_chng,
		     double dplane)
{
  static int pathological_flag=0;
  static int dplane_violation=0;
  double vol=0;
  double dplane_ab=0,dplane_ac=0,dplane_bc=0;
  struct matrix L_old,L_new;

  /* make the changes to structure cell here */
  L_old = L_e3(cell);
  (*cell).a += lat_chng.x;
  (*cell).b += lat_chng.y;
  (*cell).c += lat_chng.z;
  (*cell).alph += ang_chng.x;
  (*cell).beta += ang_chng.y;
  (*cell).gamm += ang_chng.z;
  L_new = L_e3(cell);

  /* reject if the distance between any two planes of the cell becomes less than DPLANE*d */
  vol=cell_volume(cell); dplane_violation=0;
  dplane_ab = vol / ( (*cell).a * (*cell).b );
  dplane_ac = vol / ( (*cell).a * (*cell).c );
  dplane_bc = vol / ( (*cell).b * (*cell).c );
  if ( dplane_ab <= dplane ) dplane_violation=2;
  if ( dplane_ac <= dplane ) dplane_violation=2;
  if ( dplane_bc <= dplane ) dplane_violation=2;

  /* check for sqrt( neg number ) in cell.bas.C.z */
  pathological_flag = calc_cell_basis(cell);

  /* undo the changes to structure cell here */
  L_old = L_e3(cell);
  (*cell).a -= lat_chng.x;
  (*cell).b -= lat_chng.y;
  (*cell).c -= lat_chng.z;
  (*cell).alph -= ang_chng.x;
  (*cell).beta -= ang_chng.y;
  (*cell).gamm -= ang_chng.z;
  L_new = L_e3(cell);

  if ( dplane_violation || pathological_flag ) return (dplane_violation + pathological_flag);
  else return 0;

}
