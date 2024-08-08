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


extern int numat;

extern struct obj *object;

/**************************************/
/*                                    */
/* Re-scale after a lat prm change    */
/*                                    */
/**************************************/
void rescale_objects(struct cellprm *cell,
		     struct vector lat_chng,
		     struct vector ang_chng)
{
  int i,j,nv;

  struct vector ctrans;
  struct matrix L_old,L_new;

  L_old = L_e3(cell);

  (*cell).a += lat_chng.x;
  (*cell).b += lat_chng.y;
  (*cell).c += lat_chng.z;
  (*cell).alph += ang_chng.x;
  (*cell).beta += ang_chng.y;
  (*cell).gamm += ang_chng.z;

  L_new = L_e3(cell);
  /* must do this before *return*ing */

  /* rescale should be called for any lattice param change */
  /* flag all atoms 'moved' for energy calc */
  mark_all_objects_moved();

  for (i=0; i<MAX_OBJS && object[i].used==1; i++){
    /* fract coord rep stays the same for anion centers */
    /* now go over j=1 to all verts and rescale */
    for (j=1; j<MAX_VERTS && object[i].vused[j]==1; j++) {
      ctrans = Ctrans( &L_old, &L_new, &(*(object+i)).v[0], &(*(object+i)).v[0] );
      (*(object+i)).v[j] = vertex_new( &L_old, &L_new, &ctrans, &(*(object+i)).v[j] );
    }
    
  }

  return;
}

