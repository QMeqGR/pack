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


extern int non_per_ewald;
extern int debug;
extern int numat;
extern double rep_epsilon;

extern struct cellprm cell;
extern struct atom *p;
extern struct atom *s;

void center_np()
{

  int i;
  int cushion_violation=0;
  static int rejection_count=0;
  struct vector temp,temp2;
  struct Energy Etmp={0};

  temp = geometric_center();
  temp2 = center_of_mass();

  if ( debug > 3 ){
    print_poscar_out(0);
    Etmp = (*Total_Energy)();
    printf("* center_np: (before) geom cnt = %7.3f%7.3f%7.3f"
	   "  cm = %7.3f%7.3f%7.3f",
	   temp.x,temp.y,temp.z,temp2.x,temp2.y,temp2.z);
    printf(" E.tot= %f\n",Etmp.tot);
    fflush(stdout);
  }

  temp = vsub( center_of_mass(), makevec(0.5,0.5,0.5) );

  /* check to see if centering the cluster will move
     any atoms outside the wall cushion */
  for(i=0;i<numat;i++){
    if ( 0==trans_accept( -1, vadd(makevec((*(s+i)).x,(*(s+i)).y,(*(s+i)).z),temp), non_per_ewald ) ) cushion_violation=1;
  }
  if ( cushion_violation ){
    rejection_count++;
    if ( debug > -2 ) printf("* WARN: center_np: cushion_violation  -- reject centering.  rejection_count= %d.\n",rejection_count);
  }
  
  if ( !cushion_violation ) trans_all( temp );
  
  if ( debug > 3 ){
    temp = geometric_center();
    temp2 = center_of_mass();

    print_poscar_out(0);
    Etmp = (*Total_Energy)();
    printf("* center_np: (after) geom cnt = %7.3f%7.3f%7.3f"
	   "  cm = %7.3f%7.3f%7.3f",
	   temp.x,temp.y,temp.z,temp2.x,temp2.y,temp2.z);
    printf(" E.tot= %f\n",Etmp.tot);
  }
  
  return;
}


void recenter()
{
  struct Energy Etmp={0},Etmp2={0};

  /* RECENTER */
  if ( debug > 3 ) printf("* --> recentering cluster <--\n");
  Etmp = (*Total_Energy)();
  center_np();
  Etmp2 = (*Total_Energy)();
  if ( debug > -2 ){
    if ( !tol(Etmp.tot,Etmp2.tot,1e-5) ){
      printf("* recenter: WARN: Re-centering of nanoparticle changed the energy!!\n");
      printf("* recenter: WARN: --> before tot= %f, after tot= %f\n",Etmp.tot,Etmp2.tot);
      printf("* recenter: WARN: Probable wall violation - it may remove itself.\n");
      printf("* recenter: WARN: Checing for wall violation\n");
    }
  }
  if ( check_wall_violation() > 0 ){
    printf("* ERROR: wall violation!\n");
    printf("* recenter: check_wall_violation = %d. 1=cation,  10=anion \n",check_wall_violation() );
  }

  return;  
}
