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

void np_force_mv(void)
{

  extern int debug;
  int iii,moves=CEN_MOVES;
  struct vector temp,temp2;

  if ( debug>3 ){
    printf("****************************************************************\n");
    printf("*BEG*F-S********************************************************\n");
    printf("****************************************************************\n");
    printf("* Entering Initial Forces Routine with %d centering moves.\n",moves);
  }
  
  for(iii=0;iii<moves;iii++){
    
    if ( debug>3 && iii%5 == 0 ) {
      temp = geometric_center();
      temp2 = center_of_mass();
      printf("* geom cnt [%4d] = %10.4f%10.4f%10.4f"
	     "  cm = %10.4f%10.4f%10.4f\n",
	     iii,temp.x,temp.y,temp.z,temp2.x,temp2.y,temp2.z);
      fflush(stdout);
    }
    
    force_step();
    temp = geometric_center();
    temp2 = center_of_mass();
    temp = vsub( center_of_mass(), makevec(0.5,0.5,0.5) );
    trans_all( temp );
    
  }
  if ( debug>3 ){
    printf("****************************************************************\n");
    printf("*END*F-S********************************************************\n");
    printf("****************************************************************\n");
  }
  
  return;
}
