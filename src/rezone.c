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
/* Re-zone an atom                    */
/*                                    */
/**************************************/
struct vector rezone(struct vector v)
{

  if ( v.x > 1.0 ) v.x -= 1.0;
  if ( v.x < 0.0 ) v.x += 1.0;
  if ( v.x < -1.0 || v.x > 2.0 ) {
    printf("rezone out of range! v.x = %f\n",v.x);
    kill(getpid(),SIGINT);
  }

  if ( v.y > 1.0 ) v.y -= 1.0;
  if ( v.y < 0.0 ) v.y += 1.0;
  if ( v.y < -1.0 || v.y > 2.0 ) {
    printf("rezone out of range! v.y = %f\n",v.y);
    kill(getpid(),SIGINT);
  }

  if ( v.z > 1.0 ) v.z -= 1.0;
  if ( v.z < 0.0 ) v.z += 1.0;
  if ( v.z < -1.0 || v.z > 2.0 ) {
    printf("rezone out of range! v.z = %f\n",v.z);
    kill(getpid(),SIGINT);
  }

  return(v);
}
