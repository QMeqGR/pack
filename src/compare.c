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

extern char tempfile[25];
extern struct Energy E_best;
extern FILE *restartinput;

void best_compare(void)
{
  struct Energy E_post_runs;
  double ecc_diff,eng_diff,mag_ecc_diff,mag_eng_diff;

  /* check to see if the last metro or simplex structure is better than the 'best' stored */
  E_post_runs = Total_Energy();
  ecc_diff = E_best.ecc - E_post_runs.ecc;
  eng_diff = E_best.eng - E_post_runs.eng;
  mag_ecc_diff = fabs( ecc_diff );
  mag_eng_diff = fabs( eng_diff );

  printf("*\n****************************************************************\n*\n");
  printf("* %20s%12.4e\n","best ecc  =",E_best.ecc);
  printf("* %20s%12.4e\n","best eng  =",E_best.eng);
  printf("* %20s%12.4e\n","ecc+eng  =",E_best.ecc_eng);
  printf("* %20s%12.4e\n","best vol  =",E_best.vol);
  printf("* %20s%12.4e\n","best pf  =",E_best.pf);
  printf("*\n");
  printf("* %20s%12.4e\n","post_runs ecc  =",E_post_runs.ecc);
  printf("* %20s%12.4e\n","post_runs eng  =",E_post_runs.eng);
  printf("* %20s%12.4e\n","ecc+eng  =",E_post_runs.ecc_eng);
  printf("* %20s%12.4e\n","post_runs vol  =",E_post_runs.vol);
  printf("* %20s%12.4e\n","post_runs pf  =",E_post_runs.pf);
  printf("*\n");
  printf("* %20s%12.4e\n","dV / V_runs  =",fabs( E_best.vol - E_post_runs.vol ) / E_post_runs.vol);
  printf("*\n");
  printf("* %20s%12.4e\n","ecc_diff  =", ecc_diff );
  printf("* %20s%12.4e\n","eng_diff  =", eng_diff );
  printf("*\n");
  
  /* We want to avoid the situation where we swap the structures when
   * the total energy of the 'best' structure is lower than the last
   * metro output, but only because of some difference in the repulsive energy.
   * The repulsive energy is quickly taken care of in the simplex routine.
   * Therefore, we want the difference in the electrostatic energy to be
   * larger than the difference in the repulsive energy, and that this cause
   * the total energy to be lower.  Then we are picking the structure with
   * the lowest electrostatic energy, which is what we want.
   */

  /* best ecc+eng is lower, best and post are not the same
   * and magnitude of the ecc diff is bigger than mag of eng diff,
   * and the volume change is within BEST_VOL.
   * (or all the above and the best volume is simply lower)
   * Then swap!
   */
  if (

      ( E_best.ecc_eng < E_post_runs.ecc_eng &&
	!tol(E_best.ecc_eng,E_post_runs.ecc_eng,1e-6) &&
	mag_ecc_diff > mag_eng_diff  &&
	( fabs( E_best.vol - E_post_runs.vol ) / E_post_runs.vol < BEST_VOL ) )
      
      ||

      ( E_best.ecc_eng < E_post_runs.ecc_eng &&
	!tol(E_best.ecc_eng,E_post_runs.ecc_eng,1e-6) &&
	mag_ecc_diff > mag_eng_diff  &&
	E_best.vol < E_post_runs.vol )
	
	) {

    printf("*\n****************************************************************\n*\n");
    printf("*  mag_ecc_diff > mag_eng_diff AND within BEST_VOL AND\n");
    printf("* 'best' has lower energy than post runs! --> swapping.\n");
    
    /* Open restart file and get initialization conditions */
    restartinput = fopen(tempfile,"r");
    if ( restartinput == NULL )
      {
	printf("* Error opening 'best' file\n");
	printf("* %20s%35s\n","tmpfile =",tempfile);
	exit_pack(0);
      }
    get_init_parms(restartinput,0);    
    fclose(restartinput);
    printf("*\n****************************************************************\n*\n");
    printf("*\n*          - Parameters at Post 'best-swap' -\n");
    print_info_block();
    
  } else if ( E_best.ecc_eng < E_post_runs.ecc_eng &&
	      !tol(E_best.ecc_eng,E_post_runs.ecc_eng,1e-6) &&
	      mag_ecc_diff < mag_eng_diff &&
	      E_best.eng > E_post_runs.eng) {
    printf("* 'best' has lower energy than post runs!\n");
    printf("*  BUT mag_ecc_diff < mag_eng_diff  ---> not swapping.\n");
    printf("*\n****************************************************************\n*\n");
    
  } else if ( tol(E_best.ecc_eng,E_post_runs.ecc_eng,1e-6) ) {

    printf("* post runs has same energy as 'best' --> not swapping.\n");
    printf("*\n****************************************************************\n*\n");
    
  } else {
    
    printf("* volume difference is unsatisfactory OR\n");
    printf("* post runs has lower energy than 'best' --> not swapping.\n");
    printf("*\n****************************************************************\n*\n");
    
  }


  return;
}
