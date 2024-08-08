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
extern int done_flag,get_estat,debug;
extern int n_simp,ecol;

extern double d,rtol;
extern double simp_errlim,simp_lamb;
extern double **sim_pp,*sim_p,*sim_y;

extern struct cellprm cell;

  /********************************************************************/
  /*                                                                  */
  /*                    Simplex Relaxation Routine                    */
  /*                                                                  */
  /********************************************************************/
void Simplex(int Nrestarts, int Nmoves)
{

  int i;

  /*
    Note:  This routine works with the code as follows.

    - place cell parms, cations, and anion center positions in array fn
    - set orientation variables to zero for all anions.  They still retain
      their orientations from the end of the metro runs, but this orientation
      will be the starting orientation for the simplex relaxation.
    - fill vertex points of the simplex (create matrix sim_p and calculate
      the energy of each vertex configuration for the array sim_y)
    - start the simplex routine
    - for each move of the simplex, it will move vertex points, corresponding
      to lattie parameter changes and moves of the cations or anions, and
      also rotations of the anions.  The rotations go from (0 0 0) to
      (phi the psi).
    - update the structs tetr, octa, at, cell to reflect these changes
    - don't update fn (it holds all current positions, therefore no need)
    - when doing subsequent rotations in the simplex routine, must first
      rotate back to phi=0 the=0 psi=0 (the starting config) and then to
      the new set.  Rotations go from (phi the psi) to (0 0 0) to
      (new_phi new_the new_psi).
  */

    /* variables to set before beginning simplex run */

    done_flag=0;  /* reset the done_flag from the end of the Metro runs, use this to exit with SIGHUP */
    get_estat = 1; /* still want to calculate the electrostatic energy if -E flag was set */

    if ( debug > -1 ){
      fflush(stdout);
      printf("****************************************************************\n*\n");
      printf("*\n*        - Simplex Minimiziation Routine -\n");
      printf("*\n");
    }
    if ( debug > -1 ) {
      printf("*\n");
      printf("* %20s%10d%25s\n","n_simp  =",n_simp,"D.O.F.");
      printf("* %20s%10d%25s\n","simp_nmax  =",Nmoves,"max moves");
      printf("* %20s%10.2e\n","simp_errlim  =",simp_errlim);
      printf("* %20s%10.2e\n","simp_lamb  =",simp_lamb);
      printf("*\n");
    }
    for ( i=0; i<Nrestarts; i++) {
      if ( done_flag ) break;
      if ( i != 0 && debug > -1 ) printf("* Simplex restart no. %d of %d\n",i,Nrestarts-1);
      simplex_init();
      
      if ( debug > -1 ) printf("*\n");
      if ( ! ecol && debug > -1 )
	printf("%2s%10s%15s%10s%10s%10s%18s%18s%18s\n","* ","nfunk","rtol","pf","ar","ortho","e_ecc","e_ss","e_tot");
      fflush(stdout);
      
      amoeba_ehm(sim_pp,sim_y,n_simp,simp_errlim,Total_Energy,&Nmoves);
      if ( debug > -1 ){
	printf("*\n");
	printf("* %20s%10d%25s\n","simp moves  =",Nmoves,"simp moves");
	printf("* %20s%10.2e%25s\n","rtol  =",rtol,"final rtol");
	printf("*\n");
      }
    }
    if ( debug > -1 ) {
      /* printf("*********************************************\n"); */
      printf("*\n*          - End Simplex -\n");
      /* print_info_block(); */
    }
    if ( debug > 3 ) print_c2v( );

  return;
}

