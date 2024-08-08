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

#define AA_PCNT_FLUC_HI 2e1
#define AA_PCNT_FLUC_LO AUTOADJUST_PCNT_LT_LO
#define CYCLE_COUNT_MAX 5
#define ALLOW_FACTOR 2.0

extern int get_estat;
extern int autoadjust_flag;
extern int autoadjust_t_runs_pcnt;
extern int debug;
extern int autoadjusting_hi_temp;
extern int autoadjusting_lo_temp;

extern double autoadjust_pcnt;
extern double pCnt_global;

/****************************************/
/*                                      */
/*       Autoadjust Temperature         */
/*                                      */
/****************************************/
double autoadjust_temp_pCnt(double temp_init,
			    double pcnt_lo,
			    double pcnt_hi,
			    double cell_min,
			    double cell_max,
			    double trans_frac,
			    double volum_frac,
			    int trn_chng,
			    int rot_chng,
			    int lat_chng,
			    int swp_chng)
{
  
  int aa_going_up=0;
  int aa_going_dn=0;
  int aa_overshot=0;
  int autoadjust_count=0;
  int autoadjust_runs=AUTOADJUST_RUNS;
  int cycle_count=0;
  
  double autoadjust_temp_init;
  double autoadjust_temp_finl;
  double autoadjust_temp_step=AUTOADJUST_TEMP_STEP;
  double pCnt_last= -1;
  double pCnt_fluc=0;
  double pCnt_fluc_allowed=0;
  
  autoadjust_temp_init = temp_init;
  autoadjust_temp_finl = autoadjust_temp_init - 1e-9;
  autoadjust_flag=1;

  if ( autoadjusting_hi_temp ) pCnt_fluc_allowed = AA_PCNT_FLUC_HI;
  if ( autoadjusting_lo_temp ) pCnt_fluc_allowed = AA_PCNT_FLUC_LO;
  
  do {
    
    autoadjust_count++;
    
    printf("*\n**********************************************************************\n*\n");
    printf("*                 -- Autoadjust Metropolis Run --\n*\n");
    printf("* %25s%12.4e\n","autoadjust_temp_init  =",autoadjust_temp_init);
    printf("* %25s%12.4e\n","autoadjust_temp_step  =",autoadjust_temp_step);
    printf("* %25s%12d%25s\n","autoadjust_runs  =",autoadjust_runs,"configs per temp");
    printf("* %25s%12d\n","autoadjust_count  =",autoadjust_count);
    printf("*\n");
    
    get_estat=1;
    pCnt_last= -1;
    cycle_count=0;
    do {
      Metro(autoadjust_runs,autoadjust_t_runs_pcnt,trn_chng,rot_chng,lat_chng,swp_chng,
	    cell_min,cell_max,trans_frac,volum_frac,autoadjust_temp_init,autoadjust_temp_finl);

      if ( pCnt_last == -1 ) {
	pCnt_fluc = -1;
	pCnt_last = pCnt_global;
      } else {
	/* what is the difference in pCnt between the last run and this run */
	pCnt_fluc = fabs( pCnt_last - pCnt_global );
	pCnt_last = pCnt_global;
      }
      
      printf("*\n*---------------------------------------------------------------------------------\n");
      printf("*    PCNT = %10.4e       PCNT_FLUC = %+10.4e  (allowed range = %10.3e)\n",
	     pCnt_global, pCnt_fluc, pCnt_fluc_allowed);
      if ( pCnt_fluc < 0 || pCnt_fluc > pCnt_fluc_allowed )
	printf("*    fluc > allowed, cycle again at same temperature  (%d)\n",++cycle_count);
      if ( cycle_count == CYCLE_COUNT_MAX ) {
	cycle_count=0;
	printf("*    cycle_count %d exceeded, mult fluct allowance by %.3f\n",CYCLE_COUNT_MAX,ALLOW_FACTOR);
	pCnt_fluc_allowed *= ALLOW_FACTOR;
      }
      printf("*---------------------------------------------------------------------------------\n");

    } while ( pCnt_fluc < 0 || pCnt_fluc > pCnt_fluc_allowed );
    
    printf("*\n**********************************************************************\n*\n");
    
    if ( debug > 1 ) {
      printf("*\n*          - Parameters at Post AUTOADJUST [%d] -\n",autoadjust_count);
      print_info_block();
    }
    
    printf("*                 -- AUTOADJUST [%d] OUT --\n*\n",autoadjust_count);
    printf("* %18s%12.3e%12s%12.3e\n","target: lo  =", pcnt_lo, "hi  =" , pcnt_hi );
    printf("* %18s%12.3e\n","Autoadjust pCnt =", autoadjust_pcnt );
    printf("*\n");
    
    
    /* if we are going up or down and we overshot the target percentage */
    if ( ( aa_going_dn && ( autoadjust_pcnt > pcnt_hi ) ) || 
	 ( aa_going_up && ( autoadjust_pcnt < pcnt_lo ) )     ) aa_overshot = 1;
    if ( (aa_overshot && aa_going_dn) || (aa_overshot && aa_going_up) ) {
      autoadjust_temp_step *= AUTOADJUST_TEMP_STEP_MULT;
      aa_overshot =0;
      printf("*     Autoadjust: overshot! drop temp step to: %.3e [a.u.]\n",autoadjust_temp_step);
    }
    
    /* 
     * check if we are above or below the target percentage and change
     * the temperature accordingly.  special case is if going down in
     * temp will take us below zero temp.  We have to have a special
     * test in the going down case.
     */
    if ( autoadjust_pcnt > pcnt_hi ) {
      aa_going_up = 1; aa_going_dn = 0;
      autoadjust_temp_init += autoadjust_temp_step;
      printf("*     Autoadjust: Raising the temperature %.3e [a.u.]\n",autoadjust_temp_step);
    }
    else if ( autoadjust_pcnt < pcnt_lo ) {
      aa_going_up = 0; aa_going_dn = 1;
      /* special case: if (T-Tstep) < 0 */
      if ( autoadjust_temp_init - autoadjust_temp_step <= 0.0 ) {
	do{
	  printf("*     Autoadjust: special case: (T-Tstep) = %.3e [a.u.]\n",autoadjust_temp_init-autoadjust_temp_step);
	  autoadjust_temp_step *= AUTOADJUST_TEMP_STEP_MULT;
	  printf("*     Autoadjust: special case: ===> drop temp step to: %.3e [a.u.]\n",autoadjust_temp_step);
	} while ( autoadjust_temp_init - autoadjust_temp_step < 0.0 );
      }
      autoadjust_temp_init -= autoadjust_temp_step;
      printf("*     Autoadjust: Lowering the temperature %.3e [a.u.]\n",autoadjust_temp_step);
      /* exit if we are so low in temperature changes that there are problems */
      if ( autoadjust_temp_step < AUTOADJUST_TEMP_STEP_MIN ){
	printf("*     WARN: Autoadjust: T step very low. Algorithm having problems. Exiting autoadjust.\n");
	printf("*\n**********************************************************************\n*\n");
	autoadjust_flag=0;
      }
    }
    else {
      aa_going_up = 0; aa_going_dn = 0;
      printf("*     Autoadjust: temp = %.3e  is good!  Exiting autoadjust.\n",autoadjust_temp_init);
      printf("*\n**********************************************************************\n*\n");
      autoadjust_flag=0;
    }
    autoadjust_temp_finl = autoadjust_temp_init * 0.999;
    if ( autoadjust_runs < AUTOADJUST_RUNS_MAX )
      autoadjust_runs = AUTOADJUST_RUNS + AUTOADJUST_RUNS_STEP*autoadjust_count;
    
  } while ( autoadjust_flag );
  
  
  return(autoadjust_temp_init);
}
