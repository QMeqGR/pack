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

extern int pf_anneal_flag,get_estat;
extern double cvd_1,pack_frac_min;

/****************************************/
/*                                      */
/*           PF anneal                  */
/*                                      */
/****************************************/
void pf_anneal(int runs,int t_runs,double temp_init, double temp_finl,
	       int trn_chng,int rot_chng,int lat_chng,int swp_chng)
{

  double cell_min = cvd_1;
  double cell_max = 100;
  double trans_frac = TRANS_FRAC_MAX;
  double volum_frac = VOLUM_FRAC_MAX;

  struct Energy E;

  printf("*\n**********************************************************************\n*\n");
  printf("*           -- Get Pack Frac Routine --\n*\n");
  pf_anneal_flag = 1;
  
  pack_frac_min = 1.0;

  get_estat=1;
  
  /* set lat change to zero, set by hand in this algorithm */
  Metro(runs,t_runs,trn_chng,rot_chng,lat_chng,swp_chng,
	cell_min,cell_max,trans_frac,volum_frac,temp_init,temp_finl);
  
  Simplex();
  
  E = Total_Energy();
  
  printf("*\n");
  printf("* E.eng = %lf\n",E.eng);
  printf("* E.pf = %lf\n",E.pf);
  printf("* Setting pack_frac_min to: %lf\n",E.pf);
  
  printf("*\n");
  
  pack_frac_min = E.pf;
  pf_anneal_flag = 0;
  printf("*\n**********************************************************************\n*\n");

  return;
}

