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

extern int efunc,numat,num_ob;
extern int done_flag;
extern int anneal_sched;
extern int ecol,mvout;
extern int ortho;
extern int debug;
extern int halt_feature;
extern int autoadjust_flag;
extern int autoadjust_chng;
extern int non_per_ewald;
extern int *SWPtab;

extern double dplane;
extern double cvd_max,c2vtol;
extern double hold_acc;
extern double max_ang_chng,aspect_max,sqsh_vol_pct;
extern double autoadjust_pct,autoadjust_pcnt,pCnt_global;

extern struct cellprm cell;
extern struct atom *s;
extern struct obj *object;
extern struct Energy E_best;

/******************************************************************/
/*                                                                */
/*                    Metropolis algorithm                        */
/*                                                                */
/******************************************************************/
void Metro(int runs,
	   int t_runs,
	   int trn_chng,
	   int rot_chng,
	   int lat_chng,
	   int swp_chng,
	   double cell_min,
	   double cell_max,
	   double trans_frac_max,
	   double volum_frac_max,
	   double temp_init,
	   double temp_finl)
{

  char K; /* for printing out: T in temp autoadjust, M in Metro */

  int i,order_loop;
  int trans_frac_const=1;
  int obj_number,swp_number=0;
  int n_steps=0;
  int tot_trn_attempts=0,tot_lat_attempts=0;
  int tot_rot_attempts=0,tot_swp_attempts=0;
  int nnze=0,tw=0;
  int conf_changes=0;
  int T_runs=t_runs,temp_iteration=0;
  int cell_chng_rand,order_rand_num,op1=0,op2=0,op3=0,op4=0,OP=0;
  int reject_trn=0,reject_lat=0,reject_rot=0,reject_swp=0;
  int trn_halt=0,lat_halt=0,rot_halt=0,swp_halt=0;
  int hold_it=0;
  int hold_type=floor(hold_acc);
  int ok,lat_not_ok_loops=0;
  int dplane_violation=0;
  int *w,current_config=0,last_accptd_config=0;
  int exceed_trans_attempts_REJECT=0;
  int n_trans_attempts=0;
  int tf_sign=1;
  int vi;

  double T=0,T_0=temp_init,T_L=temp_finl,T_alpha=0;
  double trans_frac=trans_frac_max,trans_alpha=0;
  double volum_frac=volum_frac_max;
  double *EE,Eave,SigEsq,SigE,E2ave,Eave2,pCnt;
  double *PF,PFave,SigPFsq,SigPF,PF2ave,PFave2;
  double dx,dy,dz;
  double dsx,dsy,dsz;
  double da,db,dc,dalph=0.0,dbeta=0.0,dgamm=0.0;
  double da_frac,db_frac,dc_frac;
  double vol,dvol;
  double phi,the,psi;
  double trn_pct=0,swp_pct=0,lat_pct=0,rot_pct=0;
  double hold_pct=(hold_acc - floor(hold_acc));
  double hold_rej_pct=0;
  static double min_cell_volume=0;

  struct Energy En,Elast,Etmp;
  struct vector trans_vec;
  struct vector s1={0},s2={0},swap_to_vec,swap_from_vec;
  struct vector cell_chng,ang_chng;
  struct vector pressure;
  struct vector TP={0};
  struct matrix R,L;
  struct obj obj_orig={0};
  
  static int been_called=0;

  if ( !been_called ) {

    min_cell_volume = SPHERE_VOL * ( num_ob * pow(cvd_max,3.0) );
    if ( debug > 2 ) printf("* Minimum cell volume: %10.3f\n",min_cell_volume);    

    En.ecc=0; En.eng=0;
    been_called=1;

  }

  /************************************************************************/
  /* Calculate the value for T_alpha */
  /* there are two extra runs, one at high temp, one at low temp */
  T_alpha = pow( (double)(T_runs-1), anneal_sched ) / log ( T_0 / T_L );

  /* calculate the value for trans_alpha */
  /* technical NOTE: this should still satisfy detailed balance within each temperature */
  if ( trans_frac_max < 0 ) {
    trans_frac_max = fabs(trans_frac_max);
    trans_alpha = (double)(T_runs) / log( trans_frac_max / TRANS_FRAC_MIN );
    trans_frac_const = 0;
  }
  else trans_frac_max = fabs(trans_frac_max);

  if ( debug > -1 ){
    printf("**********************************************************************\n");
    printf("*\n* Metropolis: %d temperature runs with these parameters:\n*\n",T_runs);
    printf("* %20s%10d%25s\n","runs  =",runs,"loops per T");
    printf("* %20s%10.3e%25s\n","T_0  =",T_0,"initial temp");
    printf("* %20s%10.3e%25s\n","T_L  =",T_L,"final temp");
    printf("* %20s%10.3e%25s\n","T_alpha  =",T_alpha,"annealing alpha");
    printf("*\n");
    if ( anneal_sched >0 ) printf("*            T = T_0 * exp( -run_num**%d / T_alpha )\n*\n",anneal_sched);
    
    if ( !trans_frac_const ) {
      printf("* %20s%10.3e%25s\n","trans_frac_max  =",trans_frac_max,"annealed");
      printf("* %20s%10.3e\n","trans_frac_min  =",TRANS_FRAC_MIN);
      printf("* %20s%10.3e\n","trans_alpha  =",trans_alpha);
      printf("*\n* trans_frac = trans_frac_max * exp( -run_num / trans_alpha )\n*\n");
    }
    printf("**********************************************************************\n");
  }
  /* Allocate memory for energy fluctuations calculation */
  EE = (double *)malloc( 4 * runs * sizeof(double) );
  if ( !EE ) {
      printf("%s\n","Not enough memory. Metro loc 1.\n");
      exit_pack(0);
  }
  /* Allocate memory for packing fraction fluctuations calculation */
  PF = (double *)malloc( 4 * runs * sizeof(double) );
  if ( !PF ) {
      printf("%s\n","Not enough memory. Metro loc 2.\n");
      exit_pack(0);
  }
  /* Allocate memory for weighting calculation */
  w = (int *)malloc( 4 * runs * sizeof(int) );
  if ( !w ) {
      printf("%s\n","Not enough memory. Metro loc 3.\n");
      exit_pack(0);
  }
  /* initialize the EE , PF, and w structures */
  for(i=0;i<4*runs;i++) { EE[i]=0.0; PF[i]=0.0; w[i]=0; }

  /************************************************************************/

  if ( ecol==0 && debug > -1 ) {
    
    if ( efunc==0 ) printf("%-3s%3s%10s%12s%12s%10s%12s%9s%10s%12s%28s\n",
			   "*H","run","Temp","E_ion","E_ss","cel vol","SigPF","ar","ortho","pCnt","% Rejected Changes");
    if ( efunc==1 ) printf("%-3s%3s%10s%12s%12s%10s%12s%9s%10s%12s%28s\n",
			   "*H","run","Temp","E_ion","E_lj","cel vol","SigPF","ar","ortho","pCnt","% Rejected Changes");
    printf("%-3s%90s%9s%9s%9s%9s\n","*H", " ","trn","rot","lat","swp");
    
    fflush(NULL);
  }

  /* initialize the seed for random number generation */
  /* srand48(1234); removed this line in version 4.6.1.0 */

  /*
   * The structure of this algorithm is:
   *
   *     MAIN TEMP LOOP (variable temp_iteration)
   *             set temperature
   *             set trans_max
   *
   *             MAIN CONFIG LOOP (variable n_steps)
   *
   *                     choose an object ( anion, cation )
   *                     choose order of (trans, rot, lat, swp)
   *
   *                     SUB CONFIG LOOP (4 loops)
   *
   *                           (1 trans, 2 rot, 3 lat, 4 swp) change
   *                            note: cations will NOT be rotated
   *
   *                     END SUB CONFIG LOOP
   *
   *             END MAIN CONFIG LOOP
   *
   *     END MAIN TEMP LOOP
   *
   */

  /***************************************/
  /* This is the main (temperature) loop */
  /***************************************/
  temp_iteration=0;
  done_flag=0;
  do {
    /* set the temperature for this run */
    if ( anneal_sched >0 ) T = T_0 * exp( -pow(temp_iteration,anneal_sched) / T_alpha );
    if ( anneal_sched==0 ) T = T_0 - temp_iteration * (T_0-T_L)/(T_runs);
    if ( !trans_frac_const ) trans_frac = trans_frac_max * exp( -temp_iteration / trans_alpha );
    temp_iteration++;
    da_frac = db_frac = dc_frac = trans_frac;
    if ( debug > 3 ) printf("* trans_frac = %.3f\n",trans_frac);

    
    /* initialize variables for temperature run */
    n_steps=0;

    tot_trn_attempts=0;
    tot_rot_attempts=0;
    tot_lat_attempts=0;
    tot_swp_attempts=0;
    reject_trn=0;
    reject_rot=0;
    reject_lat=0;
    reject_swp=0;

    Eave = 0;
    E2ave = 0;
    Eave2 = 0;
    SigEsq = 0;

    PFave = 0;
    PF2ave = 0;
    PFave2 = 0;
    SigPFsq = 0;
    
    pCnt = 0;

    /* DEBUG */
    if ( debug > 3 ){
      printf("* Top of temperature loop, temp_iteration= %d\n",temp_iteration);
      debug_block_metro( );
    }

    /* find initial system energy */
    En = (*Total_Energy)();
    Elast = En;
    if ( debug > 3 ) {
      printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
      printf("* init: E: %10.6e     temp_iteration = %d\n",Elast.tot,temp_iteration);
    }
    

    /*********************************/
    /* MAIN CONFIG LOOP (RUNS loops) */
    /*********************************/
    do {
      
      /* choose any random object */
      obj_number = (int)( drand48() * num_ob );
      if ( debug > 3 ) printf("* ************* METRO: CHOOSE OBJECT: no. %2d\n",obj_number);
      
      /*****************************
       * Below are the configuration changes:
       * 
       * 1) translational movement of an object
       *
       * 2) rotation of a polygon
       *
       * 3) changes in the lattice parameters
       *
       * 4) swap a cation with an anion
       *
       *****************************/

      /*
       * To maintain detailed-balance, we need to perform these operations,
       * not in sequence, but in random order.  We will therefore, choose a
       * random number from 1 to 24 and set the order of operations as:
       *
       *    0=trn      1=rot      2=lat     3=swp
       *
       * 0123    0132     0213    0231   0312    0321
       * 1023    1032     1203    1230   1302    1320
       * 2013    2031     2103    2130   2301    2310
       * 3012    3021     3102    3120   3201    3210
       *
       * we will then run through the set of configs three times, each time
       * setting the appropriate operation number 0,1,2,3.
       *
       */
      order_rand_num = (int)( drand48() * 24.0 );
      if ( debug > 4 ) printf("* order_rand_num = %d\n",order_rand_num);
      switch ( order_rand_num ) {
      case  0: op1=0; op2=1; op3=2; op4=3; break;
      case  1: op1=0; op2=1; op3=3; op4=2; break;
      case  2: op1=0; op2=2; op3=1; op4=3; break;
      case  3: op1=0; op2=2; op3=3; op4=1; break;
      case  4: op1=0; op2=3; op3=1; op4=2; break;
      case  5: op1=0; op2=3; op3=2; op4=1; break;

      case  6: op1=1; op2=0; op3=2; op4=3; break;
      case  7: op1=1; op2=0; op3=3; op4=2; break;
      case  8: op1=1; op2=2; op3=0; op4=3; break;
      case  9: op1=1; op2=2; op3=3; op4=0; break;
      case 10: op1=1; op2=3; op3=0; op4=2; break;
      case 11: op1=1; op2=3; op3=2; op4=0; break;

      case 12: op1=2; op2=0; op3=1; op4=3; break;
      case 13: op1=2; op2=0; op3=3; op4=1; break;
      case 14: op1=2; op2=1; op3=0; op4=3; break;
      case 15: op1=2; op2=1; op3=3; op4=0; break;
      case 16: op1=2; op2=3; op3=0; op4=1; break;
      case 17: op1=2; op2=3; op3=1; op4=0; break;

      case 18: op1=3; op2=0; op3=1; op4=2; break;
      case 19: op1=3; op2=0; op3=2; op4=1; break;
      case 20: op1=3; op2=1; op3=0; op4=2; break;
      case 21: op1=3; op2=1; op3=2; op4=0; break;
      case 22: op1=3; op2=2; op3=0; op4=1; break;
      case 23: op1=3; op2=2; op3=1; op4=0; break;
      }
      if ( debug > 3 ) {
	printf("* op1 op2 op3 op4 = %d %d %d %d\t\t0=trn  1=rot  2=lat  3=swp\n",op1,op2,op3,op4);
      }


      /*********************************/
      /* SUB CONFIG  LOOP  (3 loops)   */
      /*********************************/      
      for(order_loop=0; order_loop<4; order_loop++) {

	if ( order_loop == 0 ) OP = op1;
	if ( order_loop == 1 ) OP = op2;
	if ( order_loop == 2 ) OP = op3;
	if ( order_loop == 3 ) OP = op4;
	if ( debug > 3 ) {
	  printf("* order_loop = %d     OP = %d   (0=trn 1=rot 2=lat 3=swp)\n",order_loop,OP);
	}

	/*****************************************************************************************/
	/*                                 **  TRANSLATION  **                                   */
	/*****************************************************************************************/

	/* 1. translational movement of an object */
	if ( trn_chng && OP==0 && object[obj_number].vfxd[0]==0 ) {

	  tot_trn_attempts += 1;

	  /* cut here --------------- */
	  
	  /* get the current location of the object */
	  TP = object[obj_number].v[0];

	  /* randomly pick the translation numbers */
	  /* dx, dy, dz should be between -1...+1  */
	  dx = drand48() * da_frac * rsign();
	  dy = drand48() * db_frac * rsign();
	  dz = drand48() * dc_frac * rsign();

	  trans_vec = vsub( rezone( makevec( TP.x+dx, TP.y+dy, TP.z+dz ) ) , TP );

	  /* make sure translations don't take objects too close to the wall */
	  n_trans_attempts=0;
	  exceed_trans_attempts_REJECT=0;
	  /* check constraints and non_per_ewald */
	  while ( 0==trans_accept( obj_number, vadd(TP,trans_vec), non_per_ewald ) && exceed_trans_attempts_REJECT==0 ) {
	    dx = drand48() * da_frac * rsign();
	    dy = drand48() * db_frac * rsign();
	    dz = drand48() * dc_frac * rsign();
	    trans_vec = vsub( rezone( makevec( TP.x+dx, TP.y+dy, TP.z+dz ) ) , TP );
	    
	    if ( debug > 4 ){
	      printf("* obj[%2d] transvec[ %7d ] = %9.4f%9.4f%9.4f  TP= %9.4f%9.4f%9.4f  New= %9.4f%9.4f%9.4f\n",
		     obj_number,
		     n_trans_attempts,trans_vec.x,trans_vec.y,trans_vec.z,
		     TP.x,TP.y,TP.z,
		     makevec( TP.x+dx, TP.y+dy, TP.z+dz ).x,
		     makevec( TP.x+dx, TP.y+dy, TP.z+dz ).y,
		     makevec( TP.x+dx, TP.y+dy, TP.z+dz ).z);
	    }
	    
	    if ( n_trans_attempts > (int)MAX_TRANS_ATTEMPTS ) {
	      if ( debug > -1 ){
		trans_frac = trans_frac_max + trans_frac_max*TRANS_FRAC_INCREASE_FACTOR*tf_sign; tf_sign *= -1;
		printf("* WL: TRANS: WARN: --------->> CHANGING trans_frac by %2d * %.2f to %10.5e\n",tf_sign,TRANS_FRAC_INCREASE_FACTOR,trans_frac);
		printf("* MT: TRANS: WARN: Exceeded max number of attempts at translation (%d).\n",(int)MAX_TRANS_ATTEMPTS);
		printf("* MT: TRANS: WARN: obj_number = %d\n",obj_number);
		printf("* obj[%2d] transvec[ %7d ] = %9.4f%9.4f%9.4f  TP= %9.4f%9.4f%9.4f  New= %9.4f%9.4f%9.4f\n",
		       obj_number,
		       n_trans_attempts,trans_vec.x,trans_vec.y,trans_vec.z,
		       TP.x,TP.y,TP.z,
		       makevec( TP.x+dx, TP.y+dy, TP.z+dz ).x,
		       makevec( TP.x+dx, TP.y+dy, TP.z+dz ).y,
		       makevec( TP.x+dx, TP.y+dy, TP.z+dz ).z);
	      }
	      
	      exceed_trans_attempts_REJECT=1;
	    } /* if n_trans_attmempts > ... */	      
	    
	  } /* end while : check constraints and non_per_ewald */
	  
	  if ( debug > 3 ) {
	    printf("* MT: n_trans_attempts= %d\n",n_trans_attempts);
	    if ( exceed_trans_attempts_REJECT==1 ){
	      printf("* REJECT CHANGE: TRANS - rejected because exceeded trans_attempts max\n");
	    }
	  }
	  
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: TRANS obj_number=[%d]\n",obj_number);
	    printf("* dx dy dz = %15.9f%15.9f%15.9f\n",dx,dy,dz);
	    printf("* TP       = %15.10f%15.10f%15.10f\n",TP.x,TP.y,TP.z);
	    if ( TP.x>1.0 || TP.x<0.0 || TP.y>1.0 || TP.y<0.0 || TP.z>1.0 || TP.z<0.0 ){
	      printf("* Metro: atom is out of bounds! Exiting.\n"); exit_pack(0);
	    }
	    printf("* transvec = %15.10f%15.10f%15.10f\n",trans_vec.x,trans_vec.y,trans_vec.z);
	    printf("* exceed_trans_attempts_REJECT = %d\n",exceed_trans_attempts_REJECT);
	  }

	  if ( exceed_trans_attempts_REJECT==0 ) {
	    /* if we are here, we have to recalculate L, the transformation
	       matrix, becuase we may have changed the lattice parameters
	       in the lat parm change below. */
	    L = L_e3(&cell);
	    trans_object( obj_number, &L, trans_vec );
	    if ( debug > 3 ) print_c2v( );
	  }

	  /* what are the energy consequences */
	  En = (*Total_Energy)();


	  if ( debug > 3 ) {
	    printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
	    printf("* trans: Elast.tot = %10.6e\tEn.tot = %10.6e\n",Elast.tot,En.tot);
	  }

	  /* reject the change if the temperature is low and En.tot > Elast.tot */
	  if ( ( En.tot > Elast.tot &&  (drand48() > exp( -(double)(En.tot - Elast.tot)/T )) ) ||
	       exceed_trans_attempts_REJECT==1 ) {
	    reject_trn++;	
	    if ( debug > 3 ) printf("* REJECT CHANGE: TRANS\n");
	    
	    trans_vec = makevec( -trans_vec.x , -trans_vec.y , -trans_vec.z );
	    
	    if ( exceed_trans_attempts_REJECT==0 ) {	    
	      trans_object( obj_number, &L, trans_vec );
	    }
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
	    w[last_accptd_config] += 1;
	  } else {
	    /* if move was accepted, then weight the energy */
	    current_config = (4 * n_steps + 0);
	    w[current_config] += 1;
	    last_accptd_config = current_config;
	    EE[ current_config ] = En.tot;
	    PF[ current_config ] = En.pf;
	    if ( debug > 5 ) {
	      printf("* EE[%d]=%e\n",current_config,EE[current_config]);
	      printf("* PF[%d]=%e\n",current_config,PF[current_config]);
	      printf("* w[%d]=%d\n",current_config,w[current_config]);
	    }
	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}
	

	/*****************************************************************************************/
	/*                                   **  ROTATION  **                                    */
	/*****************************************************************************************/

	/* 2. Re-orient an object */
	if ( rot_chng  && 
	     object[obj_number].num_vert>1 &&
	     OP==1 && 
	     object[obj_number].allowrot==1 ) {

	  tot_rot_attempts += 1;

	  /* recalculate L */
	  L = L_e3(&cell);

	  /* initial system energy is Elast.tot */

	  /* choose random rotation angles */
	  phi = 2*PI * drand48();
	  the = PI/2.0 + asin( 2 * drand48() - 1.0 );
	  psi = 2*PI * drand48();

	  R = R_ptp(phi,the,psi);
	  
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: ROT [%d]\n",obj_number);
	    printf("* rot: phi the psi = %8.4f%8.4f%8.4f\n",
		   phi,the,psi);
	  }

	  obj_orig = object[obj_number];
	  rot_object( obj_number, &L, &cell, R );
	  
	  /* what are the energy consequences */
	  En = (*Total_Energy)();


	  if ( debug > 3 ) {
	    printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
	    printf("* rotation: Elast.tot = %10.6e\tEn.tot = %10.6e\n",Elast.tot,En.tot);
	  }
	  
	  /* reject the change if the temperature is low and En.tot > Elast.tot */
	  if ( En.tot > Elast.tot &&  (drand48() > exp( -(double)(En.tot-Elast.tot)/T )) ) {
	    reject_rot++;	
	    if ( debug > 3 ) printf("* REJECT CHANGE: ROT\n");
	    
	    /* replace the rotated object with the original */
	    object[obj_number] = obj_orig;
	    for(vi=0; vi<MAX_VERTS && object[obj_number].vused[vi]==1; vi++){
	      object[obj_number].vmoved[vi] = 1; /* flag for energy calc */
	    }
	    
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
	    w[last_accptd_config] += 1;
	  } else {
	    /* if move was accepted, then weight the energy */
	    current_config = (4 * n_steps + 1);
	    w[current_config] += 1;
	    last_accptd_config = current_config;
	    EE[ current_config ] = En.tot;
	    PF[ current_config ] = En.pf;
	    if ( debug > 5 ) {
	      printf("* EE[%d]=%e\n",current_config,EE[current_config]);
	      printf("* PF[%d]=%e\n",current_config,PF[current_config]);
	      printf("* w[%d]=%d\n",current_config,w[current_config]);
	    }
	    
	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}
	

	/*****************************************************************************************/
	/*                       **  Lattice Parameter Change  **                                */
	/*****************************************************************************************/

	/* 3. Change the lattice parameters */
	
	/* Notes, this code:
	   (1) should not scale the bond lengths within the anions.
	   (2) will change the lattice vectors lengths and angles.
	   (3) will prevent aspect ratio from going bigger than aspect_max.
	*/
	
	if ( lat_chng && OP==2 ) {
	  /* initial system energy is Elast.tot */

	  tot_lat_attempts += 1;

	  lat_not_ok_loops=0;
	  do {

	    da = 0; db = 0; dc = 0; dalph =0; dbeta = 0; dgamm = 0;
	    ok=1;

	    /* randomly choose the volume to change up to VOLUM_FRAC_MAX */
	    vol = cell_volume(&cell);
	    dvol = drand48() * volum_frac * vol;

	    /* randomly pick which cell value to change: 0,1,2 lat vec length */
	    /*                                           3,4,5 angle          */
	    /*                                           6,7,8 scale change   */
	    cell_chng_rand = (int)( drand48() * 9 );
	    if ( debug > 4 ) printf("* cell_chng_rand = %d  0-2 lat vec mag, 3-5 angle, 6-8 scale\n",cell_chng_rand);

	    if ( cell_chng_rand == 0 ) { da = rsign() * dvol * cell.a / vol; }
	    if ( cell_chng_rand == 1 ) { db = rsign() * dvol * cell.b / vol; }
	    if ( cell_chng_rand == 2 ) { dc = rsign() * dvol * cell.c / vol; }

	    if ( cell_chng_rand == 3 ) { dalph = rsign() * drand48() * max_ang_chng * (1-fabs(PIo2-cell.alph)/PIo2); }
	    if ( cell_chng_rand == 4 ) { dbeta = rsign() * drand48() * max_ang_chng * (1-fabs(PIo2-cell.beta)/PIo2); }
	    if ( cell_chng_rand == 5 ) { dgamm = rsign() * drand48() * max_ang_chng * (1-fabs(PIo2-cell.gamm)/PIo2); }

	    if ( cell_chng_rand > 5 ) {
	      if ( rsign() > 0 ) {
		da =  dvol * cell.a / vol;  db =  dvol * cell.b / vol;  dc =  dvol * cell.c / vol;
	      }
	      else {
		da = -dvol * cell.a / vol;  db = -dvol * cell.b / vol;  dc = -dvol * cell.c / vol;
	      }
	    }

	    if ( ortho ) { dalph = 0.0; dbeta = 0.0; dgamm = 0.0; }

	    /* if the cell changes are too much, set the cell max or min size and
	       recalculate the change parameters.
	    */

	    if ( cell_chng_rand == 0 ) {
	      if ( fmax(cell.a+da,cell.b)/fmin(cell.a+da,cell.b) > aspect_max ) ok=0;
	      if ( fmax(cell.a+da,cell.c)/fmin(cell.a+da,cell.c) > aspect_max ) ok=0;
	      if ( cell.a + da > cell_max ) ok=0;
	      if ( cell.a + da < cell_min ) ok=0;
	    }
	    if ( cell_chng_rand == 1 ) {
	      if ( fmax(cell.b+db,cell.a)/fmin(cell.b+db,cell.a) > aspect_max ) ok=0;
	      if ( fmax(cell.b+db,cell.c)/fmin(cell.b+db,cell.c) > aspect_max ) ok=0;
	      if ( cell.b + db > cell_max ) ok=0;
	      if ( cell.b + db < cell_min ) ok=0;
	    }
	    if ( cell_chng_rand == 2 ) {
	      if ( fmax(cell.c+dc,cell.a)/fmin(cell.c+dc,cell.a) > aspect_max ) ok=0;
	      if ( fmax(cell.c+dc,cell.b)/fmin(cell.c+dc,cell.b) > aspect_max ) ok=0;
	      if ( cell.c + dc > cell_max ) ok=0;
	      if ( cell.c + dc < cell_min ) ok=0;
	    }
	    if ( cell_chng_rand == 3 ) {
	      if ( (cell.alph + dalph) < ANG_ALPH_MIN || (cell.alph + dalph) > ANG_ALPH_MAX ) ok=0;
	    }
	    if ( cell_chng_rand == 4 ) {
	      if ( (cell.beta + dbeta) < ANG_BETA_MIN || (cell.beta + dbeta) > ANG_BETA_MAX ) ok=0;
	    }
	    if ( cell_chng_rand == 5 ) {
	      if ( (cell.gamm + dgamm) < ANG_GAMM_MIN || (cell.gamm + dgamm) > ANG_GAMM_MAX ) ok=0;
	    }


	    /****************************************************************/

	    cell_chng = makevec(da,db,dc);
	    ang_chng = makevec(dalph,dbeta,dgamm);

	    /* check for dplane violations and lat change pathology sqrt(neg number ) */
	    dplane_violation = lat_pathological( &cell, cell_chng, ang_chng, dplane );
	    if ( dplane_violation ) ok=0;
	    if ( debug > 4 && dplane_violation ) {
	      printf("* dplane_violation = %d (2=dplane 1=sqrt(-num) 3=both)\n",dplane_violation);
	    }

	    if ( lat_not_ok_loops++ >= 5e6 ) {
	      printf("* FATAL ERROR: lat change code in metro routine:\n");
	      printf("* dplane_violation = %d (2=dplane 1=sqrt(-num) 3=both)\n",dplane_violation);
	      printf("* Can't find an OK lat change after %d tries! Cell conformation is bad.\n",lat_not_ok_loops);
	      kill(getpid(),SIGINT);
	    }

	  } while ( !ok );

	  
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: LAT\n");
	    if ( debug>4 ){
	      printf("* latvec a = %12.6f%12.6f%12.6f\n",cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
	      printf("* latvec b = %12.6f%12.6f%12.6f\n",cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
	      printf("* latvec c = %12.6f%12.6f%12.6f\n",cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
	    }
	    printf("*  a  b  c = %12.6f%12.6f%12.6f\n",cell.a,cell.b,cell.c);
	    printf("* da db dc = %12.6f%12.6f%12.6f\n",da,db,dc);
	    printf("*  angles  = %12.6f%12.6f%12.6f\n",cell.alph,cell.beta,cell.gamm);
	    printf("* dangles  = %12.6f%12.6f%12.6f\n",dalph,dbeta,dgamm);
	    printf("* dvol = %12.6f    cell_vol = %12.6f     ratio = %12.6f\n",
		   dvol,cell_volume(&cell),dvol/cell_volume(&cell));
	  }

	  /* rescale the object locations */
	  /* the rescaling needs the pre-changed cell parameters and the change params */
	  rescale_objects( &cell, cell_chng, ang_chng);
	  if ( debug>3 ) print_c2v();

	  /* what are the energy consequences */
	  mark_all_objects_moved();
	  En = (*Total_Energy)();

	  if ( debug > 3 ) {
	    printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
	    printf("* cell chng: Elast.tot = %10.6e  En.tot = %10.6e\n",Elast.tot,En.tot);

	  }

	  /* reject the change if the temperature is low and E3 > E2 */
	  /* Also reject if the packing fraction is unreasonably high (cell squashing) */
	  /* Also reject if the dplane criterion is violated (cell flattening) */
	  if ( ( En.tot > Elast.tot &&  (drand48() > exp( -(double)(En.tot-Elast.tot)/T ) ) ) ||
	       ( cell_volume(&cell) < min_cell_volume ) ||
	       ( En.orth < sqsh_vol_pct ) ) {
	    reject_lat++;
	    if ( debug > 3 ) printf("* REJECT CHANGE: LAT PARM\n");
	    
	    cell_chng = makevec(-da,-db,-dc);
	    ang_chng = makevec(-dalph,-dbeta,-dgamm);
	    
	    if ( debug > 3 ) {
	      printf("*  a  b  c = %15.9f%15.9f%15.9f\n",cell.a,cell.b,cell.c);
	      printf("* da db dc = %15.9f%15.9f%15.9f\n",-da,-db,-dc);
	      printf("*  angles  = %12.6f%12.6f%12.6f\n",cell.alph,cell.beta,cell.gamm);
	      printf("* dangles  = %12.6f%12.6f%12.6f\n",-dalph,-dbeta,-dgamm);
	    }
	    
	    /* rescale the objects */
	    rescale_objects( &cell, cell_chng, ang_chng);
	    
	    /* reset energy */
	    mark_all_objects_moved();
	    En = Elast;
	    w[last_accptd_config] += 1;
	  } else {
	    /* if move was accepted, then weight the energy */
	    current_config = (4 * n_steps + 2);
	    w[current_config] += 1;
	    last_accptd_config = current_config;
	    EE[ current_config ] = En.tot;
	    PF[ current_config ] = En.pf;
	    if ( debug > 5 ) {
	      printf("* EE[%d]=%e\n",current_config,EE[current_config]);
	      printf("* PF[%d]=%e\n",current_config,PF[current_config]);
	      printf("* w[%d]=%d\n",current_config,w[current_config]);
	    }

	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}
	

	/*****************************************************************************************/
	/*                                  **  SWAPPING  **                                     */
	/*****************************************************************************************/	
	
	/* 4. swapping of a cation with an anion */
	if ( swp_chng && OP==3 &&
	     num_ob>1 && object[obj_number].vfxd[0]!=1 &&
	     object[obj_number].swappable==1 ) {

	  tot_swp_attempts += 1;

	  /* recalculate L */
	  L = L_e3(&cell);

	  /* get object coordinates */
	  s1 = object[obj_number].v[0];

	  /* pick a different object */
	  do {
	    swp_number = (int)( drand48() * num_ob );
	  } while ( swp_number == obj_number
		    || 0==SWPtab[ object[obj_number].typ * SWPTAB + object[swp_number].typ ]
		   );
	  s2 = object[swp_number].v[0];

	  dsx = s2.x - s1.x;
	  dsy = s2.y - s1.y;
	  dsz = s2.z - s1.z;
	  swap_to_vec   = makevec(  dsx,  dsy,  dsz );
	  swap_from_vec = makevec( -dsx, -dsy, -dsz );
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: SWAP\n");
	    printf("* swp: s1 %9.4f%9.4f%9.4f   (no. %2d) with:\n"
		   "       s2 %9.4f%9.4f%9.4f   (no. %2d)\n",
		   s1.x,s1.y,s1.z,obj_number,
		   s2.x,s2.y,s2.z,swp_number);
	    printf("* dsx dsy dsz = %15.9f%15.9f%15.9f (swap_to_vec)\n",dsx,dsy,dsz);
	  }
	  
	  /* swap the objects */
	  trans_object( obj_number, &L, swap_to_vec );
	  trans_object( swp_number, &L, swap_from_vec );
	  
	  /* what are the energy consequences */
	  En = (*Total_Energy)();

	  if ( debug > 3 ) {
	    printf("* swap: ecc = %10.6e        eng = %10.6e    vol = %10.6f\n",
		   En.ecc,En.eng,En.vol);
	    printf("*  Elast.cc = %10.6e  Elast.eng = %10.6e\n",
		   Elast.ecc,Elast.eng);
	  }

	  /* reject the change if the temperature is low and En.tot > Elast.tot */
	  if ( En.tot > Elast.tot &&  (drand48() > exp( -(double)(En.tot - Elast.tot)/T )) ) {
	    reject_swp++;	
	    if ( debug > 3 ) printf("* REJECT CHANGE: SWAP\n");
	    
	    swap_to_vec   = makevec( -dsx, -dsy, -dsz );
	    swap_from_vec = makevec(  dsx,  dsy,  dsz );	    
	    
	    /* swap the objects back */
	    trans_object( obj_number, &L, swap_to_vec );
	    trans_object( swp_number, &L, swap_from_vec );		

	    if ( debug > 3 ) {
	      printf("* Swap back check:\n");
	      Etmp = (*Total_Energy)();
	      printf("* swp back chk: tot= %10.6e  ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",
		     Etmp.tot,Etmp.ecc,Etmp.eng,Etmp.vol);
	    }
	    
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
	    w[last_accptd_config] += 1;
	  } else {
	    current_config = (4 * n_steps + 3);
	    w[current_config] += 1;
	    last_accptd_config = current_config;
	    EE[ current_config ] = En.tot;
	    PF[ current_config ] = En.pf;
	    if ( debug > 5 ) {
	      printf("* EE[%d]=%e\n",current_config,EE[current_config]);
	      printf("* PF[%d]=%e\n",current_config,PF[current_config]);
	      printf("* w[%d]=%d\n",current_config,w[current_config]);
	    }
	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}

    

	
	/* DEBUG */
	if ( debug > 3 ){
	  printf("* Post config change\n");
	  debug_block_metro( );
	  /* find and print the minimum distances */
	  get_min_distances( );
	}

	if ( mvout ) print_xbs_mvframe();

      } /* end order loop */

      n_steps++;
    } while ( n_steps < runs );
    

    /************ METROPOLIS TEMP LOOP FINISHED ********************/

    /* Energy fluctuations calculations */
    /* count the number of non-zero energies for the average.
       We do this in the case that some of the config changes
       above were skipped over, such as the case with no
       rotations when doing spherical ions.  Some of the
       energies in EE[] will be zero and we don't want to
       count those energies.  No system should have zero
       total energy, so this is a fine way of weeding
       out the unwanted zeros. */
    for (i=0,nnze=0; i<4*runs; i++) if ( !tol(EE[i],0.0,1e-6) ) nnze++;
    /* the new code for weighting removes the need for nnze variable,
     but we'll calculate it anyway, just for debugging purposes      */

    /* tw is total weight */
    for (i=0,tw=0; i<4*runs; i++) tw += w[i];
    for (i=0; i<4*runs; i++) {
      /* if ( w[i] >0 ) printf("EE[%d] = %20.10e    w = %10d\n",i,EE[i],w[i]); */
      Eave   += EE[i] * w[i];
      E2ave  += EE[i] * EE[i] * w[i];
      PFave  += PF[i];
      PF2ave += PFave * PFave;
    }
    Eave   /= tw;
    E2ave  /= tw;
    PFave  /= tw;
    PF2ave /= tw;
    
    for (i=0; i<4*runs; i++) {
      SigEsq += (EE[i] - Eave)*(EE[i] - Eave)*w[i];
      SigPFsq += (PF[i] - PFave)*(PF[i] - PFave)*w[i];
    }
    SigEsq /= tw;
    SigPFsq /= tw;

    SigE = sqrt(SigEsq);
    SigPF = sqrt(SigPFsq);
    
    Eave2 = Eave*Eave;

    /* pseudo Heat Cap stuff */
    pCnt = fabs( E2ave - Eave2 ) / (T * T);
    pCnt_global = pCnt;
    autoadjust_pcnt = pCnt;

    if ( debug>3 ) printf("* tw = %d  Eave= %f   E2ave= %f  Eave2= %f\n",tw,Eave,E2ave,Eave2);
    PFave2 = PFave*PFave;



    /* reset the variables to zero for the next temperature loop */
    for(i=0;i<4*runs;i++) { EE[i]=0.0; PF[i]=0.0; w[i]=0; }

    conf_changes = runs;

    if ( debug>3 ) { 
      printf("* conf_changes=%d     4*runs=%d\n",conf_changes,4*runs);
      printf("* nnze=%d (number of non-zero energies)\n",nnze);
      printf("* tw = %d (total weight)\n",tw);
      printf("* tot_trn_attempts = %10d\n",tot_trn_attempts);
      printf("* tot_rot_attempts = %10d\n",tot_rot_attempts);
      printf("* tot_lat_attempts = %10d\n",tot_lat_attempts);
      printf("* tot_swp_attempts = %10d\n",tot_swp_attempts);
      printf("* reject_trn = %10d\n",reject_trn);
      printf("* reject_rot = %10d\n",reject_rot);
      printf("* reject_lat = %10d\n",reject_lat);
      printf("* reject_swp = %10d\n",reject_swp);
    }

    pressure = Pressure();

    if ( tot_trn_attempts == 0 ) trn_pct = 0;
    else trn_pct = (double)reject_trn/tot_trn_attempts*100;
    if ( tot_rot_attempts == 0 ) rot_pct = 0;
    else rot_pct = (double)reject_rot/tot_rot_attempts*100;
    if ( tot_lat_attempts == 0 ) lat_pct = 0;
    else lat_pct = (double)reject_lat/tot_lat_attempts*100;
    if ( tot_swp_attempts == 0 ) swp_pct = 0;
    else swp_pct = (double)reject_swp/tot_swp_attempts*100;

    switch ( autoadjust_chng ){
    case 1: autoadjust_pct = trn_pct/100.0; break;
    case 2: autoadjust_pct = rot_pct/100.0; break;
    case 3: autoadjust_pct = lat_pct/100.0; break;
    case 4: autoadjust_pct = swp_pct/100.0; break;
    }

    if ( autoadjust_flag == 1 ) K='T';
    else K='M';

    if ( ecol && debug > -1 ){
      printf("*%c------------------------------------------------------------\n",K);
      printf("*%c %20s%12d%25s\n",K,"T_iter  =",temp_iteration,"temperature iteration");
      printf("*%c %20s%12.3e%25s\n",K,"T  =",T,"temperature");
      printf("*%c %20s%12.3e%25s\n",K,"ecc  =",Elast.ecc,"electrostatic energy");
      printf("*%c %20s%12.3e%25s\n",K,"eng  =",Elast.eng,"repulsive energy");
      printf("*%c %20s%12.3e%25s\n",K,"bnd  =",Elast.bnd,"bounds");
      printf("*%c %20s%12.3f%25s\n",K,"Sigpf  =",SigPF,"pack frac sigma");
      printf("*%c %20s%12.3f%25s\n",K,"vol  =",Elast.vol,"cell volume");
      printf("*%c %20s%12.3f%25s\n",K,"lat a  =",Elast.a,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"lat b  =",Elast.b,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"lat c  =",Elast.c,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"ar  =",aspect(),"aspect ratio");
      printf("*%c %20s%12.3f%25s\n",K,"ortho  =",Elast.orth,"cv/ortho");
      printf("*%c %20s%12.3f%25s\n",K,"alph  =",Elast.alph*R2D,"angle 2--3");
      printf("*%c %20s%12.3f%25s\n",K,"beta  =",Elast.beta*R2D,"angle 1--3");
      printf("*%c %20s%12.3f%25s\n",K,"gamm  =",Elast.gamm*R2D,"angle 1--2");
      printf("*%c %20s%12.3e%25s\n",K,"pCnt  =",pCnt,"pseudo heat capacity");
      printf("*%c %20s%12.3e%25s\n",K,"SigEsq  =",SigEsq,"energy sigma");
      printf("*%c %20s%12.3e%25s\n",K,"Px  =",pressure.x,"pressure x-dir");
      printf("*%c %20s%12.3e%25s\n",K,"Py  =",pressure.y,"pressure y-dir");
      printf("*%c %20s%12.3e%25s\n",K,"Pz  =",pressure.z,"pressure z-dir");
      printf("*%c %20s%12.2f%25s\n",K,"trn  =",trn_pct,"% rej  translations");
      printf("*%c %20s%12.2f%25s\n",K,"rot  =",rot_pct,"% rej     rotations");
      printf("*%c %20s%12.2f%25s\n",K,"lat  =",lat_pct,"% rej lattice parms");
      printf("*%c %20s%12.2f%25s\n",K,"swp  =",swp_pct,"% rej   ca-an swaps");
      printf("*%c %20s%12.3e%25s\n",K,"best tot  =",E_best.tot,"may include cell vol");
      printf("*%c %20s%12.3e%25s\n",K,"best ecc  =",E_best.ecc,"electrostatic");
      printf("*%c %20s%12.3e%25s\n",K,"best eng  =",E_best.eng,"repulsive");
      printf("*%c %20s%12.3e%25s\n",K,"best vol  =",E_best.vol,"best volume");
      
    } else if ( debug > -1 ) {
      printf("*%c %3d%10.2e%+12.3e%+12.3e%10.2f%12.3e%9.3f%10.3f%12.3e%9.2f%9.2f%9.2f%9.2f\n",K,
	     temp_iteration,T,Elast.ecc,Elast.eng,Elast.vol,SigPF,aspect(),
	     Elast.orth,pCnt,trn_pct,rot_pct,lat_pct,swp_pct);
    }
    fflush(NULL);

    /*****************************************************************/
    /* The code that shuts down certain config changes if they have
     * been at 100% rejection for more than N temp steps
     */
    if ( halt_feature ) {
      /* if the percentage reaches 100, then increment
       * the count for the halt variable.
       */
      if ( tol( trn_pct, 100, 0.001) ) trn_halt++;
      if ( tol( rot_pct, 100, 0.001) ) rot_halt++;
      if ( tol( lat_pct, 100, 0.001) ) lat_halt++;
      if ( tol( swp_pct, 100, 0.001) ) swp_halt++;

      /* if the percent was at 100 and then dropped
       * down, restart the count.  We don't want to
       * cut out config changes if they are working.
       */
      if ( trn_pct>0.0 && trn_pct<100.0 && trn_halt>0 ) trn_halt=0;
      if ( rot_pct>0.0 && rot_pct<100.0 && rot_halt>0 ) rot_halt=0;
      if ( lat_pct>0.0 && lat_pct<100.0 && lat_halt>0 ) lat_halt=0;
      if ( swp_pct>0.0 && swp_pct<100.0 && swp_halt>0 ) swp_halt=0;
      
      if ( trn_halt == N_HALT_CHNG ) trn_chng=0;
      if ( rot_halt == N_HALT_CHNG ) rot_chng=0;
      if ( lat_halt == N_HALT_CHNG ) lat_chng=0;
      if ( swp_halt == N_HALT_CHNG ) swp_chng=0;
    }
    
    /*****************************************************************/
    /* the 'holding' code, which keeps the temperature the same while
       rejection percentage is in the desired range */
    /* which type of rejection percentage to look at */
    /* if number of iterations is greater than max, increment the hold percent */
    /* only hold for a 20% range, or up to 98%, then stop holding */
    if ( debug>3 ) printf("* hold_type=%d    0=trn    1=rot    2=lat   3=swp\n",hold_type);
    if ( debug>3 ) printf("* hold_acc=%f\n",hold_acc);
    if ( debug>3 ) printf("* hold_pct=%f\n",hold_pct);

    if      ( hold_type==0 ) hold_rej_pct = (double)reject_trn/conf_changes;
    else if ( hold_type==1 ) hold_rej_pct = (double)reject_rot/conf_changes;
    else if ( hold_type==2 ) hold_rej_pct = (double)reject_lat/conf_changes;
    else                     hold_rej_pct = (double)reject_swp/conf_changes;

    if ( hold_acc > 0.0  &&  hold_rej_pct > hold_pct ) {
      if ( hold_it++ >= HOLD_IT_MAX ) { hold_pct += 0.1; hold_it=0; }
      if ( ( hold_pct >= (hold_acc-floor(hold_acc)+0.2) ) ||
	   ( hold_rej_pct >= 0.98 ) ) hold_acc = -1.0; /* yes, hold_acc here; to turn it off */
      if ( debug>3 ) printf("* hold_it = %d    hold_pct = %.2f\n",hold_it,hold_pct);
      temp_iteration--;
    }
    
    if ( tol(T,T_L,1e-10) || T < T_L ) done_flag = 1;

  } while ( !done_flag );
  

  free(EE);
  free(w);
  return;
}


void debug_block_metro(void)
{
  struct matrix L;
  extern struct cellprm cell;

  L = L_e3(&cell);

  /* check whether the fractional coords are out of bounds */
  check_objects( );

  /* print the cell parms */
  printf("* cell: %15.10f%15.10f%15.10f  %15.10f%15.10f%15.10f\n",cell.a,cell.b,cell.c,
	 R2D*cell.alph,R2D*cell.beta,R2D*cell.gamm);
  
  /* print coordinates of objects */
  print_all_objs( );
  
  /* print out the center to vertex distances */
  print_c2v( );
  return;
}

