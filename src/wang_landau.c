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
extern char restfile[25];

int wl_reinit=0;
int n_steps=0;  /* This is made external for automatic resetting of the WL_INIT 
		   when a new energy goes below the lowest pre-set bin number.
		   Setting this to a large value will force an exit of the loops. */
#ifdef MPI
double *HSTprev, *DOSprev, *HSTinc, *DOSinc;
#include "mpi.h"
#endif

extern int efunc,numat,num_ob;
extern int wl_done_flag;
extern int anneal_sched;
extern int ecol,mvout;
extern int ortho;
extern int debug;
extern int halt_feature;
extern int autoadjust_flag;
extern int autoadjust_chng;
extern int simp_restarts,simp_nmax;
extern int non_per_ewald;
extern int wl_nEbins,wl_iter_print;
extern int wl_weighted_dos;
extern int wl_simp_nmax,wl_lnfmod_style;
extern int wl_hst_flat_method;
extern int wl_prt_all_poscars;
extern int *SWPtab;

extern double dplane;
extern double cvd_max,cell_vol_max;
extern double hold_acc;
extern double max_ang_chng,aspect_max,sqsh_vol_pct;
extern double autoadjust_pct,autoadjust_pcnt;
extern double wl_fmin_conv;
extern double wl_lnfmod_init,wl_emax_mult,wl_emin_mult;
extern double wl_qual_fact,wl_hst_flat_pct,wl_hst_flat_window,wl_hst_rpt_pct;
extern double wl_bin_width,wl_en_max,wl_prt_pos_tol;
extern double move_prob[MAX_MOVE_TYPE];

extern struct cellprm cell;
extern struct atom *s;
extern struct obj *object;
extern struct Energy E_best;

struct Energy E_wl_reinit={0};

extern FILE *tmpfp;
extern FILE *restartinput;

extern int numat;

/* added for faster deltaE Ewald sum */
extern int use_deltaEion;


// for MPI
extern int process_num, process_rank;


/******************************************************************/
/*                                                                */
/*                    Wang-Landau algorithm                       */
/*                                                                */
/******************************************************************/
void Wang_Landau(int runs,
		 int t_runs,
		 int trn_chng,
		 int rot_chng,
		 int lat_chng,
		 int swp_chng,
		 double cell_min,
		 double cell_max,
		 double trans_frac_max,
		 double volum_frac_max)
{

  char K; /* for printing out: T in temp autoadjust, M in Metro */

  int i,vi,order_loop;
  int obj_number,swp_number=0;
  extern int n_steps;
  int tot_trn_attempts=0,tot_lat_attempts=0;
  int tot_rot_attempts=0,tot_swp_attempts=0;
  int tot_attempts[MAX_MOVE_TYPE]={0}, rejects[MAX_MOVE_TYPE]={0};
  int nnze=0;
  int conf_changes=0;
  int T_runs=t_runs,T_iteration=0;
  int cell_chng_rand,order_rand_num,op1=0,op2=0,op3=0,op4=0,OP=0;
  int reject_trn=0,reject_lat=0,reject_rot=0,reject_swp=0;
  int ok,lat_not_ok_loops=0;
  int dplane_violation=0;
  int fflush_error=0;
  int qual_flag=0;
  int exceed_trans_attempts_REJECT=0;
  int n_trans_attempts=0,swappable=0;
  int mm_ok;

  double T=0;
  double trans_frac=trans_frac_max;
  double volum_frac=volum_frac_max;
  double Eave,SigEsq,SigE,E2ave,Eave2;
  double dx,dy,dz;
  double dsx,dsy,dsz;
  double da,db,dc,dalph=0.0,dbeta=0.0,dgamm=0.0;
  double da_frac,db_frac,dc_frac;
  double vol,dvol;
  double phi,the,psi;
  double trn_pct=0,swp_pct=0,lat_pct=0,rot_pct=0;
  double rej_pct[MAX_MOVE_TYPE]={0};
  double tot_time=0;
  static double min_cell_volume=0;

  struct Energy En={0},Elast,Eprt={0},Etmp={0};
  struct vector trans_vec;
  struct vector s1={0},s2={0},swap_to_vec,swap_from_vec;
  struct vector cell_chng,ang_chng;
  struct vector TP={0};
  struct matrix R,L, Linv;
  struct obj obj_orig={{0}};
  struct rusage r_usage;
  struct hstdat hstdata;
  
  static int been_called=0;

  /* new variables for WL */
  int nEbins=wl_nEbins;
  int jj,bin_of_minE=0,found_new_bin_of_minE=0;
  int Emax_bin=0,Emin_bin=0;
  int structure_count=1;
  int nz=0;
  int number_hst=1;
  int tf_sign=1;

  double Emax=50000,Emin= -100,Ebin_width=0;
  double old_lognE=0,DlognE=0,tot_lognE=0,lognEl=0,lognEg=0;
  double hst_flat_pct=0;
  double lnfmod=wl_lnfmod_init;
  double sum_pct=0;
  double E_low=1e6;

  struct endos *DOS={0},*HST={0};
  struct str_wl_init WL_INIT={0};

#ifdef MPI
  double global_E_best;
#endif

  if ( !been_called ) {

    min_cell_volume = SPHERE_VOL * ( num_ob * pow(cvd_max,3.0) );
    if ( debug > 2 ) printf("* Minimum cell volume: %10.3f\n",min_cell_volume);

    En.ecc=0; En.eng=0;
    been_called=1;
  }


//  if ( trans_frac_max < 0 ) trans_frac_max = fabs(trans_frac_max);

  /* initialize the WL-routine */
  WL_INIT    = wl_init(DOS,HST);
  nEbins     = WL_INIT.nEbins;
  Emin_bin   = WL_INIT.Emin_bin;
  Emax_bin   = WL_INIT.Emax_bin;
  Emin       = WL_INIT.Emin;
  Emax       = WL_INIT.Emax;
  Ebin_width = WL_INIT.Ebin_width;
  E_low      = WL_INIT.E_low;
  DOS        = WL_INIT.DOS;
  HST        = WL_INIT.HST;

  if ( debug>1 ){
    printf("DOS pointer = %p\n",DOS);
    printf("HST pointer = %p\n",HST);
    printf("nEbins = %d\n",nEbins);
    printf("Emax_bin = %d\n",Emax_bin);
    printf("Emin_bin = %d\n",Emin_bin);
    printf("Emax = %f\n",Emax);
    printf("Emin = %f\n",Emin);
    printf("Ebin_width= %f\n",Ebin_width);
    printf("E_low = %f\n",E_low);
  }
  
  if ( is_master() ) {
    if ( ecol==0 ) {
      if ( efunc==0 ) printf("%-3s%3s%10s%12s%12s%10s%10s%10s%8s%8s%30s%23s%8s%11s\n",
                             "*H","itr","lnfmod","E_ion","E_ss","minEbin","SigE","ortho","a/c","b/c","% Rejected Changes","%vst","%flat","DlognE");
      if ( efunc==1 ) printf("%-3s%3s%10s%12s%12s%10s%10s%10s%8s%8s%30s%23s%8s%11s\n",
                             "*H","itr","lnfmod","E_ion","E_lj","minEbin","SigE","ortho","a/c","b/c","% Rejected Changes","%vst","%flat","DlognE");
      printf("%-3s%83s%9s%9s%9s%9s%9s\n","*H"," ","trn","rot","lat","swp","mm");
      fflush(NULL);
    }
  }

  /* initialize the seed for random number generation */
  /* srand48(1234); removed this line in version 4.6.1.0 */

  /*
   * The structure of this algorithm is:
   *
   *     MAIN LOOP (variable iteration)
   *
   *             MAIN CONFIG LOOP (variable n_steps)
   *
   *                     choose an object 
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
   *     END MAIN LOOP
   *
   */

  /***************************************/
  /*         This is the main loop       */
  /***************************************/
  wl_done_flag=0;
  do {

    da_frac = db_frac = dc_frac = trans_frac;
    if ( debug > 3 ) printf("* trans_frac = %.3f\n",trans_frac);
    
    /* initialize variables for temperature run */
    n_steps=0;

    tot_trn_attempts=0;
    tot_rot_attempts=0;
    tot_lat_attempts=0;
    tot_swp_attempts=0;
    memset(tot_attempts, 0, sizeof(int)* MAX_MOVE_TYPE);
    reject_trn=0;
    reject_rot=0;
    reject_lat=0;
    reject_swp=0;
	memset(rejects, 0, sizeof(int)* MAX_MOVE_TYPE);

    Eave = 0;
    E2ave = 0;
    Eave2 = 0;
    SigEsq = 0;

    /* DEBUG */
    if ( debug > 3 ){
      printf("* Top of W-L loop, at T_iteration %d\n",T_iteration);
      debug_block_WL();
    }

    /* find initial system energy */
    En = (*Total_Energy)();
    Elast = En;
    if ( debug > 3 ) {
      printf("* W-L: ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
      printf("* W-L: init: E: %10.6e     T_iteration = %d\n",Elast.tot,T_iteration);
    }
    

    /*********************************/
    /* MAIN CONFIG LOOP (RUNS loops) */
    /*********************************/
    do {
      
      /* choose any random object */
      obj_number = (int)( drand48() * num_ob );
      if ( debug > 3 ) printf("* ************* W-L: CHOOSE OBJECT: no. %2d\n",obj_number);
      
      /*****************************
       * Below are the configuration changes:
       * 
       * 1) translational movement of an object
       *
       * 2) rotation of a polygon
       *
       * 3) changes in the lattice parameters
       *
       * 4) swap objects
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
		  /* 
		   choose move type by probability
		   */
	if (move_prob[2*MAX_MOVE_TYPE-1]>0) OP=RandomChoice(move_prob,MAX_MOVE_TYPE);
		  
		  
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
	  /* checking constraints and non_per_ewald */
	  while ( 0==trans_accept( obj_number, vadd(TP,trans_vec), non_per_ewald ) && exceed_trans_attempts_REJECT==0 ) { 
	    n_trans_attempts++;
	    /* if ( n_trans_attempts > 0.9*MAX_TRANS_ATTEMPTS) printf("* n_trans_attempts= %d\n",n_trans_attempts);  REMOVE ME */
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
	      /* printf("* exceeded n_trans_attempts\n"); debug=0; REMOVE ME */
	      trans_frac = trans_frac_max + trans_frac_max*TRANS_FRAC_INCREASE_FACTOR*tf_sign; tf_sign *= -1;
	      if ( debug > -1 ){
		printf("* WL: TRANS: WARN: --------->> CHANGING trans_frac by %2d * %.2f to %10.5e\n",tf_sign,TRANS_FRAC_INCREASE_FACTOR,trans_frac);
		printf("* WL: TRANS: WARN: Exceeded max number of attempts at translation (%d).\n"
		       "* WL: TRANS: WARN: Try decreasing the range of accessible high energies through\n"
		       "* WL: TRANS: WARN: wl_emax_mult.\n",(int)MAX_TRANS_ATTEMPTS);
		printf("* WL: TRANS: WARN: obj_number = %d\n",obj_number);
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
	    
	  } /* end while : checking constraints and non_per_ewald */
	  
	  if ( debug > 3 ) {
	    printf("* WL: n_trans_attempts= %d\n",n_trans_attempts);
	    if ( exceed_trans_attempts_REJECT==1 ){
	      printf("* REJECT CHANGE: TRANS - rejected because exceeded trans_attempts max\n");
	    }
	  }
	  
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: TRANS  obj_number=[%d]\n",obj_number);
	    printf("* dx dy dz = %15.9f%15.9f%15.9f\n",dx,dy,dz);
	    printf("* TP       = %15.10f%15.10f%15.10f\n",TP.x,TP.y,TP.z);
	    if ( TP.x>1.0 || TP.x<0.0 || TP.y>1.0 || TP.y<0.0 || TP.z>1.0 || TP.z<0.0 ){
	      printf("* W-L: WARN: atom is out of bounds!\n");
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
          //use_deltaEion=1;
	  En = (*Total_Energy)();
          //use_deltaEion=0;

	  if ( debug > 3 ) {
	    printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
	    printf("* trans: Elast.tot = %10.6e\tEn.tot = %10.6e\n",Elast.tot,En.tot);
	  }

	  /* Is the move accepted ? */
	  if ( !wl_accept_move( Emin_bin, bin_of_minE , Elast.tot, En.tot, Emax, Ebin_width, lnfmod, DOS, HST ) ||
	       exceed_trans_attempts_REJECT==1 ) {
	    reject_trn++;	
	    if ( debug > 3 ) printf("* REJECT CHANGE: TRANS\n");
	    
	    trans_vec = makevec( -trans_vec.x , -trans_vec.y , -trans_vec.z );

	    if ( exceed_trans_attempts_REJECT==0 ) {
	      trans_object( obj_number, &L, trans_vec );
	    }
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
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
	    printf("* CONFIG CHANGE: ROT object [%d]\n",obj_number);
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
	  
	  /* Is the move accepted ? */
	  if ( !wl_accept_move( Emin_bin, bin_of_minE, Elast.tot, En.tot, Emax, Ebin_width, lnfmod, DOS, HST ) ) {
	    reject_rot++;	
	    if ( debug > 3 ) printf("* REJECT CHANGE: ROT\n");
	    
	    /* replace the rotated object with the original */
	    object[obj_number] = obj_orig;
	    for(vi=0;vi<MAX_VERTS && object[obj_number].vused[vi]==1;vi++){
	      object[obj_number].vmoved[vi] = 1; /* flag for energy calc */
	    }
	    
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
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
	   (1) should not scale the bond lengths within the objects.
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
	      printf("* FATAL ERROR: lat change code in w-l routine:\n");
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

	  /* what are the energy consequences */
	  mark_all_objects_moved();
	  En = (*Total_Energy)();

	  if ( debug > 3 ) {
	    printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
	    printf("* cell chng: Elast.tot = %10.6e  En.tot = %10.6e\n",Elast.tot,En.tot);
	  }

	  /* Is the move accepted ? */
	  /* Also reject if the packing fraction is unreasonably high (cell squashing) */
	  /* Also reject if the dplane criterion is violated (cell flattening) */
	  /* Also reject if the cell volume is becoming too large */
	  if ( ( !wl_accept_move( Emin_bin, bin_of_minE, Elast.tot, En.tot, Emax, Ebin_width, lnfmod, DOS, HST  ) ) ||
	       ( cell_volume(&cell) < min_cell_volume ) ||
	       ( En.orth < sqsh_vol_pct ) ||
	       ( cell_volume(&cell) >= cell_vol_max ) ) {
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
	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}
	

	/*****************************************************************************************/
	/*                                  **  SWAPPING  **                                     */
	/*****************************************************************************************/	
	
	/* 4. swapping of objects */
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
	    swappable = SWPtab[ object[obj_number].typ * SWPTAB + object[swp_number].typ ];
	    /* printf("obnum=%d swpnum=%d typ1=%d typ2=%d swappable=%d\n",
	       obj_number,swp_number,object[obj_number].typ,object[swp_number].typ,swappable); */
	  } while ( swp_number == obj_number || swappable == 0 );
	  s2 = object[swp_number].v[0];

	  if ( debug > 3 ) printf("* ***loc1a******* SWAP OBJ %2d WITH OBJ %2d \n",obj_number,swp_number);

	  dsx = s2.x - s1.x;
	  dsy = s2.y - s1.y;
	  dsz = s2.z - s1.z;
	  swap_to_vec   = makevec(  dsx,  dsy,  dsz );
	  swap_from_vec = makevec( -dsx, -dsy, -dsz );
	  if ( debug > 3 ) {
	    printf("* CONFIG CHANGE: SWAP\n");
	    printf("* swp: s1 %9.4f%9.4f%9.4f   with:\n"
		   "       s2 %9.4f%9.4f%9.4f\n",
		   s1.x,s1.y,s1.z,
		   s2.x,s2.y,s2.z);
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
	    printf("* Swapped structure:\n");
	    debug_block_WL();
	  }

	  /* Is the move accepted ? */
	  if ( !wl_accept_move( Emin_bin, bin_of_minE, Elast.tot, En.tot, Emax, Ebin_width, lnfmod, DOS, HST ) ) {
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
	      printf("* swp back chk: tot= %10.6e ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",
		     Etmp.tot,Etmp.ecc,Etmp.eng,Etmp.vol);
	      printf("* Swapped back structure:\n");
	      debug_block_WL();
	    }

	    
	    /* reset energy */
	    En = Elast; /* restore the energy structure to prev. value */
	  }
	  /* if translation change was ACCEPTED set Elast struct
	     for the next configuration change */
	  Elast = En;

	}

    

		  /*****************************************************************************************/
		  /*                                 ** OP==4 Min-Map  **                                  */
		  /*****************************************************************************************/
		  
		  /* min-map */
		  if ( OP==4 ) {
            while ((object[obj_number].vfxd[0]==1) || (object[obj_number].num_vert >1)) {
              obj_number = (int)( drand48() * num_ob );
            }
            tot_attempts[OP] += 1;
            /* get the current location of the object */
            TP = object[obj_number].v[0];
            
            // FZ: ????
            /* if we are here, we have to recalculate L, the transformation
             matrix, becuase we may have changed the lattice parameters
             in the lat parm change below. */
            L = L_e3(&cell);
            Linv = Linverse(&L);
            mm_ok=MinMapMove(obj_number, &En, & L, & Linv);
			  

            if ( debug > 3 ) print_c2v( );

            En = (*Total_Energy)();
			  
            if ( debug > 3 ) {
				  printf("* ecc = %10.6e   eng = %10.6e   vol = %10.6f\n",En.ecc,En.eng,En.vol);
				  printf("* trans: Elast.tot = %10.6e\tEn.tot = %10.6e\n",Elast.tot,En.tot);
            }
			  
			  /* Is the move accepted ? */
            if ((!wl_accept_move( Emin_bin, bin_of_minE , Elast.tot, En.tot, Emax, Ebin_width, lnfmod, DOS, HST ))
                ||(!mm_ok)) {
              rejects[OP]++;	
              if ( debug > 3 ) printf("* REJECT CHANGE: min-map\n");
              trans_vec = vsub( TP, object[obj_number].v[0]);
              trans_object( obj_number, &L, trans_vec );
              /* reset energy */
              En = Elast; /* restore the energy structure to prev. value */
            }
			  Elast = En;
			  
		  }

		  
	/* DEBUG */
	if ( debug > 3 ){
	  printf("* Post config change\n");
	  debug_block_WL();
	  get_min_distances( );
	}

	if ( mvout ) print_xbs_mvframe();

      } /* end order loop */


      /* REMOVE ME  if ( check_wall_violation() > 0 ){
	printf("* !!!! WALL violation. exiting.\n");
	exit_pack(0);
	} */

      n_steps++;
    } while ( n_steps < runs );

#ifdef MPI
// communication among all processes
    if (process_num>1) { // no need to bother for single process job
//       MPI_Allreduce ( &bin_of_minE, &jj, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
//printf("debug> process %d: current Emin at %d overall at %d\n", process_rank, bin_of_minE,jj);
//       if (jj>bin_of_minE) bin_of_minE=jj;

       for(jj=0;jj<nEbins;jj++) {
        DOSinc[jj]= DOS[jj].lngE - DOSprev[jj];
        HSTinc[jj]= HST[jj].lngE - HSTprev[jj];
       }
//printf(">debug DOSinc[700]= %f HSTinc[700]=%f DOS[700=%f HST[700]=%f\n", DOSinc[700], HSTinc[700], DOS[700].lngE, HST[700].lngE);
       MPI_Allreduce ( DOSinc, DOSprev, nEbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
       MPI_Allreduce ( HSTinc, HSTprev, nEbins, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

       MPI_Allreduce ( &E_best.ecc_eng, &global_E_best, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
       E_best.ecc_eng = global_E_best;

       for(jj=0;jj<nEbins;jj++) {
            DOS[jj].lngE+= DOSprev[jj]-DOSinc[jj];
            DOSprev[jj]=DOS[jj].lngE;
            HST[jj].lngE+= HSTprev[jj]-HSTinc[jj];
            HSTprev[jj]=HST[jj].lngE;
       }
//printf(">debug DOSinc[700]= %f HSTinc[700]=%f DOS[700=%f HST[700]=%f\n", DOSinc[700], HSTinc[700], DOS[700].lngE, HST[700].lngE);
    }
#endif

    if ( wl_reinit== 1 ){
      /* re-initialize the WL-routine */

      wl_read_restart();

      WL_INIT    = wl_init(DOS,HST);
      nEbins     = WL_INIT.nEbins;
      Emin_bin   = WL_INIT.Emin_bin;
      Emax_bin   = WL_INIT.Emax_bin;
      Emin       = WL_INIT.Emin;
      Emax       = WL_INIT.Emax;
      Ebin_width = WL_INIT.Ebin_width;
      E_low      = WL_INIT.E_low;
      DOS        = WL_INIT.DOS;
      HST        = WL_INIT.HST;
      wl_reinit=0;
    }

    /************ Wang-Landau LOOP FINISHED ********************/

    /* Calculate the energy at the end of the WL loops */
    Eprt = (*Total_Energy)();
    if ( debug > 3 ) printf("* WL: end of runs loop: tot= %10.6e  ecc= %10.6e  eng=%10.6e\n",
			    Eprt.tot,Eprt.ecc,Eprt.eng);

    /* Energy fluctuations calculations */
    for (i=0,nnze=0; i<Emin_bin; i++) {
      if ( HST[i].lngE < 1 ) continue;
      Eave   += HST[i].lngE;
      E2ave  += HST[i].lngE * HST[i].lngE;
      nnze++;
    }
    Eave  /= nnze;
    E2ave /= nnze;
    
    for (i=0; i<Emin_bin; i++) {
      if ( HST[i].lngE < 1 ) continue;
      SigEsq += (HST[i].lngE - Eave)*(HST[i].lngE - Eave);
    }
    SigE = sqrt(SigEsq);
    Eave2 = Eave*Eave;

    if ( debug>3 ) printf("* Eave= %f   E2ave= %f  Eave2= %f\n",
			  Eave,E2ave,Eave2);

    lognEg=0; lognEl=0;
    old_lognE = tot_lognE;
    for(i=0; i<Emax_bin; i++) lognEg += DOS[i].lngE;
    for(i=Emax_bin; i<Emin_bin; i++) lognEl += DOS[i].lngE;
    tot_lognE = lognEl+lognEg;
    DlognE = tot_lognE - old_lognE;

    /*********************************************************/


    /*    if ( DlognE<1 && debug>1 ){
      printf("* WL: not sampling anything. DlognE= %e ??  Exiting.\n",DlognE);
      printf("*PP Begin POSCAR N %d  ecc=%20.10e  eng=%20.10e  etot=%20.10e\n",
	     structure_count,Eprt.ecc,Eprt.eng,Eprt.tot);
      print_poscar_out(0);
      exit_pack(0);
      }*/
    
    conf_changes = runs;

    if ( debug>3 ) {
      printf("* conf_changes=%d     4*runs=%d\n",conf_changes,4*runs);
      printf("* nnze=%d (number of non-zero energies)\n",nnze);
      printf("* tot_trn_attempts = %10d\n",tot_trn_attempts);
      printf("* tot_rot_attempts = %10d\n",tot_rot_attempts);
      printf("* tot_lat_attempts = %10d\n",tot_lat_attempts);
      printf("* tot_swp_attempts = %10d\n",tot_swp_attempts);
	  printf("* tot_attempts  = %10d %10d %10d %10d %10d\n",
			 tot_attempts[0],tot_attempts[1],tot_attempts[2],tot_attempts[3],tot_attempts[4]);
      printf("* reject_trn = %10d\n",reject_trn);
      printf("* reject_rot = %10d\n",reject_rot);
      printf("* reject_lat = %10d\n",reject_lat);
      printf("* reject_swp = %10d\n",reject_swp);
	  printf("* reject  = %10d %10d %10d %10d %10d\n",
			   rejects[0],rejects[1],rejects[2],rejects[3],rejects[4]);
    }

    if ( tot_trn_attempts == 0 ) trn_pct = 0;
    else trn_pct = (double)reject_trn/tot_trn_attempts*100;
    if ( tot_rot_attempts == 0 ) rot_pct = 0;
    else rot_pct = (double)reject_rot/tot_rot_attempts*100;
    if ( tot_lat_attempts == 0 ) lat_pct = 0;
    else lat_pct = (double)reject_lat/tot_lat_attempts*100;
    if ( tot_swp_attempts == 0 ) swp_pct = 0;
    else swp_pct = (double)reject_swp/tot_swp_attempts*100;
    for (i=0; i<MAX_MOVE_TYPE; i++) {
		if (tot_attempts[i]<=0) rej_pct[i]=0;
		else rej_pct[i]=(double)rejects[i]/tot_attempts[i]*100;
	}

    switch ( autoadjust_chng ){
    case 1: autoadjust_pct = trn_pct/100.0; break;
    case 2: autoadjust_pct = rot_pct/100.0; break;
    case 3: autoadjust_pct = lat_pct/100.0; break;
    case 4: autoadjust_pct = swp_pct/100.0; break;
    }


    T_iteration++;

    /* if ( T_iteration== 1849 ) debug=5;    REMOVE ME */
    hstdata = wl_hst_flat( HST , DOS , Emin_bin );
    hst_flat_pct = hstdata.flt;
    for(jj=1;jj<Emin_bin;jj++){
      if ( DOS[jj].lngE > 1 && jj>bin_of_minE ) {
        bin_of_minE = jj;
        found_new_bin_of_minE=1;
      }
    }

// resetting lnfmod to init value to fasten convergence
#ifdef MPI
    MPI_Allreduce ( &found_new_bin_of_minE, &jj, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
    found_new_bin_of_minE=jj;
#endif
//    if (found_new_bin_of_minE && (2>bin_of_minE*hst_flat_pct)) {
    if (found_new_bin_of_minE && (lnfmod< wl_lnfmod_init)) {
      lnfmod=wl_lnfmod_init;
      found_new_bin_of_minE=0;
      if (is_master()) printf("* Found new minE bin %d; resetting lnf to original %f\n", bin_of_minE,lnfmod);
    }


    /*****************************************************/
    /*       Calculate the quality factor                */
    /*****************************************************/

    /* Print the structure if the rej pcts are all high */
    nz=0;
    /* count the number that are turned off */
    if ( tol( trn_pct , 0.0 , 0.001 ) ) nz++;
    if ( tol( rot_pct , 0.0 , 0.001 ) ) nz++;
    if ( tol( lat_pct , 0.0 , 0.001 ) ) nz++;
    if ( tol( swp_pct , 0.0 , 0.001 ) ) nz++;
    if ( tol( rej_pct[4] , 0.0 , 0.001 ) ) nz++;

    /* If the sum percent is say 395/400 or 295/300, then we
       are in a low energy region and we should print out
       the structure! */
    sum_pct = trn_pct + rot_pct + lat_pct + swp_pct + rej_pct[4];

    if ( debug > 0 )
      printf("* WL:  qual_fact = %f   sum_pct= %f (need sum_pct > qual_fact)\n",
	     (400-nz*100)*wl_qual_fact, sum_pct);

    if ( ( sum_pct > (100*MAX_MOVE_TYPE-nz*100)*wl_qual_fact ) || debug>0 ) {
      if ( debug > -1 ) printf("* WL: setting qual_flag to 1\n");
      qual_flag=1;
    }
    else {
      qual_flag=0;
    }


    /*********************************************************************/
    /*        RECENTER THE CLUSTER and SHORT (OPTIONAL) SIMPLEX          */
    /*********************************************************************/
    
    if ( non_per_ewald == 1 ) { /* ONLY RECENTER of non_per_ewald == 1 */

      /*  If we are in region where the rej% are high and energy
          is close to E_low, then simplex and recenter the np ... OR...
          If the iter_print has passed, it's probably time to recenter  */
      
      if ( ( qual_flag==1 && ( Eprt.tot < E_low*WL_SIMP_WINDOW ) ) ||
	   ( T_iteration%wl_iter_print==1 ) ){
	
	if ( debug > 2 ) printf("* (optional) Centering and (optional) Simplex section of WL\n");
	
	/* recenter();	 RECENTER the nanoparticle */

	/* VERY SHORT simplex relaxation here !!! */
	if ( wl_simp_nmax>0 ){
	  do {
	    Simplex( WL_SIMP_RESTARTS , wl_simp_nmax );
	    Etmp = (*Total_Energy)();
	  } while ( Etmp.tot > 0.0 );
	}
	
	/* this call to Total_Energy is in case of energy change
	   when re-centering or simplexing the cluster above */
	Eprt = (*Total_Energy)();
	
      }
      
    }

    /***************************************************************************/
    /*                  PRINT OUT STRUCTURE IF ENERGY IS LOW                   */
    /***************************************************************************/
    /* If the structure has the lowest energy so far, then print it
       to output! */

    if ( 
	wl_prt_all_poscars == 1 || 
	( qual_flag==1 && wl_prt_pos_tol>0.0 && Eprt.tot<(E_low-fabs(E_low*wl_prt_pos_tol))  ) || 
	( qual_flag==1 && wl_prt_pos_tol<0.0 && tol( Eprt.tot, E_low, wl_prt_pos_tol ) && Eprt.tot<E_low )  )
      {

      if ( debug > -1 ) {
	printf("* WL: E_low= %f   Eprt= %f   Elow-Eprt= %f\n",
	       E_low,Eprt.tot,E_low-Eprt.tot);
      }


      /* set all 'moved' flags to one */
      mark_all_objects_moved();
      Eprt = (*Total_Energy)();
      if ( debug > 0 ) {
	printf("* WL: --------> qual_flag == 1 <---------\n");
	printf("* WL: Eprt.tot= %.20f  E_low= %.20f  diff= %.20f\n",Eprt.tot,E_low,E_low-Eprt.tot);
	printf("* WL: copying struct to restfile\n");
	printf("* WL: wl_simp_nmax= %d\n",wl_simp_nmax);
      }

      /* SAVE the structure for further minimization BEFORE the
      simplex relaxation:
 
      This is necessary in case the simplex just pushes
      the structure back into local icosahedral symmetry
      (for LJ minimization).  It is also the procedurally
      correct method... We don't want to bias the structure
      with simplex. */

      print_restart();

      /* perform VERY LONG simplex relaxation for output structures !!! */
      if ( wl_simp_nmax > 0 ){
	do {
	  Simplex( WL_SIMP_RESTARTS_OUTPUT , WL_SIMP_NMAX_OUTPUT );
	  Etmp = (*Total_Energy)();
	} while ( Etmp.tot > 0.0 );
      }
      
      /* get energy after simplex */
      Eprt = (*Total_Energy)();
      if ( getrusage(RUSAGE_SELF,&r_usage) == -1 )
	printf("* Error getting timer information.\n");

      tot_time=((double)r_usage.ru_utime.tv_sec/60);
      printf("*PP Begin POSCAR N %10d  ecc= %12.5e  tot= %20.10e  pf=%10.2e  b/a=%6.3f  c/a=%6.3f  orth=%6.3f  time=%.2fm\n",
	     structure_count++,Eprt.ecc,Eprt.tot,Eprt.pf,Eprt.b/Eprt.a,Eprt.c/Eprt.a,Eprt.orth,tot_time);
      print_poscar_out(0);
      /* set new E_low if Simplex didn't muck anything up */
      if ( Eprt.tot < E_low || tol( Eprt.tot, E_low, 1e-9 )==1 ) {
	E_low = Eprt.tot;
	printf("* --->  New E_low =  %.20f\n",E_low);
      } else if ( E_low < Eprt.tot ) {
	if ( wl_simp_nmax>0 ) {
	  printf("* WL: Simplex raised energy of particle!! (overlap?)  Will not set new E_low!\n");
	  printf("* WL: E_low= %.20f    simplex= %.20f    low-simplex= %.20f\n",E_low,Eprt.tot,E_low-Eprt.tot);
	}
      }

      /* Replace the simplexed structure with the one we had before to
	 continue the W-L algo. */

      wl_read_restart();
      
      /* set all 'moved' flags to one */
      mark_all_objects_moved();
      Elast = (*Total_Energy)();
      if ( debug > 0 ) printf("* WL: qual_flag==1, replaced struct from restart file. Etot= %.20f\n",Elast.tot);
      
    }


    /********************************************************************/
    /*               PRINTING WL output line in output file             */
    /********************************************************************/
    /* if ( tol( Elast.eng, 44.06, 0.01) ) exit_pack(0); REMOVE ME */
    if ( autoadjust_flag == 1 ) K='T';
    else K='W';

    if ( (ecol && qual_flag && debug > -2) ||
	 (ecol && T_iteration%wl_iter_print==1) ){
      printf("*%c------------------------------------------------------------\n",K);
      printf("*%c %20s%12d%25s\n",K,"T_iter  =",T_iteration,"T_iteration");
      printf("*%c %20s%12.3e%25s\n",K,"T  =",T,"temperature");
      printf("*%c %20s%12.3e%25s\n",K,"ecc  =",Elast.ecc,"electrostatic energy");
      printf("*%c %20s%12.3e%25s\n",K,"eng  =",Elast.eng,"repulsive energy");
      printf("*%c %20s%12.3e%25s\n",K,"bnd  =",Elast.bnd,"bounds");
      printf("*%c %20s%12.3f%25s\n",K,"vol  =",Elast.vol,"cell volume");
      printf("*%c %20s%12.3f%25s\n",K,"lat a  =",Elast.a,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"lat b  =",Elast.b,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"lat c  =",Elast.c,"lattice vector");
      printf("*%c %20s%12.3f%25s\n",K,"ar  =",aspect(),"aspect ratio");
      printf("*%c %20s%12.3f%25s\n",K,"ortho  =",Elast.orth,"cv/ortho");
      printf("*%c %20s%12.3f%25s\n",K,"alph  =",Elast.alph*R2D,"angle 2--3");
      printf("*%c %20s%12.3f%25s\n",K,"beta  =",Elast.beta*R2D,"angle 1--3");
      printf("*%c %20s%12.3f%25s\n",K,"gamm  =",Elast.gamm*R2D,"angle 1--2");
      printf("*%c %20s%12.3e%25s\n",K,"SigEsq  =",SigEsq,"energy sigma");
      printf("*%c %20s%12.2f%25s\n",K,"trn  =",trn_pct,"% rej  translations");
      printf("*%c %20s%12.2f%25s\n",K,"rot  =",rot_pct,"% rej     rotations");
      printf("*%c %20s%12.2f%25s\n",K,"lat  =",lat_pct,"% rej lattice parms");
      printf("*%c %20s%12.2f%25s\n",K,"swp  =",swp_pct,"% rej   ca-an swaps");
      printf("*%c %20s%12.3e%25s\n",K,"best tot  =",E_best.tot,"may include cell vol");
      printf("*%c %20s%12.3e%25s\n",K,"best ecc  =",E_best.ecc,"electrostatic");
      printf("*%c %20s%12.3e%25s\n",K,"best eng  =",E_best.eng,"repulsive");
      printf("*%c %20s%12.3e%25s\n",K,"best vol  =",E_best.vol,"best volume");
      
    } else if ( ((qual_flag && debug > -2) ||
		 (T_iteration%wl_iter_print==1)) && is_master()) {
      printf("*%c %5d%10.2e%+12.3e%+12.3e%6d%12.3e%10.3f%8.2f%8.2f%9.2f%9.2f%9.2f%9.2f%9.2f%7.0f%8.1f%12.2e\n",K,
	     T_iteration,lnfmod,Elast.ecc,Elast.eng,bin_of_minE,SigE,
	     Elast.orth,Elast.a/Elast.c,Elast.b/Elast.c,trn_pct,rot_pct,lat_pct,swp_pct,rej_pct[4],hstdata.vst,100*hst_flat_pct,DlognE);
    }

    if ( 0 != (fflush_error = fflush(NULL)) ) {
      printf("* Warn: wl: fflush returns non-zero\n");
    }

    /* END PRINTING */

    /* CHECK FOR FLAT HST CODE */
    if ( (hst_flat_pct>wl_hst_flat_pct) && (lnfmod < wl_fmin_conv) ) {
      if (is_master()){
      printf("*wl HST check: %%flat=%f  lnfmod = %e\n",hst_flat_pct,lnfmod);
      printf("*wl HST is flat! lnfmod = %20.5e. Within tol! Exiting!\n",lnfmod);
      printf("*wl HST +++ Done. +++\n");
    }
      wl_done_flag=1;
    } else if ( hst_flat_pct>wl_hst_flat_pct && lnfmod > wl_fmin_conv ) {
      if (is_master()) {
      print_dos(DOS,Emin_bin,1);
      print_hst(HST,Emin_bin,1,number_hst);
      }
      reset_HST( HST , Emin_bin );  hst_flat_pct=0; number_hst++;
      lnfmod = modify_lnfmod(lnfmod, wl_lnfmod_style);
      T_iteration=0;
      if ( (debug >= -2)&&is_master()  )
        printf("*wl HST is flat!  New lnfmod = %20.5e\n",lnfmod);
    } else if ( hst_flat_pct<wl_hst_flat_pct &&
		lnfmod > wl_fmin_conv &&
		T_iteration >= T_runs ) {
      if (is_master()) {
      print_dos(DOS,Emin_bin,1);
      print_hst(HST,Emin_bin,1,number_hst);
      }
      reset_HST( HST , Emin_bin );  hst_flat_pct=0; number_hst++;
      lnfmod = modify_lnfmod(lnfmod, wl_lnfmod_style);
      T_iteration=0;
      if ( (debug >= -2)&& is_master()) {
        printf("*wl HST is not flat!\n");
        printf("*wl Exceeded max number of passes at lnmodf value.\n"
               "*wl Decreasing lnfmod even though HST is not flat.\n"
               "*wl New lnfmod = %20.5e\n",lnfmod);
      }
    } else if ( hst_flat_pct<wl_hst_flat_pct &&
		lnfmod > wl_fmin_conv &&
		T_iteration < T_runs ) {
//      print_dos(DOS,Emin_bin,1);
//      print_hst(HST,Emin_bin,1,number_hst);
/* modified by Fei */
if ( (T_iteration%wl_iter_print==1)&&(is_master()) ) {
      print_dos(DOS,Emin_bin,1);
      print_hst(HST,Emin_bin,1,number_hst);
}
      if ( (debug > 0)&&(is_master()) ) printf("*wl HST is not flat! Continuing. T_iteration=%d\n",
			      T_iteration);
    } else if ( hst_flat_pct<wl_hst_flat_pct &&
		lnfmod < wl_fmin_conv &&
		T_iteration >= T_runs ) {
      if (is_master()) printf("*wl Exceeded max number of passes at lnmodf value.\n"
	     "*wl HST is not flat, and lnfmod is too low. Exiting.\n");
      wl_done_flag=1;
    }


  } while ( !wl_done_flag );
  
  if(is_master() ){
  print_dos(DOS,Emin_bin,1);
  print_hst(HST,Emin_bin,1,number_hst);
  }

  free(DOS);
  free(HST);

#ifdef MPI
  free(HSTprev);
  free(DOSprev);
  free(HSTinc);
  free(DOSinc);
#endif

  return;
}



/***************************************/
/***************************************/
struct str_wl_init wl_init(struct endos *DOS, struct endos *HST)
{
  
  int i;
  int nEbins=wl_nEbins;
  int Emax_bin,Emin_bin;
  extern int wl_reinit;

  double Emin= -100,Emax=50000,E_low=1e6,Ebin_width=0;

  static struct str_wl_init WL_INIT={0};
  struct Energy En;
  
  printf("*\n**********************************************************************\n*\n");
  En = (*Total_Energy)();
  if ( wl_reinit==0 ) printf("*          - Wang-Landau Init -\n*\n");
  if ( wl_reinit==1 ) printf("*          - Wang-Landau RE-Init -\n*\n");
  printf("*          Total Energy = %15.5e ecc=%15.5e\teng=%15.5e\n*\n",En.tot,En.ecc,En.eng);
  if ( En.eng > 0.0 && efunc==0 ) {
    printf("* WARN! WL-init: Overlap still exists before start.\n"
	   "*                simp_lamb too big or too small, or not enough simplex steps,\n"
	   "*                or there is non-removable overlap on the objects themselves.\n");
  }

  Emin =  En.tot - fabs(En.tot) * wl_emin_mult;
  if ( !tol( wl_en_max, WL_EN_MAX, 1e-2 ) ){

    Emax = wl_en_max;

  } else {

    Emax = En.tot + fabs(En.tot) * wl_emax_mult;

  }

  E_low = En.tot;

  /* if the bin width was not set in the input file, use default number
     of energy bins */
  if ( tol( wl_bin_width, WL_BIN_WIDTH, 1e-6 ) ) {
    Ebin_width = (Emax-Emin) / (double)nEbins;
  } else {
    Ebin_width=wl_bin_width;
    nEbins = (int)( (Emax-Emin)/Ebin_width );
  }
  if ( nEbins > WL_MAX_EN_BINS ){
    printf("* WARN: exceeded %d energy bins! Consider increasing bin width!\n",WL_MAX_EN_BINS);
  }
  
  printf("**********************************************************************\n");
  printf("*\n* Wang-Landau MC with the following parameters:\n*\n");
  printf("* %20s%15.3e%25s\n","Emax  =",Emax,"Energy maximum");
  printf("* %20s%15.3e%25s\n","Emin  =",Emin,"Energy minimum");
  if ( Emax < 0.0 ) {
    printf("*\n*           WARN: WL-MC init routine. Emax < 0 !\n*\n");
  }

  if ( !tol( wl_en_max, WL_EN_MAX, 1e-2 ) ){
    printf("* %20s%15.3e%25s\n","wl_en_max  =",wl_en_max,"wl_en_max");
  } else {
    printf("* %20s%15.3e%25s\n","wl_emax_mult  =",wl_emax_mult,"emax multiplier");
  }
  printf("* %20s%15.3e%25s\n","wl_emin_mult  =",wl_emin_mult,"emin multiplier");
  printf("* %20s%15.3e%25s\n","Ebin_width  =",Ebin_width,"Energy bin width");
  printf("* %20s%15d%25s\n","nEbins  =",nEbins,"Number energy bins");
  Emax_bin = get_bin_number(Emax,Emax,Ebin_width);
  Emin_bin = get_bin_number(Emin,Emax,Ebin_width);
  printf("* %20s%15d%25s\n","Emax bin  =",Emax_bin,"");
  printf("* %20s%15d%25s\n","starting bin  =",get_bin_number(En.tot,Emax,Ebin_width),"");
  printf("* %20s%15d%25s\n","Emin bin  =",Emin_bin,"");
  printf("*\n");
  printf("* %20s%15.3e%25s\n","wl_hst_flat_pct  =",wl_hst_flat_pct,"flat pct criterion");
  printf("* %20s%15.3e%25s\n","wl_qual_fact  =",wl_qual_fact,"rej pct criterion");
  
  printf("* %20s%15.3e%25s\n","wl_fmin_conv  =",wl_fmin_conv,"fmin converg criterion");
  printf("*\n");
  
  printf("**********************************************************************\n");
  
  
  /* Allocate memory for DOS */
  if ( DOS ) free ( DOS );
  DOS = (struct endos *)malloc( nEbins * sizeof(struct endos) );
  if ( !DOS ) {
    printf("%s\n","Not enough memory. Wang-Landau loc DOS.\n");
    exit_pack(0);
  }
  
  /* Allocate memory for HISTOGRAM */
  if ( HST ) free ( HST );
  HST = (struct endos *)malloc( nEbins * sizeof(struct endos) );
  if ( !HST ) {
    printf("%s\n","Not enough memory. Wang-Landau loc HST.\n");
    exit_pack(0);
  }

#ifdef MPI
  if ( DOSprev ) free ( DOSprev );
  if ( DOSinc ) free ( DOSinc );
  if ( HSTprev ) free ( HSTprev );
  if ( HSTinc ) free ( HSTinc );

  DOSprev = (double *)malloc( nEbins * sizeof(double) );
  DOSinc = (double *)malloc( nEbins * sizeof(double) );
  HSTprev = (double *)malloc( nEbins * sizeof(double) );
  HSTinc = (double *)malloc( nEbins * sizeof(double) );
  if ( (!DOSprev)||(!DOSinc)||(!HSTprev)||(!HSTinc) ) {
    printf("%s\n","Not enough memory. Wang-Landau loc.\n");
    exit_pack(0);
  }
#endif

 
  /************************************************************************/
  
  /* initialize the DOS and HST */
  for(i=0; i<Emin_bin; i++) {
    DOS[i].lngE= 1;
    DOS[i].refen = Emax - i*Ebin_width;
    
    HST[i].lngE= 0;
    HST[i].refen = Emax - i*Ebin_width;

#ifdef MPI
    HSTinc[i]=HSTprev[i]=DOSinc[i]=0;
    DOSprev[i]=1; // see above: init DOS is 1 everywhere!
#endif
  }
  /* added by Fei to mark the starting bin */
  /*  DOS[get_bin_number(En.tot,Emax,Ebin_width)].lngE= 2;
      HST[get_bin_number(En.tot,Emax,Ebin_width)].lngE= 1; */

  
  /************************************************************************/
  
  WL_INIT.nEbins     = nEbins;
  WL_INIT.Emax_bin   = Emax_bin;
  WL_INIT.Emin_bin   = Emin_bin;
  WL_INIT.Emin       = Emin;
  WL_INIT.Emax       = Emax;
  WL_INIT.Ebin_width = Ebin_width;
  WL_INIT.E_low      = E_low;
  WL_INIT.DOS        = DOS;
  WL_INIT.HST        = HST;

  return(WL_INIT);
}

/**************************************/
/*  Accept or reject WL move          */
/**************************************/
int wl_accept_move(int Emin_bin,
		   int bin_of_minE,
		   double en,
		   double en_try,
		   double Enmax,
		   double Ebin_width,
		   double lnfmod,
		   struct endos *DOS,
		   struct endos *HST)
{

  int dbg=0; /* INSANE amount of output */
  int bin,bin_try,bin_en_ret,bin_tr_ret;
  static int reinit_structure_count=0;
  double lnDOS_diff;
  double prob_trans=0;

  if ( debug > 4 ) printf("* wl_accept_move: en=%f en_try=%f Enmax=%f\n",en,en_try,Enmax);
  
  bin     = get_bin_number(en,    Enmax,Ebin_width);
  bin_try = get_bin_number(en_try,Enmax,Ebin_width);

  if ( bin<0 ){
    printf("* FATAL ERROR: en= %10.3e  bin= %d   en_try= %10.3e  bin_try= %d\n",
	   en,bin,en_try,bin_try);
    print_info_block();
    exit_pack(0);
  }

  //printf("Emin_bin = %d   bin= %d    bin_try= %d\n",Emin_bin,bin,bin_try);
  if ( bin_try > Emin_bin ) {
    E_wl_reinit = (*Total_Energy)();
    printf("* Holy shit. The energy is lower than your lowest bin. Re-initializing the WL routine...\n");
    printf("* Saving this configuration to the restart file and printing the POSCAR.\n");
    printf("*RR restart output RN %10d\n",reinit_structure_count);
    print_restart();
    print_restart_inline();

    printf("*PP Begin POSCAR RN %10d  ecc= %12.5e  tot= %20.10e  pf=%10.2e  b/a=%6.3f  c/a=%6.3f  orth=%6.3f\n",
	     reinit_structure_count++,E_wl_reinit.ecc,E_wl_reinit.tot,E_wl_reinit.pf,E_wl_reinit.b/E_wl_reinit.a,E_wl_reinit.c/E_wl_reinit.a,E_wl_reinit.orth);
    print_poscar_out(0);

    wl_reinit=1;
    n_steps=WL_REINIT; /* force the current loops to finish */
    fflush(stdout);
    return(0);
  }

  if (( bin_try > 0 ) && ( en_try <= Enmax )){

    if ( dbg>0 ) {
      printf("*\n");
      printf("* wl_accept_move: bin crt= %13d  bin try= %13d\n",bin,bin_try);
    }

    /* criterion is min[ g1/g2 , 1 ]
       Now ln(g1/g2) = lng1-lng2, and ln(1)=0,
       so if lng1-lng2>0 then return 1 */
    lnDOS_diff = DOS[bin].lngE - DOS[bin_try].lngE;

    if (lnDOS_diff<0) prob_trans = exp( lnDOS_diff );

    if ( dbg>0 ) {
      printf("* wl_accept_move: e_crt= %15.3e  e_try= %15.3e\n",en,en_try);
      printf("* wl_accept_move: lnDOS= %15.3e  lnDOS= %15.3e\n",DOS[bin].lngE,DOS[bin_try].lngE);
      printf("* wl_accept_move: lnDOS_diff= %15.3e  exp[lnDOS_diff]= %10.4e  ",lnDOS_diff,exp(lnDOS_diff));
    }

  } else {

    lnDOS_diff = -999;

  }

  if ( (lnDOS_diff >= 0) || (prob_trans>drand48()) ) {

    bin_tr_ret = bin_energy_E(Emin_bin,bin_of_minE,en_try,Enmax,Ebin_width, lnfmod, DOS, HST);
    if ( bin_tr_ret == -1 ){
      printf("* !!!! WARN: en_try= %10.3e, bin_try<0\n",en_try);
      fflush(stdout);
    }

    bin_energy_hst( bin_try , HST );
    if ( dbg>0 ) printf("* gtry<gcrt: returning 1\n");
    return(1);
    
  } else {
    
    bin_en_ret = bin_energy_E(Emin_bin,bin_of_minE,en,    Enmax,Ebin_width, lnfmod, DOS, HST);
    if ( bin_en_ret == -1 ){
      printf("* !!!! WARN: en= %10.3e, bin<0\n",en);
      fflush(stdout);
    }

    bin_energy_hst( bin , HST );
    if ( dbg>0 ) printf("* gtry>gcrt: returning 0\n");
    return(0);
    
  }
}

/* 
 * Energy Binning Code
*/

int get_bin_number(double en,
		   double Emax,
		   double Ebin_width)

{
  int bin;
  bin = (int)( (Emax-en)/Ebin_width );
  //printf("------> bin number is: %d\n",bin);
  /* REMOVE ME */
  /* if ( bin== -1050101 ) {
    printf("********************* Emax=%12.4e, en= %12.4e, width= %f  bin_num= %d\n",
	   Emax, en, Ebin_width, bin);
    print_poscar_out(0);
    print_restart_inline();
    exit_pack(0);
    } */
  return(bin);
}

int bin_energy_E(int Emin_bin,
		 int bin_of_minE,
		 double en,
		 double Emax,
		 double Ebin_width,
		 double lnfmod,
		 struct endos *DOS,
		 struct endos *HST)
{
  int bin;
  extern int n_steps,wl_reinit;
  double weighted_contr;
  
  bin = get_bin_number(en, Emax, Ebin_width);
  
  if ( bin < 0 ) { /* energy is higher than highest bin */

    DOS[0].lngE += lnfmod;
    return (-1);

  } else {         /* energy is lower than lowest bin */

    if ( bin > Emin_bin ) {
      printf("* WARN: Wang-Landau\n");
      printf("* WARN: Error: Max bin number exceeded. Energy very low - re-setting HST and DOS!!\n");
      printf("* WARN: en=%e  Emax=%e  Ebin_width=%e\n",en,Emax,Ebin_width);
      printf("* WARN: Emin_bin=%d   bin=%d\n",Emin_bin,bin);
      printf("* WARN: Try setting wl_emin_mult to larger magnitude value to avoid this.\n");
      fflush(stdout);
      exit_pack(0);
      
#ifdef MPI
      //      MPI_Allreduce ( &wl_reinit, &tmp, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
      // simply exit
      printf("FATAL ERROR: debug> process %d: abort due to very low E_min_bin=%d\nConsider a lower bound", process_rank, bin);
      MPI_Abort(MPI_COMM_WORLD,-1);
#endif      
      
    }
    
    if ( wl_weighted_dos==1 && bin_of_minE>0 ) {

      weighted_contr = ( lnfmod - (lnfmod-wl_fmin_conv)*bin/(bin_of_minE));
      /* printf("weighted_contr = %f   bin_of_minE=%d\n",weighted_contr,bin_of_minE); */
      DOS[bin].lngE += weighted_contr;

    } else {

      DOS[bin].lngE += lnfmod;

    }

    return (1);
  }

}

/**************************************/
void bin_energy_hst(int bin, struct endos *HST)
{
  HST[bin].lngE += 1;
}

/****************************************/
struct hstdat wl_hst_flat(struct endos *HST, struct endos *DOS, int nmax)
{
  int i,bin_minE=0;
  int num_flat=0;
  int non_zero_bins=0;
  double hst_ave=0,hst_sig=0;
  double pct_flat=0;
/* added by Fei */
  int n_occupied_bins=0;
  
  struct hstdat dat;

  dat.vst=0; dat.flt=0;

  /* bin_minE is here the number of histogram bars to be checked */

  for(i=1;i<nmax;i++){
    if ( DOS[i].lngE > 1 && i>bin_minE ) { 
       bin_minE = i;
/* modified by Fei 
   only occupied sites are accounted for, since there might be gaps in the DOS
*/
       n_occupied_bins++; }
  }
  /* calculate average and std dev of HST */
// Fei: should be i<=bin_minE, not i<bin_minE
  for(i=1;i<=bin_minE;i++){
    if ( HST[i].lngE>0 ) non_zero_bins++;
    hst_ave += HST[i].lngE;
  }
//  hst_ave /= bin_minE;
  if (!n_occupied_bins) {
    printf("* WARN: no occupied bin found!\n");
    return(dat);
  }
  hst_ave /= n_occupied_bins;
  for(i=1;i<=bin_minE;i++){
    hst_sig += ( hst_ave - HST[i].lngE )*( hst_ave - HST[i].lngE );
  }
//  hst_sig = sqrt(hst_sig/bin_minE);
  hst_sig = sqrt(hst_sig/n_occupied_bins);

  if ( debug > 3 ) printf("* hst_ave= %e  hst_sig= %e  bin_minE = %d  nnzrbins= %d\n",
			  hst_ave, hst_sig, bin_minE,non_zero_bins);

  if ( wl_hst_flat_method == 0 ){
    for(i=1;i<=bin_minE;i++){
      if ( HST[i].lngE > wl_hst_flat_window * hst_ave ) num_flat++;
    }
  } else if ( wl_hst_flat_method == 1 ) {
    for(i=1;i<=bin_minE;i++){
      if ( HST[i].lngE > wl_hst_flat_window * hst_ave - hst_sig &&
	   HST[i].lngE < wl_hst_flat_window * hst_ave + hst_sig    ) num_flat++;
    }
  }


//  pct_flat = (double)num_flat/bin_minE;
  dat.flt = (double)num_flat/n_occupied_bins;

  /* if there are too few populated bins in the histogram */
//  if ( ( (double)non_zero_bins/bin_minE ) < wl_hst_rpt_pct ) return(0);
  dat.vst = 100*non_zero_bins/n_occupied_bins;
  return(dat);

}


/****************************************/
void reset_HST(struct endos *HST, int nmax)
{
  int i;
#ifdef MPI
  if (is_master()){
#endif
  printf("*wl Reseting HST.\n");
#ifdef MPI
  }
#endif
  for(i=1;i<nmax;i++) HST[i].lngE = 0;
#ifdef MPI
  for(i=1;i<nmax;i++) {
    HSTprev[i]=0;
  }
#endif
  return;
}

/****************************************/
void debug_block_WL(void)
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

void wl_read_restart(void)
{
  if (is_master())
    {
      restartinput = fopen(restfile,"r");
      if ( restartinput == NULL )
        {
          printf("* WL:  Error opening 'rest' file!\n");
          printf("* %20s%35s\n","restfile =",restfile);
          exit_pack(0);
        }
      get_init_parms(restartinput,0);    
      fclose(restartinput);
    }
}
