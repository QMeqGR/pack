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

extern int autoadjust;
extern int best_compare_off;
extern int debug;
extern int restart;
extern int runs,t_runs;
extern int wl_runs,wl_t_runs;
extern int efunc,forces;
extern int non_per_ewald;
extern int retain_temp_files;
extern int simp_nmax,simp_restarts;
extern int use_wlmc,wl_nEbins,wl_iter_print;
extern int wl_weighted_dos;
extern int wl_prt_dos,wl_prt_hst,wl_prt_all_poscars;
extern int wl_simp_nmax,wl_lnfmod_style,wl_hst_flat_method;
extern int wl_short_init,wl_simp_restarts_init,wl_simp_nmax_init;
extern int *SWPtab;

extern double cvd_max,cell_vol_max,cell_vol_max_fact;
extern double c2vtol,dplane;
extern double simp_errlim,simp_lamb;
extern double errlim;
extern double ext_press,aspect_max;
extern double basin_hop_scale_factor_hi;
extern double trans_frac;
extern double sqsh_vol_pct,rep_epsilon;
extern double wl_fmin_conv,wl_lnfmod_init;
extern double wl_emax_mult,wl_emin_mult,wl_qual_fact;
extern double wl_hst_rpt_pct,wl_hst_flat_window,wl_hst_flat_pct;
extern double wl_bin_width,wl_en_max,wl_prt_pos_tol,wl_cell_init_dx;
extern double wl_simp_ediff_init;
extern double LJeps[MAX_TABLE_SQR];
extern double LJrer[MAX_TABLE_SQR];

// added by Fei for min-map move type
extern double move_prob[MAX_MOVE_TYPE];
extern double MM_transfrac[2];

extern struct cellprm cell;
extern struct Energy E_best;
extern struct obj *init_obj,*object,*object_orig;
extern struct plane surf_plane;

extern FILE *input;


/***************************************/
void get_init_parms(FILE *getfile,int exearly)
{
  char *line=NULL,*token,*copy=NULL,*tag=NULL,*value=NULL;
  char *gaF=" ";

  int flag_cell_a=1,flag_cell_b=1,flag_cell_c=1;
  int flag_cell_alph=1,flag_cell_beta=1,flag_cell_gamm=1;
  int flag_obj=1;
  int flag_count=0;
  int cao,cav,gao,gav,gonum,govrt=0,gorep=0;
  int objnum=0;
  int swp_typ1,swp_typ2,swp_tab;
  int rest_onum,rest_nvrt,rest_nrep;
  static int been_here=0;
  

  int errval;
  int i,j;
  int Z1,Z2,gaZ,gafxd,allowrot;
  
  double cvd=0,feps,frer;
  double gax,gay,gaz,gaR,gachrg, tot_prob;
  double caxmin,caxmax,caymin,caymax,cazmin,cazmax;
  static double SGMA = 0.8908987184; /* kittel R/sigma minimum */

  struct vector center,vec;

  //  printf(" Debug level at 1 = %12d\n",debug);

  if ( been_here == 0 ) {

    /* allocate space for object structures */
    init_obj = (struct obj *)malloc( MAX_OBJS * sizeof(struct obj) );
    if ( !init_obj ) {
      printf("%s\n","Not enough memory for init_obj (init_parms.c)\n");
      exit_pack(0);
    }

    /* allocate space for general object structures */
    object = (struct obj *)malloc( MAX_OBJS * sizeof(struct obj) );
    if ( !object ) {
      printf("%s\n","Not enough memory for object (init_parms.c)\n");
      exit_pack(0);
    }

    /* allocate space for general object structures */
    object_orig = (struct obj *)malloc( MAX_OBJS * sizeof(struct obj) );
    if ( !object_orig ) {
      printf("%s\n","Not enough memory for object_orig (init_parms.c)\n");
      exit_pack(0);
    }

    //    printf(" Debug level at 2 = %12d\n",debug);

    /* initialize the numbr of repeats to zero */
    for(i=0;i<MAX_OBJS;i++){
      object[i].used = 0; /* set all to 0 for init */
      object[i].num_rept = 0;
      for(j=0;j<MAX_VERTS;j++){
	object[i].vused[j] = 0;
	object[i].Z[j] = init_obj[i].Z[j] = -1;
      }
    }
    for(i=0;i<MAX_INIT_OBJS;i++){
      init_obj[i].num_rept = 0;
    }

    /* allocate space for SWPtab matrix */
    SWPtab = (int *)malloc( MAX_TABLE_SWP * sizeof(int) );
    if ( !SWPtab ) {
      printf("%s\n","Not enough memory for SWPtab array.\n");
      exit_pack(0);
    }
    /* initialize the swap objects allow matrix */
    for(i=0;i<MAX_TABLE_SWP;i++){
      SWPtab[i] = 0;
    }

    /* initialize the Lennard-Jones matrices */
    for(i=0;i<MAX_TABLE;i++){
      for(j=0;j<MAX_TABLE;j++){
	LJeps[ i * MAX_TABLE + j ] = EPS_DEFAULT; /* 4.0 is the standard for homogeneous systems */
	LJrer[ i * MAX_TABLE + j ] = RER_DEFAULT; /* 0.0 forces div by zero error if undefined */
      }
    }
    been_here=1;
  }
  if ( exearly==1 ) {

    /* read the input file by parsing around '=' and '!' characters */
    /* the tag values lines should look like:  tag = value ! comment stuff... */
    line  = (char *)malloc( 100 * sizeof(char) );
    copy  = (char *)malloc( 100 * sizeof(char) );
    tag   = (char *)malloc( 100 * sizeof(char) );
    value = (char *)malloc( 100 * sizeof(char) );

    while( fgets(line,100,input) != NULL ) {

      if ( debug > 2 ) {
	printf("line string is: %s", line);
      }
      if ( NULL==strpbrk(line, "=") ) continue;
      if ( NULL!=strstr(line,"#") ) continue;

      copy = strcpy(copy,line);
      
      token = strtok(copy, "=!");
      if ( debug > 2 ) printf("tag = %s\n",token);
      tag = strcpy(tag,token);
      
      token = strtok(NULL, "=!"); 
      if ( debug > 2 ) printf("value = %s\n",token);
      value = strcpy(value,token);
      
      token = strtok(NULL, ""); 
      if ( debug > 2 ) printf("comment = %s\n",token);
      


      /*------------------------------------------------------------------
	CELL VARIABLES
	------------------------------------------------------------------*/

      if ( NULL!=strstr(tag,"cella") ){
	cell.a = atof(value); flag_cell_a=0;
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.a",cell.a);
      }
      if ( NULL!=strstr(tag,"cellb") ){
	cell.b = atof(value); flag_cell_b=0;
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.b",cell.b);
      }
      if ( NULL!=strstr(tag,"cellc") ){
	cell.c = atof(value); flag_cell_c=0;
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.c",cell.c);
      }
      if ( NULL!=strstr(tag,"alph") ){
	cell.alph = atof(value); flag_cell_alph=0;
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.alph",cell.alph);
      }
      if ( NULL!=strstr(tag,"beta") ){
	cell.beta = atof(value); flag_cell_beta=0;
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.beta",cell.beta);
      }
      if ( NULL!=strstr(tag,"gamm") ){
	cell.gamm = atof(value); flag_cell_gamm=0;
      if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell.gamm",cell.gamm);
      }


      /*------------------------------------------------------------------
	Surface vectors for optimized Ewald (surface calcs)
	------------------------------------------------------------------*/
      if ( NULL!=strstr(tag,"surf") ){

	surf_plane.u.x = atof( strtok(value," ") );
	surf_plane.u.y = atof( strtok(NULL," ") );
	surf_plane.u.z = atof( strtok(NULL," ") );

	surf_plane.v.x = atof( strtok(NULL," ") );
	surf_plane.v.y = atof( strtok(NULL," ") );
	surf_plane.v.z = atof( strtok(NULL," ") );
	
	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","surf");
	  printf("*                u.x    = %12f\n",surf_plane.u.x);
	  printf("*                u.y    = %12f\n",surf_plane.u.y);
	  printf("*                u.x    = %12f\n",surf_plane.u.z);
	  printf("*                v.x    = %12f\n",surf_plane.v.x);
	  printf("*                v.y    = %12f\n",surf_plane.v.y);
	  printf("*                v.x    = %12f\n",surf_plane.v.z);
	}
      }

      /*------------------------------------------------------------------
	Generalized object variables
	------------------------------------------------------------------*/

      if ( NULL!=strstr(tag,"gen_obj") ){

	objnum += gorep; /* this only works if the input file is in correct order */

	gonum = atoi( strtok(value," ") );
	govrt = atoi( strtok(NULL," ") );
	gorep = atof( strtok(NULL," ") );
	
	init_obj[gonum].typ      = gonum;
	init_obj[gonum].num_vert = govrt;
	init_obj[gonum].num_rept = gorep;
	flag_obj=0;

	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","gen_obj");
	  printf("*            typ/num    = %12d\n",gonum);
	  printf("*                vrt    = %12d\n",govrt);
	  printf("*                rep    = %12d\n",gorep);
	}
      }

      gafxd = 0; allowrot=1; /* identifies if there are fixed objects or vertices here */
      if ( NULL!=strstr(tag,"def_gen") ){
	gao = atoi( strtok(value," ") );
	gav = atoi( strtok(NULL," ") );
	gax = atof( strtok(NULL," ") );
	gay = atof( strtok(NULL," ") );
	gaz = atof( strtok(NULL," ") );
	gaR = atof( strtok(NULL," ") );
	gaZ = atof( strtok(NULL," ") );
	gachrg = atof( strtok(NULL," ") );
	gaF = strtok(NULL," ");
	if ( debug>1 && NULL!=gaF ) printf("Fixed char: %s gaf[]=%c\n",gaF,gaF[0]);
	if ( gaF!=NULL && gaF[0]=='F' ) gafxd = 1;
	if ( gafxd==1 && gav>0 ) allowrot=0;

	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","def_gen");
	  printf("*                a    = %12d\n",gao);
	  printf("*              typ    = %12d\n",gao);
	  printf("*                v    = %12d\n",gav);
	  printf("*                x    = %12.4f\n",gax);
	  printf("*                y    = %12.4f\n",gay);
	  printf("*                z    = %12.4f\n",gaz);
	  printf("*                R    = %12.4f\n",gaR);
	  printf("*                Z    = %12d\n",gaZ);
	  printf("*             chrg    = %12.4f\n",gachrg);
	  printf("*              fxd    = %12d\n",gafxd);
	  printf("*         allowrot    = %12d\n",allowrot);
	}

	for ( i=objnum; i<( objnum + init_obj[gao].num_rept ); i++) {
	  if ( debug>1 ) {
	    printf("* %5s Input: creating object %d, vertex %d"," ",i,gav);
	  }
	  fflush(stdout);
	  /* set object_orig here for when cell_init() is called more than once
	     and the objects are remade in make_obj.c */
	  object_orig[i].used        = object[i].used        = 1;
	  object_orig[i].typ         = object[i].typ         = init_obj[gao].typ;
	  object_orig[i].vused[gav]  = object[i].vused[gav]  = 1;
	  object_orig[i].vmoved[gav] = object[i].vmoved[gav] = 1;
	  object_orig[i].vfxd[gav]   = object[i].vfxd[gav]   = gafxd;
	  object_orig[i].allowrot    = object[i].allowrot    = allowrot;
	  object_orig[i].v[gav].x    = object[i].v[gav].x    = gax;
	  object_orig[i].v[gav].y    = object[i].v[gav].y    = gay;
	  object_orig[i].v[gav].z    = object[i].v[gav].z    = gaz;
	  object_orig[i].R[gav]      = object[i].R[gav]      = gaR;
	  object_orig[i].Z[gav]      = object[i].Z[gav]      = gaZ;
	  object_orig[i].vchrg[gav]  = object[i].vchrg[gav]  = gachrg;

	  /* set defaul constraints (none) */
	  object_orig[i].vcnstr[gav].xmin = object[i].vcnstr[gav].xmin = CONSTR_XMIN;
	  object_orig[i].vcnstr[gav].xmax = object[i].vcnstr[gav].xmax = CONSTR_XMAX;
	  object_orig[i].vcnstr[gav].ymin = object[i].vcnstr[gav].ymin = CONSTR_YMIN;
	  object_orig[i].vcnstr[gav].ymax = object[i].vcnstr[gav].ymax = CONSTR_YMAX;
	  object_orig[i].vcnstr[gav].zmin = object[i].vcnstr[gav].zmin = CONSTR_ZMIN;
	  object_orig[i].vcnstr[gav].zmax = object[i].vcnstr[gav].zmax = CONSTR_ZMAX;

	  if ( gav == 0 ) center = makevec(gax,gay,gaz);
	  vec = makevec(gax,gay,gaz);

	  cvd = vmag( vsub(vec,center) );
	  object[i].cvd[gav]   = cvd;
	  object[i].cvdsq[gav] = cvd*cvd;
	  if ( debug>1 ) printf("  cvd %5.3f\n",cvd);
	  if ( object[i].vfxd[gav] != 1 ) cvd_max = fmax( fmax( cvd_max, cvd ), gaR);

	  object[i].num_vert = govrt;
	  object[i].num_rept = gorep;
	}

      }

      if ( NULL!=strstr(tag,"swp_objs") ){
	swp_typ1 = atoi( strtok(value," ") );
	swp_typ2 = atoi( strtok(NULL," ") );
	swp_tab  = atoi( strtok(NULL," ") );

	SWPtab[ swp_typ1 * SWPTAB + swp_typ2 ] = swp_tab;
	SWPtab[ swp_typ2 * SWPTAB + swp_typ1 ] = swp_tab;

	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","swp_objs");
	  printf("*                swp_obj1    = %12d\n",swp_typ1);
	  printf("*                swp_obj2    = %12d\n",swp_typ2);
	  printf("*                swp_tab     = %12d\n",swp_tab);
	}
      }

      /*************** constraints *********************/
      if ( NULL!=strstr(tag,"at_constr") ){
	cao    = atoi( strtok(value," ") );
	cav    = atoi( strtok(NULL," ") );
	caxmin = atof( strtok(NULL," ") );
	caxmax = atof( strtok(NULL," ") );
	caymin = atof( strtok(NULL," ") );
	caymax = atof( strtok(NULL," ") );
	cazmin = atof( strtok(NULL," ") );
	cazmax = atof( strtok(NULL," ") );

	printf(" Debug level = %12d\n",debug);
	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","at_constr");
	  printf("*               ob    = %12d\n",cao);
	  printf("*                v    = %12d\n",cav);
	  printf("*             xmin    = %12.3f\n",caxmin);
	  printf("*             xmax    = %12.3f\n",caxmax);
	  printf("*             ymin    = %12.3f\n",caymin);
	  printf("*             ymax    = %12.3f\n",caymax);
	  printf("*             zmin    = %12.3f\n",cazmin);
	  printf("*             zmax    = %12.3f\n",cazmax);
	}

	for ( i=objnum; i<( objnum + init_obj[cao].num_rept ); i++) {
	  if ( debug>1 ) {
	    printf("* %5s Input: setting constraints object %d, vertex %d\n"," ",i,cav);
	  }
	  fflush(stdout);
	  /* set object_orig here for when cell_init() is called more than once
	     and the objects are remade in make_obj.c */
	  object_orig[i].vcnstr[cav].constrained = object[i].vcnstr[cav].constrained = 1;
	  object_orig[i].vcnstr[cav].xmin  = object[i].vcnstr[cav].xmin        = caxmin;
	  object_orig[i].vcnstr[cav].xmax  = object[i].vcnstr[cav].xmax        = caxmax;
	  object_orig[i].vcnstr[cav].ymin  = object[i].vcnstr[cav].ymin        = caymin;
	  object_orig[i].vcnstr[cav].ymax  = object[i].vcnstr[cav].ymax        = caymax;
	  object_orig[i].vcnstr[cav].zmin  = object[i].vcnstr[cav].zmin        = cazmin;
	  object_orig[i].vcnstr[cav].zmax  = object[i].vcnstr[cav].zmax        = cazmax;
	}

      }



      /*------------------------------------------------------------------
	OTHER RUN VARIABLES
	------------------------------------------------------------------*/

      if ( NULL!=strstr(tag,"c2vtol") ){
	c2vtol = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","c2vtol",c2vtol);
      }
      if ( NULL!=strstr(tag,"cell_vol_max_fact") ){
	cell_vol_max_fact = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","cell_vol_max_fact",cell_vol_max_fact);
      }
      if ( NULL!=strstr(tag,"dplane") ){
	dplane = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","dplane",dplane);
      }
      if ( NULL!=strstr(tag,"basin_hop_scale_factor_hi") ){
	basin_hop_scale_factor_hi = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10f\n","basin_hop_scale_factor_hi",basin_hop_scale_factor_hi);
      }
      if ( NULL!=strstr(tag,"non_per_ewald") ){
	non_per_ewald = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","non_per_ewald",non_per_ewald);
      }
      if ( NULL!=strstr(tag,"forces") ){
	forces = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","forces",forces);
      }
      if ( NULL!=strstr(tag,"debug") ){
	debug = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","debug",debug);
      }
      if ( NULL!=strstr(tag,"retain_temp_files") ){
	retain_temp_files = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","retain_temp_files",retain_temp_files);
      }
      if ( NULL!=strstr(tag,"runs")
	   && NULL==strstr(tag,"t_runs")
	   && NULL==strstr(tag,"wl_runs")
	   && NULL==strstr(tag,"wl_t_runs") ){
	runs = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","runs",runs);
      }
      if ( NULL!=strstr(tag,"t_runs")
	   && NULL==strstr(tag,"wl_t_runs") ){
	t_runs = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","t_runs",t_runs);
      }
      if ( NULL!=strstr(tag,"wl_runs") && NULL==strstr(tag,"wl_t_runs") ){
	wl_runs = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_runs",wl_runs);
      }
      if ( NULL!=strstr(tag,"wl_t_runs") ){
	wl_t_runs = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_t_runs",wl_t_runs);
      }
      if ( NULL!=strstr(tag,"errlim") 
	   && NULL==strstr(tag,"simp_errlim") ){
	errlim = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","errlim",errlim);
      }
      if ( NULL!=strstr(tag,"simp_errlim") ){
	simp_errlim = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","simp_errlim",simp_errlim);
      }
      if ( NULL!=strstr(tag,"simp_lamb") ){
	simp_lamb = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","simp_lamb",simp_lamb);
      }
      if ( NULL!=strstr(tag,"simp_nmax") 
	   && NULL==strstr(tag,"wl_simp_nmax")
	   && NULL==strstr(tag,"wl_simp_nmax_init") ){
	simp_nmax = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","simp_nmax",simp_nmax);
      }
      if ( NULL!=strstr(tag,"simp_restarts")
 	   && NULL==strstr(tag,"wl_simp_restarts_init") ){
	simp_restarts = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","simp_restarts",simp_restarts);
      }
      if ( NULL!=strstr(tag,"restart")
	   && NULL==strstr(tag,"wl_simp_restarts_init")
	   && NULL==strstr(tag,"simp_restarts") ){
	restart = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","restart",restart);
      }
      if ( NULL!=strstr(tag,"sqsh_vol_pct") ){
	sqsh_vol_pct = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","sqsh_vol_pct",sqsh_vol_pct);
      }
      if ( NULL!=strstr(tag,"rep_epsilon") ){
	rep_epsilon = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","rep_epsilon",rep_epsilon);
      }
      if ( NULL!=strstr(tag,"autoadjust") ){
	autoadjust = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","autoadjust",autoadjust);
      }
      if ( NULL!=strstr(tag,"ext_press") ){
	ext_press = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","ext_press",ext_press);
      }
      if ( NULL!=strstr(tag,"aspect_max") ){
	aspect_max = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","aspect_max",aspect_max);
      }
      if ( NULL!=strstr(tag,"trans_frac") ){
	trans_frac = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","trans_frac",trans_frac);
      }
      if ( NULL!=strstr(tag,"wl_fmin_conv") ){
	wl_fmin_conv = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_fmin_conv",wl_fmin_conv);
      }
      if ( NULL!=strstr(tag,"wl_en_max") ){
	wl_en_max = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_en_max",wl_en_max);
      }
      if ( NULL!=strstr(tag,"wl_qual_fact") ){
	wl_qual_fact = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_qual_fact",wl_qual_fact);
      }
      if ( NULL!=strstr(tag,"wl_hst_flat_pct") ){
	wl_hst_flat_pct = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_hst_flat_pct",wl_hst_flat_pct);
      }
      if ( NULL!=strstr(tag,"wl_cell_init_dx") ){
	wl_cell_init_dx = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_cell_init_dx",wl_cell_init_dx);
      }
      if ( NULL!=strstr(tag,"wl_prt_pos_tol") ){
	wl_prt_pos_tol = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_prt_pos_tol",wl_prt_pos_tol);
      }
      if ( NULL!=strstr(tag,"wl_hst_rpt_pct") ){
	wl_hst_rpt_pct = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_hst_rpt_pct",wl_hst_rpt_pct);
      }
      if ( NULL!=strstr(tag,"wl_hst_flat_window") ){
	wl_hst_flat_window = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_hst_flat_window",wl_hst_flat_window);
      }
      if ( NULL!=strstr(tag,"wl_emax_mult") ){
	wl_emax_mult = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_emax_mult",wl_emax_mult);
      }
      if ( NULL!=strstr(tag,"wl_emin_mult") ){
	wl_emin_mult = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_emin_mult",wl_emin_mult);
      }
      if ( NULL!=strstr(tag,"wl_lnfmod_init") ){
	wl_lnfmod_init = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_lnfmod_init",wl_lnfmod_init);
      }
      if ( NULL!=strstr(tag,"best_compare_off") ){
	best_compare_off = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","best_compare_off",best_compare_off);
      }
      if ( NULL!=strstr(tag,"wl_nEbins") ){
	wl_nEbins = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_nEbins",wl_nEbins);
      }
      if ( NULL!=strstr(tag,"wl_simp_nmax")
	   && NULL==strstr(tag,"wl_simp_nmax_init") ){
	wl_simp_nmax = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_simp_nmax",wl_simp_nmax);
      }
      if ( NULL!=strstr(tag,"wl_simp_nmax_init") ){
	wl_simp_nmax_init = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_simp_nmax_init",wl_simp_nmax_init);
      }
      if ( NULL!=strstr(tag,"wl_simp_restarts_init") ){
	wl_simp_restarts_init = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_simp_restarts_init",wl_simp_restarts_init);
      }
      if ( NULL!=strstr(tag,"wl_simp_ediff_init") ){
	wl_simp_ediff_init = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_simp_ediff_init",wl_simp_ediff_init);
      }
      if ( NULL!=strstr(tag,"wl_iter_print") ){
	wl_iter_print = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_iter_print",wl_iter_print);
      }
      if ( NULL!=strstr(tag,"wl_hst_flat_method") ){
	wl_hst_flat_method = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_hst_flat_method",wl_hst_flat_method);
      }
      if ( NULL!=strstr(tag,"wl_lnfmod_style") ){
	wl_lnfmod_style = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_lnfmod_style",wl_lnfmod_style);
	if ( wl_lnfmod_style > 1 || wl_lnfmod_style < 0 ){
	  printf("Fatal error in input file: Unknown wl_lnfmod_style = %d\n",wl_lnfmod_style);
	  exit_pack(0);
	}
      }
      if ( NULL!=strstr(tag,"wl_prt_dos") ){
	wl_prt_dos = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_prt_dos",wl_prt_dos);
      }
      if ( NULL!=strstr(tag,"wl_prt_hst") ){
	wl_prt_hst = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_prt_hst",wl_prt_hst);
      }
      if ( NULL!=strstr(tag,"wl_prt_all_poscars") ){
	wl_prt_all_poscars = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_prt_all_poscars",wl_prt_all_poscars);
      }
      if ( NULL!=strstr(tag,"wl_weighted_dos") ){
	wl_weighted_dos = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_weighted_dos",wl_weighted_dos);
      }
      if ( NULL!=strstr(tag,"wl_short_init") ){
	wl_short_init = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","wl_short_init",wl_short_init);
      }
      if ( NULL!=strstr(tag,"wl_bin_width") ){
	wl_bin_width = atof(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20.10e\n","wl_bin_width",wl_bin_width);
      }
      if ( NULL!=strstr(tag,"use_wlmc") ){
	use_wlmc = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","use_wlmc",use_wlmc);
      }
      if ( NULL!=strstr(tag,"efunc") ){
	efunc = atoi(value);
	if ( debug > 1 ) printf("* Input READING: %20s  = %20d\n","efunc",efunc);
      }
      //      printf(" Debug level at 3 = %12d\n",debug);
      if ( NULL!=strstr(tag,"move_prob") ){
      	  move_prob[0]=  atof( strtok(value," ") );
        for(i=1; i<MAX_MOVE_TYPE;i++){
        	 move_prob[i]= atof( strtok(NULL," ") );
			if(move_prob[i]<0) move_prob[i]=0;
        }
        tot_prob=0;
		  for(i=0; i<MAX_MOVE_TYPE;i++) {
			  tot_prob+= move_prob[i];
			  move_prob[MAX_MOVE_TYPE+i]= tot_prob;
		  }
		  if (tot_prob<=0) {
			  printf("ERROR: sum of move_prob <=0");
			  exit_pack(-1);
		  }
		  for(i=0; i<2*MAX_MOVE_TYPE;i++) move_prob[i]/= tot_prob;
          move_prob[2*MAX_MOVE_TYPE-1]=1.0; 
		  
      	}
      //      printf(" Debug level at 4 = %12d\n",debug);
      if ( NULL!=strstr(tag,"MM_transfrac") ){
        MM_transfrac[0]=  atof( strtok(value," ") );
        for(i=1; i<2;i++){
          MM_transfrac[i]= atof( strtok(NULL," ") );
        }
      }
      
      //      printf(" Debug level at 5 = %12d\n",debug);


      /*------------------------------------------------------------------
	Lennard-Jones VARIABLES
	------------------------------------------------------------------*/

      if ( NULL!=strstr(tag,"def_eps") ){
	Z1 =   atoi( strtok(value," ") );
	Z2 =   atoi( strtok(NULL," ") );
	feps = atof( strtok(NULL," ") );

	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","def_eps");
	  printf("*                Z1   = %d\n",Z1);
	  printf("*                Z2   = %d\n",Z2);
	  printf("*                Z1 * MAX_TABLE + Z2 = %d\n",Z1*MAX_TABLE+Z2);
	  printf("*                Z2 * MAX_TABLE + Z1 = %d\n",Z2*MAX_TABLE+Z1);
	  printf("*                feps = %f\n",feps);
	}

	LJeps[ Z1 * MAX_TABLE + Z2 ] = feps;
	LJeps[ Z2 * MAX_TABLE + Z1 ] = feps;

      }


      if ( NULL!=strstr(tag,"def_rer") ){
	Z1 =   atoi( strtok(value," ") );
	Z2 =   atoi( strtok(NULL," ") );
	frer = atof( strtok(NULL," ") );

	if ( debug > 1 ) {
	  printf("* Input READING: %20s\n","def_rer");
	  printf("*                Z1   = %d\n",Z1);
	  printf("*                Z2   = %d\n",Z2);
	  printf("*                Z1 * MAX_TABLE + Z2 = %d\n",Z1*MAX_TABLE+Z2);
	  printf("*                Z2 * MAX_TABLE + Z1 = %d\n",Z2*MAX_TABLE+Z1);
	  printf("*                frer = %f\n",frer);
	}

	LJrer[ Z1 * MAX_TABLE + Z2 ] = frer * SGMA;
	LJrer[ Z2 * MAX_TABLE + Z1 ] = frer * SGMA;

      }

      
      if ( debug > 2 ) printf("*------------------------------------------------------------\n");
    }

    /* convert common values to Radians, if the user put in degrees */
    
    if ( tol( cell.alph, 30.0, 1e-5) ) cell.alph = 0.5235987755982988;
    if ( tol( cell.beta, 30.0, 1e-5) ) cell.beta = 0.5235987755982988;
    if ( tol( cell.gamm, 30.0, 1e-5) ) cell.gamm = 0.5235987755982988;
    
    if ( tol( cell.alph, 45.0, 1e-5) ) cell.alph = 0.7853981633974483;
    if ( tol( cell.beta, 45.0, 1e-5) ) cell.beta = 0.7853981633974483;
    if ( tol( cell.gamm, 45.0, 1e-5) ) cell.gamm = 0.7853981633974483;
    
    if ( tol( cell.alph, 90.0, 1e-5) ) cell.alph = 1.5707963267948966;
    if ( tol( cell.beta, 90.0, 1e-5) ) cell.beta = 1.5707963267948966;
    if ( tol( cell.gamm, 90.0, 1e-5) ) cell.gamm = 1.5707963267948966;
    
    if ( tol( cell.alph, 120.0, 1e-5) ) cell.alph = 2.094395102393195;
    if ( tol( cell.beta, 120.0, 1e-5) ) cell.beta = 2.094395102393195;
    if ( tol( cell.gamm, 120.0, 1e-5) ) cell.gamm = 2.094395102393195;

    flag_count = ( flag_cell_a +flag_cell_b +flag_cell_c +
		   flag_cell_alph + flag_cell_beta + flag_cell_gamm + 
		   flag_obj );

    if (  0 != flag_count ) {
      printf("Error reading input file. flag_count = %d    Check variable names.\n",flag_count);
      printf("Exiting.\n");
      exit(0);
    }

    cell_vol_max = cell_vol_max_fact*cell.a*cell.b*cell.c;

    /* set swappable flag */
    if ( debug > 1 ) printf("* Setting swappable flags.\n");
    for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
      for( j=0; j<MAX_OBJS && object[j].used==1; j++ ){
	if ( SWPtab[ object[i].typ * SWPTAB + object[j].typ ] == 1 ) {
	  if ( debug>2 ) printf("* Swappable: typ= %d    typ= %d\n",object[i].typ,object[j].typ);
	  object[i].swappable=1; object[j].swappable=1;
	}
      }
    }

    return; /* end here if this is not a restart calculation */
  }

  /********************************************************/
  /* If this is a restart calculation, then this reads the
     restart.dat file, which has a different format */
  /********************************************************/

  printf("* Reading from restart file\n");
  
  /* get cell parameters */
  errval = fscanf(getfile," %lf %lf %lf %lf %lf %lf ",
		  &cell.a,&cell.b,&cell.c,&cell.alph,&cell.beta,&cell.gamm);
  if ( errval != 6 ) { printf("Error reading restart file at 1 errval=%d.\n",errval); exit_pack(0); }
  cell_vol_max = cell_vol_max_fact*cell.a*cell.b*cell.c;

  /* begin reading in the objects */

  while ( ( errval = fscanf(getfile," %d %d %d",&rest_onum,&rest_nvrt,&rest_nrep) ) == 3 )
    {
      object[rest_onum].num_vert        = rest_nvrt;
      object[rest_onum].num_rept        = rest_nrep;
      /* get object vertex information */
      for ( i=0; i<rest_nrep; i++ ){
	for ( j=0; j<rest_nvrt; j++ ){
	  errval = fscanf(getfile," %d %d %lf %lf %lf %lf %d %lf %d",
			  &gao,&gav,&gax,&gay,&gaz,&gaR,&gaZ,&gachrg,&gafxd);
	  if ( errval != 9 ) { printf("Error reading restart file at 3 errval=%d.\n",errval); exit_pack(0); }
	  object[gao].used        = 1;
	  object[gao].vused[gav]  = 1;
	  object[gao].vmoved[gav] = 1;
	  object[gao].v[gav].x    = gax;
	  object[gao].v[gav].y    = gay;
	  object[gao].v[gav].z    = gaz;
	  object[gao].R[gav]      = gaR;
	  object[gao].Z[gav]      = gaZ;
	  object[gao].vchrg[gav]  = gachrg;
	  object[gao].vfxd[gav]   = gafxd;
	}
      }
    }
  
  if ( debug > 1 ){
    printf("* get parms: filepointer: %p\n",getfile);
    printf("* get parms: cell: %15.10f%15.10f%15.10f\n",cell.a,cell.b,cell.c);
    printf("* get parms: cell: %15.10f%15.10f%15.10f\n",cell.alph,cell.beta,cell.gamm);
    printf("* get parms: exearly     = %d (1=don't get atom pos)\n",exearly);
  }
  
  return;
}

