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


/*******************************/
/* External (global) Variables */
/*******************************/

#ifndef __GLOBAL_VARS_H__
#define __GLOBAL_VARS_H__

char tmpname[L_tmpnam];
char pid_string[15];
char tempfile[25]="PR",tmp_string[9]="-best.tmp";
char restfile[25]="PR",rst_string[9]="-rest.tmp";

int anneal_sched=ANNEAL_SCHED;
int autoadjust=AUTOADJUST;
int autoadjusting_hi_temp=0;
int autoadjusting_lo_temp=0;
int autoadjust_using=AUTOADJUST_USING; /* 1=rejection percentages, 2=pseudo heat cap */
int autoadjust_chng=AUTOADJUST_CHNG;   /* adjust on change type: 1=trn 2=rot 3=lat 4=swp */
int autoadjust_chng_hi=AUTOADJUST_CHNG_HI;
int autoadjust_chng_lo=AUTOADJUST_CHNG_LO;
int autoadjust_flag=0;
int autoadjust_t_runs=AUTOADJUST_T_RUNS;
int autoadjust_t_runs_pcnt=AUTOADJUST_T_RUNS_PCNT;
int basin_hop=0,dsm=0;
int best_compare_off=0; /* for testing only */
int debug=0;
int done_flag=0,wl_done_flag=0;
int efunc=0;
int ecol=0;
int forces=0;
int get_estat=0;
int halt_feature=1;
int pf_anneal_flag = 0;
int mvout=0;
int non_per_ewald=0;
int num_ob=0;
int numat=0;
int *numZ;
int n_simp=0,simp_runs=0;
int ortho=ORTHO;
int restart=0;
int retain_temp_files=0;
int runs=RUNS,t_runs=TEMP_RUNS;
int wl_runs=WL_RUNS,wl_t_runs=WL_T_RUNS;
int scale_cell_do=0;
int simp_nmax=SIMP_NMAX;
int simp_restarts=SIMP_RESTARTS;
int use_wlmc=0;
int wl_nEbins=WL_NUM_EBINS;
int wl_iter_print=WL_ITER_PRINT;
int wl_weighted_dos=0;
int wl_prt_dos=0,wl_prt_hst=0;
int wl_prt_all_poscars=0;
int wl_simp_restarts_init=WL_SIMP_RESTARTS_INIT;
int wl_simp_nmax=WL_SIMP_NMAX;
int wl_simp_nmax_init=WL_SIMP_NMAX_INIT;
int wl_lnfmod_style=WL_LNFMOD_STYLE;
int wl_hst_flat_method=WL_HST_FLAT_METHOD;
int wl_short_init=0;
int lat_parm_chng=1;
int *SWPtab;

double autoadjust_pct_ht_hi=AUTOADJUST_PCT_HT_HI;
double autoadjust_pct_ht_lo=AUTOADJUST_PCT_HT_LO;
double autoadjust_pct_lt_hi=AUTOADJUST_PCT_LT_HI;
double autoadjust_pct_lt_lo=AUTOADJUST_PCT_LT_LO;
double autoadjust_pct;
double autoadjust_pcnt_ht_hi=AUTOADJUST_PCNT_HT_HI;
double autoadjust_pcnt_ht_lo=AUTOADJUST_PCNT_HT_LO;
double autoadjust_pcnt_lt_hi=AUTOADJUST_PCNT_LT_HI;
double autoadjust_pcnt_lt_lo=AUTOADJUST_PCNT_LT_LO;
double autoadjust_pcnt;
double aspect_max=ASPECT_MAX;
double cvd_max=0;
double cell_vol_max=0,cell_vol_max_fact=CELL_VOL_MAX_FACT;
double basin_hop_temp_alpha=0.0;
double basin_hop_scale_factor_hi=BASIN_HOP_SCALE_FACTOR_HI;
double c2vtol=C2VTOL;
double dplane= -1;
double errlim=ERRLIM;
double hold_acc=HOLD_ACC;
double max_ang_chng=MAX_ANG_CHNG;
double ext_press=EXT_PRESS;
double pCnt_global=0;
double rtol=0;
double rep_epsilon=REP_EPSILON;
double scale_factor=0.0;
double sqsh_vol_pct=SQSH_VOL_PCT;
double simp_lamb=SIMP_LAMB;
double simp_errlim=SIMP_ERRLIM;
double trans_frac=TRANS_FRAC_MAX;
double trans_frac_min=TRANS_FRAC_MIN;
/* fn holds (in order) cellparms, cations, anions (centers),
   euler angles for each anion */
double *fn=NULL,*rot_var=NULL;
double **sim_pp,*sim_p,*sim_y;
double *psum,*ptry;
double *etable,*ecctable;
double wl_cell_init_dx=WL_CELL_INIT_DX;
double wl_fmin_conv=WL_FMIN_CONV;
double wl_lnfmod_init=WL_LNFMOD_INIT;
double wl_en_max=WL_EN_MAX;
double wl_emax_mult=WL_EMAX_MULT;
double wl_emin_mult=WL_EMIN_MULT;
double wl_qual_fact=WL_QUAL_FACTOR;
double wl_bin_width=WL_BIN_WIDTH;
double wl_hst_rpt_pct=WL_HST_RPT_PCT;
double wl_hst_flat_window=WL_HST_FLAT_WINDOW;
double wl_hst_flat_pct=WL_HST_FLAT_PCT;
double wl_prt_pos_tol=WL_PRT_POS_TOL;
double wl_simp_ediff_init=WL_SIMP_EDIFF_INIT;
double LJeps[MAX_TABLE_SQR];
double LJrer[MAX_TABLE_SQR];
double *SShrd,*SShrd_2,*SShrd_6,*SShrd12,*SSmin;
double move_prob[2*MAX_MOVE_TYPE]={0.0};
double MM_transfrac[2];

pid_t processid;

struct cellprm cell; /* lattice parameters */
struct atom *p=NULL;
struct atom *s=NULL;
struct cation *at=NULL;
struct anion *an=NULL;
struct Energy E_best;
struct obj *init_obj=NULL; /* for initialization of complex anions v 4.15.x.x */
struct obj *object=NULL; /* for initialization of complex anions v 4.15.x.x */
struct obj *object_orig=NULL; /* for initialization of complex anions v 4.15.x.x */
struct plane surf_plane; /* for defining the surface plane for optimized Ewald routine */

extern char *optarg;
extern int optind,opterr,optopt;

FILE *input=NULL;
FILE *outfile=NULL;
FILE *restartfile=NULL;  /* prints 'restart' file, also used for WL, reads restart file */
FILE *restartinput=NULL; /* used in 'best' comparison, reads tmpfp */
FILE *tmpfp=NULL;
FILE *dosfp=NULL;
FILE *hstfp=NULL;

/* for MPI */
extern int process_num, process_rank;
#endif // __GLOBAL_VARS_H__
