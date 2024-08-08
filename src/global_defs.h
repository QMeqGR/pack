#ifndef __GLOBAL_DEFS_H__
#define __GLOBAL_DEFS_H__

#define INIT_DATA_1 20
#define INIT_DATA_2 24
#define SEED 1234

#define MAX_INIT_OBJS  12 /* max number of object types */
#define MAX_OBJS  200 /* max number of objects */
#define MAX_VERTS 100 /* max number of atoms on an object */
#define MAX_TABLE 150
#define MAX_TABLE_SQR 22500
#define SWPTAB 12   /* see next entry */
#define MAX_TABLE_SWP 144  /* max 12x12 */
#define ZMAX 150 /* max Z to allow */
#define CELL_VOL_MAX_FACT 2 /* timex the cell volume in the input file */
#define EPS_DEFAULT (4.0)
#define RER_DEFAULT (0.0) /* so undefined will force div by zero */
#define MAX_TRANS_ATTEMPTS (1000) /* 1e5 */
#define TRANS_FRAC_INCREASE_FACTOR 0.10  /* increase trans_frac by x% to avoid
					    problems with exceeding trans attempts */
#define TRANS_FRAC_NPE_FACTOR 0.10    /* use 10% of dplane for npe calculations */

#define TYP_CA_1 100  /* for categorization and setting radii etc */
#define TYP_CA_2 200
#define TYP_AN_1 1000
#define TYP_AN_1_V 1100
#define TYP_AN_2 2000
#define TYP_AN_2_V 2100

#define TINY (1e-12)
#define SMALL (1e-5)
#define C2VTOL (1e-6)
#define Bohr2Ang 0.529177249
#define PI (3.14159265358979323846264338327950288419716939937510)
#define PIo2 (1.5707963267949)
#define R2D ( 180.0 / PI )
#define D2R ( PI / 180.0 )
#define SPHERE_VOL 4.18879 /* 4*pi/3 */
#define ERRLIM 1e-6
#define REP_EPSILON 1.0  /* repulsion epsilon */
#define LJSGMA 0.8908987184
#define SELF_ENERGY_TOLERANCE (1e-7) /* to eliminate self energy of overlapping */
                                     /* vertices on the same anion */

#define PRESS_DVOL 1e-4  /* dvol for pressure calculation */
#define BEST_VOL 0.1
#define VOLUM_FRAC_MAX 0.07  /* max cell volume change per lat step */
#define ASPECT_MAX 20.0
#define ORTHO 0             /* orthorhomic restriction */ 
#define MAX_ANG_CHNG 0.174  /*  in RADS (latvec angles) */
#define SQSH_VOL_PCT 0.2

// Cell squashing parameters (dplane) and keeping things from the walls in nanoparticle
// conformation calcultions:
#define DPLANE 0.1     /* min distance between planes of the cell. units: Angstrom (try c2v dist) */
/* Set the following to 1.0 in version 5.2.0.0. Set dplane in the input file instead of using these. */
#define VIOLAT_NPE_DPLANE_WALL_CUSHION 1.0 /* non-periodic ewald dplane cushion, units Angstrom*/
#define ACCEPT_NPE_DPLANE_WALL_CUSHION 1.0 /* for trans_accept, units Angstrom */

#define CELL_INIT_MAX_TRY 10000000

#define CEN_MOVES 5
#define FORCE_STEPS 20
#define FORCE_STEP_SIZE_AN 0.01
#define FORCE_STEP_SIZE_CA 0.01

/* PACK_FRAC_MIN is set to -ext_press when this is input as a neg number */
#define EXT_PRESS 0.0 /* default is not to have external pressure */
#define PACK_FRAC_MIN 1e-4   /* 1.0 ==> effectively uses cell volume all the time */

#define AUTOADJUST_RUNS 1200
#define AUTOADJUST_RUNS_STEP 10
#define AUTOADJUST_RUNS_MAX 1500
#define AUTOADJUST_T_RUNS 3
#define AUTOADJUST_T_RUNS_PCNT 2

#define AUTOADJUST_TEMP_STEP 0.5
#define AUTOADJUST_TEMP_STEP_MULT 0.45
#define AUTOADJUST_TEMP_STEP_MIN 1e-6

#define AUTOADJUST 1        /* to use or not to use autoadjusting */
#define AUTOADJUST_USING 1  /* 1 = use rej percentages, 2 = use psedo heat capacity */

#define AUTOADJUST_CHNG 3       /* initial setting for rej pct change on lat */
#define AUTOADJUST_CHNG_HI 3    /* adjust on change type: 1=trn 2=rot 3=lat 4=swp */
#define AUTOADJUST_CHNG_LO 3    /* adjust on change type: 1=trn 2=rot 3=lat 4=swp */

#define AUTOADJUST_TEMP_HT_INIT 1.0
#define AUTOADJUST_PCT_HT_HI 0.650
#define AUTOADJUST_PCT_HT_LO 0.600
#define AUTOADJUST_PCNT_HT_HI 9e1
#define AUTOADJUST_PCNT_HT_LO 5e0

#define AUTOADJUST_TEMP_LT_INIT 0.02
#define AUTOADJUST_PCT_LT_HI 0.999
#define AUTOADJUST_PCT_LT_LO 0.950
#define AUTOADJUST_PCNT_LT_HI 1e8
#define AUTOADJUST_PCNT_LT_LO 1e7

#define AUTOADJUST_LO_T_SAFETY 3


#define RUNS 2000         /* to be multiplied by the number of dof */
#define TEMP_RUNS 15
#define TEMP_INIT 1.0
#define TEMP_FINL 1e-3
#define ANNEAL_SCHED 3
#define TRANS_FRAC_MAX 0.5
#define TRANS_FRAC_MIN 0.001 /* was 0.02 up until v 4.14.1.10 */
#define HOLD_ACC -1.0


#define SIMP_ERRLIM 1e-8
#define SIMP_NMAX 4500     /* max number of simplex iterations */
#define SIMP_LAMB 0.01     /* percent change in starting vertices (1.0+SIMP_LAMB) */
#define SIMP_DEG 1.0       /* degree of rotation for anion starting vertices */
#define SIMP_RESTARTS 11   /* number of times to restart the simplex routine */
#define SIMP_RESTARTS_BASIN_HOP 21

#define BASIN_HOP 0       /* number of basin hops */
#define BASIN_HOP_SCALE_FACTOR_HI 1.25

#define ANG_ALPH_MIN 0.174  /* about 10 degrees */
#define ANG_BETA_MIN 0.174  /* about 10 degrees */
#define ANG_GAMM_MIN 0.174  /* about 10 degrees */
#define ANG_ALPH_MAX 2.200  /* just over 120 degrees */
#define ANG_BETA_MAX 2.200  /* just over 120 degrees */
#define ANG_GAMM_MAX 2.200  /* just over 120 degrees */

#define HOLD_IT_MAX 5
#define N_HALT_CHNG 15      /* stops config changes if they're at 100% for N temp steps */

/**********************************/
/*     Wang-Landau Variables      */
/**********************************/

#define WL_SIMP_RESTARTS_INIT 10
#define WL_SIMP_NMAX_INIT 5000
#define WL_SIMP_RESTARTS 1
#define WL_SIMP_NMAX 10
#define WL_SIMP_RESTARTS_OUTPUT 15
#define WL_SIMP_NMAX_OUTPUT 5000

#define WL_CELL_INIT_DX 0.25 /* locations in 0.25 x 0.25 x 0.25 frac box */
#define WL_SIMP_EDIFF_INIT (5e-3)
#define WL_SIMP_WINDOW 0.95 /* be within 5% to start a simplex */

#define WL_RUNS 200     /* controls the loops per d.o.f. */
#define WL_T_RUNS 20000

#define WL_POSITIVE_STARTING_ENERGY_FACTOR 1.2 /* 20% */
#define WL_EN_MAX 1e10 /* in case user wants to set this expicitly */
#define WL_BIN_WIDTH 0.0 /* can set either nEbins or width in input file */
#define WL_NUM_EBINS 200 /* can set either nEbins or width in input file */
#define WL_MAX_EN_BINS 1500   /* maximum number of energy bins default */

#define WL_HST_RPT_PCT (0.90) /* number of bins filled before calc of flatness */
#define WL_HST_FLAT_METHOD 0 /* 0=x%  1=sigma_window */
#define WL_HST_FLAT_WINDOW (0.5) /* within how much of 1 st dev (method==1)*/
#define WL_HST_FLAT_PCT (0.95) /* flatness requirement default for HST */

#define WL_LNFMOD_INIT 1.0 /* 1==> f0=e */
#define WL_FMIN_CONV (1e-7)
#define WL_LNFMOD_STYLE 0 /* 0=sqrt(lnfmod), 1=lnfmod^(1/n) */
#define WL_QUAL_FACTOR (0.95)
#define WL_EMAX_MULT (0.1) /* max allowed energy multiplier */
#define WL_EMIN_MULT (0.1) /* min allowed energy multiplier */

#define WL_ITER_PRINT 100 /* print info regardless of quality or energy */
#define WL_PRT_POS_TOL 1e-3 /* print poscar if within tol of E_low or lower */

#define WL_REINIT 999999

#define MAX_MOVE_TYPE 5

/**********************************/
/*     Position Constraints       */
/**********************************/
#define CONSTR_XMIN 0.0
#define CONSTR_XMAX 1.0
#define CONSTR_YMIN 0.0
#define CONSTR_YMAX 1.0
#define CONSTR_ZMIN 0.0
#define CONSTR_ZMAX 1.0

#endif // __GLOBAL_DEFS_H__
