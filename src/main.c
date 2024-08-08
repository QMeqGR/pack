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


/**************************************************************/
/*                                                            */
/*  pack.c                                                    */
/*                                                            */
/*  This program calculates the minimum energy and finds the  */
/*  location of tetrahedral or octahedral anions in a cation  */
/*  matrix. The calculation uses a Metropolis algorithm.      */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  file format for input                                     */
/*                                                            */
/*  See header and help region in program                     */
/*                                                            */
/*                                                            */
/**************************************************************/
/*                                                            */
/*  Eric Majzoub                                              */
/*  Sandia National Laboratories                              */
/*  Began programming on                                      */
/*  11 May 2005                                               */
/*                                                            */
/**************************************************************/
/*                                                            */
/*                                                            */
/**************************************************************/

#define PROGRAM_NAME "pack"
#define VERSION_MAJOR 5
#define VERSION_MINOR 2
#define VERSION_BUGFX 3
#define VERSION_TRIVL 1


/*  version notes:

  5.2.3.1    Thu May 10 12:00:20 CDT 2012
  -printing issues with new fix.

  5.2.3.0    Wed May  9 16:35:53 CDT 2012
  -fix WL crashing after finding very low energy bin.
  -small bug fixes for printing.

  5.2.2.0    Thu Apr 19 14:23:17 CDT 2012
  -fixed bug for same-atom in neighboring cell distance.

  5.2.1.1    Wed Apr 18 17:47:30 CDT 2012
  -small bug fix for large Z (LJ tables)

  5.2.1.0    Tue Apr 17 17:43:45 CDT 2012
  -fix constraint bug associated with non_per_ewald setting

  5.2.0.4    Thu Mar  1 10:01:26 CST 2012
  -add input file tag for retaining temp files

  5.2.0.3    Wed Feb  15
  -fix initialization issues with WL routine, minor printing fixes.

  5.2.0.2    Tue Feb  7 12:44:12 CST 2012
  -add wl_simp_nmax_init and wl_simp_restarts_init variables for input file.
   This allows for skipping all init routins in WL.

  5.2.0.1    Wed Dec 21 10:45:08 CST 2011
  -trivial printing issues.

  5.2.0.0    Fri Nov 11 10:20:38 CST 2011
  -add constraints on atom movement by confining their
   fractional coordinates in the cell.

  5.1.2.0    Wed May 11 11:51:53 PDT 2011
  -introduction of Min-map

  5.1.1.2    Wed Jan 26 12:24:59 CST 2011
  -fix swapping hang

  5.1.1.1    Tue Jan 18 14:35:19 CST 2011
  -triv fix for obj choice of fixed objects for rotations
  and allow surface definition in input file for optimized
  Ewald routine for surface calcs

  5.1.0.0    Sat Jan 15 06:34:45 CST 2011
  -add ability to fix atoms

  5.0.1.6    Thu Jan 13 18:34:14 CST 2011
  -change wl_emin_mult and wl_emax_mult behavior

  5.0.1.5    Thu Jan 13 18:34:14 CST 2011
  -printing fixes

  5.0.1.3    Wed Jan 12 08:15:25 CST 2011
  -fix for wl_prt_pos_tol behavior

  5.0.1.2    Mon Jan 10 13:52:51 CST 2011
  -fixed WL init issue for objects with self-overlap

  5.0.1.0    Tue Jan  4 ~11pm
  -fixed cell_init() bug for gen objs.

  5.0.0.2    Tue Jan  4 later
  -testing rotations fix

  5.0.0.0    Tue Jan  4 12:04:07 CST 2011
  -first working version with gen objects

  4.15.0.0   Thu Dec  9 12:04:27 CST 2010
  -generalize the anion definition

See pack_changelog.txt for earlier version changes.

Don't forget to update this file when making changes!!!!

*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "global_defs.h"
#include "global_vars.h"
#include "packlib.h"

#ifdef MPI
#include "mpi.h"
#endif

extern int autoadjust,autoadjust_using;
extern int autoadjusting_hi_temp,autoadjusting_lo_temp;
extern int get_estat;
extern int efunc; /* 0=SS 1=LJ */
extern int forces; /* calculate and use forces */
extern int debug;
extern int anneal_sched;
extern int restart;
extern int runs,t_runs;
extern int mvout;
extern int non_per_ewald;
extern int num_ob,numat;
extern int halt_feature;
extern int ortho;
extern int use_wlmc,wl_nEbins,wl_iter_print;
extern int wl_simp_nmax,lat_parm_chng;
extern int retain_temp_files;
extern int *numZ;
extern int *SWPtab;

extern double errlim,simp_errlim;
extern double rep_epsilon;
extern double hold_acc;
extern double aspect_max;
extern double *fn;
extern double pCnt_global;
extern double wl_fmin_conv;
extern double trans_frac;
extern double *etable,*ecctable;

extern struct cellprm cell;

extern struct cation *at;
extern struct anion *an;
extern struct obj *object;
/* 
   The atom structure will hold all the atoms regardless of whether they are
   in cations or anions.  This structure will be used for the electrostatic
   energy (Ewald) sums.
*/
extern struct atom *p;
extern struct atom *s;


/* Simplex algorithm variables */
extern int n_simp;
extern int simp_nmax;
extern double **sim_pp,*sim_p,*sim_y;


/* Basin hopping variables */
extern int basin_hop,dsm;

int process_num=1, process_rank=0;

/* added for faster deltaE Ewald sum */
int use_deltaEion=0;

/**********************************/
/*                                */
/*          MAIN                  */
/*                                */
/**********************************/
int main(int argc, char *argv[])
{
  char *inputfile="data.dat";
  char tempc;

  int i,ii,j;
  int optnum;
  int ewald_flag=1;
  int do_simplex=1;
  int withQ;
  int dof=0,ncc=0; /* degree of freedom, num config changes */
  int total_basin_hops=0;
  int basin_hop_print=0;
  int temp_flag=0;
  int swpflag=0,num_swappable=0;

  double basin_hop_T0=0;
  double volum_frac=VOLUM_FRAC_MAX;
  double temp_init=TEMP_INIT;
  double temp_finl=TEMP_FINL;

  double obj_charge=0,total_charge=0;


  /* Autoadjust variables */
  double autoadjust_temp_init=0;
  double autoadjust_temp_finl=0;
  
  int trn_chng=1;
  int lat_chng=1;
  int rot_chng=1;
  int swp_chng=1;

  double cell_min,cell_max;

  /* Basin hopping variables */

  double basin_hop_scale_alpha=0.0;
  struct Energy Ebasin={0};

  /* Wang-Landau variables */
  double simp_diff=1e10, simp_old=0, simp_new=0;
  double ewl_eng=0,ewl_eng_old=0;
  struct Energy Ewl={0};

#ifdef MPI
  int ierr;
  char stdoutMPI[25];
#endif
  int seed;

  /**********************************************/
  /*                                            */
  /*              Begin the program             */
  /*                                            */
  /**********************************************/


#ifdef MPI
/*
  Initialize MPI.
*/
  ierr = MPI_Init ( &argc, &argv );

  if ( ierr != 0 )
  {
    printf ( "\n" );
    printf ( "BUFFON_LAPLACE: Warning!\n" );
    printf ( "  MPI_INIT returns IERR = %d\n", ierr );
    ierr = MPI_Finalize ( );
    exit ( 1 );
  }
/*
  Get the number of processes.
*/
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &process_num );
/*
  Get the rank of this process.
*/
  ierr = MPI_Comm_rank ( MPI_COMM_WORLD, &process_rank );
/*
  The master process prints a message.
*/
  if ( is_master() ) 
  {
    printf ( "* MPI mode (WL only) with %d processes\n", process_num );
  }

/* redirecting stdout to different files
*/
  if (process_num>1) {
    sprintf(stdoutMPI,"pegs.out-%d",process_rank);
    freopen( stdoutMPI, "w", stdout );
  }
#endif

  /* initialize the seed for random number generation */
  seed=1234 + process_rank * 100;
  srand48(seed);

  /* be ready to catch a SIGINT signal.  this will cause the program
   * to exit and print the xbs file.
   */
  (void) signal( SIGHUP, exit_and_print );
  (void) signal( SIGINT, exit_immediately );
  (void) signal( SIGUSR1, print_status );
  (void) signal( SIGUSR2, increase_debug );


  if ( argc == 1 ) {
      printf("\nTry '%s -h' for help\n\n",PROGRAM_NAME);
      exit(0);
  }
  while ( (optnum = getopt(argc, argv,
			   "a:A:b:B:cC:d:D:e:Ef:Fg:G:hI:J:K:LMp:q:Qr:Rs:St:TUvVWXY:Z:8")) != -1 ){
    switch(optnum) {
    case 'h':
      printf("#################################################\n"
	     "\t%s\t(version %d.%d.%d.%d)\n"
	     "\tEric Majzoub, %s\n"
	     "\tSandia National Laboratories\n"
             "#################################################\n\n",
	     PROGRAM_NAME,
	     VERSION_MAJOR,VERSION_MINOR,VERSION_BUGFX,VERSION_TRIVL
	     ,__DATE__);
      printf("Source file: " __FILE__ "\n");
      printf("Compile date: " __DATE__ "\n");
      printf("Compile time: " __TIME__ "\n");
      printf("\nCommand line options (* indicates required)\n"
	     "                     (x indicates not implemented yet)\n"
	     "  -h --- help\n"
	     "  -v --- print version number and exit\n"
	     "  -F --- list file formats\n"

	     "\n         Input options:\n\n"

	     "  -f -*- input filename\n"
	     "  -g --- variable SEED, random number generator (default %d)\n"
	     "  -C --- variable TRANS_FRAC_MAX (default %.3f)\n"
	     "           normally constant valued, if TRANS_FRAC_MAX < 0\n"
	     "           it will change with T from |TRANS_FRAC_MAX| to %.3e\n"
	     "  -D --- variable VOLUM_FRAC_MAX (default %.3e)\n"
	     "  -I --- variable TEMP_INIT (default %.2e)\n"
	     "  -J --- variable TEMP_FINL (default %.2e)\n"
	     "  -K --- variable REP_EPSILON (default %.2e)\n"
	     "  -e --- variable ERRLIM (default %.3e)\n"
	     "  -s --- variable SIMP_ERRLIM (default %.2e)\n"
	     "  -p --- variable SIMP_LAMB (%% change in start vertices, default %.3e)\n"
	     "  -q --- variable SIMP_NMAX (max simplex iterations, default %d)\n"
	     "  -r --- variable RUNS (default = %d)\n"
	     "           Num config changes = RUNS * D.O.F\n"
	     "           DOF = 3*numat\n"
	     "  -t --- variable TEMP_RUNS (temp steps, default %d)\n"
	     "  -a --- variable anneal_sched (default %d)\n"
	     "           0: linear\n"
	     "           1: T_0 exp(-run_num**1 / alpha)\n"
	     "           n: T_0 exp(-run_num**n / alpha)\n"
	     "  -b --- variable HOLD_ACC (default %.3e)\n"
	     "           holds temp if (rej %% of n) > HOLD_ACC\n"
	     "           0=trn      1=rot      2=lat     3=swp\n"
	     "           Holds occurs HOLD_IT_MAX times.\n"
	     "           HOLD_IT_MAX is %d, then HOLD_ACC += 0.1\n"
	     "           ex: -B 2.7, holds on lat changes at 70%%\n"
	     "  -G --- varialbe EXT_PRESS (default %.3e)\n"
	     "           turns on pV term in Energy functional.\n"
	     "  -Y --- variable ASPECT_MAX (default %.3e)\n"
	     "\n"
	     "  -T --- don't do translations\n"
	     "  -R --- don't do rotations\n"
	     "  -L --- don't change the lattice parameters\n"
	     "  -S --- don't swap cations and anions\n"
	     "  -c --- confine to orthorhombic symmetry\n"
	     "  -Q --- turn off 'halting' feature. Keeps doing\n"
	     "         config changes even though rej %% is 100.\n"
	     "\n"
	     "  -U --- use LJ potential (default SS with n=12)\n"
	     "  -B --- variable basin_hop, number of hops (int, default 0)\n"
	     "         int > 0 will use distance scaling and no cell reset\n"
	     "         int < 0 will NOT use distance scaling and WILL reset cell\n"
	     "  -V --- don't AUTOADJUST the temperature\n"
	     "  -A --- temperature autoadjust options:\n"
	     "         i)  if A == 1 (use default) or A>0 rejection percentages\n"
	     "         example A= 600650950999 (600 650 950 999)\n"
	     "         HT_LO=0.600  HT_HI=0.650  LT_LO=0.950  LT_HI=0.999\n"
	     "         ii) if A == -1 (use default) or A<0 pCnt acceptance ranges\n"
	     "         example A= -100500202103 (100 500 202 103) (a.b x 10^c)\n"
	     "         HT_LO= 1.0e0  HT_HI= 5.0e0  LT_LO= 2.0e2  LT_HI= 1.0e3\n"
	     "  -Z --- autoadjust HI,LO on (1=trn 2=rot 3=lat 4=swp) (option -A,i)\n"
	     "         44= swp,swp    43=swp,lat, etc... default( %d%d )\n"
	     "  -E --- don't do electrostatic Metropolis algorithm\n"
	     "  -W --- don't do simplex relaxation\n"
	     "  -X --- do a restart calculation (with file 'restart.dat')\n"
	     "\n"
	     "\n         Output options:\n\n"

	     "  -8 --- 80 column output\n"
	     "  -M --- make xbs mv movie file (in.mv)\n"
	     "  -d --- debug level (integer) (default 0..9=everything)\n"
	     "\n\n",SEED,TRANS_FRAC_MAX,TRANS_FRAC_MIN,VOLUM_FRAC_MAX,TEMP_INIT,TEMP_FINL,
	     REP_EPSILON,ERRLIM,SIMP_ERRLIM,SIMP_LAMB,SIMP_NMAX,RUNS,TEMP_RUNS,
	     anneal_sched,HOLD_ACC,HOLD_IT_MAX,EXT_PRESS,ASPECT_MAX,
	     AUTOADJUST_CHNG_HI,AUTOADJUST_CHNG_LO);
      exit(0);
      
    case ':':
      printf("option needs a value\n");
      exit(0);
      
    case 'v':
      printf("#################################################\n"
	     "\t%s\t(version %d.%d.%d.%d)\n"
	     "\tEric Majzoub, %s\n"
	     "\tSandia National Laboratories\n"
             "#################################################\n\n",
	     PROGRAM_NAME,
	     VERSION_MAJOR,VERSION_MINOR,VERSION_BUGFX,VERSION_TRIVL
	     ,__DATE__);
      exit(0);
      break;

    case 'f':
      inputfile  = optarg;
      break;

    case 'F':
      printf("# Input file format                          distance units [ang]\n\n"
	     "# The file format is: tag = value   ! comments after exclamation\n\n"
	     "cella = 10            !   lattice vector magnitues in angstrom\n"
	     "cellb = 10            !   lattice vector magnitues in angstrom\n"
	     "cellc = 10            !   lattice vector magnitues in angstrom\n"
	     "alph = 1.570       !   in RADIANS!!!  90 deg = 1.5707963267949\n"
	     "beta = 1.570       !   in RADIANS!!!  90 deg = 1.5707963267949\n"
	     "gamm = 1.570       !   in RADIANS!!!  90 deg = 1.5707963267949\n"
	     "basin_hop_scale_factor_hi = 1.25            !  default is 1.25\n"
	     "trans_frac = 0.5           ! max frac of cell for translations\n"
	     "sqsh_vol_pct = 0.5    ! min orth val (prevents cell squashing)\n"
	     "dplane = 1.2     ! sets cell_min (autoset by code, suggest you\n"
	     "                 ! leave this alone.                          \n"
	     "restart = 1         ! (default 0) 1 = do a restart calculation\n"
	     "                                                              \n"
	     "                                                              \n"
	     "# Using the min-map routine: (WL only)                        \n"
	     "#              move probabitilities of trn rot lat swp min-map\n"
	     "move_prob = 0 0 0 0 0 ! MUST SET and REPLACING move   switches\n"
	     "MM_transfrac = min max      ! MUST SET for min-map translation\n"
	     "                                                              \n"
	     "                                                              \n"
	     "use_wlmc = 1                       ! use Wang-Landau algorithm\n"
	     "non_per_ewald = 0        ! set to 1 for ALL non-periodic calcs\n"
	     "efunc = 0                    ! 0=soft sphere   1=Lennard-Jones\n"
	     "rep_epsilon = 0.01     ! LJ and SS NOTE: use 1.0 for true LJ!!\n"
	     "forces    = 0                    ! uses forces to move objects\n"
	     "runs = 1000                                        ! dsmc~2000\n"
	     "t_runs = 15                                          ! dsmc~15\n"
	     "wl_runs = 200                                     ! wl~200-400\n"
	     "wl_t_runs = 5000                               ! wl~5000-50000\n"
	     "                                                              \n"
	     "simp_errlim = 1e-5    ! default 1e-8, relax for faster simplex\n"
	     "simp_lamb = 0.01  ! larger sets simplex vertices further apart\n"
	     "simp_nmax = 5000   ! max no. of simplex moves for each restart\n"
	     "simp_restarts = 5   ! no. of restarts each time Simplex called\n"
	     "                                                              \n"
	     "wl_nEbins = 200                        ! number of energy bins\n"
	     "wl_bin_width = 0.1     ! optional to nEbins:  energy bin width\n"
	     "wl_fmin_conv = 2.5e-6      ! conv criterion for W-L f variable\n"
	     "wl_lnfmod_init = 1            ! inital value of W-L f variable\n"
	     "wl_lnfmod_style = 0                          ! 0=sqrt 1=^(1/n)\n"
	     "wl_en_max = -12.5 ! Set max en explicitly (overrides max mult)\n"
	     "wl_emin_mult = 0.2 ! 20pct    Emin = E-fabs(E)*mult in wl_init\n"
	     "wl_emax_mult = 0.1 ! 10pct    Emax = E+fabs(E)*mult in wl_init\n"
	     "#                   set wl_emax_mult to -0.9 or lower for PEGS\n"
	     "#                   Hamiltonian or wall_violations will result\n"
	     "wl_hst_rpt_pct = 0.99    ! pct of HST bins that must be filled\n"
	     "wl_hst_flat_method = 1        ! 0=x%% (WL paper) 1=sigma window\n"
	     "wl_hst_flat_window = 0.80     !  bins not < 80%% ave (method=0)\n"
	     "wl_hst_flat_window = 0.6     ! within +/- 0.6*sigma (method=1)\n"
	     "wl_hst_flat_pct = 0.95      ! pct of HST bins fullfilling reqs\n"
	     "wl_qual_fact = 0.985                    ! rejection pct factor\n"
	     "wl_iter_print = 200        ! lower num check struct more often\n"
	     "wl_weighed_dos = 1         ! weight DOS to sample lower energy\n"
	     "                                                              \n"
	     "#   set WL simp nmax and restarts to zero to skip simp in init\n"
	     "                                                              \n"
	     "wl_simp_nmax_init = 5000               ! simp moves in WL init\n"
	     "wl_simp_restarts_init = 10          ! simp restarts in WL init\n"
	     "wl_simp_nmax = 1                    ! simplex max moves in W-L\n"
	     "wl_cell_init_dx = 0.25  ! init positions in 0.25x0.25x0.25 box\n"
	     "                        !  try larger value for quicker starup\n"
	     "wl_short_init = 0               !       1=only removes overlap\n"
	     "wl_simp_ediff_init = 0.005  ! efunc==0     simplex dE stopping\n"
	     "                 ! try 75 for efunc==1     simplex dE stopping\n"
	     "wl_prt_dos = 1                  ! print the DOS to output file\n"
	     "wl_prt_hst = 1                  ! print the HST to output file\n"
	     "wl_prt_all_poscars = 1              ! generates LOTS of output\n"
	     "wl_prt_pos_tol = 1e-3 ! prints if lower than E_low by this tol\n"
	     "               <0     ! prints if within tol of E_low or lower\n"
	     "                                                              \n"
	     "def_eps =  Z1 Z2 eps_Z1_Z2   ! Lennard-Jones epsilon for Z1 Z2\n"
	     "def_eps =  Z1 Z1 eps_Z1_Z1   ! Lennard-Jones epsilon for Z1 Z1\n"
	     "def_eps = ...                                as many as needed\n"
	     "                                                              \n"
	     "def_rer =  Z1 Z2 rer_Z1_Z2      ! Lennard-Jones dist for Z1 Z2\n"
	     "def_rer =  Z1 Z1 rer_Z1_Z1      ! Lennard-Jones dist for Z1 Z1\n"
	     "def_rer =  ...                               as many as needed\n"
	     "                                                              \n"
	     "                                                              \n"
	     "# The order matters, you must have gen_obj before each def_gen\n"
	     "# Start object count from 0, it also defines the type         \n"
	     "#    A restart file has to be started with the corresponding  \n"
	     "#    input file so that the types will match for swapping.    \n"
	     "                                                              \n"
	     "gen_obj = 0 7 2   ! obj num/typ, num atoms, num repeat objects\n"
	     "#         ob   v         x      y      z     R   Z  chrg  fxd \n"
	     "def_gen =  0   0    0.0000 0.0000 0.0000   0.3   13   7       \n"
	     "def_gen =  0   1    1.0000 0.0000 0.0000   0.2   1   -1       \n"
	     "def_gen =  0   2    0.0000 1.0000 0.0000   0.3   1   -1       \n"
	     "def_gen =  0   3    0.0000 0.0000 1.0000   0.4   1   -1       \n"
	     "def_gen =  0   4    -1.000 0.0000 0.0000   0.5   1   -1       \n"
	     "def_gen =  0   5    0.0000 -1.000 0.0000   0.6   1   -1       \n"
	     "def_gen =  0   6    0.0000 0.0000 -1.000   0.7   1   -1       \n"
	     "                                                              \n"
	     "gen_obj =  1 5 3  ! obj num/typ, num atoms, num repeat objects\n"
	     "#         ob   v         x      y      z     R   Z   chrg fxd \n"
	     "def_gen =  1   0    0.0000 0.0000 0.0000   0.3   13    6   F  \n"
	     "def_gen =  1   1    1.0000 0.0000 0.0000   0.2   1    -1   F  \n"
	     "def_gen =  1   2    0.0000 1.0000 0.0000   0.3   1    -1   F  \n"
	     "def_gen =  1   3    0.0000 0.0000 1.0000   0.4   1    -1   F  \n"
	     "def_gen =  1   4    -1.000 0.0000 0.0000   0.5   1    -1   F  \n"
	     "                                                              \n"
	     "# Do not swap constrained objects with unconstrained!!        \n"
	     "swp_objs = 0 1 1        !! typ1 typ2 1=yes 0=no               \n"
	     "swp_objs = 0 0 0        !! max 12x12 matrix here !!!          \n"
	     "swp_objs = 1 1 0        !! default is zero !!                 \n"
	     "                                                              \n"
	     "# atomic constraints on allowed positions within the cell,    \n"
	     "# input in fractional coordinates                             \n"
	     "                                                              \n"
	     "            ob  v  xmin xmax ymin ymax zmin zmax              \n"
	     "at_constr =  2  1  0.0  1.0  0.0  1.0  0.0  0.2               \n"
	     "at_constr =  2  2  0.0  0.5  0.0  0.5  0.0  0.2               \n"
	     "                                                              \n"
	     "                                                              \n"
	     "#       u.x u.y u.z     v.x v.y v.z !  surface vectors        \n"
	     "surf =     1  0   0       0   1   0                           \n"
	     "                                                              \n"
	     "debug = -2                     ! for WL calcs use -2 or -1    \n"
	     "retain_temp_files = 1          ! retains tempfile and restfile\n"
	     "                                                              \n"
	     "#c2vtol = 0.1  ! WARN def= 1e-8, change only to relax CONTCARs\n"
	     "\n\n\n\n"
	     "# No spaces before variable names allowed.\n\n");
      exit(0);
      break;

    case 'd':
      debug  = atoi(optarg);
      break;

    case 'g':
      srand48( atoi(optarg) );
      break;

    case 'c':
      ortho = 1;
      sqsh_vol_pct = 0.999;
      break;

    case 'e':
      errlim  = atof(optarg);
      break;

    case 'p':
      simp_lamb  = atof(optarg);
      break;

    case 's':
      simp_errlim  = atof(optarg);
      break;

    case 'r':
      runs  = atoi(optarg);
      break;

    case 'q':
      simp_nmax  = atoi(optarg);
      break;

    case 'Q':
      halt_feature = 0;
      break;

    case 'M':
      mvout = 1;
      break;

    case 'U':
      efunc  = 1;
      break;

    case 'T':
      trn_chng = 0;
      break;

    case 'R':
      rot_chng = 0;
      break;

    case 'L':
      lat_chng = 0;
      lat_parm_chng = 0; /* global variable to turn off in simplex */
      break;

    case 'S':
      swp_chng = 0;
      break;

    case 'E':
      ewald_flag = 0;
      break;

    case 'V':
      autoadjust = 0;
      break;

    case 'B':
      basin_hop= atoi(optarg);
      if ( basin_hop > 0 ) dsm=1;
      debug = -1;
      break;

    case 'b':
      hold_acc = atof(optarg);
      break;

    case 'G':
      ext_press = atof(optarg);
      break;

    case 'Y':
      aspect_max = atof(optarg);
      break;

    case 'a':
      anneal_sched = atoi(optarg);
      break;

    case 'Z':
      autoadjust_chng_hi = (int)(atof(optarg)/10.0);
      autoadjust_chng_lo = (int)(atoi(optarg)%10);
      break;
      
    case 'A':
      if ( ( strlen(optarg) != 12 && atoi(optarg) != 1  ) &&
	   ( strlen(optarg) != 13 && atoi(optarg) != -1  ) ){
	printf("String length for option -A is %d characters."
	       "Needs to be 12 or 13(-). Exiting.\n",(int)strlen(optarg));
	printf("A = %d\n",atoi(optarg));
	exit(0);
      }

      if ( atoi(optarg) == 1 ){
	autoadjust_using=1;
	autoadjust_pct_ht_hi = AUTOADJUST_PCT_HT_HI;
	autoadjust_pct_ht_lo = AUTOADJUST_PCT_HT_LO;
	autoadjust_pct_lt_hi = AUTOADJUST_PCT_LT_HI;
	autoadjust_pct_lt_lo = AUTOADJUST_PCT_LT_LO;
      }
      if ( atoi(optarg) > 0 && atoi(optarg) != 1 ){
	autoadjust_using=1;
	for (i=0;i<12;i++){
	  tempc = optarg[i];

	  if ( i==0 ) { autoadjust_pct_ht_lo  = atoi(&tempc) * 1e-1;  }
	  if ( i==1 ) { autoadjust_pct_ht_lo += atoi(&tempc) * 1e-2;  }
	  if ( i==2 ) { autoadjust_pct_ht_lo += atoi(&tempc) * 1e-3;  }
	  
	  if ( i==3 ) { autoadjust_pct_ht_hi  = atoi(&tempc) * 1e-1;  }
	  if ( i==4 ) { autoadjust_pct_ht_hi += atoi(&tempc) * 1e-2;  }
	  if ( i==5 ) { autoadjust_pct_ht_hi += atoi(&tempc) * 1e-3;  }
	  
	  if ( i==6 ) { autoadjust_pct_lt_lo  = atoi(&tempc) * 1e-1;  }
	  if ( i==7 ) { autoadjust_pct_lt_lo += atoi(&tempc) * 1e-2;  }
	  if ( i==8 ) { autoadjust_pct_lt_lo += atoi(&tempc) * 1e-3;  }
	  
	  if ( i==9 )  { autoadjust_pct_lt_hi  = atoi(&tempc) * 1e-1;  }
	  if ( i==10 ) { autoadjust_pct_lt_hi += atoi(&tempc) * 1e-2;  }
	  if ( i==11 ) { autoadjust_pct_lt_hi += atoi(&tempc) * 1e-3;  }
	}
      }

      if ( atoi(optarg) == -1 ){
	autoadjust_using=2;
	autoadjust_pcnt_ht_hi = AUTOADJUST_PCNT_HT_HI;
	autoadjust_pcnt_ht_lo = AUTOADJUST_PCNT_HT_LO;
	autoadjust_pcnt_lt_hi = AUTOADJUST_PCNT_LT_HI;
	autoadjust_pcnt_lt_lo = AUTOADJUST_PCNT_LT_LO;
      }
      if ( atoi(optarg) < 0 && atoi(optarg) != -1 ) {
	autoadjust_using=2;
	for (i=1;i<13;i++){ /* note these are shifted by 1 to make room for the - sign */
	  tempc = optarg[i];
	  
	  if ( i==1 )  { autoadjust_pcnt_ht_lo  = atoi(&tempc) * 1e0;  }
	  if ( i==2 )  { autoadjust_pcnt_ht_lo += atoi(&tempc) * 1e-1;  }
	  if ( i==3 )  { autoadjust_pcnt_ht_lo *= pow( 10, atoi(&tempc) ) ; }
	  
	  if ( i==4 )  { autoadjust_pcnt_ht_hi  = atoi(&tempc) * 1e0;  }
	  if ( i==5 )  { autoadjust_pcnt_ht_hi += atoi(&tempc) * 1e-1;  }
	  if ( i==6 )  { autoadjust_pcnt_ht_hi *= pow( 10, atoi(&tempc) );  }
	  
	  if ( i==7 )  { autoadjust_pcnt_lt_lo  = atoi(&tempc) * 1e0;  }
	  if ( i==8 )  { autoadjust_pcnt_lt_lo += atoi(&tempc) * 1e-1;  }
	  if ( i==9 )  { autoadjust_pcnt_lt_lo *= pow( 10, atoi(&tempc) );  }
	  
	  if ( i==10 ) { autoadjust_pcnt_lt_hi  = atoi(&tempc) * 1e0;  }
	  if ( i==11 ) { autoadjust_pcnt_lt_hi += atoi(&tempc) * 1e-1;  }
	  if ( i==12 ) { autoadjust_pcnt_lt_hi *= pow( 10, atoi(&tempc) );  }
	}
      }

      break;

    case 'C':
      trans_frac  = atof(optarg);
      break;

    case 'D':
      volum_frac  = atof(optarg);
      break;

    case 'I':
      temp_init  = atof(optarg);
      temp_flag=1;
      break;

    case 'J':
      temp_finl  = atof(optarg);
      temp_flag=1;
      break;

    case 'K':
      rep_epsilon  = atof(optarg);
      break;

    case 't':
      t_runs  = atoi(optarg);
      break;

    case 'X':
      restart = 1;
      break;

    case 'W':
      do_simplex = 0;
      break;

    case '8':
      ecol = 1;
      break;

    case '?':
      printf("unknown option: %c\n",optopt);
      exit(0);
      break;
    }
  }

  /*******************************************************************/



  /**********************/
  /*  Output Header     */
  /**********************/
  printf("*************************************************\n"
	 "*\t%s\t(version %d.%d.%d.%d)\n"
	 "*\tEric Majzoub\n"
	 "*\tSandia National Laboratories\n"
	 "*************************************************\n",
	 PROGRAM_NAME,
	 VERSION_MAJOR,VERSION_MINOR,VERSION_BUGFX,VERSION_TRIVL);
  printf("*  Source file: " __FILE__ "\n");
  printf("*  Compile date: " __DATE__ "\n");
  printf("*  Compile time: " __TIME__ "\n");
  printf("*  commandline = ");
  for(i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n");
  printf("*  input file = %s\n",inputfile);
  if (debug) printf("*  Debug level is %d\n",debug);
  printf("**********************************************************************\n");

  printf("*\n*   -  Reading init parameters from command line  -\n*\n");


  /* exit if the number of temperature steps is less than 2 */
  if ( t_runs < 2 ) {
    printf("Must have at least two temperature steps!  Correct the -t flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if anneal_sched is out of bounds */
  if ( anneal_sched < 0 ) {
    printf("Annealing schedule is out of bounds!  Correct the -A flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* warn if one of the temp flags was set and don't autoadjust */
  if ( temp_init != TEMP_INIT || temp_finl != TEMP_FINL || temp_flag==1 ){
    printf("* WARN: -I  AND  -J switch set.  Turning off autoadjust!\n");
    autoadjust = 0;
  }

  /* exit if the initial temperature is given as negative number */
  if ( temp_init < 0 ) {
    printf("Temperature cannot be negative here!  Correct the -I flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if the final temperature is greater than the initial temp */
  if ( temp_finl > temp_init ) {
    printf("Final temperature greater than initial temperature!  Correct the -I or -J flag(s).\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if the final temperature is less than zero */
  if ( temp_finl < 0.0 ) {
    printf("Final temperature less than zero!  Correct the -J flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if the initial temperature is less than zero */
  if ( temp_init < 0.0 ) {
    printf("Initial temperature less than zero!  Correct the -I flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if the trans_frac is out of bounds */
  if ( trans_frac < -0.999 || trans_frac > 0.999 ) {
    printf("TRANS_FRAC_MAX must be > -0.999 and  < 0.999,  Correct the -C flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }
  /* warn if trans_frac is less than zero */
  if ( trans_frac < 0.0 ){
    printf("*  ********* WARN: TRANS_FRAC_MAX < 0.0, will scale with temp.\n"
	   "*  ********* Autoadjusting temp will not give target rejection %%!!!\n");
  }

  /* exit if the volum_frac is out of bounds */
  if ( volum_frac < 0.0 || volum_frac > 1.0 ) {
    printf("VOLUM_FRAC_MAX must be > 0.0  and  < 1.0.  Correct the -D flag.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if wl_cell_init_dx is larger than 1.0 */
  if ( wl_cell_init_dx < 0.0 || wl_cell_init_dx > 1.0 ) {
    printf("wl_cell_init_dx must be > 0.0  and  < 1.0.  Correct the parameter value.\n");
    printf("Exiting program.\n");
    exit_pack(0);
  }

  /* exit if the aspect_max is out of bounds */
  if ( aspect_max <= 1.0 ) {
    printf("ASPECT_MAX must be > 1.0.  Correct the -Y flag.\n");
    printf("Exiting program.\n");
    exit(0);
  }

  /* exit if autoadjust values don't make sense */
  if ( autoadjust_using == 1 ){
    if ( autoadjust_pct_lt_hi <= autoadjust_pct_lt_lo || 
	 autoadjust_pct_ht_hi <= autoadjust_pct_ht_lo  ) {
      printf("* ERROR: autoadjust pct value: hi < low.\n");
      printf("%10.4f%10.4f%10.4f%10.4f\n",
	     autoadjust_pct_ht_lo,
	     autoadjust_pct_ht_hi,
	     autoadjust_pct_lt_lo,
	     autoadjust_pct_lt_hi);
      printf("* Exiting\n");
      exit(0);
    }
  }
  if ( autoadjust_using == 2 ){
    if ( autoadjust_pcnt_lt_hi <= autoadjust_pcnt_lt_lo || 
	 autoadjust_pcnt_ht_hi <= autoadjust_pcnt_ht_lo  ) {
      printf("* ERROR: autoadjust pcnt value: hi < low.\n");
      printf("%12.3e%12.3e%12.3e%12.3e\n",
	     autoadjust_pcnt_ht_lo,
	     autoadjust_pcnt_ht_hi,
	     autoadjust_pcnt_lt_lo,
	     autoadjust_pcnt_lt_hi);
      printf("* Exiting\n");
      exit(0);
    }
  }

  /* warn if flags are set such that no MC run can be performed  */
  if ( ewald_flag==0 && autoadjust==0 ) {
    printf("* WARN: As flags are set, no MC runs will be done!  Check -E and -V flags.\n");
  }

  printf("**********************************************************************\n");

  /**************************************************/
  /*                                                */
  /*               READ INPUT FILE                  */
  /*                                                */
  /**************************************************/

  print_input_out(inputfile);
  printf("**********************************************************************\n");
  printf("*\n*   -  Reading init parameters from input file %s -\n*\n",inputfile);
  /* Open input file and get initialization conditions */
  input = fopen(inputfile,"r");
  if ( input == NULL )
    {
      printf("Error opening input file\n");
      exit_pack(0);
    }
  get_init_parms(input,1);
  fclose(input);


  /* warn if c2vtol is set  */
  if ( c2vtol > C2VTOL ) {
    printf("* WARN: C2VTOL is %f : strict checking has been relaxed!\n",c2vtol);
    printf("*       Are you running a restart calculation with non-ideal anions?\n");
    printf("*       Only ideal, analytical vr-vr overlaps are removed in E_ss.\n");
    printf("*       You may get anion-anion overlap contributions in E_ss.\n");
  }
  /* warn if best_compare_off is set  */
  if ( best_compare_off==1 ) printf("* WARN: best_compare() is turned off!\n");

  /* turn off lat changes if non_per_ewald is set */
  if ( non_per_ewald && lat_chng==1 ){
    printf("* WARN: Found non_per_ewald=1 in input file and lat_chng=1!\n");
    printf("*       Turning lattice parameter changes off.\n");
    lat_chng=0;
  }

  /* warn if one use_wlmc flag is set and turn off autoadjust */
  if ( use_wlmc==1 && non_per_ewald==1 && autoadjust==1 ) {
    printf("* WARN: Wang-Landu and non_per_ewald flags set.  Turning off autoadjust!\n");
    autoadjust = 0;
  }

  /* Warn the user if rep_epsilon was redefined in the input file: the LJ energies will not be correct */
  if ( !tol(rep_epsilon , REP_EPSILON, 1e-8) ){
    printf("* WARN: rep_epsilon = %10.5e ... default is 1.0---> L-J energies will not be correct!!!\n",rep_epsilon);
  }


  /*******************************/
  /*                             */
  /*        dplane setting       */
  /*                             */
  /*******************************/
  /* set this for later comparisons when the cell size changes */
  /* if dplane = -1 then it was not set in the input file */
  if ( dplane < 0 ) {
    cell_min = cvd_max;
    cell_max = 100 * cvd_max;
    dplane = DPLANE * cvd_max;
  } else {
    /* dplane was set in input file */
    cell_min = dplane;
    cell_max = 100*dplane;
  }

  /* If doing cluster calc, set trans_frac automatically */
  if ( tol(trans_frac, TRANS_FRAC_MAX, 1e-8) && non_per_ewald==1 ){
    if ( 0 == is_cubic() ) {
      printf("* init: FATAL ERROR:  Cell does not look cubic\n");
      printf("* Check settings in input file.\n");
      exit_pack(0);
    }
    trans_frac = dplane * TRANS_FRAC_NPE_FACTOR / cell.a;
  }

  /* if trans frac min is greater than trans frac max, then lower trans frac min */
  /* The variable trans_frac_max is only defined in Metropolis and Wang-Landau source
     files.  Main passes trans_frac to Metro and WL. */
  if ( trans_frac_min > trans_frac ){
    printf("* WARN: found trans_frac_min > trans_frac.  Setting trans_frac_min to 0.1*max.\n");
    trans_frac_min = 0.10 * trans_frac;
  }


  printf("*\n");
  printf("**********************************************************************\n");

  /*******************************/
  /*                             */
  /*    allocate memory          */
  /*                             */
  /*******************************/

  /* total number of atoms */
  for (i=0; i<MAX_OBJS && object[i].used==1 ;i++){
    num_ob++;
    numat  += object[i].num_vert;
  }
  if ( debug>1 ){
    printf("* Total number of objects is: %d\n",num_ob);
    printf("* Total number of atoms is: %d\n",numat);
    fflush(stdout);
  }

  /* exit if the number of objects is zero */
  if ( num_ob == 0 ){
    printf("Found num_ob == 0. Check input file.\n");
    exit(0);
  }

  /* allocate space for the SS hard distance matrix */
  SShrd = (double *)malloc( numat * numat * sizeof(double) );
  if ( !SShrd ) {
    printf("%s\n","Not enough memory for SShrd.\n");
    exit_pack(0);
  }
  /* allocate space for the SS hard distance matrix */
  SShrd_2 = (double *)malloc( numat * numat * sizeof(double) );
  if ( !SShrd_2 ) {
    printf("%s\n","Not enough memory for SShrd_2.\n");
    exit_pack(0);
  }
  /* allocate space for the SS hard distance matrix */
  SShrd_6 = (double *)malloc( numat * numat * sizeof(double) );
  if ( !SShrd_6 ) {
    printf("%s\n","Not enough memory for SShrd_6.\n");
    exit_pack(0);
  }
  /* allocate space for the SS hard distance matrix */
  SShrd12 = (double *)malloc( numat * numat * sizeof(double) );
  if ( !SShrd12 ) {
    printf("%s\n","Not enough memory for SShrd12.\n");
    exit_pack(0);
  }

  /* allocate space for the SS min distance matrix */
  SSmin = (double *)malloc( numat * numat * sizeof(double) );
  if ( !SSmin ) {
    printf("%s\n","Not enough memory for SSmin.\n");
    exit_pack(0);
  }


  /* allocate space for atoms for ss and lj */
  s = (struct atom *)malloc( numat * sizeof(struct atom) );
  if ( !s ) {
    printf("%s\n","Not enough memory for s atom struct.\n");
    exit_pack(0);
  }
  /* allocate space for all charged atoms for electrostatic sums */
  p = (struct atom *)malloc( numat * sizeof(struct atom) );
  if ( !p ) {
    printf("%s\n","Not enough memory for all atoms.\n");
    exit_pack(0);
  }
  /* allocate space for numZ array for printing */
  numZ = (int *)malloc( ZMAX * sizeof(int) );
  if ( !numZ ) {
    printf("%s\n","Not enough memory for numZ array.\n");
    exit_pack(0);
  }


  /* allocate space for energy table for optimized ss and lj calculations */
  /* n^2 may be a waste of space, but indexing the upper triangular portion
     takes computation (see Knuth!!) , so we will settle for wasted space */
  etable = (double *)malloc( numat * numat * sizeof(double) );
  if ( !etable ) {
    printf("%s\n","Not enough memory for etable.\n");
    exit_pack(0);
  }
  ecctable = (double *)malloc( numat * numat * sizeof(double) );
  if ( !ecctable ) {
    printf("%s\n","Not enough memory for ecctable.\n");
    exit_pack(0);
  }
  /* initialize etable */
  for(i=0;i<numat;i++){
      for(ii=0;ii<numat;ii++){
	etable[ numat * i + ii ]=0;
	ecctable[ numat * i + ii ]=0;
      }
  }

  /* allocate space for simplex array, the number is for 6 cell parameters,
     then (x,y,z) for all cation centers and all anion centers, and 
     3 degrees of rotation freedom for each anion */
  n_simp = 6 + 3*num_ob + 3*num_ob;
  fn = (double *)malloc( n_simp * sizeof(double) );
  if ( !fn ) {
    printf("%s\n","Not enough memory for simplex array fn.\n");
    exit_pack(0);
  }
  /* follow NR and make sim_pp an array of pointers to the rows of sim_p */
  sim_pp = (double **)malloc( ( n_simp+1 ) * sizeof(double *) );
  if ( !sim_pp ) {
    printf("%s\n","Not enough memory for simplex pointer array sim_pp.\n");
    exit_pack(0);
  }
  sim_p = (double *)malloc( ( n_simp+1 ) * n_simp * sizeof(double) );
  if ( !sim_p ) {
    printf("%s\n","Not enough memory for simplex array sim_p.\n");
    exit_pack(0);
  }
  for(i=0; i<(n_simp+1); i++) sim_pp[i] = sim_p + i*n_simp;
  sim_y = (double *)malloc( ( n_simp+1 ) * sizeof(double) );
  if ( !sim_y ) {
    printf("%s\n","Not enough memory for simplex array sim_y.\n");
    exit_pack(0);
  }

  /* space for rotation variables for the simplex routine */
  /* needed to keep track of each anion rotation for each new
     call to simplex_trans */
  rot_var = (double *)malloc( 3* num_ob * sizeof(double) );
  if ( !rot_var ) {
    printf("%s\n","Not enough memory for simplex rotation variables.\n");
    exit_pack(0);
  }

  /* allocate vectors for psum and ptry, for the routines in
     amotry_ehm and amoeba_ehm.  this should save overhead
     for the system calls and speed up the app. */
  psum = (double *)malloc( n_simp * sizeof(double) );
  if ( !psum ) {
    printf("%s\n","Not enough memory for simplex vector psum.\n");
    exit_pack(0);
  }
  ptry = (double *)malloc( n_simp * sizeof(double) );
  if ( !ptry ) {
    printf("%s\n","Not enough memory for simplex vector psum.\n");
    exit_pack(0);
  }

  if (is_master()){
  /* open output mv file if needed */
    if ( mvout ) {
      outfile = fopen("in.mv","w");
      if ( outfile==NULL ){
	printf("Can't open in.mv file.\n");
	exit_pack(0);
      }
    }
  }

  /* get temporary filename for storing the 'best structure so far' */
  processid = getpid();
  if ( debug>0 ) printf("* processid=%d\n",processid);
  snprintf(pid_string, 15, "%d", processid );
  strncat(tempfile, pid_string, 15);
  strncat(tempfile, tmp_string, 9);

  /* get restart filename for storing 'restart' file (used in WL too!) */
  strncat(restfile, pid_string, 15);
  strncat(restfile, rst_string, 9);

  snprintf(pid_string, 15, "%d", process_rank );
  strncat(tempfile, pid_string, 15);
  strncat(restfile, pid_string, 15);

  
  /*******************************************************************/
  
  /* turn off swapping if all objects are fixed */
  for (i=0; i<MAX_OBJS && object[i].used==1; i++){
    if ( object[i].vfxd[0] == 0 ) num_swappable++;
  }
  if ( num_swappable<2 ) {
    printf("* WARN: less than 2 swappable items. Turning off swapping.\n");
    swp_chng=0;
  }
  
  /* calculate the charges */
  for (i=0; i<MAX_OBJS && object[i].used==1; i++){
    for (j=0, obj_charge=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      total_charge += object[i].vchrg[j];
      obj_charge   += object[i].vchrg[j];
    }
    if ( debug>1) printf("* Object %3d   object_charge %10.5f\n",i,obj_charge);
  }

  /* WARN if total charge is not zero. */
  if ( !tol( total_charge , 0.0 , 1e-6 ) ) {
    printf("* WARN: total system charge is non-zero: %.20f\n",total_charge);
  }


  /* 
   * initialize the cell structure
   * and initialize the random locations for the objects
   * in the starting cell, Z values, radii, etc.
   */

  printf("*  Calling cell_init ... \n");
  cell_init();
  printf("*  +++ cell init done +++\n");
  
  
  
  /**************************************************/
  /*                                                */
  /*                RESTART STUFF                   */
  /*                                                */
  /**************************************************/
  /*
   *  The program starts as usual, but gets
   *  interrupted here and just replaces the
   *  values that were initialized above.
   *
   */
  if ( restart ) {
    /* Open restart file and get initialization conditions */
    restartinput = fopen("restart.dat","r");
    if ( restartinput == NULL )
      {
	printf("Error opening restart file\n");
	exit_pack(0);
      }
    get_init_parms(restartinput,0);    
    fclose(restartinput);
  }
  check_constraints();


  /***************************************************************************/
  /***************************************************************************/

  /* check cell parameters  */
  if ( !restart && ( fabs(cell.a) < cell_min || fabs(cell.a) > cell_max ||
		     fabs(cell.b) < cell_min || fabs(cell.b) > cell_max ||
		     fabs(cell.c) < cell_min || fabs(cell.c) > cell_max ) ) {
    printf("* !!!!!!  cell parameters out of range, check the parameter d\n");
    printf("* !!!!!!  cell_min = %f   cell_max = %f [ang]\n",cell_min,cell_max);
    printf("* Exiting program\n");
    exit_pack(0);
  }

  /* check the swap table for swappable items */
  if ( debug > 2 ) {
    printf("* Swap table ( SWPtab )\n*");
    for(i=0;i<SWPTAB;i++) printf("%3d",i);
    printf("\n*\n");
  }
  for(i=0,swpflag=0;i<MAX_TABLE_SWP;i++){
    if ( debug > 2 ) {
      if ( (i+1)%SWPTAB==1 ) printf("*");
      printf("%3d",SWPtab[i]);
      if ( (i+1)%SWPTAB==0 ) printf("\n");
    }
    if ( SWPtab[i]==1 ) swpflag=1;
  }
  if ( swpflag==0 ) swp_chng=0; /* no objects allow swapping, turn off swaps! */

  printf("**********************************************************************\n");
  printf("* \n");
  printf("* %30s\n","- General run comments -");
  printf("* \n");
  if ( ext_press > 1e-12 ) printf("* Using external pressure, and H = E + pV as functional.\n");
  if ( trn_chng == 0 ) printf("* WARN: no Translations (-T flag).\n");
  if ( rot_chng == 0 ) printf("* WARN: no Rotations (-R flag).\n");
  if ( lat_chng == 0 ) printf("* WARN: no Lat Param Changes (-L flag).\n");
  if ( swp_chng == 0 ) printf("* WARN: no Swaps (-S flag).\n");
  if ( ortho == 1 ) printf("* WARN: run restricted to orthorhombic symmetry.\n");
  printf("*\n");
  printf("**********************************************************************\n");

  /* calculate the number of degrees of freedom and the appropriate number of runs */
  dof = 3*numat;
  runs *= dof;
  ncc = 4*runs;

  /*****************************/
  /* print of input parameters */
  /*****************************/
  printf("* %40s\n","- Input parameters -");
  printf("*\n");
  printf("* %30s%20d%30s\n","debug  =",debug,"debug level");
  printf("* %30s%20d%30s\n","restart  =",restart,"1 = restart calc");
  printf("* %30s%20d%30s\n","autoadjust  =",autoadjust,"automatically set temp");
  if ( autoadjust ) printf("* %30s%20d%30s\n","aa_using  =",autoadjust_using,"1=rej pct, 2=pCnt");
  printf("* %30s%20d%30s\n","basin hop  =",basin_hop,"number of hops");
  printf("* %30s%20d%30s\n","non_per_ewald  =",non_per_ewald,"1=nanopart calc");
  printf("* %30s%20d%30s\n","forces  =",forces,"1=use forces");
  printf("* %30s%20d%30s\n","use_wlmc  =",use_wlmc,"1=use Wang-Landau MC");
  if ( use_wlmc == 1 ){
    printf("*\n");
    printf("* %30s%20d%30s\n","wl_runs  =",wl_runs,"loops per T");
    printf("* %30s%20d%30s\n","wl_t_runs  =",wl_t_runs,"num t steps");
    if ( tol( wl_bin_width, WL_BIN_WIDTH, 1e-6 ) ) {
      printf("* %30s%20d%30s\n","wl_nEbins  =",wl_nEbins,"no. of energy bins");
    } else {
      printf("* %30s%20.5e%30s\n","wl_bin_width  =",wl_bin_width,"energy bin width");
    }
    printf("* %30s%20.5e%30s\n","wl_emin_mult  =",wl_emin_mult,"wl emin_mult");
    if ( !tol( wl_en_max, WL_EN_MAX, 1e-2 ) ){
      printf("* %30s%20.5e%30s\n","wl_en_max  =",wl_en_max,"wl_en_max");
    } else {
      printf("* %30s%20.5e%30s\n","wl_emax_mult  =",wl_emax_mult,"wl emax_mult");
    }
    printf("* %30s%20.5e%30s\n","wl_lnfmod_init  =",wl_lnfmod_init,"wl lnfmod init");
    printf("* %30s%20d%30s\n","wl_lnfmod_style  =",wl_lnfmod_style,"0=sqrt 1=^1/n");
    printf("* %30s%20.5e%30s\n","wl_fmin_conv  =",wl_fmin_conv,"wl conv factor");
    printf("* %30s%20d%30s\n","wl_iter_print  =",wl_iter_print,"iter print modulus");
    printf("* %30s%20d%30s\n","wl_weighted_dos  =",wl_weighted_dos,"use weighted DOS");
    printf("* %30s%20d%30s\n","wl_simp_nmax  =",wl_simp_nmax,"max moves in WL");
    printf("* %30s%20d%30s\n","wl_simp_nmax_init  =",wl_simp_nmax_init,"max moves in WL init");
    printf("* %30s%20d%30s\n","wl_simp_restarts_init  =",wl_simp_restarts_init,"max restarts in WL init");
    printf("* %30s%20.5f%30s\n","wl_cell_init_dx  =",wl_cell_init_dx,"WL init cell positions");
    printf("* %30s%20.5f%30s\n","wl_simp_ediff_init  =",wl_simp_ediff_init,"WL init simplex dE stopping");
    printf("* %30s%20.5f%30s\n","wl_hst_rpt_pct  =",wl_hst_rpt_pct,"pct filled");
    printf("* %30s%20d%30s\n","wl_hst_flat_method  =",wl_hst_flat_method,"0=(WL paper) 1=sigma window");
    printf("* %30s%20.5f%30s\n","wl_hst_flat_window  =",wl_hst_flat_window,">w*ave (mth=0) ,sgma (mth=1)");
    printf("* %30s%20.5f%30s\n","wl_hst_flat_pct  =",wl_hst_flat_pct,"pct bins req flat");
    printf("* %30s%20d%30s\n","wl_prt_hst  =",wl_prt_hst,"print histogram");
    printf("* %30s%20d%30s\n","wl_prt_dos  =",wl_prt_dos,"print d.o.s.");
    printf("* %30s%20d%30s\n","wl_prt_all_poscars  =",wl_prt_all_poscars,"print ALL poscars");
    printf("* %30s%20.5e%30s\n","wl_prt_pos_tol  =",wl_prt_pos_tol,"within Elow or lower");
    printf("* %30s%20d%30s\n","wl_short_init  =",wl_short_init,"1=only removes overlap");
  }
  if ( efunc == 1 ){
    printf("*\n");
    printf("* %10s%5s%5s%10s%20s\n"," ","Z1","Z2","  ","eps");
    for(i=0;i<MAX_TABLE;i++){
      for(ii=0;ii<MAX_TABLE;ii++){
	if ( !tol( EPS_DEFAULT, LJeps[i*MAX_TABLE+ii], 1e-6) )
	printf("* %10s%5d%5d%10s%20.5e%30s\n"," ",i,ii,"  =",LJeps[i*MAX_TABLE+ii],"Lennard-Jones eps");
      }
    }
    printf("*\n");
    printf("* %10s%5s%5s%10s%20s\n"," ","Z1","Z2","  ","rer");
    for(i=0;i<MAX_TABLE;i++){
      for(ii=0;ii<MAX_TABLE;ii++){
	if ( !tol( RER_DEFAULT, LJrer[i*MAX_TABLE+ii], 1e-6) )
	printf("* %10s%5d%5d%10s%20.5e%30s\n"," ",i,ii,"  =",LJrer[i*MAX_TABLE+ii]/LJSGMA,"Lennard-Jones rer");
      }
    }
  }
  printf("*\n");
  printf("* %30s%20.5f\n","lat a  =",cell.a);
  printf("* %30s%20.5f\n","lat b  =",cell.b);
  printf("* %30s%20.5f\n","lat c  =",cell.c);
  printf("* %30s%20.5f\n","alpha  =",cell.alph);
  printf("* %30s%20.5f\n","beta   =",cell.beta);
  printf("* %30s%20.5f\n","gamma  =",cell.gamm);
  printf("* %30s%20.5g\n","cell vol  =",cell_volume(&cell));
  printf("* %30s%20.5f\n","sqsh vol pct  =",sqsh_vol_pct);
  printf("* \n");
  printf("* %30s%20.5f%30s\n","cell_min  =",cell_min,"min lat vec");
  printf("* %30s%20.5f%30s\n","cell_max  =",cell_max,"max lat vec");
  printf("* \n");
  printf("* %30s%20d%30s\n","num_ob  =",num_ob,"tot num objects");
  printf("* %30s%20d%30s\n","numat  =",numat,"tot num atoms");
  printf("* %30s%20.5e%30s\n","total charge  =",total_charge,"total system charge");
  printf("* \n");
  printf("* %30s%20.5e%30s\n","C2VTOL  =",c2vtol,"c2v tolerance");
  printf("* %30s%20.5e%30s\n","errlim  =",errlim,"ewald error setting");
  if ( efunc == 0 )
    printf("* %30s%20.5e%30s\n","rep_epsilon  =",rep_epsilon,"soft sphere repul");
  if ( efunc == 1 )
    printf("* %30s%20.5e%30s\n","rep_epsilon  =",rep_epsilon,"epsilon repul");
  if ( trans_frac > 0 )
    printf("* %30s%20.5f%30s\n","trans_frac_max  =",trans_frac,"fixed fract of cell parm");
  if ( trans_frac < 0 )
    printf("* %30s%20.5f%30s\n","trans_frac_max  =",trans_frac,"annealed");
  printf("* %30s%20.5f\n","trans_frac_min  =",trans_frac_min);
  printf("* %30s%20.5f%30s\n","MAX_ANG_CHNG  =",MAX_ANG_CHNG*R2D,"degrees");
  printf("* %30s%20.5f\n","volum_frac_max  =",volum_frac);
  printf("* %30s%20.5e%30s\n","ext_press  =",ext_press,"external pressure");
  printf("* %30s%20.5f%30s\n","sqsh_vol_pct  =",sqsh_vol_pct,"min orth val");
  printf("* %30s%20.5f%30s\n","aspect_max  =",aspect_max,"cell parm");
  printf("* %30s%20.5f%30s\n","dplane  =",dplane,"min plane sep");
  if ( non_per_ewald==1 ){
    printf("* %30s%20.5f%30s\n","wall cushion  =",ACCEPT_NPE_DPLANE_WALL_CUSHION*dplane/cell.a,"NPE wall cushion");
    printf("* %30s%20.5f%30s\n","wall violation  =",VIOLAT_NPE_DPLANE_WALL_CUSHION*dplane/cell.a,"NPE wall violation");
  }
  if ( hold_acc > 0 )
    printf("* %30s%20.5f%30s\n","HOLD_ACC  =",hold_acc,"max rej%%");
  printf("* %30s%20.5f%30s\n","SCALE_FACTOR  =",basin_hop_scale_factor_hi,"bhop max scaling");
  printf("* \n");
  printf("* %30s%20d%30s\n","d.o.f.  =",dof,"degrees of freedom");
  printf("* %30s%20d%30s\n","runs/dof  =",runs/dof,"loops per d.o.f.");
  printf("* %30s%20d%30s\n","runs  =",runs,"loops per T");
  printf("* %30s%20d%30s\n","ncc  =",ncc,"config changes per T");
  printf("* %30s%20d%30s\n","t_runs  =",t_runs,"num temp steps");
  printf("* %30s%20d%30s\n","efunc  =",efunc,"0=SS  1=LJ");
  printf("* %30s%20d%30s\n","N_HALT_CHNG  =",N_HALT_CHNG,"halts after N 100%s");
  printf("*\n");
  printf("* %30s%20d%30s\n","anneal_sched  =",anneal_sched,"0=lin or exp(-run_n^n)");
  printf("* %30s%20.5e%30s\n","temp_init  =",temp_init,"initial temp");
  printf("* %30s%20.5e%30s\n","temp_finl  =",temp_finl,"final temp");
  printf("* \n");
  printf("* %30s%20d%30s\n","trn_chng  =",trn_chng,"0=off  1=on");
  printf("* %30s%20d%30s\n","rot_chng  =",rot_chng,"0=off  1=on");
  printf("* %30s%20d%30s\n","lat_chng  =",lat_chng,"0=off  1=on");
  printf("* %30s%20d%30s\n","swp_chng  =",swp_chng,"0=off  1=on");
  printf("* \n");
  printf("* %30s%20.5f%30s\n","move_prob (trn)  =",move_prob[0],"translation");
  printf("* %30s%20.5f%30s\n","move_prob (rot)  =",move_prob[1],"rotation");
  printf("* %30s%20.5f%30s\n","move_prob (lat)  =",move_prob[2],"lat parm");
  printf("* %30s%20.5f%30s\n","move_prob (swp)  =",move_prob[3],"swaps");
  printf("* %30s%20.5f%30s\n","move_prob (m-m)  =",move_prob[4],"min-map");
  printf("* %30s%20.5f%30s\n","MM_transfrac (min)  =",MM_transfrac[0],"min trans");
  printf("* %30s%20.5f%30s\n","MM_transfrac (max)  =",MM_transfrac[1],"max trans");
  printf("* \n");
  printf("* %30s%20d%30s\n","simp_restarts  =",simp_restarts,"no. simplex restarts");
  printf("* %30s%20d%30s\n","SIMP_RESTARTS BH  =",SIMP_RESTARTS_BASIN_HOP,"basin hop simp restarts");
  printf("* %30s%20d%30s\n","simp_nmax  =",simp_nmax,"max simplex iters");
  printf("* %30s%20.5f%30s\n","SIMP_DEG  =",SIMP_DEG,"anion start deg");
  printf("* %30s%20.5f%30s\n","simp_lamb  =",simp_lamb,"vertex start frac");
  printf("* %30s%20.5e%30s\n","simp_errlim  =",simp_errlim,"error limit");
  printf("* \n");
  printf("* %30s%20d%30s\n","aa runs  =",AUTOADJUST_RUNS,"config changes");
  printf("* %30s%20d%30s\n","aa runs_step  =",AUTOADJUST_RUNS_STEP,"config change increase");
  printf("* %30s%20d%30s\n","aa runs_max  =",AUTOADJUST_RUNS_MAX,"max auto config changes");
  printf("* %30s%20d%30s\n","aa t_runs  =",autoadjust_t_runs,"num temp runs");
  printf("* %30s%20.5f%30s\n","aa temp_step  =",AUTOADJUST_TEMP_STEP,"temp step");
  printf("* \n");
  if ( autoadjust_using == 1 ){
    printf("* %30s%20d%30s\n","aa chng hi  =",autoadjust_chng_hi,"1=trn 2=rot 3=lat 4=swp");
    printf("* %30s%20d%30s\n","aa chng lo  =",autoadjust_chng_lo,"1=trn 2=rot 3=lat 4=swp");
    printf("* \n");
    printf("* %30s%20.5f%30s\n","aa ht_pct_hi  =",autoadjust_pct_ht_hi,"high rej pct");
    printf("* %30s%20.5f%30s\n","aa ht_pct_lo  =",autoadjust_pct_ht_lo,"low  rej pct");
    printf("* \n");
    printf("* %30s%20.5f%30s\n","aa lt_pct_hi  =",autoadjust_pct_lt_hi,"high rej pct");
    printf("* %30s%20.5f%30s\n","aa lt_pct_lo  =",autoadjust_pct_lt_lo,"low  rej pct");
  }
  if ( autoadjust_using == 2 ){
    printf("* %30s%20.5e%30s\n","aa ht_pcnt_hi  =",autoadjust_pcnt_ht_hi,"high pCnt");
    printf("* %30s%20.5e%30s\n","aa ht_pcnt_lo  =",autoadjust_pcnt_ht_lo,"low  pCnt");
    printf("* \n");
    printf("* %30s%20.5e%30s\n","aa lt_pcnt_hi  =",autoadjust_pcnt_lt_hi,"high pCnt");
    printf("* %30s%20.5e%30s\n","aa lt_pcnt_lo  =",autoadjust_pcnt_lt_lo,"low  pCnt");
  }
  printf("* \n");
  printf("* %30s%20d%30s\n","filetag  =",processid,"process ID");
  printf("* %30s%20s\n","temporary file  =",tempfile);
  printf("* %30s%20s\n","restart file  =",restfile);
  printf("* \n");
  

  /******************************************************/
  /*                                                    */
  /*   Initialize the structures for estat and energ    */
  /*                                                    */
  /******************************************************/

  /* initialize the p structure for ewald sums */
  estat_init( );
  /* initialize the s structure for ss and lj  */
  energ_init( );

  /* print out the objects for debugging */
  if ( debug>1 ) {
    printf("* STARTING POSITIONS of input objects (xbs)\n");
    print_xbs( );
  }

  if ( debug>1 ) {
    printf("* STARTING POSITIONS of input objects (fract coords)\n");
    print_all_objs( );

  /* print out the center to vertex distances */
    print_c2v( );
  }

  /* check whether the fractional coords are out of bounds */
  if ( debug > -1 ){
    printf("* Checking for out of bounds fractional coords.\n*\n");
    check_objects( );
  }


  /* initialize the 'moved' flags and calculate the initial energy */
  mark_all_objects_moved();

  printf("**********************************************************************\n");

  /* print info block before starting metro runs */
  printf("*\n*          - Parameters at START -\n");
  print_info_block();
  
  /******************************************************************/
  /*                                                                */
  /*                    Metropolis algorithms                       */
  /*                                                                */
  /******************************************************************/
  /*
    There may be several metropolis runs.  The startup runs are
    for the autoadjust process to find the proper annealing temperature
    defined by the #define variables in the header.

    11 aug 2005, fixes to version 2.2.5 need to be made:
    During the Metro runs, every time the energy functional is called,
    the energy will be compared to a 'best-so-far'.  If the new structure
    has a lower energy, it will become the 'best' and it's structure
    will be saved to a temporary file.  At the end of the Metro
    runs, if this structure is still the lowest energy, it will
    replace the structure at the end of the metro runs.

  */


  /********************************************************************/
  /*                                                                  */
  /*              Temperature AUTOADJUST Metro runs                   */
  /*                                                                  */
  /********************************************************************/

  if ( autoadjust ) {
    
    printf("*\n**********************************************************************\n*\n");

    if ( temp_init == TEMP_INIT ) {
      printf("* ------------------------------------------------------------\n");
      printf("*              Autoadjust high temperatrure:\n");
      printf("* ------------------------------------------------------------\n");
      autoadjusting_hi_temp=1;
      autoadjusting_lo_temp=0;

      if ( autoadjust_using == 1 ) {
	printf("*              Autoadjusting on %d (1=trn,2=rot,3=lat,4=swp)\n",autoadjust_chng=autoadjust_chng_hi);
	autoadjust_temp_init = autoadjust_temp_rejpct(AUTOADJUST_TEMP_HT_INIT,
						      autoadjust_pct_ht_lo,autoadjust_pct_ht_hi,
						      cell_min,cell_max,trans_frac,volum_frac,
						      trn_chng,rot_chng,lat_chng,swp_chng);
      }
      
      if ( autoadjust_using == 2 )
	autoadjust_temp_init = autoadjust_temp_pCnt(AUTOADJUST_TEMP_HT_INIT,
						    autoadjust_pcnt_ht_lo,autoadjust_pcnt_ht_hi,
						    cell_min,cell_max,trans_frac,volum_frac,
						    trn_chng,rot_chng,lat_chng,swp_chng);
      
      temp_init = autoadjust_temp_init;
    }


    if ( temp_finl == TEMP_FINL ) {
      printf("* ------------------------------------------------------------\n");
      printf("*              Autoadjust low temperatrure:\n");
      printf("* ------------------------------------------------------------\n");
      autoadjusting_hi_temp=0;
      autoadjusting_lo_temp=1;

      if ( autoadjust_using == 1 ) {
	printf("*              Autoadjusting on %d (1=trn,2=rot,3=lat,4=swp)\n",autoadjust_chng=autoadjust_chng_lo);
	autoadjust_temp_finl = autoadjust_temp_rejpct(AUTOADJUST_TEMP_LT_INIT,
						      autoadjust_pct_lt_lo,autoadjust_pct_lt_hi,
						      cell_min,cell_max,trans_frac,volum_frac,
						      trn_chng,rot_chng,lat_chng,swp_chng);
      }

      if ( autoadjust_using == 2 )
	autoadjust_temp_finl = autoadjust_temp_pCnt(AUTOADJUST_TEMP_LT_INIT,
						    autoadjust_pcnt_lt_lo,autoadjust_pcnt_lt_hi,
						    cell_min,cell_max,trans_frac,volum_frac,
						    trn_chng,rot_chng,lat_chng,swp_chng);

      temp_finl = autoadjust_temp_finl;
      autoadjust_temp_finl /= AUTOADJUST_LO_T_SAFETY;
    }
    
    if ( temp_finl > temp_init ) {
      printf("*     WARN: temp_finl > temp_init.\n");
      while ( temp_finl > temp_init ) temp_finl /= 10;
      printf("*     Resetting temp_finl = %e\n*\n",temp_finl);
    }
    
    /*   re-initialize cell contents to remove autoadjust artifacts    */
    if ( debug > -1 ){
      printf("*\n*          -- Initializing cell contents --\n*\n");
    }
    cell_init();

  } /* end autoadjust */  
  


  /********************************************************************/
  /*                                                                  */
  /*                         MAIN ALGORITHM                           */
  /*                                                                  */
  /********************************************************************/

  if ( ewald_flag ) {

    /*********************************************/
    /*                                           */
    /*  Distance scaling and basin hopping code  */
    /*                                           */
    /*********************************************/
    if ( basin_hop > 0 ){

      /* if greater than zero, then use distance scaling, and do NOT reset the cell contents
       * at every hop
       */

      /* set the temperature parameters */
      total_basin_hops=basin_hop;
      basin_hop_T0 = temp_init;
      basin_hop_temp_alpha  = pow( (double)(total_basin_hops-1), anneal_sched ) / log ( basin_hop_T0 / temp_finl );
      basin_hop_scale_alpha = pow( (double)(total_basin_hops-1), anneal_sched ) / log ( basin_hop_scale_factor_hi / 1.0 );

      printf("*\n**********************************************************************\n*\n");
      if ( basin_hop > 0 ) printf("*          - Initializing DSM Routine -\n*\n");
      if ( basin_hop < 0 ) printf("*          - Initializing Basin Hopping Routine -\n*\n");
      printf("* %20s%10.2e%25s\n","temp_init  =",temp_init,"from aa or cmd line");
      printf("* %20s%10.2e%25s\n","temp_finl  =",temp_finl,"from aa or cmd line");
      printf("* %20s%10.2e%25s\n","temp_alph  =",basin_hop_temp_alpha,"for bhop t_init");
      printf("* %20s%10.2e%25s\n","scale_alph  =",basin_hop_scale_alpha,"for bhop scale");

      /* set the initial scaling factor */
      scale_cell_do = 1;  /* turn on the cell scaling */
      scale_factor = basin_hop_scale_factor_hi;

    }

    if ( basin_hop < 0 ){

      /*  if basin_hop was set less than zero, don't use scaling, and also reset the cell contents
       *  at every hop.
       */

      /* set the initial scaling factor */
      scale_cell_do = 0;  /* turn off the cell scaling */
      scale_factor = 1.0;

    }

    /* print output if basin hopping is used */
    if ( basin_hop > 0 || basin_hop < 0 ){

      printf("*\n**********************************************************************\n*\n");
      if ( basin_hop > 0 ) printf("*           -- DSM Routine --\n*\n");
      if ( basin_hop < 0 ) printf("*           -- Basin Hopping Routine --\n*\n");
      if ( ecol == 0 )
      printf("%-3s%6s%15s%15s%10s%10s%10s%8s%8s%8s%10s%10s%10s%10s\n",
	     "*H","step","ecc","ess","a","b/a","c/a","alph","beta","gamm","orth","pCnt","scale","t_init");

    }

    /************************************************************************/
    /*                  Standard Metropolis or Wang-Landau                  */
    /************************************************************************/


    do {

      
      if ( debug >= -2 && use_wlmc==1 ){
	printf("*\n");
	printf("*\n**********************************************************************\n*\n");
	printf("*            -- Wang-Landau calculation --\n*\n");
      }
      else if ( debug > -1 && use_wlmc==0 ) {
	printf("*            -- Starting Metropolis calculation --\n*\n");
	printf("* %25s%10.3e\n","temp_init  =",temp_init);
      }
      fflush(stdout);
      
      get_estat=1;
      if ( use_wlmc==1 ) {
	
//	if ( non_per_ewald == 1 ){
        if ( (non_per_ewald == 1) && (!restart) ){ /* not a bulk structure */

	  printf("*            -- WL Nanoparticle Run,  Initial Simplex To Find Low Energy Bound --\n");
	  if ( wl_short_init ) printf("*          >>>> running short init, only remove overlap <<<<\n");
	  printf("*                   dE stopping criteria = %.3e\n",
		 wl_simp_ediff_init);
	  Ewl = Total_Energy();
	  printf("*\n");
	  printf("%1s%15s%15s%15s%15s\n","*","ecc","eng","tot","dE (tot)");
	  printf("* --------------------------------------------------------------\n");
	  fflush(stdout);
	  while ( simp_diff>wl_simp_ediff_init ) {

	    Simplex( wl_simp_restarts_init , wl_simp_nmax_init );

	    ewl_eng_old = ewl_eng;
	    Ewl = Total_Energy();
	    ewl_eng = Ewl.eng;

	    simp_old = simp_new;
	    simp_new = Ewl.tot;
	    simp_diff = fabs(simp_old-simp_new);

	    if ( wl_short_init ){
	      if ( efunc==0 && tol(Ewl.eng,0.0,1e-4) ) simp_diff=1e-10;
	      if ( efunc==1 && (simp_diff<wl_simp_ediff_init) ) simp_diff=1e-10;
	    }

	    if ( non_per_ewald==1 ) {
	      center_np();
	      if ( 0 != check_wall_violation() ){
		printf("* WARN: Wall violation! return code= %d  1=ca   10=an\n",check_wall_violation());
	      }
	    }
	    printf("*%15.6e%15.6e%15.6e%15.6e",Ewl.ecc,Ewl.eng,simp_new,simp_diff); 	fflush(stdout);
	    if ( (ewl_eng - ewl_eng_old)>errlim ) {
	      printf(" ---> overlap : re-initialize cell\n");
	      cell_init(); simp_diff=1e10;
	    } else printf("\n");
	  }
	  printf("*     Simplex convergence reached.  Starting structure:\n");
	} else if ( !restart ) { /* if we are doing a bulk structure */
	  
	  printf("*         -- WL Bulk Init Metro To Find Low Energy Bound --\n");
	  do {
	    
	    Metro(runs,t_runs,trn_chng,rot_chng,lat_chng,swp_chng,
		  cell_min,cell_max,trans_frac,volum_frac,temp_init,temp_finl);
	    Simplex( wl_simp_restarts_init , wl_simp_nmax_init );

	    ewl_eng_old = ewl_eng;
	    Ewl = Total_Energy();
	    ewl_eng = Ewl.eng;
	    /* the second condition is if overlap is on the anions themselves */
	    if ( debug > -2 ) {
	      printf("* --------------------------------------------------\n");
	      printf("*\n* WL Bulk Init Metro and Simplex: ediff routine:\n");
	      printf("* Stops on wl_simp_ediff_init= %f\n",wl_simp_ediff_init);
	      printf("* Check this vale in input if stuck in loop.\n");
	      printf("* ewl_eng= %10.4e    ewl_eng_old= %10.4e   fabsdiff= %10.4e\n*\n",
				     ewl_eng,ewl_eng_old,fabs(ewl_eng - ewl_eng_old));
	      printf("* --------------------------------------------------\n");
	    }

	  } while ( efunc==0 && (fabs(ewl_eng - ewl_eng_old)>wl_simp_ediff_init) );

	}
	
	printf("*PP Begin POSCAR N 0 ecc=0 tot=0\n"); /* dummy header for script */
	print_poscar_out(0);
	printf("*rr restart output\n");
	printf("*rr (use: | awk '($1==\"*r\"){print $0}' | sed 's/*r //' > temp.out )\n");
	print_restart_fileptr(stdout,'r');

	printf("*            -- Starting Wang-Landau calculation --\n*\n");
	Wang_Landau(wl_runs,wl_t_runs,trn_chng,rot_chng,lat_chng,swp_chng,
		    cell_min,cell_max,trans_frac,volum_frac);

      } else { /* just run straight Metro */

	Metro(runs,t_runs,trn_chng,rot_chng,lat_chng,swp_chng,
	      cell_min,cell_max,trans_frac,volum_frac,temp_init,temp_finl);
	
      }
      
      if ( debug > -1 ){
	printf("*\n**********************************************************************\n*\n");
	printf("*\n*          - Parameters at Post Ewald-MC -\n");
	print_info_block();
      }
      
      if ( basin_hop > 0 ) {
	
	basin_hop_print=1;
	scale_cell_do=0; /* turn off scaling for actual energy value */
	if ( forces == 1 ) {
	  np_force_mv();
	}
	if ( do_simplex ) Simplex( SIMP_RESTARTS_BASIN_HOP , simp_nmax );
	
	
	Ebasin = Total_Energy();
	scale_cell_do=1; /* turn scaling back on for next basin hop */
	
	scale_factor = basin_hop_scale_factor_hi * exp( -pow( total_basin_hops-basin_hop , anneal_sched ) / basin_hop_scale_alpha );
	temp_init    = basin_hop_T0 * exp( -pow( total_basin_hops-basin_hop , anneal_sched ) / basin_hop_temp_alpha );
	basin_hop--;
	
      }
      if ( basin_hop < 0 ) {
	
	basin_hop_print=1;
	scale_cell_do=0; /* turn off */
	if ( do_simplex ) Simplex( SIMP_RESTARTS_BASIN_HOP , simp_nmax  );
	Ebasin = Total_Energy();

	scale_factor = 1.0;
	basin_hop++;
	if ( basin_hop != 0 ) cell_init();

      }

      /* the print statement for basin hop output */
      if ( basin_hop_print == 1 ){

	if ( ecol == 0 )
	  /* number of colums is 143 here */  
	  printf("*B %6d%15.5e%15.5e%10.3f%10.4f%10.4f%8.1f%8.1f%8.1f%10.4f%10.2e%10.4f%10.2e\n",
		 basin_hop,
		 Ebasin.ecc,Ebasin.eng,
		 Ebasin.a,Ebasin.b/Ebasin.a,Ebasin.c/Ebasin.a,
		 R2D*Ebasin.alph,R2D*Ebasin.beta,R2D*Ebasin.gamm,
		 Ebasin.orth,pCnt_global,
		 scale_factor,temp_init);
	if ( ecol == 1 ){
	  printf("*%c------------------------------------------------------------\n",'B');
	  printf("*B %20s%12d%25s\n","basin  =",basin_hop,"hop number");
	  printf("*B %20s%12.3f%25s\n","Ebasin  =",Ebasin.ecc,"ecc actual");
	  printf("*B %20s%12.3f%25s\n","Ebasin  =",Ebasin.eng,"eng actual");
	  printf("*B %20s%12.3f%25s\n","a  =",Ebasin.a,"latt a");
	  printf("*B %20s%12.3f%25s\n","b  =",Ebasin.b,"latt b");
	  printf("*B %20s%12.3f%25s\n","c  =",Ebasin.c,"latt c");
	  printf("*B %20s%12.3f%25s\n","alph  =",R2D*Ebasin.alph,"alpha");
	  printf("*B %20s%12.3f%25s\n","beta  =",R2D*Ebasin.beta,"beta");
	  printf("*B %20s%12.3f%25s\n","gamm  =",R2D*Ebasin.gamm,"gamma");
	  printf("*B %20s%12.3f%25s\n","pf  =",Ebasin.pf,"packing frac");
	  printf("*B %20s%12.3f%25s\n","orth  =",Ebasin.orth,"max 1.0");
	  printf("*B %20s%12.3e%25s\n","pCnt  =",pCnt_global,"pseudo heat cap");
	  printf("*B %20s%12.3f%25s\n","scale  =",scale_factor,"min 1.0");
	  printf("*B %20s%12.3e%25s\n","t_init  =",temp_init,"new t_init");
	  timing_info( );
	}

      }

      
    } while ( basin_hop > 0 || basin_hop < 0 );
    /* BASIN HOPPING CODE */
    
    
    scale_cell_do=0;  /* turn off scaling just to be sure */
    
  } /* if ewald_flag */
  
  
  /********************************************************************/
  /*                                                                  */
  /*         Comparison of 'BEST' and post metro structures           */
  /*                                                                  */
  /********************************************************************/

  /* only do this if DSM basin hopping was NOT used */
  if ( dsm==0 && best_compare_off==0 ) {
    printf("*\n**********************************************************************\n*\n");
    printf("*         - Comparing with 'best' structure -\n");
    fflush(stdout);
    best_compare();
  }

  /********************************************************************/
  /*                                                                  */
  /*                    Simplex Relaxation Routine                    */
  /*                                                                  */
  /********************************************************************/
  if ( do_simplex ) {
    printf("*\n**********************************************************************\n*\n");
    printf("*         - Final Simplex relaxation with restarts -\n");
    printf("*\n**********************************************************************\n*\n");
    printf("*\n");
    simp_restarts=SIMP_RESTARTS;
    Simplex( simp_restarts, simp_nmax );
    print_info_block();
  }

  /******************************/
  /*                            */
  /*       PRINTING Out         */
  /*                            */
  /******************************/

  /* Print symsearch, findsym, and diffraction input files  */
  print_symsearch_out();
  print_findsym_out();
  print_powder_out();
  printf("*PP Begin POSCAR N 1000000 ecc= %.20f tot= %.20f\n",Total_Energy().ecc,Total_Energy().tot); /* dummy header for script */
  print_poscar_out(withQ=0);
  /* print_poscar_out(withQ=1); */
  print_restart_inline();
  
  /* XBS output */
  print_xbs();

  /* print info block last time for processing scripts to grab info */
  printf("*-------------------------------------------------------------\n*\n");
  printf("*             - Structure parameters at end of run -\n");
  print_info_block();
  printf("*-------------------------------------------------------------\n");

  /* print timing information */
  timing_info( );
  
  /* close open files */
  if ( mvout ) fclose( outfile );

  if ( tempfile != NULL && retain_temp_files == 0 ) remove_tempfile();
  if ( restfile != NULL && retain_temp_files == 0 ) remove_restfile();

  exit_pack(0);
#ifdef MPI
  MPI_Finalize();
#endif
  return(0); 
}
/*****************************************/
/*                                       */
/*           END of Program              */
/*                                       */
/*****************************************/









/*****************************************/
/*                                       */
/*           TIMING INFO                 */
/*                                       */
/*****************************************/
void timing_info(void)
{
  struct rusage r_usage;

  printf("***********************************************************\n");
  printf("*        - Process and timing information -\n");
  printf("*\n");
  if ( getrusage(RUSAGE_SELF,&r_usage) == -1 )
    printf("* Error getting timer information.\n");
  
  printf("* user time: %3ld h %02ld m %02ld.%06ld s\n",
	 r_usage.ru_utime.tv_sec/3600,
	 (r_usage.ru_utime.tv_sec%3600)/60,
	 r_usage.ru_utime.tv_sec%60,
	 r_usage.ru_utime.tv_usec);
  printf("* syst time: %3ld h %02ld m %02ld.%06ld s\n",
	 r_usage.ru_stime.tv_sec/3600,
	 (r_usage.ru_stime.tv_sec%3600)/60,
	 r_usage.ru_stime.tv_sec%60,
	 r_usage.ru_stime.tv_usec);
  printf("* signals: %ld\n",r_usage.ru_nsignals);
  printf("* voluntary context switches: %ld\n",r_usage.ru_nvcsw);
  printf("* involuntary context switches: %ld\n",r_usage.ru_nivcsw);
  
  printf("*\n");
  printf("***********************************************************\n");
  
  return;
}


/*****************************************/
/*                                       */
/*        Signal Catching Code           */
/*                                       */
/*****************************************/
void increase_debug(int sig)
{
  extern int debug;

  printf("*########################################################################\n");
  printf("* Caught signal %d !  Increasing debug value by 1!\n",sig);
  printf("*########################################################################\n");
  debug = debug+1;
  return;
}

void exit_and_print(int sig)
{
  static int call_number=0;
  extern int done_flag,wl_done_flag;
  extern int basin_hop;
  extern int simp_nmax,simp_runs;
  extern double hold_acc;

  printf("*########################################################################\n");
  printf("* Caught SIGNAL    %10d\n",sig);
  printf("* call number   =  %10d\n",call_number);
  printf("* done_flag     =  %10d   (for Metro routine)\n",done_flag);
  printf("* wl_done_flag  =  %10d   (for W-L MC routine)\n",wl_done_flag);
  printf("* basin_hop     =  %10d   (basin hoppings left)\n",basin_hop);
  printf("* simp_runs     =  %10d   (number of simp runs)\n",simp_runs);
  printf("* hold_acc      =  %10.4f\n",hold_acc);
  printf("*\n");
  
  
  if ( call_number==0 && (hold_acc<0.0) && simp_runs==0 ) {
    printf("* Caught signal %d, call_number %d. Will exit program at next loop!\n",
	   sig,call_number);
    done_flag=1;
    wl_done_flag=1;
    basin_hop=0;
  }
  
  if ( call_number==0 && (hold_acc>0.0) && simp_runs==0 ) {
    printf("* Caught signal %d, call_number %d. Will set hold_acc to -1.\n",
	   sig,call_number);
    hold_acc=-1.0;
    basin_hop=0;
  }
  
  if ( call_number==1 && (hold_acc<0.0) && simp_runs==0 ) {
    printf("* Caught signal %d, call_number %d. Will exit Metro at next loop!\n",
	   sig,call_number);
    done_flag=1;
    wl_done_flag=1;
    basin_hop=0;
  }
  
  if ( simp_runs > 0 ) {
    simp_nmax=0;
    done_flag=1;
    basin_hop=0;
    printf("* simp_runs > 0 ---> Will exit at next opportunity!\n");
  }
  
  printf("*########################################################################\n");
  
  call_number++;
  return;
}

void exit_immediately(int sig)
{
  double vol=0;
  extern struct cellprm cell;
  vol = cell_volume(&cell);

  printf("*########################################################################\n");
  printf("* Balls! Caught signal %d !\n",sig);
  printf("* Exit Immediately Error\n");
  printf("*########################################################################\n");
  print_restart_inline();
  printf("*PP Begin POSCAR N 999999 ecc= %.20f tot= %.20f\n",Total_Energy().ecc,Total_Energy().tot); /* dummy header for script */
  print_poscar_out(0);
  signal_print_status();
  timing_info( );
  printf("* SIGHUP will exit gracefully (if this was an control-c exit).\n");
  printf("* Cleaning up temporary file.\n");
  printf("*########################################################################\n");

  exit_pack(0);
}

void print_status(int sig)
{
  double vol=0;
  extern struct cellprm cell;
  vol = cell_volume(&cell);

  printf("*########################################################################\n");
  printf("* Caught signal %d !\n",sig);
  printf("*########################################################################\n");
  signal_print_status();
  printf("*PP Begin POSCAR N 999999 ecc= %.20f tot= %.20f\n",Total_Energy().ecc,Total_Energy().tot); /* dummy header for script */
  print_poscar_out(0);
  printf("*########################################################################\n");
  return;
}

void signal_print_status(void)
{
  double vol=0;
  extern double dplane;
  extern struct cellprm cell;

  vol = cell_volume(&cell);

  printf("*########################################################################\n");
  printf("*\n* Lattice parameters:\n*\n");
  printf("* %20s%20.4f : %10.4f%10.4f%10.4f\n","lat a  =",cell.a,cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
  printf("* %20s%20.4f : %10.4f%10.4f%10.4f\n","lat b  =",cell.b,cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
  printf("* %20s%20.4f : %10.4f%10.4f%10.4f\n","lat c  =",cell.c,cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
  printf("* %20s%20.4f\n","alpha  =",cell.alph * R2D);
  printf("* %20s%20.4f\n","beta   =",cell.beta * R2D);
  printf("* %20s%20.4f\n","gamma  =",cell.gamm * R2D);
  printf("* %20s%20.4f\n","cell vol  =",vol);
  printf("* %20s%20.4f\n","ortho  =",vol/(cell.a*cell.b*cell.c));
  printf("* %20s%20.4f\n","sa  =",surface_area());
  printf("* %20s%20.4f\n","pf  =",pack_frac());
  printf("* %20s%20.4f\n","dplane  =",dplane);
  printf("* %20s%20.4f\n","dab  =",vol / ( cell.a * cell.b ) );
  printf("* %20s%20.4f\n","dac  =",vol / ( cell.a * cell.c ) );
  printf("* %20s%20.4f\n","dbc  =",vol / ( cell.b * cell.c ) );
  printf("*\n");
  printf("*########################################################################\n");
  fflush(stdout);

  return;
}


void remove_tempfile(void)
{

  if ( -1 == unlink(tempfile) ){
    printf("* Error removing temporary file: %s\n",tempfile);
  }

  return;
}

void remove_restfile(void)
{

  if ( -1 == unlink(restfile) ){
    printf("* Error removing restart file: %s ... file may not exist.\n",restfile);
  }

  return;
}

void exit_pack(int n)
{
  
  free_memory();

  if ( tempfile != NULL && retain_temp_files == 0 ) remove_tempfile();
  if ( restfile != NULL && retain_temp_files == 0 ) remove_restfile();

  exit(n);

  return;
}

void free_memory(void)
{
  
  if ( at )      free( at );
  if ( fn )      free( fn );
  if ( sim_p )   free( sim_p );
  if ( sim_y )   free( sim_y );
  if ( rot_var ) free( rot_var );
  if ( psum )    free( psum );
  if ( ptry )    free( ptry );

  /*  Problems freeing these ??

    if ( s )       free( s );
    if ( p )       free( p );
    if ( sim_pp ) free( sim_pp );

  */
  
  return;
}
