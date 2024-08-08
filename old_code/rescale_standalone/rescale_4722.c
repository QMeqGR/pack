/* the process is:

1. read in the 'mixed' new restart file
2. read the 'old' cell parms
3. do the rescale
4. print out the new restart file

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "packlib.h"

#define PROGRAM_NAME "pack"
#define VERSION_MAJOR 4
#define VERSION_MINOR 6
#define VERSION_BUGFX 0
#define VERSION_TRIVL 3

/* note bugfix version numbers include small enhancements such as changes in
 * starting parameters below, or changes in output printing styles, etc...
 */

#define INIT_DATA_1 20
#define INIT_DATA_2 24

#define TINY (1e-12)
#define SMALL (1e-5)
#define Bohr2Ang 0.529177249
#define PI (3.14159265358979323846264338327950288419716939937510)
#define R2D ( 180.0 / PI )
#define D2R ( PI / 180.0 )
#define SPHERE_VOL 4.18879 /* 4*pi/3 */
#define ERRLIM 1e-4
#define REP_EPSILON 1.0  /* repulsion epsilon */

#define PRESS_DVOL 1e-4  /* dvol for pressure calculation */
#define BEST_VOL 0.1
#define VOLUM_FRAC_MAX 0.07  /* max cell volume change per lat step */
#define ASPECT_MAX 20.0
#define ORTHO 0             /* orthorhomic restriction */ 
#define MAX_ANG_CHNG 0.003   /* max angle change in RADIANS for lat vector angles */
#define ANG_ALPH_MIN 0.174  /* about 10 degrees */
#define ANG_BETA_MIN 0.174  /* about 10 degrees */
#define ANG_GAMM_MIN 0.174  /* about 10 degrees */
#define ANG_ALPH_MAX 2.094  /* about 120 degrees */
#define ANG_BETA_MAX 2.094  /* about 120 degrees */
#define ANG_GAMM_MAX 2.094  /* about 120 degrees */
#define SQSH_VOL_PCT 0.85

/*******************************/
/* External (global) Variables */
/*******************************/

int cac=0,debug=0;
int an_coord=0;
int num_an=0;
int n_type_ca=0;
int ca_chrg_1=0;
int ca_chrg_2=0;
int num_ca=0;
int num_ca_1=0;
int num_ca_2=0;
int Zcat_1=0,Zcat_2=0,Zcent=0,Zvert=0;

double Rcat_1=0,Rcat_2=0,Rcent=0,Rvert=0;
double an_chrg_c=0;
double an_chrg_v=0;
double d=0;

/* Structures  (see packlib.h for other definitions) */
struct hrd_dist{
  /*  double a_a_min,a_c_min,a_v_min,c_c_min,c_v_min,v_v_min; */
  double an_an_min,an_c1_min,an_c2_min,an_vr_min;
  double c1_c1_min,c1_c2_min,c1_vr_min;
  double c2_c2_min,c2_vr_min;
  double vr_vr_min;
};

struct cellprm cell; /* lattice parameters */
struct atom *p=NULL;
struct atom *s=NULL;
struct tetrahedron *tetr=NULL;
struct octahedron *octa=NULL;
struct dimer *dimr=NULL;
struct cation *at=NULL;
struct Energy E_best;

extern char *optarg;
extern int optind,opterr,optopt;

FILE *mixed_restart=NULL;
FILE *old_header=NULL;
FILE *new_restart=NULL;
FILE *restartfile=NULL;


void print_restart(void);
void print_restart_fileptr(FILE *,char);
void get_init_parms(FILE *,
		    int exearly,
		    int *Zcat_1,
		    int *Zcat_2,
		    int *Zcent,
		    int *Zvert,
		    double *Rcat_1,
		    double *Rcat_2,
		    double *Rcent,
		    double *Rvert);


/**********************************/
/*                                */
/*          MAIN                  */
/*                                */
/**********************************/
int main(int argc, char *argv[])
{
  double a_old=0,b_old=0,c_old=0,a_new=0,b_new=0,c_new=0;
  double alph_old=0,beta_old=0,gamm_old=0,alph_new=0,beta_new=0,gamm_new=0;
  double da=0,db=0,dc=0;
  double dalph=0,dbeta=0,dgamm=0;
  double dvol=0;
  
  struct vector cell_chng,ang_chng;
  extern struct cellprm cell;

  if ( debug>0 ) printf("Starting the rescale program!\n");

  old_header = fopen("old_header.dat","r");               /* header with old latt parms */
  if ( old_header == NULL ){
    printf("File error on old header.  Exiting.\n");
    exit(0);
  }

  mixed_restart = fopen("mixed_restart.dat","r");               /* header with old new parms */
  if ( mixed_restart == NULL ){
    printf("File error on mixed restart.  Exiting.\n");
    exit(0);
  }

  if ( debug>0 ) printf("Getting init_parms from mixed restart file for memory allocation.\n");
  get_init_parms(mixed_restart,1,&Zcat_1,&Zcat_2,&Zcent,&Zvert,&Rcat_1,&Rcat_2,&Rcent,&Rvert);


  if ( debug>0 ) printf("Allocating Memory.\n");
  /*******************************/
  /*                             */
  /*    allocate memory          */
  /*                             */
  /*******************************/

  /* allocate space for cation structures */
  at = (struct cation *)malloc( num_ca * sizeof(struct cation) );
  if ( !at ) {
    printf("%s\n","Not enough memory for atoms.\n");
    exit(0);
  }
  
  if ( an_coord==0 ) {
    /* allocate space for tetrahedra structures */
    tetr = (struct tetrahedron *)malloc( num_an * sizeof(struct tetrahedron) );
    if ( !tetr )  {
      printf("%s\n","Not enough memory for tetrahedra.\n");
      exit(0);
    }
  }
  
  if ( an_coord==1 ) {
    /* allocate space for octahedra structures */
    octa = (struct octahedron *)malloc( num_an * sizeof(struct octahedron) );
    if ( !octa ) {
      printf("%s\n","Not enough memory for octahedra.\n");
      exit(0);
    }
  }

  if ( an_coord==2 ) {
    /* allocate space for dimer structures */
    dimr = (struct dimer *)malloc( num_an * sizeof(struct dimer) );
    if ( !dimr ) {
      printf("%s\n","Not enough memory for dimers.\n");
      exit(0);
    }
  }

  if ( debug>0 ) printf("Getting atom positions from mixed restart file.\n");
  get_init_parms(mixed_restart,0,&Zcat_1,&Zcat_2,&Zcent,&Zvert,&Rcat_1,&Rcat_2,&Rcent,&Rvert);

  if ( debug>0 ) printf("Getting cell parms from old header file.\n");
  /* read in the old lat parms */
  fscanf(old_header," %lf %lf %lf ",&a_old,&b_old,&c_old);
  fscanf(old_header," %lf %lf %lf ",&alph_old,&beta_old,&gamm_old);

  if ( debug>0 ){
    printf("Old Cell parms\n");
    printf("%20.10f%20.10f%20.10f\n",a_old,b_old,c_old);
    printf("%20.10f%20.10f%20.10f\n",alph_old,beta_old,gamm_old);
  }
  cell.a = a_old;
  cell.b = b_old;
  cell.c = c_old;
  cell.alph = alph_old;
  cell.beta = beta_old;
  cell.gamm = gamm_old;

  if ( debug>0 ) printf("Getting cell parms from mixed restart (new header) file.\n");
  /* read in the new lat parms */
  fscanf(mixed_restart," %lf %lf %lf ",&a_new,&b_new,&c_new);
  fscanf(mixed_restart," %lf %lf %lf ",&alph_new,&beta_new,&gamm_new);
  if ( debug>0 ){
    printf("New Cell parms\n");
    printf("%20.10f%20.10f%20.10f\n",a_new,b_new,c_new);
    printf("%20.10f%20.10f%20.10f\n",alph_new,beta_new,gamm_new);
  }

  fclose(old_header);
  fclose(mixed_restart);

  da = a_new - a_old;
  db = b_new - b_old;
  dc = c_new - c_old;

  dalph = alph_new - alph_old;
  dbeta = beta_new - beta_old;
  dgamm = gamm_new - gamm_old;

  cell_chng = makevec(da,db,dc);
  ang_chng = makevec(dalph,dbeta,dgamm);
  
  if ( debug ) {
    printf("* CONFIG CHANGE: LAT\n");
    printf("*  a  b  c = %12.6f%12.6f%12.6f\n",cell.a,cell.b,cell.c);
    printf("* da db dc = %12.6f%12.6f%12.6f\n",da,db,dc);
    printf("*  angles  = %12.6f%12.6f%12.6f\n",cell.alph,cell.beta,cell.gamm);
    printf("* dangles  = %12.6f%12.6f%12.6f\n",dalph,dbeta,dgamm);
    printf("* dvol = %12.6f    cell_vol = %12.6f     ratio = %12.6f\n",
	   dvol,cell_volume(&cell),dvol/cell_volume(&cell));
  }
  
  if( debug ) printf("Rescaling the thing!\n");
  /* rescale the object locations */
  /* the rescaling needs the pre-changed cell parameters and the change params */
  rescale_anions( an_coord, tetr, octa, dimr, num_an, &cell, cell_chng, ang_chng);

  if( debug ) printf("Printing restart file\n");
  print_restart();

  return (0);
}

/* copied and modified from pack.c */
void get_init_parms(
		    FILE *getfile,
		    int exearly,
		    int *Zcat_1,
		    int *Zcat_2,
		    int *Zcent,
		    int *Zvert,
		    double *Rcat_1,
		    double *Rcat_2,
		    double *Rcent,
		    double *Rvert
		    )
{
  int errval,num_vert=0;
  int i,j;

  extern struct cellprm cell;

  if ( an_coord==0 ) num_vert=5;
  if ( an_coord==1 ) num_vert=7;
  if ( an_coord==2 ) num_vert=2;

  /* get cell parameter */
  errval = fscanf(getfile," %lf %lf %lf %lf %lf %lf ",
		  &cell.a,&cell.b,&cell.c,&cell.alph,&cell.beta,&cell.gamm);
  if ( errval != 6 ) { printf("Error reading input file: loc. 1 errval=%d.\n",errval); exit(0); }

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


  /* get anion information */
  errval = fscanf(getfile," %d %d %lf %lf %lf %lf %lf %d %d ",
		  &an_coord,&num_an,
		  &an_chrg_c,&an_chrg_v,
		  &d,Rcent,Rvert,Zcent,Zvert);
  if ( errval != 9 ) { printf("Error reading input file: loc. 2 errval=%d.\n",errval); exit(0); }

  /* get number of cations */
  errval = fscanf(getfile," %d ",&n_type_ca);
  if ( errval != 1 ) { printf("Error reading input file: loc. 3 errval=%d.\n",errval); exit(0); }

  /* get cation information */
  if ( n_type_ca == 1 ) {
    errval = fscanf(getfile," %d %d %lf %d ",&num_ca_1,&ca_chrg_1,Rcat_1,Zcat_1);
    if ( errval != 4 ) { printf("Error reading input file: loc. 4 errval=%d.\n",errval); exit(0); }
    num_ca = num_ca_1;
  } else {
    errval = fscanf(getfile," %d %d %lf %d %d %d %lf %d ",
		    &num_ca_1,&ca_chrg_1,Rcat_1,Zcat_1,
		    &num_ca_2,&ca_chrg_2,Rcat_2,Zcat_2);
    if ( errval != 8 ) { printf("Error reading input file: loc. 5 errval=%d.\n",errval); exit(0); }
    num_ca = num_ca_1 + num_ca_2;
  }

  if ( debug > 3 ){
    printf("* get parms: filepointer: %p\n",getfile);
    printf("* get parms: cell: %15.10f%15.10f%15.10f\n",cell.a,cell.b,cell.c);
    printf("* get parms: cell: %15.10f%15.10f%15.10f\n",cell.alph,cell.beta,cell.gamm);
    printf("* get parms: an_coord    = %d\n",an_coord);
    printf("* get parms: num_an      = %d\n",num_an);
    printf("* get parms: an_chrg_c   = %f\n",an_chrg_c);
    printf("* get parms: an_chrg_v   = %f\n",an_chrg_v);
    printf("* get parms: d           = %f\n",d);
    printf("* get parms: Rcent       = %f\n",*Rcent);
    printf("* get parms: Rvert       = %f\n",*Rvert);
    printf("* get parms: Zcent       = %d\n",*Zcent);
    printf("* get parms: Zvert       = %d\n",*Zvert);
    printf("* get parms: n_type_ca   = %d\n",n_type_ca);
    printf("* get parms: num_ca_1    = %d\n",num_ca_1);
    printf("* get parms: ca_chrg_1   = %d (type is int, 1.0, 2.0, etc. will cause err reading subsequent values)\n",ca_chrg_1);
    printf("* get parms: Rcat_1      = %f\n",*Rcat_1);
    printf("* get parms: Zcat_1      = %d\n",*Zcat_1);
    printf("* get parms: num_ca_2    = %d\n",num_ca_2);
    printf("* get parms: ca_chrg_2   = %d\n",ca_chrg_2);
    printf("* get parms: Rcat_2      = %f\n",*Rcat_2);
    printf("* get parms: Zcat_2      = %d\n",*Zcat_2);
    printf("* get parms: exearly     = %d (1=don't get atom pos)\n",exearly);
  }

  if (exearly==1) { rewind(getfile); return; }

  if ( debug>3 ) printf("Getting cations.\n");
  for(i=0;i<num_ca;i++) {
    if ( debug>3 ) printf("Getting cation %d.\n",i);
    fscanf(getfile," %lf %lf %lf ",
	   &at[i].center.x,
	   &at[i].center.y,
	   &at[i].center.z);
  }
  
  if ( debug>3 ) printf("Getting anions.\n");
  for(i=0;i<num_an;i++) {
    for(j=0;j<num_vert;j++){
      if ( an_coord==0 ) fscanf(getfile," %lf %lf %lf ",
				&tetr[i].v[j].x,
				&tetr[i].v[j].y,
				&tetr[i].v[j].z);
      if ( an_coord==1 ) fscanf(getfile," %lf %lf %lf ",
				&octa[i].v[j].x,
				&octa[i].v[j].y,
				&octa[i].v[j].z);
      if ( an_coord==2 ) fscanf(getfile," %lf %lf %lf ",
				&dimr[i].v[j].x,
				&dimr[i].v[j].y,
				&dimr[i].v[j].z);
    }
  }

  rewind(getfile);
  return;
}

void print_restart(void)
{
  restartfile = fopen("restart.dat","w");
  if ( restartfile==NULL ){
    printf("Can't open restart.dat file.\n");
    exit(0);
  }

  /* print out restart parameters */
  print_restart_fileptr(restartfile,' ');
  
  fclose(restartfile);
  return;
}

void print_restart_fileptr(FILE *fp,char K)
{

  char lead;

  int i,j,num_vert=0;
  
  if ( an_coord==0 ) num_vert=5;
  if ( an_coord==1 ) num_vert=7;
  if ( an_coord==2 ) num_vert=2;

  if ( K == ' ' ) lead = ' ';
  else lead = '*';


  /* print out restart parameters */
  fprintf(fp,
	  "%c%c %18.12f%18.12f%18.12f\n"
	  "%c%c %18.12f%18.12f%18.12f\n",
	  lead,K,cell.a,cell.b,cell.c,
	  lead,K,cell.alph,cell.beta,cell.gamm);
  
  if ( n_type_ca == 1 )
    fprintf(fp,
	    "%c%c %6d%6d\n"
	    "%c%c %10.4f%10.4f\n"
	    "%c%c %18.12f%18.12f%18.12f%6d%6d\n"
	    "%c%c %6d\n"
	    "%c%c %6d%6d\n"
	    "%c%c %18.12f%6d\n",
	    lead,K,an_coord,num_an,
	    lead,K,an_chrg_c,an_chrg_v,
	    lead,K,d,Rcent,Rvert,Zcent,Zvert,
	    lead,K,n_type_ca,
	    lead,K,num_ca_1,ca_chrg_1,
	    lead,K,Rcat_1,Zcat_1);
  else
    fprintf(fp,
	    "%c%c %6d%6d\n"
	    "%c%c %10.4f%10.4f\n"
	    "%c%c %18.12f%18.12f%18.12f%6d%6d\n"
	    "%c%c %6d\n"
	    "%c%c %6d%6d\n"
	    "%c%c %18.12f%6d\n"
	    "%c%c %6d%6d\n"
	    "%c%c %18.12f%6d\n",
	    lead,K,an_coord,num_an,
	    lead,K,an_chrg_c,an_chrg_v,
	    lead,K,d,Rcent,Rvert,Zcent,Zvert,
	    lead,K,n_type_ca,
	    lead,K,num_ca_1,ca_chrg_1,
	    lead,K,Rcat_1,Zcat_1,
	    lead,K,num_ca_2,ca_chrg_2,
	    lead,K,Rcat_2,Zcat_2);
  
  for(i=0; i<num_ca; i++) {
    fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
	    lead,K,
	    at[i].center.x,
	    at[i].center.y,
	    at[i].center.z);
  }
  for(i=0; i<num_an; i++) {
    for(j=0; j<num_vert; j++) {
      if ( an_coord==0 ) {
	if ( j==0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    tetr[i].v[j].x,
			    tetr[i].v[j].y,
			    tetr[i].v[j].z);
	if ( j!=0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    tetr[i].v[j].x,
			    tetr[i].v[j].y,
			    tetr[i].v[j].z);
      }
      if ( an_coord==1 ) {
	if ( j==0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    octa[i].v[j].x,
			    octa[i].v[j].y,
			    octa[i].v[j].z);
	if ( j!=0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    octa[i].v[j].x,
			    octa[i].v[j].y,
			    octa[i].v[j].z);
      }
      if ( an_coord==2 ) {
	if ( j==0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    dimr[i].v[j].x,
			    dimr[i].v[j].y,
			    dimr[i].v[j].z);
	if ( j!=0 ) fprintf(fp,"%c%c %18.12f%18.12f%18.12f\n",
			    lead,K,
			    dimr[i].v[j].x,
			    dimr[i].v[j].y,
			    dimr[i].v[j].z);
      }
    }
  }
  

}
