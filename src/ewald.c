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


/***************************************************************/
/*                                                             */
/*  ewald.c                                                    */
/*  Calculates the ewald energy of a structure given a modified*/
/*  CONTCAR file.                                              */
/*                                                             */
/***************************************************************/
/*                                                             */
/*                                                             */
/*                                                             */
/***************************************************************/
/*                                                             */
/*  adapted from contcar_pdf.c                                 */
/*  02 Feb 2006                                                */
/*                                                             */
/***************************************************************/


/*
 * Compile:
 *
 * gcc ~/src/pack/ewald.c ~/src/pack/packlib.c -lm -Wall -o ~/bin/ewald
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global_defs.h"
#include "global_vars.h"
#include "packlib.h"

#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_BUGFX 1

/*

 1.0.1
 Tue Jan 16 14:30:40 PST 2007
 -fixed input reading of first line...


*/

#define Pi 3.14159265358979
#define MAX_ATOMS 1000
#define NUM_SKIP 11
#define ERRORLIMIT 1e-10

/**************************/
/*   global variables     */
/**************************/
FILE *filepointer = NULL;
extern char *optarg;
extern int optind,opterr,optopt;

/*************************/
/*  function prototypes  */
/*************************/
int getopt(int argc, char * const argv[], const char *optstring);
int main(int,char **);

/*************************/
/*     Main              */
/*************************/
int main(int argc,char *argv[])
{
  char *fname="QPOSCAR";
  char tmp,tmp2[2];
  
  int atcount=0;
  int wordcount=0;
  int on_word=0;
  int on_colon=0;
  int on_char=0;
  int on_space=0;
  int past_colon=0;
  int i;
  extern int debug;
  int errval;
  int optnum;
  int count=0;
  int ntypat=0;
  int nat=0;
  int *typat;
  int *ntyp;

  double scalefact=0;
  double *q;
  double eta=0,ecc=0;
  double errlim=ERRORLIMIT;
  
  struct atom *at;
  struct vector r1,r2,r3; /* lattice vectors */
  struct vector T;
  struct cellprm cell;
  struct matrix L;
  

  printf("#\n");
  printf("###############################\n");
  printf("#  ewald  version %d.%d.%d\n",
	 VERSION_MAJOR,VERSION_MINOR,VERSION_BUGFX);
  printf("###############################\n");
  printf("# Source file: " __FILE__ "\n");
  printf("# Compile date: " __DATE__ "\n");
  printf("# Compile time: " __TIME__ "\n");
  printf("# Commandline: ");
  for(i=0; i<argc; i++) printf("%s ",argv[i]); printf("\n#\n");

  /* get command line options */
  while ( (optnum = getopt(argc, argv,
			   "hFf:d:e:E:")) != -1 ) {
    switch(optnum) {
    case 'h':
      printf("\nCommand line options (* indicates required):\n"
	     "  -h --- help\n"
	     "\n         Input parameters:\n\n"
	     "  -f -*- input filename (default QPOSCAR)\n"
	     "  -F --- file format\n"
	     "  -e --- input errlim (default %.2e)\n"
	     "  -E --- input eta\n"
	     "  -d --- debug level (integer) 0...4\n"
	     "\n\n",ERRLIM);
      exit(0);
      
    case ':':
      printf("option needs a value\n");
      exit(0);
      
    case 'f':
      fname  = optarg;
      break;

    case 'F':
      printf("\n\nAssumed file format is QPOSCAR output from\n"
	     "a pack run. (A modified CONTCAR or POSCAR file.)\n\n"
	     "Z: z1 z2 z3 z4 \n"
	     "q1 q2 q3 q4 \n"
	     "1.0 [no text following scale factor!!!]\n"
	     "rest of file normally\n\n");
      exit(0);
      break;

    case 'd':
      debug  = atoi(optarg);
      break;

    case 'e':
      errlim  = atof(optarg);
      break;

    case 'E':
      eta  = atof(optarg);
      break;
      
    case '?':
      printf("unknown option: %c\n",optopt);
      exit(0);
      break;
    }
  }


  /**************************************************************************/
  /**************************************************************************/

  /* open file: will be CONTCAR if not given with -f switch */
  if ( NULL==(filepointer=fopen(fname,"r")) ) {
    printf("Error opening file %s\n",fname);
    printf("Try: contcar_pdf -h\n");
    exit(0);
  }


  /* Read in the header data.
   * Count the number of words
   * on the line and subtract one.
   * This is for a line that looks like
   * Z: 1 2 3 4 5
   */
  while ( (tmp=fgetc(filepointer)) != '\n' ){
    if ( debug>4) printf("tmp=%c\n",tmp);
    if ( tmp == ' ' && on_colon==1 ) past_colon=1;
    if ( tmp == ' ' ) { on_space=1; on_char=0; }
    if ( tmp != ' ' ) { on_space=0; on_char=1; }
    if ( tmp == ' ' && on_word==1 ) on_word=0;
    if ( tmp == ':' ) on_colon=1;
    if ( past_colon && on_char && on_word==0){
      on_word=1;
      on_char=0;
      if ( debug>4) printf("wordcount=%d, incrementing to %d\n",wordcount,wordcount+1);
      wordcount++;
    }
  }
  ntypat = wordcount;
  if ( debug > 0 ) printf("# ntypat = %d\n",ntypat);

  /* allocate space for atom types and number of each */
  typat = (int *)malloc( ntypat * sizeof( int ) );
  if ( !typat )
    {
      printf("%s\n","Error allocating memory for num atom types.");
      exit(0);
    }
  ntyp = (int *)malloc( ntypat * sizeof( int ) );
  if ( !ntyp )
    {
      printf("%s\n","Error allocating memory for atom type numbers.");
      exit(0);
    }
  q = (double *)malloc( ntypat * sizeof( double ) );
  if ( !q )
    {
      printf("%s\n","Error allocating memory for atom charges.");
      exit(0);
    }

  /* next line is the MODIFIED line with the charges for each Z */
  for(i=0; i<ntypat; i++){
    errval = fscanf(filepointer, " %lf ",&q[i]);
    if ( debug>0 ) printf("# q[%d] = %9.6f\n",i,q[i]);
    if ( errval != 1 ){
      printf("Error reading charges. errval=%d\n",errval);
      exit(0);
    }
  }
  
  /* next line is the scale factor */
  errval=fscanf(filepointer," %lf ", &scalefact);
  if ( errval != 1 ){
    printf("Error reading scale factor.\n");
    exit(0);
  }
  if ( debug > 0 ) printf("# scale = %f\n",scalefact);

  /* next three lines are the lattice vectors */
  errval=fscanf(filepointer," %lf %lf %lf ",&r1.x,&r1.y,&r1.z);
  if ( errval != 3 ){
    printf("Error reading lattice vector 1.\n");
    exit(0);
  }
  if ( debug > 0 ) printf("# r1: %20.15f%20.15f%20.15f\n",r1.x,r1.y,r1.z);

  errval=fscanf(filepointer," %lf %lf %lf ",&r2.x,&r2.y,&r2.z);
  if ( errval != 3 ){
    printf("Error reading lattice vector 2.\n");
    exit(0);
  }
  if ( debug > 0 ) printf("# r2: %20.15f%20.15f%20.15f\n",r2.x,r2.y,r2.z);

  errval=fscanf(filepointer," %lf %lf %lf ",&r3.x,&r3.y,&r3.z);
  if ( errval != 3 ){
    printf("Error reading lattice vector 3.\n");
    exit(0);
  }
  if ( debug > 0 ) printf("# r3: %20.15f%20.15f%20.15f\n",r3.x,r3.y,r3.z);

  /* read in the number of each type of atom */
  for(i=0; i<ntypat; i++) {
    errval=fscanf(filepointer," %d ", &ntyp[i] );
    nat += ntyp[i];
    if ( errval != 1 ){
      printf("Error filling typat.\n");
      free(ntyp);
      exit(0);
    }
    if ( debug > 0 ) printf("# ntyp[%d] = %d\n",i,ntyp[i]);
  }

  if ( debug > 0 ) printf("# nat = %d\n",nat);

  /* allocate space for data points */
  at = (struct atom *)malloc( nat * sizeof( struct atom ) );
  if ( !at )
    {
      printf("%s\n","Error allocating memory for atoms.");
      free(typat);
      exit(0);
    }

  /* read in the data from the file */
  rewind(filepointer);
  fscanf(filepointer," %s ", tmp2);
  for (i=0; i<ntypat; i++) {
    errval=fscanf(filepointer," %d ",&typat[i]);
    if (errval != 1) {
      printf("Bad file type. Reading Z values.\n");
      printf("errval = %d\n",errval);
      exit(0);
    }
    if ( debug > 0 ) printf("# typat[%d] = %d\n",i,typat[i]);
  }

  /* skip the next NUM_SKIP pieces of data, they were already read above */
  for(i=0;i<(NUM_SKIP+2*ntypat);i++){
    errval = fscanf(filepointer," %*s ");
  }
  
  /* read in atom positions */
  for (i=0; i<nat; i++) {
    errval = fscanf(filepointer," %lf %lf %lf ", &at[i].x, &at[i].y, &at[i].z );
    if ( debug > 2 ) printf("# %20.15lf%20.15lf%20.15lf\n",at[i].x,at[i].y,at[i].z);
    if (errval != 3) {
      printf("Bad file type.  Reading atom positions.\n");
      printf("[%d] errval = %d\n",i,errval);
      free(typat);
      free(at);
      exit(0);
    }

  }

  /* assign the atom Zs and Qs */
  atcount=1;
  count=0;
  for(i=0;i<nat;i++){
    if ( debug > 4 ) printf("atcount = %d   count = %d    ntyp[] = %d\n",
			    atcount,count,ntyp[count]);
    if ( atcount <= ntyp[count] ) {
      at[i].Z = typat[count];
      at[i].chrg = q[count];
      atcount++;
    }
    if ( atcount-1 == ntyp[count] ) { count++; atcount=1; }
  }

  if ( debug > 0 )
  for(i=0;i<nat;i++){
    printf("# atom [%3d]: %15.10f%15.10f%15.10f%5d%10.5f\n",i,
	   at[i].x,at[i].y,at[i].z,at[i].Z,at[i].chrg);
  }  


  /***********************************************************************/
  /***********************************************************************/

  cell.bas.A.x = r1.x;
  cell.bas.A.y = r1.y;
  cell.bas.A.z = r1.z;

  cell.bas.B.x = r2.x;
  cell.bas.B.y = r2.y;
  cell.bas.B.z = r2.z;

  cell.bas.C.x = r3.x;
  cell.bas.C.y = r3.y;
  cell.bas.C.z = r3.z;

  cell.a = vmag( r1 );
  cell.b = vmag( r2 );
  cell.c = vmag( r3 );

  cell.alph = acos( vdotprod(cell.bas.B,cell.bas.C)/(cell.b*cell.c) );
  cell.beta = acos( vdotprod(cell.bas.A,cell.bas.C)/(cell.a*cell.c) );
  cell.gamm = acos( vdotprod(cell.bas.B,cell.bas.A)/(cell.b*cell.a) );

  if ( debug > 0 ){
    printf("cell: %15.10f%15.10f%15.10f\n",cell.a,cell.b,cell.c);
    printf("angl: %15.10f%15.10f%15.10f\n",R2D*cell.alph,R2D*cell.beta,R2D*cell.gamm);
  }


  /* refill the atom positions with cartesian coordinates */
  L = L_e3(&cell);

  if ( debug > 0 ){
    printf("L:\n");
    printf("%15.10f%15.10f%15.10f\n",L.r11,L.r12,L.r13);
    printf("%15.10f%15.10f%15.10f\n",L.r21,L.r22,L.r23);
    printf("%15.10f%15.10f%15.10f\n",L.r31,L.r32,L.r33);
  }

  if ( debug > 2 ) printf("> %15s%15s%15s%15s\n","X","Y","Z","chrg");
  for(i=0; i<nat; i++) {
    T = makevec( at[i].x, at[i].y, at[i].z );
    T = rezone( T );
    T = cart(&L,&T);
    at[i].x = T.x;
    at[i].y = T.y;
    at[i].z = T.z;
    if ( debug > 2 ) printf("> %15.10f%15.10f%15.10f%15.10f\n",at[i].x,at[i].y,at[i].z,at[i].chrg);
  }


  if ( tol(eta,0.0,1e-9) ) eta = find_Ewald_eta(cell.bas,errlim,nat,at);
  ecc = Eion(cell,at,nat,errlim,eta);

  if ( debug > 0 ) printf("Ewald eta: %15.10f\n",eta);
  printf("Ewald energy for this structure: %15.10f [Ha]\n",ecc);

  /*  END of Program */

  if (typat) free(typat);
  if (at) free(at);
  
  return(0);
}


