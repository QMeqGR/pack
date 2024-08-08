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

extern int debug,n_simp,num_ob;

extern double **sim_pp,*sim_y,*fn,d,simp_lamb;

extern struct cellprm cell;
extern struct obj *object;


/* The simplex_init function fills the array fn with initial values
   of the cell parameters, the cation atoms, the anion centers only,
   and an initial orientation of (0 0 0) for the anions.  It also
   calls the energy function to get the energy for each starting
   position of the simplex and fills the sim_pp and sim_y arrays.*/
void simplex_init(void)
{
  int i,j,loc;
  double lambda=0;
  double *orig=NULL;
  struct Energy E;

  orig = (double *)malloc( n_simp * sizeof( double ) );
  if ( !orig ) {
    printf("* Memory allocation error for orig.\n");
    exit_pack(0);
  }

  if ( debug > -1 ) printf("*         - Initializing simplex arrays -\n");

  fn[0] = cell.a;
  fn[1] = cell.b;
  fn[2] = cell.c;
  fn[3] = cell.alph;
  fn[4] = cell.beta;
  fn[5] = cell.gamm;

  /* the object centers */
  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    loc = 6 + 3*i;
    fn[loc+0] = (*(object+i)).v[0].x;
    fn[loc+1] = (*(object+i)).v[0].y;
    fn[loc+2] = (*(object+i)).v[0].z;
  }

  /* for the anion orientations, intialize all as zero */
  for(i=0; i<num_ob; i++) {
    loc = (6 + 3*num_ob) + 3*i;
    fn[loc+0] = 0;
    fn[loc+1] = 0;
    fn[loc+2] = 0;
  }

  if ( debug > 4 ) for(i=0;i<n_simp;i++) printf("* fn[%3d] = %12.8f\n",i,fn[i]);

  /* fill the simplex starting vertices */
  for(i=0; i<(n_simp+1); i++){
    for(j=0; j<n_simp; j++){
      lambda=0.0;
      if (i==j+1) lambda=simp_lamb;
      sim_pp[i][j] = fn[j]*(1+lambda);
      if (i==j+1 && i>(n_simp-3*num_ob)) sim_pp[i][j]= SIMP_DEG * 0.0174533;  /* no. degrees */
      orig[j] = fn[j];
    }
  }

  /* find the energy for the starting vertices */
  for(i=0; i<(n_simp+1); i++){

    if ( debug > 4 ) {
      printf("simplex_init: before[%d]:\n",i);
      for(j=0;j<n_simp;j++) printf("%.9f ",sim_pp[i][j]); printf("\n");
      debug_block_metro( );
    }

    simplex_trans(sim_pp[i]);

    if ( debug > 4 ) {
      printf("simplex_init: after[%d]:\n",i);
      debug_block_metro( );
    }

    E = Total_Energy();

    simplex_trans(orig);
    if ( debug > 4 ) {
      printf("simplex_init: after return [%d]:\n",i);
      debug_block_metro( );
    }

    sim_y[i] = E.tot;
  }

  if ( debug > 2 ) {
    printf("* Initial Simplex array:\n*");
    for(i=0; i<(n_simp+1); i++){
      for(j=0; j<n_simp; j++){
	printf(" %10.6f ",sim_pp[i][j]);
	if (j==n_simp-1) printf("\n*");
      }
    }
    printf("\n");
    printf("* Energies:\n*");
    for(i=0; i<(n_simp+1); i++){
      printf(" %10.6e ",sim_y[i]);
    }
    printf("\n");
  }

  if (debug>2) printf("*\n* ------------ Simplex Init Finished -------------------- *\n");

  free(orig);
  return;
}
