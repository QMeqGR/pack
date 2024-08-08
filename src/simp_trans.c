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

extern int numat,num_ob;
extern int debug,n_simp,mvout;
extern int non_per_ewald,lat_parm_chng;

extern double d,c2vtol;
extern double *rot_var;

extern struct cellprm cell;
extern struct atom *s;
extern struct obj *object;

/*
 * This function calculates the translations that the
 * minimization routine wants to do to change the energy
 * function value.
 *
 */
void simplex_trans(double *fn)
{
  int i,loc=0;
  static int frame_reduce=0;
  static int been_called=0;
  static int nat=0;
  
  double da,db,dc,dalph,dbeta,dgamm,dx,dy,dz;
  static double phi=0,the=0,psi=0;

  struct vector TP={0};
  struct vector v_lat_chng,ang_chng,trans_vec;
  struct matrix R,L;

  if ( !been_called ) {

    /* initialize rot_var */
    for (i=0;i<3*num_ob;i++) rot_var[i]=0.0;

    nat = numat;
    been_called=1;
  }

  if ( debug > 4 ) {

    printf("* simplex_trans: (beginning)\n");
    debug_block_metro( );
    get_min_distances( );

    printf("\n* simplex_trans: fn:\n");
    for(i=0;i<n_simp;i++) {
      printf(" %f",fn[i]);
    }
    printf("\n");

  }

  /* these are loc 0-5 */
  da = fn[0] - cell.a;
  db = fn[1] - cell.b;
  dc = fn[2] - cell.c;
  dalph = fn[3] - cell.alph;
  dbeta = fn[4] - cell.beta;
  dgamm = fn[5] - cell.gamm;
  if ( debug > 4 ) {
    printf("* simplex_trans: cell: %15.10f%15.10f%15.10f%15.10f%15.10f%15.10f\n",
	   cell.a,cell.b,cell.c,cell.alph,cell.beta,cell.gamm);
    printf("* simplex_trans: cell_chng = %e %e %e %e %e %e\n",
	   da,db,dc,dalph,dbeta,dgamm);
  }
  v_lat_chng = makevec(da   ,db   ,dc);
  ang_chng = makevec(dalph,dbeta,dgamm);
  /* Rescale the object locations. */
  /* The rescaling needs the pre-changed cell parameters and the change params. */
  /* Note: rescale_anions will call L_e3 */
  if ( non_per_ewald==0 && lat_parm_chng==1 )
    rescale_objects( &cell, v_lat_chng, ang_chng );
  /* simply don't update the lat parm changes if doing nanoparticle */
  /* or if lat parm changes are turned off */

  if ( debug > 4 ) {
    printf("* simplex_trans: post rescale:\n");
    debug_block_metro( );
    get_min_distances( );
    printf("* simplex_trans: begin object translations:\n");
  }

  /* before doing object translations, must recalc L */
  L = L_e3(&cell);

  /* object translations */
  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    if ( object[i].vfxd[0]==1 ) continue; /* skip over fixed objects */
    loc = 6 + 3*i;
    TP = object[i].v[0];
    dx = fn[loc+0] - object[i].v[0].x;
    dy = fn[loc+1] - object[i].v[0].y;
    dz = fn[loc+2] - object[i].v[0].z;
    trans_vec = vsub( rezone( makevec( TP.x+dx, TP.y+dy, TP.z+dz ) ) , TP );
    /* don't move the object if doing a cluster calc and it gets too close to wall */
    if ( 1==trans_accept( i, vadd(TP,trans_vec), non_per_ewald ) ){
      trans_object( i, &L, trans_vec );
    }
  }

  if ( debug > 4 ) {
    printf("* simplex_trans: post object translation:\n");
    debug_block_metro( );
    get_min_distances( );
    printf("* simplex_trans: begin object rotations:\n");
  }

  /* object rotations */
  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    if ( object[i].allowrot==0 ) continue;
    loc = (6 + 3*num_ob) + 3*i;

    if (debug>4) printf("* simplex_trans: objec rotations: accessing loc[%d,%d,%d]\n",loc+0,loc+1,loc+2);

    /* undo the previous rotation */
    phi = rot_var[3*i+0];
    the = rot_var[3*i+1];
    psi = rot_var[3*i+2];
    
    R = R_ptp_inv( R_ptp(phi,the,psi) );
    if (debug>4){
      print_mat( R );
      printf("determinant R = %f\n",det3x3(R));
    }

    rot_object( i, &L, &cell, R );

    if (debug>4) printf("* simplex_trans: rotation stage 2:\n");

    /* set the new parms */
    phi = fn[loc+0];
    the = fn[loc+1];
    psi = fn[loc+2];
    /* save the parms for the next 'undo' */
    rot_var[3*i+0] = phi;
    rot_var[3*i+1] = the;
    rot_var[3*i+2] = psi;

    R = R_ptp(phi,the,psi);
    if (debug>4){
      print_mat( R );
      printf("determinant R = %f\n",det3x3(R));
    }

    rot_object( i, &L, &cell, R );
  }  

  if ( debug > 4 ) {
    printf("* simplex_trans: post objec rotation:\n");
    debug_block_metro( );
    get_min_distances( );
  }

  if ( mvout && frame_reduce== 0 ) {
    print_xbs_mvframe();
  }
  frame_reduce++;
  if (frame_reduce == 50) frame_reduce=0;
  
  return;
}
