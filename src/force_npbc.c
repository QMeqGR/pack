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


extern int debug;
extern int numat,efunc;
extern double rep_epsilon;

extern struct cellprm cell;
extern struct atom *p;
extern struct atom *s;
extern struct obj *object;

/* DEBUG - set to 4 for normal operation */
int db= 0;

/********************************************************

  Calculate force on an atom from electrostatics.
  
  estat_init() should be called before these functions
  to fill the p structure with information on location
  from the object structures.

  The forces are then calculated using non-periodic
  boundary conditions for simplicity.  This obviously
  only works for the nano-particle calculations since
  we will ignore interactions from atoms in neighbor cells.

  Once a force is calculated for the object, it
  will be used to move the object.  If the move is rejected,
  due to increasing overlap, then the variable reject_move
  will be set for the object.  If this variable is set,
  then the next time the forces are calculated, the force
  will be removed along the direction of the increasing
  overlap ( Fprime = F - (F . r) r ).

********************************************************/

int calc_forces(void)
{
  int i,k,err=0;
  struct vector ftmp={0};

  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    ftmp = force_on_obj(i);
    if ( debug>db ) printf("calc_forces: ob %2d: %15.4e%15.4e%15.4e\n",i,ftmp.x,ftmp.y,ftmp.z);
  }

  return err;
}

struct vector force_on_atom(int i, struct vector Frej)
{
  int iat1;
  int been_called=0;
  static int nat = 0;

  double Rsq=0;
  struct vector t={0},tperp={0},tperpn,F={0};
  struct vector Fn={0},dtau={0},dtaun={0};

  if (!been_called) {
    nat = numat;
    been_called = 1;
  }

  for(iat1=0; iat1<nat; iat1++){
    if ( i == iat1 ) continue;
    if ( p[i].chrg == 0 || p[iat1].chrg == 0 ) continue;
      dtau = makevec( p[iat1].x-p[i].x,
		      p[iat1].y-p[i].y,
		      p[iat1].z-p[i].z );

      Rsq = vmagsq( dtau );
      t = vsmult( -p[iat1].chrg * p[i].chrg / Rsq , dtau );
      tperp = vadd( t , vsmult( -vdotprod( t , Frej ) , Frej ) );
      F = vadd(tperp,F);

      if ( debug > db ){

	printf("> force_on_atom i due to iat1:\n");
	printf("> i   =%2d\t%15.4f%15.4f%15.4f\tq=%8.4f\n"
	       "> iat1=%2d\t%15.4f%15.4f%15.4f\tq=%8.4f\tR=%10.4f\n",
	       i,p[i].x,p[i].y,p[i].z,p[i].chrg,
	       iat1,p[iat1].x,p[iat1].y,p[iat1].z,p[iat1].chrg,sqrt(Rsq));

	dtaun = normvec(dtau);
	printf("> dtaun  =  \t%15.4e%15.4e%15.4e   (normalized cart)\n",
	       dtaun.x,dtaun.y,dtaun.z);

	printf("> F   =  \t%15.4e%15.4e%15.4e   (tot force cart)\n",F.x,F.y,F.z);

	Fn = normvec(F);
	printf("> Fnorm  =  \t%15.4e%15.4e%15.4e   (normalized cart)\n",
	       Fn.x,Fn.y,Fn.z);

	tperpn = normvec( tperp );
	printf("> tperp = \t%15.4e%15.4e%15.4e   (normalized cart)\n",
	       tperpn.x,tperpn.y,tperpn.z);
	printf("> Frej = \t%15.4e%15.4e%15.4e\n",Frej.x,Frej.y,Frej.z);
	printf("> Frej  . tperp =    %f\n",vdotprod(Frej,tperp));
	printf("> Fnorm . Frej  =    %f\n",vdotprod(Fn,Frej));
	printf(">\n");

      }
  }

  return (F);
}


struct vector force_on_obj(int obj_number)
{

  int j,nv=0;

  struct vector F={0},Fn={0},Frej={0};

  if ( object[obj_number].reject_move == 1 ) {
    if ( debug > db ) printf("* force_on_obj[%d]: reject_move=1 \n",obj_number);
    Frej = normvec( object[obj_number].force );
  }

  /* assume object CM is at center of anion.
     This won't work for all objects, but ...
     Sum the forces for all atoms and return
     the force */

  nv = object[obj_number].num_vert;

  for(j=0; j<nv; j++){
    F = vadd( F , force_on_atom( object[obj_number].p_num[j] , Frej ) );
  }
  
  F = normvec(F);
  object[obj_number].force = F;

  if ( debug > db ) {
    Fn = normvec(F);
    printf("force_on_an: force on an %d: %15.5f%15.5f%15.5f\n",
	   obj_number,Fn.x,Fn.y,Fn.z);
  }

  return (F);
}


void force_step(void)
{
  static int been_called=0;

  int i,j;

  struct vector trans_vec={0}, TP={0};
  struct matrix L;
  struct Energy E0={0},E={0};

  if ( been_called == 0 ){
    for(j=0; j<MAX_OBJS && object[j].used==1; j++) object[j].reject_move = 0;
    been_called = 1;
  }

  E0 = (*Total_Energy)();

  if ( debug > db ) {
    printf("force_step: E0: ecc= %10.4e   eng= %10.4e\n",E0.ecc,E0.eng);
  }
  
  for(j=0; j<FORCE_STEPS; j++){

    if ( debug > db ) printf("+++ forces +++\n");
    estat_init();
    calc_forces();

    /* LOOP over OBJECTS */
    for(i=0; i<MAX_OBJS && object[i].used==1; i++){

      L = L_e3(&cell);

      trans_vec = vsmult( FORCE_STEP_SIZE_AN , object[i].force );
      if ( debug > db ) {
	printf("force_step: moving object  %d: %10.4f%10.4f%10.4f\n",i,
			    trans_vec.x, trans_vec.y, trans_vec.z);
      }

      /* rezone the trans vec if it takes object out of bounds */
      TP = object[i].v[0];
      trans_vec = vsub( rezone( vadd(TP,trans_vec) ) , TP );
      trans_object( i , &L , trans_vec );
      E = (*Total_Energy)();

      if ( debug > db ) {
	printf("force_step: E: ecc= %10.4e   eng= %10.4e\n",E.ecc,E.eng);
      }
      if ( E.eng > E0.eng ) {
	if ( debug > db ) printf("force_step: Overlap increasing, reversing object move\n");
	trans_object( i , &L , vsmult( -1.0 , trans_vec ) );
	object[i].reject_move = 1;
      } else {
	object[i].reject_move = 0;
	E0 = E;
      }

    }

  }

  return;
}
