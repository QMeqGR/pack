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

extern int numat;
extern int basin_hop,scale_cell_do;
extern int debug,get_estat;
extern int non_per_ewald;

extern double scale_factor,ext_press,errlim;
extern double sqsh_vol_pct,cell_vol_max;

extern struct cellprm cell;
extern struct atom *s;
extern struct atom *p;
extern struct obj *object;
extern struct Energy E_best;

extern FILE *tmpfp;

extern char tempfile[25];


/* added for faster deltaE Ewald sum */
extern int use_deltaEion;

/********************************************/

struct Energy Total_Energy(void)
{


  static int nat=0;
  static int been_called=0;

  int i,j;
  double eta=0;

  static struct Energy E;

// for deltaEion
  static double Eion_last;
  int nmoved;
  static int * ixmoved;
  static struct atom * at_old;

  if ( !been_called ) {
    E.ecc=0.0; E.eng=0.0; E.pf=0.0; E.bnd=0.0;
    
    E_best.ecc = 1e10;
    E_best.eng = 1e10;
    E_best.pf = 0;
    E_best.ecc_eng = 1e10;
    E_best.tot = 1e10;

    nat = numat;

// fist time, force full Ewald evaluation
    use_deltaEion= 0;
    ixmoved= malloc(nat*sizeof(int));
    if (! at_old)   at_old = (struct atom *)malloc( nat * sizeof( struct atom ) );
    if ( !ixmoved || !at_old)
    {
      printf("%s\n","Error allocating memory for ixmoved or at_old.");
      exit(0);
    }
    memcpy(at_old, p, nat * sizeof( struct atom ));

  been_called=1;
  }

  /* calculate the packing fraction and cell volume */
  E.a = cell.a; E.b = cell.b; E.c = cell.c;
  E.alph = cell.alph; E.beta = cell.beta; E.gamm = cell.gamm;
  E.vol = cell_volume(&cell);
  E.pf = pack_frac();
  E.orth = E.vol / (cell.a*cell.b*cell.c);


  /* BASIN HOPPING CODE */
  /* if we are basin hopping, scale the cell upward before the energy calc */
  if ( basin_hop>0 ) scale_cell( scale_factor, +1, scale_cell_do );

  energ_init( );
  E.eng = Energy_rep(E.orth);

  /**************************************************************/
  /* Restrictive boundary settings:  These are rather arbitrary */
  /* add in 'out of bounds' energy if the angles get too large  */
  E.bnd = 0.0;
  if ( cell.alph < ANG_ALPH_MIN || cell.alph > ANG_ALPH_MAX ) E.bnd += 1e3;
  if ( cell.beta < ANG_BETA_MIN || cell.beta > ANG_BETA_MAX ) E.bnd += 1e3;
  if ( cell.gamm < ANG_GAMM_MIN || cell.gamm > ANG_GAMM_MAX ) E.bnd += 1e3;
  if ( E.orth < sqsh_vol_pct ) E.bnd +=1e5;
  if ( E.vol >= cell_vol_max ) E.bnd +=1e5;

  /* get the electrostatic energy */
  if ( get_estat ){

    estat_init( );
    if ( non_per_ewald ){
      E.ecc = Eion_nopbc(cell,p,nat,errlim,eta);
    } else {
      eta = find_Ewald_eta(cell.bas,errlim,nat,p);
// prepare to call deltaEion
        nmoved=0;
        for (i=0; i<nat; i++) {
          if (p[i].moved) {
            ixmoved[nmoved++]=i;
          }
        }
// auto decide whether to use fast deltaE
      if ((nmoved<8)&&(nmoved<nat/2)&&(use_deltaEion)) {
//      if (use_deltaEion) {
        E.ecc = Eion_last +  deltaEion(nmoved, ixmoved, cell, at_old, p, nat, errlim, eta);
      }
      else {
        E.ecc = Eion(cell,p,nat,errlim,eta);
        use_deltaEion=1;
      }
// save coordinates and Eion of last visit
      Eion_last = E.ecc;
      memcpy(at_old, p, nat * sizeof( struct atom ));

    }

    if ( debug > 4 ) {
      printf("* Total_Energy: E_ss   = %10.5e\n",E.eng);
      printf("* Total_Energy: eta    = %10.5e\n",eta);
      printf("* Total_Energy: E_cc   = %10.5e\n",E.ecc);
      printf("* Total_Energy: E_bnd  = %10.5e\n",E.bnd);
      printf("* Total_Energy: E_orth = %10.5e\n",E.orth);
      printf("* Total_Energy: E_vol  = %10.5e\n",E.vol);
    }
    
  }

  /* BASIN HOPPING CODE */
  /* if we are basin hopping, scale the cell downward after the energy calc */
  if ( basin_hop>0 ) scale_cell( scale_factor, -1, scale_cell_do );

  /* calculat this for the 'best' structure stuff */
  E.ecc_eng = E.ecc + E.eng;


  /* THE ENERGY */
  E.tot = E.ecc + E.eng + E.bnd;
  if ( ext_press > 1e-12 ) E.tot += ext_press * E.vol;

  /* Test against 'best' structure.  This is a more simple test than is
   * applied at the end of the metro run.  This is because of the
   * complication of having the cell volume go in and out of the
   * energy functional above (depending on the pack frac).
   * We want the best structure to have lower ecc+eng, not worrying
   * about the cell volume so much.
   */
  if ( (E.ecc_eng < E_best.ecc_eng) && E.orth>sqsh_vol_pct ) {
    if ( debug > 2 ) printf("* writing to %s E_low=%f\n", tempfile, E.ecc_eng);
    print_tmp();
    E_best = E;
  }


  /* set 'moved' flags to zero */
  for(i=0;i<MAX_OBJS && object[i].used==1;i++){
    for(j=0;j<MAX_VERTS && object[i].vused[j]==1;j++){
      object[i].vmoved[j]=0;
    }
  }

  return(E);
}

