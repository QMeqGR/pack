/*
 
 This file contains subroutines for the min-map Monte-Carlo moves,
 which can be conbined with any Markov Chain Monte-Carlo algorithm,
 in order to enhance acceptance ratio of "significant" change in the
 phase space, e.g. large or "magnified" translation
 
 Note that we are still restricting ourselves to the simplist case
 of k=1 (only one attempted move) and r=1 (forward and backford attempt
 probabilities are equal)
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


#define DRAND drand48
#define MAX_ACT 40
// #ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
// #endif

// #ifnef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
 // #endif

 #define DB -1

 extern int efunc,numat,num_ob;
 extern struct obj *object;
 extern struct cellprm cell;
 extern double LJeps[MAX_TABLE_SQR];
 extern double LJrer[MAX_TABLE_SQR];
 extern int efunc;


 /*
  locally minimize a few atoms
  supported force type
  1 = LJ (6-12)
  2 = 1/r^12 only
  supported method
  1 = Newton
  2 = Conjugated gradient
  */
 void local_min(int N_act, int * list_act, int force_type, int method, double conv) {
   double force[MAX_ACT][3];
   int i;
   double dE;

   switch (method) {
     case 1:
       // Newton
       //minByForce(N_act, list_act, conv);
       break;
     case 2:
       // CG
       printf("CG algorithm not implemented yet\n");
       exit(0);
       break;
     default:
       printf("Unknown local min algorithm %d\n", method);
       exit(0);
       break;
   }
 }


 // random vector in N dimensional hypersphere
 void hyperspherepoint(int Ndim, double veclen, double *pick) {
   //    USE MKL_VSL_TYPE
   //     USE MKL_VSL
   //      include 'mkl_vsl.f77'

   double R, x1, x2;
   int i;

   // Box-Mueler
   for(i=0; i<(Ndim+1)/2; i++) {
     x1=DRAND();
     x2=DRAND();
     pick[2*i  ]=sqrt(-2*log(x1))*cos(2*PI*x2);
     pick[2*i+1]=sqrt(-2*log(x1))*sin(2*PI*x2);
   }

   R=0;
   for (i=0; i<Ndim; i++) R+= pick[i]*pick[i];
   R=sqrt(R);
   for (i=0; i<Ndim; i++){
     pick[i]*=veclen/R;
   }
 }

double EstimateInvHessianFromForce(double f0, double rmin, double maxb, double r0) { 
   double r,r14;
   const double LJRcut=1.12246;

   // no force, not isolated atom
   if(tol(f0, 0.0, 0.00001)&& (rmin/r0<3)) return(1.0);

   if (efunc==1) {
     r = pow(48/(r0*f0), 1/13.0)*r0;
     if(DB>0) printf("   hs> f0=%f rmin=%f r0=%f r=%f\n", f0, rmin, r0, r);
     if(rmin/r0 > 1.2) { if(DB>0) printf("   hs>rmin >> r0\n");
       return((rmin - LJRcut*r0)/1.1/(f0));
     } else {
       if(f0*r0 > 30) {if(DB>0) printf("   hs> f0*r0>30\n");
	 return(fabs(( LJRcut*r0 - rmin)/1.3)/(f0));
       } else {
	 r14=pow(rmin/r0,14)/624*r0*r0;if(DB>0) printf("   hs>r14=%f\n",r14);
	 return( min(r14,maxb) );
       }
     }
   } else if (efunc==0) {
     // now force > 0 means collision, force is usually large
     return(fabs(( r0 - rmin)/1.2)/(f0));
   }

 }

#define ORTH_MINIMUM 0.7
struct vector cart_dvec(int iat1, int iat2, struct matrix *L,double orth){
  int a1,a2,a3;
  int n_a1=1,n_a2=1,n_a3=1;
  extern int non_per_ewald;
  extern struct cellprm cell;

  double d,dmin=1e12;

  struct vector ta1,ta2,Rcel,T, dtau, dtaumin;
  struct vector A,B,C;

  dtau = makevec(object[iat1].v[0].x-object[iat2].v[0].x,
		   object[iat1].v[0].y-object[iat2].v[0].y,
		   object[iat1].v[0].z-object[iat2].v[0].z);
  dtau=cart(L, &dtau);
  if (non_per_ewald==1) {
    dtaumin=dtau;
  } else if (non_per_ewald==0) { //periodic

    A = cell.bas.A;
    B = cell.bas.B;
    C = cell.bas.C;

    if ( orth < ORTH_MINIMUM ) { n_a1=2; n_a2=2; n_a3=2; }

    for (a1=-n_a1; a1<(n_a1+1); a1++) {
      ta1 = vsmult(a1,A);
      for (a2=-n_a2; a2<(n_a2+1); a2++) {
	ta2 = vsmult(a2,B);
	for (a3=-n_a3; a3<(n_a3+1); a3++) {
	  Rcel = vadd( ta1 , vadd( ta2 , vsmult(a3,C) ) );
	  T = vadd( Rcel , dtau );
	  d = vmagsq( T );
	  if ( d < dmin ) {
	    dmin=d;
	    dtaumin= T;
	  }

      }
    }
  }
  }

  return dtaumin;
}


struct vector force_on_atom_mm_tot(int i, int inc_Coulomb, struct matrix * L, double orth, double *rmin,double *rscale) {
   int iat1, Z1, Z2;
   int been_called=0;

   double Rsq, r, R6, ljrer, ljrersq;
   struct vector t={0},tperp={0},tperpn,F={0}, Fzero={0};
   struct vector Fn={0},dtau={0},dtaun={0};

   *rmin=9999;
   *rscale=0;
   Z1 = object[i].Z[0];
   for(iat1=0; iat1<numat; iat1++){
     if ( i == iat1 ) continue;
     // IF SAME OBJ or group, continue;
     dtau = cart_dvec(iat1,i,L,orth);
     Rsq = vmagsq( dtau );
     r = sqrt(Rsq);

     if (efunc==1) {
       if (*rmin> r) *rmin=r;
       Z2 = object[iat1].Z[0];
       ljrer = LJrer[ Z1 * MAX_TABLE + Z2 ];
       *rscale=ljrer;
       ljrersq = ljrer*ljrer;
       Rsq /= ljrersq;
       R6= Rsq*Rsq*Rsq;
       F = vsmult( -(48/(R6*R6) - 24/R6) /Rsq , dtau );
       // DO NOT understand the following two lines
    //    tperp = vadd( t , vsmult( -vdotprod( t , Frej ) , Frej ) );
    //    F = vadd(tperp,F);
     } else if (efunc==0) {
       Z2 = object[iat1].Z[0];
       ljrer = LJrer[ Z1 * MAX_TABLE + Z2 ];
       *rscale=ljrer;

  //  min_val = SShrd_2[ idx ];
  //  if ( (distsq+TINY)<min_val ) temp = epsilon * ( SShrd12[ idx ] * D12 );
       ljrersq = ljrer*ljrer;
       if ( Rsq + TINY < ljrersq) {
	 Rsq /= ljrersq;
	 R6= Rsq*Rsq*Rsq;
	 F = vsmult( -(12/(R6*R6) ) /Rsq , dtau );
	 // only keep the nearest distance WITH collision
	 if (*rmin> r) *rmin=r;
       }
     }
  }
  if (inc_Coulomb) F= vadd(F, force_on_atom(i, Fzero));
  return F;
}


void minByForce(int nAct,int Act[],double Etolerance, double potim,
                int nMax, struct matrix *L, struct matrix * Linv, double orth, int evaluatetot){
  //  double de0, de1, dE, ;
  const int DIM =3;
  int i, j;
  struct vector force={0}, trans_vec, TP;
  double rmin, fmag, fmagtot, temp, rscale, movesize;
  struct Energy en;
  
  for (j=0; j<nMax; j++) { 
    fmagtot = 0;
    if(DB>0) printf(" by F> act=%d round=%d en=%f\n", Act[0], j, (*Total_Energy)().tot);
    
    for (i=0; i< nAct; i++) {
      force = force_on_atom_mm_tot(Act[i],0, L,orth, &rmin, &rscale); 
      fmag= vmag(force);
      fmagtot+= fmag;
      TP=object[Act[i]].v[0];
      if(DB>0) printf(" by F > atom=%d xyz= %f %f %f\n", Act[i],TP.x,TP.y,TP.z);
      if(DB>0) printf(" by F > fmag=%f rmin=%f force=%f %f %f\n", fmag, rmin, force.x,force.y,force.z);
      movesize=EstimateInvHessianFromForce(fmag, rmin, potim, rscale)*potim;
      trans_vec= vsmult(movesize, force);
      if(DB>0) printf(" by F > hessian=%f trans vec=%f %f %f\n", movesize,trans_vec.x,trans_vec.y,trans_vec.z);
      trans_vec = vsub( rezone( makevec( TP.x+trans_vec.x, TP.y+trans_vec.y, TP.z+trans_vec.z ) ) , TP );
      trans_object( Act[i], L, trans_vec );
    }
    
    //    de0 =0;// EupdateN[pos, Na, R, fr] - de0;
    //    dE += de0;
    if( fmagtot < Etolerance) break;
  } 
  //  {If[evaluatetot, Etot[pos, N, R], 0], dE, pos}
}



/* 
 basin hopping 
 */

void BasinHoppingMove(double ** p, int N, double R, double movesize, double Enow, double * Etry) {
  double **pos;
  double dx, tmp, dE=0, Enow1, Etry1;
  int ia;
  ia = (int) DRAND()*N;
  dx = 0;//randMove[0.0, movesize, DIM];
  
  Enow1 =0;// Eupdate1[pos, Na, R, ia];
  //pos[[ia]] += dx;
  Etry1 =0;// Eupdate1[pos, Na, R, ia];
  
  //  {tmp, dE, pos} = 
  //minByForce[ pos, R, fr,  1.0*10^-5, 0.003 , 40, False];
  
  *Etry = Enow + dE + Etry1 - Enow1;
  //{pos,  Etry}
  
}

/*
 find active objs, which can be defined as those that collide with
 the specified one
 */
int find_active(int i, double Rcut, struct matrix * L, double orth, int Act[], struct vector vsave[]) {
  int j, nact;
  int same_atom=0;
  struct vector v1,v2;
  
  
  nact=0;
  v1 = makevec( object[i].v[0].x, object[i].v[0].y, object[i].v[0].z );
  for (j=0; j<numat; j++) {
    v2 = makevec( object[j].v[0].x, object[j].v[0].y, object[j].v[0].z );
    // !!!!!!!!!!! ASSUMING only the moved object is now considered active	
    if((j==i) ||( /*distsq_ssrep(L,&v1,&v2,&cell,orth,same_atom)< Rcut*/0)) {
      //printf("%d dist=%f\n",j, distsq_ssrep(L,&v1,&v2,&cell,orth,same_atom));	  
      Act[nact]=j;
	  vsave[nact]=v2;
	  nact++;
    }
  }
  return nact;
}


/*
 Min-map move
 */
int MinMapMove(int ia, struct Energy * Etry, 
	       struct matrix * L, struct matrix * Linv) {
  int i;
  double dx, dy, dz, Rtolocalmin, MMfrac, orth;
  double  rmove[MAX_ACT*3];
  int nAct0, nAct1, Act0[MAX_ACT], Act1[MAX_ACT];
  struct vector v0[MAX_ACT],v1[MAX_ACT], trans_vec, TP;
  extern double MM_transfrac[2];

  orth=(*Etry).orth;
  MMfrac= MM_transfrac[1]- MM_transfrac[0];
  TP = object[ia].v[0];
  nAct0=find_active(ia, 3.9, L, orth, Act0, v0);
  if(DB>0) printf("MM> ia=%5d nAct0=%5d Act0[0]=%5d\n",ia, nAct0,Act0[0]);
  minByForce(nAct0,Act0, 0.03, 0.035 , 10, L,Linv,orth, 0);
  Rtolocalmin = dist_pbc(L, v0, object[ia].v);
  for (i=0; i<nAct0; i++) {
    // undo minimization. ONLY work for atom, NOT group!
    object[Act0[i]].v[0]= v0[i];
  }
  
  dx = (drand48() * MMfrac +MM_transfrac[0]) * rsign();
  dy = (drand48() * MMfrac +MM_transfrac[0]) * rsign();
  dz = (drand48() * MMfrac +MM_transfrac[0]) * rsign();
  if(DB>0) printf("MM> dx=%f dy=%f dz=%f Rtolocalmin=%f\n", dx, dy, dz, Rtolocalmin);
  trans_vec = vsub( rezone( makevec( TP.x+dx, TP.y+dy, TP.z+dz ) ) , TP );
  trans_object( ia, L, trans_vec );
  
  nAct1=find_active(ia, 3.9, L, orth, Act1, v1);
  minByForce(nAct1,Act1, 0.03, 0.035 , 10, L,Linv,orth, 0);
  
  // generate a random move of length Rtolocalmin
  hyperspherepoint(3, Rtolocalmin, rmove);
  trans_vec.x = rmove[0]; 
  trans_vec.y = rmove[1]; 
  trans_vec.z = rmove[2]; 
  trans_vec =cart( Linv, &trans_vec);
  if(DB>0) printf("MM> rmove dx=%f dy=%f dz=%f\n", trans_vec.x, trans_vec.y, trans_vec.z);
  TP = object[ia].v[0];
  trans_vec = vsub( rezone( makevec( TP.x+trans_vec.x, TP.y+trans_vec.y, TP.z+trans_vec.z ) ) , TP );
  trans_object( ia, L, trans_vec );
  
  return 1; // return 0 if min-map failed. note: NOT checked!!!!
}


