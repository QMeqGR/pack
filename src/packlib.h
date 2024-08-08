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

/*

Packlib header.

Eric H Majzoub
Sandia National Laboratories.
Livermore, CA 2006

22 jan 2006
add code for non-orth latts

modified on 21 jan 2006: added support for dimers
-noticed the struct octahedra had v[9], must have
been a mistake.  replaced with v[7].

*/

#ifndef __PACKLIB_H__
#define __PACKLIB_H__

#include "global_defs.h"

/* constants */

#define PI (3.14159265358979323846264338327950288419716939937510)
#define R2D ( 180.0 / PI )

#define SRT2 1.4142135623731
#define SRT3 1.73205080756888
#define SRT6 2.44948974278318


/* Structures */

struct cmplx{
  double real;
  double imag;
};

struct vector{
  double x;
  double y;
  double z;
};

struct matrix{
  double r11,r12,r13;
  double r21,r22,r23;
  double r31,r32,r33;
};

struct intvector{
  int n1,n2,n3;
};

struct Energy{
  double ecc,eng,ecc_eng,bnd,pf,vol,tot,orth;
  double a,b,c,alph,beta,gamm;
};

struct atom{
  int Z;
  int typ;
  int moved; /* for ss,lj sums scaling improvement */

  double x;
  double y;
  double z;
  double w; /* atomic weight */
  double R; /* radius of atom */
  double chrg; /* for Ewald sums */

};

struct hrd_dist{
  double ca1_ca1,ca1_ca2,ca1_ac1,ca1_av1,ca1_ac2,ca1_av2;
  double ca2_ca2,ca2_ac1,ca2_av1,ca2_ac2,ca2_av2;
  double ac1_ac1,ac1_av1,ac1_ac2,ac1_av2;
  double ac2_av1,ac2_ac2,ac2_av2;
  double av1_av1,av1_av2;
  double av2_av2;

  double sen_av1_nv4,sen_av2_nv4;
  double sen_av1_nv5,sen_av2_nv5;
  double sen_av1_nv7,sen_av2_nv7;
};

struct hstdat{
  double vst,flt;
};

struct endos{
  double lngE;  /* this holds log(number of states) */
  double refen; /* this is the reference energy */
};

struct plane{ /* define the vectors of a plane */
 struct vector u,v;
};

struct str_wl_init{
  int nEbins,Emax_bin,Emin_bin;
  double Emin,Emax,E_low,Ebin_width;
  struct endos *DOS,*HST;
};

/* geometrical object structures */

struct basis{
  struct vector A,B,C;
};

struct cellprm{
  struct basis bas;
  double a,b,c;
  double alph,beta,gamm;
};

struct constr{
  int constrained;
  double xmin,xmax,ymin,ymax,zmin,zmax;
};

struct cation{
  int p_num;              /* correspondence to p structure */
  int Z;
  double R;
  double chrg;
  struct vector center;
  struct vector force;
  int reject_move;
  int moved;                  /* for optimized ss,lj energy */
};

struct anion{
  int p_num[MAX_VERTS];    /* correspondence to p structure */
  int num_vert;               /* number of vertices for this anion */
  int Zc,Zv;                  /* Z of center and vertices */
  double Rc,Rv;               /* radius of center and vertices */
  double cvd;                 /* center to vertex distance */
  struct vector v[MAX_VERTS]; /* [0] is the center */
  struct vector or;           /* orientation: phi,theta,psi, byron and fuller */
  struct vector force;
  int reject_move;
  int moved;                  /* for optimized ss,lj energy */
};

/* generalized object, defined initially for more complicated anions */
struct obj{

  int reject_move;
  int typ;                    /* type def for swapping table */
  int used;                   /* is the object used: 1=yes 0=no */
  int num_vert;               /* number of vertices for this anion */
  int num_rept;               /* number of repeats */
  int p_num[MAX_VERTS];       /* correspondence to p structure */
  int vmoved[MAX_VERTS];      /* for optimized ss,lj energy */
  int vused[MAX_VERTS];       /* is the vertex used: 1=yes 0=no */
  int vfxd[MAX_VERTS];        /* is the vertex fixed in position? */
  int allowrot;               /* allow rotations if only v[0] is fixed */
  int swappable;
  int Z[MAX_VERTS];           /* Z of vertices, Z == -1 for unused verts */
  

  double R[MAX_VERTS];        /* radius of vertices */
  double cvd[MAX_VERTS];      /* center to vertex distances */
  double cvdsq[MAX_VERTS];    /* square of center to vertex distances */
  double vchrg[MAX_VERTS];
  
  struct constr vcnstr[MAX_VERTS]; /* constraints vertex atom, 0=unconstrained, 1=constrained  with at_constr tag */
  struct vector v[MAX_VERTS]; /* [0] is the 'center' vertex */
  struct vector center;
  struct vector or;           /* orientation: phi,theta,psi, byron and fuller */
  struct vector force;

};


/* Function prototypes */

void print_restart_fileptr(FILE *,char);
void get_init_parms(FILE *,int exearly);

void amoeba_ehm(double **,double y[],int,double,struct Energy (*funk)(void), int *);
void bin_energy_hst(int,struct endos*);
void check_objects(void);
void invt_matrx(int,double *);
void mark_all_objects_moved(void);
void NxNmult(int N,double *,double *,double *);
void np_force_mv(void);
void trans_all(struct vector);
void print_cell_frame_xbs(struct cellprm,double,int);
void print_object(int);
void print_object_xbs(int,struct obj *,struct cellprm);
void print_hst(struct endos*, int,int,int);
void print_input_out(char *);
void print_dos(struct endos*, int,int);
void print_spe_xbs(char c,double R,double r,double g,double b);
void print_speZ_xbs(int,int,double);
void print_bnd_xbs(char c1,char c2,double,double,double);
void print_moved_flags(void);
void rescale_objects(struct cellprm *,struct vector,struct vector);
void remove_tempfile(void);
void remove_restfile(void);
void rot_object(int,struct matrix *,struct cellprm *,struct matrix);
void simplex_init(void);
void simplex_trans(double *);
void trans_object(int,struct matrix *,struct vector);
void timing_info(void);

void pf_anneal(int,int,double, double,int,int,int,int);
void Simplex(int,int);
void Basinhop(void);
void cell_init(void);
void center_np(void);
void recenter(void);
void best_compare(void);
void simplex_init(void);
void remove_tempfile(void);
void simplex_trans(double *);
void free_memory(void);
void increase_debug(int sig);
void kr_reverse(char*);
void kr_itoa(int, char*);
void estat_init(void);
void energ_init(void);
void exit_pack(int);
void exit_and_print(int sig);
void exit_immediately(int sig);
void print_status(int sig);
void signal_print_status(void);

void print_all_objs(void);
void print_mat( struct matrix );
void print_c2v(void);
void debug_block_metro(void);
void debug_block_WL(void);
void print_powder_out(void);
void print_poscar_out(int);
void print_symsearch_out(void);
void print_findsym_out(void);
void print_xbs(void);
void print_info_block(void);
void print_tmp(void);
void print_restart(void);
void print_restart_inline(void);
void print_obj_bonds_xbs(int);
void print_xbs_mvframe(void);
void get_min_distances(void);
void get_hrd_dist(void);

void Metro(int runs,
	   int t_runs,
	   int trn_chng,
	   int rot_chng,
	   int lat_chng,
	   int swp_chng,
	   double cell_min,
	   double cell_max,
	   double trans_frac,
	   double volum_frac,
	   double temp_init,
	   double temp_finl);

void Wang_Landau(int runs,
		 int t_runs,
		 int trn_chng,
		 int rot_chng,
		 int lat_chng,
		 int swp_chng,
		 double cell_min,
		 double cell_max,
		 double trans_frac,
		 double volum_frac);

void scale_cell(double sf,int updown,int scale_cell);
void force_step(void);
void reset_HST(struct endos *, int);
void wall_check_error(int,int,double);
void ztoelm(char *, int);
void numat_typeZ(void);
void check_constraints(void);
void wl_read_restart(void);

int bin_energy_E(int,int,double,double,double,double,struct endos*,struct endos*);
int calc_forces(void);
int calc_cell_basis(struct cellprm *);
int get_bin_number(double, double, double);
int is_cubic(void);
int lat_pathological(struct cellprm *, struct vector, struct vector,double);
int rsign(void);
int tol(double,double,double);
int main(int argc, char *argv[]);
int make_object(int,struct vector,double,double,double);
int tol(double,double,double);
int getopt(int argc, char * const argv[], const char *optstring);
int wl_accept_move(int,int,double,double,double,double,double,struct endos *, struct endos *);
int trans_accept(int, struct vector,int);
int check_wall_violation(void);
int check_constr_violation(int,struct vector);
int check_all_object_constraints(void);
int getz( char * );
int wttoz( double );
int get_ndummy( void );

int RandomChoice(double * prob, int n);

double amotry_ehm(double **,double y[],double psum[],int,struct Energy (*funk)(),int,double);
double cell_volume(struct cellprm *);
double dmod(double,double);
double det3x3(struct matrix);
double dist(struct vector *,struct vector *);
double dist_pbc(struct matrix *,struct vector *,struct vector *);
double dist_NOpbc(struct matrix *,struct vector *,struct vector *);
double distsq_pbc(struct matrix *,struct vector *,struct vector *);
double distsq_ssrep(struct matrix *,struct vector *,struct vector *,struct cellprm *,double,int);
double Eion(struct cellprm, struct atom *,int,double,double);
double deltaEion(int, int *, struct cellprm, struct atom *, struct atom *,int,double,double);
double Eion_nopbc(struct cellprm, struct atom *,int,double,double);
double find_Ewald_eta(struct basis,double,int,struct atom *);
double modify_lnfmod(double, int);
double vdotprod(struct vector, struct vector);
double vtriple(struct vector, struct vector, struct vector);
double vmag(struct vector);
double vmagsq(struct vector);
double autoadjust_temp_rejpct(double temp_init,
			      double pct_lo,
			      double pct_hi,
			      double cell_min,
			      double cell_max,
			      double trans_frac,
			      double volum_frac,
			      int trn_chng,
			      int rot_chng,
			      int lat_chng,
			      int swp_chng);
double autoadjust_temp_pCnt(double temp_init,
			    double pct_lo,
			    double pct_hi,
			    double cell_min,
			    double cell_max,
			    double trans_frac,
			    double volum_frac,
			    int trn_chng,
			    int rot_chng,
			    int lat_chng,
			    int swp_chng);
double aspect(void);
double molwt(int);
double fmin(double,double);
double pack_frac(void);
double surface_area(void);
double Energy_rep(double orth);
double prep_energy(double (*)(double []));
double ztowt( int );



struct Energy Total_Energy(void);
struct vector Pressure(void);


struct vector cart(struct matrix *, struct vector *);
struct vector center_of_mass(void);
struct vector crossprod(struct vector, struct vector);
struct vector Ctrans(struct matrix *, struct matrix *, struct vector *, struct vector *);
struct vector force_on_atom(int, struct vector);
struct vector force_on_obj(int);
struct vector geometric_center(void);
struct vector get_rand_pos(int);
struct vector makevec(double,double,double);
struct vector multiply(struct matrix *,struct vector *);
struct vector normvec(struct vector);
struct vector rezone(struct vector);
struct vector UTMvmult(struct matrix *, struct vector *);
struct vector vertex_new(struct matrix *,struct matrix *,struct vector *,struct vector *);
struct vector vsub(struct vector, struct vector);
struct vector vadd(struct vector, struct vector);
struct vector vsmult(double, struct vector);
struct vector ztorgb( int );

struct intvector gbox(double, double, struct basis);

struct matrix L_e3(struct cellprm *);
struct matrix Linverse(struct matrix *);
struct matrix matrix_for(struct vector,struct vector);
struct matrix R_ptp(double,double,double);
struct matrix R_ptp_inv(struct matrix Rf);

struct face facegen(struct vector *, struct vector *, struct vector *);

struct hstdat wl_hst_flat(struct endos *,struct endos *,int);
struct str_wl_init wl_init(struct endos *, struct endos *);

/* function prototypes that don't seem to be in math.h */
double fmax(double,double);

/* MPI related */
//int is_master(void);
#define is_master() (process_rank==0)
#endif // __PACKLIB_H__
