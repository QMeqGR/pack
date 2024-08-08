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


#define CAT_1_R 1.0
#define CAT_1_G 0.0
#define CAT_1_B 0.0

#define CAT_2_R 0.5
#define CAT_2_G 0.5
#define CAT_2_B 0.5

#define ANI_1_C_R 0.0
#define ANI_1_C_G 0.0
#define ANI_1_C_B 1.0
#define ANI_1_V_R 0.0
#define ANI_1_V_G 1.0
#define ANI_1_V_B 0.0

#define ANI_2_C_R 0.0
#define ANI_2_C_G 1.0
#define ANI_2_C_B 0.0
#define ANI_2_V_R 0.1
#define ANI_2_V_G 0.1
#define ANI_2_V_B 0.9

#define BND_RAD 0.01
#define FRM_RAD 0.01

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
#include "packlib.h"

extern int process_rank, process_num;

extern char tempfile[25];
extern char restfile[25];

extern int debug,numat;
extern int get_estat,mvout;

extern int non_per_ewald;
extern int wl_prt_dos,wl_prt_hst;

extern double c2vtol;

extern struct cellprm cell;
extern struct atom *p;
extern struct atom *s;
extern struct obj *object;

extern FILE *input;
extern FILE *outfile;
extern FILE *restartfile;
extern FILE *restartinput;
extern FILE *tmpfp;
extern FILE *dosfp;
extern FILE *hstfp;

/*********************************************/
/*                                           */
/*        Printing functions                 */
/*                                           */
/*********************************************/

void print_input_out(char * inputfile)
{
  if (!is_master()) return;

  char line[200] = {0};
  
  /* Open input file and get initialization conditions */
  input = fopen(inputfile,"r");
  if ( input == NULL )
    {
      printf("print_input_out: Error opening input file.\n");
      exit_pack(0);
    }
  printf("*II verbatim input file\n");
  printf("*II ( use:  | awk '($1==\"*I\"){print $0}' | sed 's/*I //' > input.dat )\n");
  rewind(input);
  while ( !feof( input ) ) {
    if ( fgets(line, 200, input) != NULL ) printf("*I %s", line);
  }
  fclose(input);
  return;
}

void print_restart(void)
{

  restartfile = fopen(restfile,"w+");
  if ( restartfile==NULL ){
    printf("ERROR: print_restart: Can't open restart output file %s.\n",restfile);
    exit_pack(0);
  }

  rewind(restartfile);
  /* print out restart parameters */
  print_restart_fileptr(restartfile,' ');
  
  fclose(restartfile);
  return;
}

void print_restart_fileptr(FILE *fp,char K)
{
  char lead,*fxd=" ";

  int i,j;

  if ( K == ' ' ) lead = ' ';
  else lead = '*';

  /* print out restart parameters */
  fprintf(fp,
	  "%c%c %18.12f%18.12f%18.12f\n"
	  "%c%c %18.12f%18.12f%18.12f\n",
	  lead,K,cell.a,cell.b,cell.c,
	  lead,K,cell.alph,cell.beta,cell.gamm);

  /* print out object information */
  for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
    fprintf(fp,
	    "%c%c %5d%5d%5d\n",
	    lead,K,i,object[i].num_vert,1); /* each restart object is unique */
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      if ( object[i].vfxd[j]==1 ) fxd="F";
      else fxd="";
      fprintf(fp,
	      "%c%c %5d%5d%18.12f%18.12f%18.12f%8.4f%5d%10.6f %d\n",
	      lead,K,i,j,
	      object[i].v[j].x,
	      object[i].v[j].y,
	      object[i].v[j].z,
	      object[i].R[j],
	      object[i].Z[j],
	      object[i].vchrg[j],
	      object[i].vfxd[j]);
    }
  }

}

void print_restart_inline(void)
{
  printf("*RR restart output\n");
  printf("*RR (use: | awk '($1==\"*R\"){print $0}' | sed 's/*R //' > restart.dat )\n");

  print_restart_fileptr(stdout,'R');

  return;
}

void print_info_block(void)
{
  if (!is_master()) return;
  struct Energy En;
  struct vector pressure;

  printf("*\n* Lattice parameters:\n*\n");
  printf("* %20s%15.4f\n","lat a  =",cell.a);
  printf("* %20s%15.4f\n","lat b  =",cell.b);
  printf("* %20s%15.4f\n","lat c  =",cell.c);
  printf("* %20s%15.4f\n","alpha  =",cell.alph * R2D);
  printf("* %20s%15.4f\n","beta   =",cell.beta * R2D);
  printf("* %20s%15.4f\n","gamma  =",cell.gamm * R2D);
  printf("* %20s%15.4e\n","cell vol  =",cell_volume(&cell));
  printf("* %20s%15.4f\n","cv/ortho  =",cell_volume(&cell)/(cell.a*cell.b*cell.c));
  
  get_estat=1;
  En = (*Total_Energy)();
  printf("*\n* Energies:\n");
  printf("* %20s%20.10e%4s\n","ecc  =",En.ecc,"Ha");
  printf("* %20s%20.10e%4s\n","ecc  =",En.ecc * 27.211383,"eV");
  printf("* %20s%20.10e\n","eng  =",En.eng);
  printf("*\n");
  printf("*\n* Geometry:\n");
  printf("*\n");
  printf("* %20s%15.4f\n","packing frac  =", En.pf );
  printf("* %20s%15.4f\n","aspect ratio  =", aspect() );
  pressure = Pressure();
  printf("*\n* (dE/dV):\n");
  printf("* %20s%15.4e%25s\n","Px  =",pressure.x,"x-dir");
  printf("* %20s%15.4e%25s\n","Py  =",pressure.y,"y-dir");
  printf("* %20s%15.4e%25s\n","Pz  =",pressure.z,"z-dir");

  /* find and print the minimum distances */
  get_min_distances( );
  
  /* print out the center to vertex distances */
  if ( debug > 1 ) print_c2v( );
  printf("*\n");

  return;
}

void print_xbs(void)
{
  int i;

  calc_cell_basis(&cell);

  printf("* XBS output\n");

  /* objects */
  for (i=0; i<MAX_OBJS && object[i].used==1 ; i++) {
    print_object_xbs(i,object,cell);
  }

  /* print the cell frame */
  print_cell_frame_xbs(cell,FRM_RAD,mvout);
  
  return;
}


void print_xbs_mvframe(void)
{
  int i,j;

  struct vector c0,c1,c2,c3,c4,c5,c6,c7;
  struct vector tmp;
  struct matrix L;

  L = L_e3(&cell);

  fprintf(outfile,"frame\n");

  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      tmp = cart(&L,&object[i].v[j]);
      fprintf(outfile,"%12.4f%12.4f%12.4f\n",tmp.x,tmp.y,tmp.z);
    }
  }

  c0 = makevec(0.0,0.0,0.0);
  c1 = cell.bas.A;
  c2 = vadd(cell.bas.A,cell.bas.B);
  c3 = cell.bas.B;
  c4 = cell.bas.C;
  c5 = vadd(cell.bas.A,cell.bas.C);
  c6 = vadd(cell.bas.A, vadd(cell.bas.B,cell.bas.C) );
  c7 = vadd(cell.bas.B, cell.bas.C);
  
  /* print the cell frame corners */
  /* must match print_cell_frame_xbs order */
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c0.x,c0.y,c0.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c1.x,c1.y,c1.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c2.x,c2.y,c2.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c3.x,c3.y,c3.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c4.x,c4.y,c4.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c5.x,c5.y,c5.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c6.x,c6.y,c6.z);
  fprintf(outfile,"%12.4f%12.4f%12.4f\n",c7.x,c7.y,c7.z);
  
  return;
}


void print_c2v(void)
{
  int i,j;
  int violation_flag=0;
  int vnum=0;
  double cvd=0;
  double dist= -1;
  double vdist=0;
  struct matrix L;

  printf("* Center to vertex distances\n");
  printf("* cell   = %15.10f %15.10f %15.10f\n",cell.a,cell.b,cell.c);
  printf("* angles = %15.10f %15.10f %15.10f\n",
	 R2D*cell.alph,
	 R2D*cell.beta,
	 R2D*cell.gamm);

  L = L_e3(&cell);

  if ( debug > 3 ){
    printf("* L1_ %15.10f%15.10f%15.10f\n",L.r11,L.r12,L.r13);
    printf("* L2_ %15.10f%15.10f%15.10f\n",L.r21,L.r22,L.r23);
    printf("* L3_ %15.10f%15.10f%15.10f\n",L.r31,L.r32,L.r33);
  }

  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      dist = dist_NOpbc( &L, &object[i].v[0], &object[i].v[j] );
      printf("* [%d %d] : %15.9f\n",i,j,dist);
      cvd = object[i].cvd[j];
      if ( !tol(cvd,dist,c2vtol) ) {
	violation_flag=1; vdist=dist; vnum=i;
	break;
      }
    }
  }

  if ( violation_flag ) {
    printf("!! Vertex to Center distance violation in object %d.\n",vnum);
    printf("!! Expecting %20.10f, found %20.10f (tol = %10.4e)\n",cvd,vdist,c2vtol);
    printf("!! Exiting.\n");
    exit_pack(0);
  }

}

void print_all_objs(void)
{
  int i;
  for (i=0; i<MAX_OBJS && object[i].used==1 ; i++) print_object(i);
}


void print_mat( struct matrix R )
{
  printf("%15.8f%15.8f%15.8f\n",R.r11,R.r12,R.r13);
  printf("%15.8f%15.8f%15.8f\n",R.r21,R.r22,R.r23);
  printf("%15.8f%15.8f%15.8f\n",R.r31,R.r32,R.r33);
  printf("\n");
  return;
}

void print_tmp(void)
{

  tmpfp = fopen(tempfile,"w+");
  if ( tmpfp==NULL ){
    printf("* Unable to open temporary file.\n");
    exit_pack(0);
  }

  rewind(tmpfp);
  /* print out best parameters */
  print_restart_fileptr(tmpfp,' ');

  fclose(tmpfp);
  return;
}

void print_powder_out(void)
{
  int i,j;
  int ndum=0;
  
  ndum = get_ndummy();

  printf("*DD powder output\n");
  printf("*DD (use: | awk '($1==\"*D\"){print $0}' | sed 's/*D //' > powder.in )\n");
  printf("*D %d\n",numat-ndum);
  printf("*D %12.5f%12.5f%12.5f\n",cell.a,cell.b,cell.c);
  printf("*D %12.5f%12.5f%12.5f\n",cell.alph*R2D,cell.beta*R2D,cell.gamm*R2D);
  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      if ( !tol(object[i].R[j],0,1e-6) ){
	printf("*D %12.5f%12.5f%12.5f%10d\n",
	       object[i].v[j].x,
	       object[i].v[j].y,
	       object[i].v[j].z,
	       object[i].Z[j]);
      }
    }
  }

  return;
}

void numat_typeZ(void)
{
  int i,j,Z;

  static int been_called=0;
  
  extern int *numZ;
  extern struct obj *object;

  if ( been_called == 0 ){
    /* initialize */
    for(i=0;i<ZMAX;i++){
      numZ[i]=0;
    }

    for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
      for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
	for ( Z=1; Z<ZMAX; Z++ ){
	  if ( Z==object[i].Z[j] && !tol(object[i].R[j],0,1e-6) ) (*(numZ+Z))++;
	}
      }
    }

  been_called=1;
  }

  return;

}

void print_poscar_out(int withQ)
{
  char charout;

  int i,j,Z;
  extern int *numZ;
  struct matrix L;

  L = L_e3(&cell);

  if ( withQ ) charout='Q';
  else charout='P';

  numat_typeZ();

  printf("*%c%c poscar output\n",charout,charout);
  printf("*%c%c (use: | awk '($1==\"*%c\"){print $0}' | sed 's/*%c //' > POSCAR )\n",charout,charout,charout,charout);
  if ( withQ==0 ){
    printf("*%c Z: ",charout);
    for(i=0;i<ZMAX;i++){
      if ( numZ[i]>0 ) printf("%d ",i);
    }
    printf("\n");
  }
  if ( withQ==1 ){
    printf("*%c ",charout);
    printf(" not available in pack v5+\n");
  }  
  printf("*%c 1.0 \n",charout);
  printf("*%c %15.8f%15.8f%15.8f\n",charout,cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
  printf("*%c %15.8f%15.8f%15.8f\n",charout,cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
  printf("*%c %15.8f%15.8f%15.8f\n",charout,cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
  
  printf("*%c ",charout);
  for(Z=0;Z<ZMAX;Z++){
    if ( numZ[Z]>0 ) printf("%d ",numZ[Z]);
  }
  printf("\n");

  printf("*%c Direct\n",charout);

  for(Z=0;Z<ZMAX;Z++){
    for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
      for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
	if ( Z==object[i].Z[j] && !tol(object[i].R[j],0,1e-6) ){
	  printf("*%c %20.12f%20.12f%20.12f\n",charout,
		 object[i].v[j].x,
		 object[i].v[j].y,
		 object[i].v[j].z);
	}
      }
    }
  }
  
  return;
}

void print_symsearch_out(void)
{
  int i,j,ndum;

  struct matrix L;

  ndum = get_ndummy();
  L = L_e3(&cell);

  printf("*SS symsearch output\n");
  printf("*SS ( use:  | awk '($1==\"*S\"){print $0}' | sed 's/*S //' > symsearch.in )\n");
  printf("*S Title\n");
  printf("*S 1                 search type\n");
  printf("*S 5.0 0.00001       displacement rms max min\n");
  printf("*S 0.1 0.00001       strain rms max min\n");
  printf("*S 1                 lattice input type\n");
  printf("*S %12.5f%12.5f%12.5f\n",cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
  printf("*S %12.5f%12.5f%12.5f\n",cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
  printf("*S %12.5f%12.5f%12.5f\n",cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
  printf("*S 1 0 0\n");
  printf("*S 0 1 0\n");
  printf("*S 0 0 1\n");
  printf("*S %d                 number of atoms\n",numat-ndum);
  printf("*S");
  for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      if ( !tol(object[i].R[j],0,1e-6) ) printf(" %d",object[i].Z[j]);
    }
  }
  printf("\n");

  for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      if ( !tol(object[i].R[j],0,1e-6) ) printf("*S %12.5f%12.5f%12.5f\n",
						object[i].v[j].x,
						object[i].v[j].y,
						object[i].v[j].z);
    }
  }

  return;
}

void print_findsym_out(void)
{
  int i,j,ndum;

  struct matrix L;

  ndum = get_ndummy();
  L = L_e3(&cell);

  printf("*FF findsym output\n");
  printf("*FF ( use:  | awk '($1==\"*F\"){print $0}' | sed 's/*F //' > findsym.in )\n");
  printf("*F Title\n");
  printf("*F 0.001       displacement rms max min\n");
  printf("*F 1                 lattice input type\n");
  printf("*F %12.5f%12.5f%12.5f\n",cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
  printf("*F %12.5f%12.5f%12.5f\n",cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
  printf("*F %12.5f%12.5f%12.5f\n",cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
  printf("*F 1               basis vector in type\n");
  printf("*F 1 0 0\n");
  printf("*F 0 1 0\n");
  printf("*F 0 0 1\n");
  printf("*F %d                 number of atoms\n",numat);
  printf("*F");
  for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      if ( !tol(object[i].R[j],0,1e-6) ) printf(" %d",object[i].Z[j]);
    }
  }
  printf("\n");

  for( i=0; i<MAX_OBJS && object[i].used==1; i++ ){
    for( j=0; j<MAX_VERTS && object[i].vused[j]==1; j++ ){
      if ( !tol(object[i].R[j],0,1e-6) ) printf("*S %12.5f%12.5f%12.5f\n",
						object[i].v[j].x,
						object[i].v[j].y,
						object[i].v[j].z);
    }
  }

  return;
}


/**************************************/
/*                                    */
/*    print cell frame xbs            */
/*                                    */
/**************************************/
void print_cell_frame_xbs(struct cellprm cell,double frm_rad,int mvout)
{

  double HI;
  double LO;
  double CORNER_RAD=0.05;

  struct vector c0,c1,c2,c3,c4,c5,c6,c7;

  c0 = makevec(0.0,0.0,0.0);
  c1 = cell.bas.A;
  c2 = vadd(cell.bas.A,cell.bas.B);
  c3 = cell.bas.B;
  c4 = cell.bas.C;
  c5 = vadd(cell.bas.A,cell.bas.C);
  c6 = vadd(cell.bas.A, vadd(cell.bas.B,cell.bas.C) );
  c7 = vadd(cell.bas.B, cell.bas.C);

  if ( mvout ) {
    HI = 100.0;
    LO = 0.00001;
  } else {
    HI = 1.001;
    LO = 0.999;
  }

  printf("* CELL PARAMETERS AND BORDER\n");
  printf("atom 0 %10.5f%10.5f%10.5f\n",c0.x,c0.y,c0.z);
  printf("atom 1 %10.5f%10.5f%10.5f\n",c1.x,c1.y,c1.z);
  printf("atom 2 %10.5f%10.5f%10.5f\n",c2.x,c2.y,c2.z);
  printf("atom 3 %10.5f%10.5f%10.5f\n",c3.x,c3.y,c3.z);

  printf("atom 4 %10.5f%10.5f%10.5f\n",c4.x,c4.y,c4.z);
  printf("atom 5 %10.5f%10.5f%10.5f\n",c5.x,c5.y,c5.z);
  printf("atom 6 %10.5f%10.5f%10.5f\n",c6.x,c6.y,c6.z);
  printf("atom 7 %10.5f%10.5f%10.5f\n",c7.x,c7.y,c7.z);

  print_spe_xbs('0',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('1',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('2',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('3',CORNER_RAD,0.0,0.0,0.0);

  print_spe_xbs('4',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('5',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('6',CORNER_RAD,0.0,0.0,0.0);
  print_spe_xbs('7',CORNER_RAD,0.0,0.0,0.0);

  print_bnd_xbs('0','1',cell.a*LO,cell.a*HI,frm_rad);
  print_bnd_xbs('3','2',cell.a*LO,cell.a*HI,frm_rad);
  print_bnd_xbs('4','5',cell.a*LO,cell.a*HI,frm_rad);
  print_bnd_xbs('6','7',cell.a*LO,cell.a*HI,frm_rad);

  print_bnd_xbs('1','2',cell.b*LO,cell.b*HI,frm_rad);
  print_bnd_xbs('0','3',cell.b*LO,cell.b*HI,frm_rad);
  print_bnd_xbs('5','6',cell.b*LO,cell.b*HI,frm_rad);
  print_bnd_xbs('4','7',cell.b*LO,cell.b*HI,frm_rad);

  print_bnd_xbs('0','4',cell.c*LO,cell.c*HI,frm_rad);
  print_bnd_xbs('1','5',cell.c*LO,cell.c*HI,frm_rad);
  print_bnd_xbs('2','6',cell.c*LO,cell.c*HI,frm_rad);
  print_bnd_xbs('3','7',cell.c*LO,cell.c*HI,frm_rad);

  printf("inc %f\n",5.0);
  if ( !mvout && non_per_ewald==0 ) {
    printf("dup %10.5f%10.5f%10.5f\n",cell.bas.A.x,cell.bas.A.y,cell.bas.A.z);
    printf("dup %10.5f%10.5f%10.5f\n",cell.bas.B.x,cell.bas.B.y,cell.bas.B.z);
    printf("dup %10.5f%10.5f%10.5f\n",cell.bas.C.x,cell.bas.C.y,cell.bas.C.z);
  }

  return;
}


/**************************************/
/*                                    */
/*    print object                    */
/*                                    */
/**************************************/
void print_object(int N)
{
  int i;
  struct vector v={0};

  printf("# OBJECT:\n");
  for (i=0; i<object[N].num_vert && object[N].vused[i]==1; i++) {
    v = object[N].v[i];
    printf("# [%d] %15.9lf%15.9lf%15.9lf\n",N,v.x,v.y,v.z);
  }
  return;
}

/**************************************/
/*                                    */
/*    print object xbs                */
/*                                    */
/**************************************/
void print_object_xbs(
		      int obj_num,
		      struct obj *object,
		      struct cellprm cell
		     )
{

  char c[2];
  int i,nv=0,Z,p;
  double R;

  struct vector tmp;
  struct matrix L;
  
  L = L_e3(&cell);

  nv = object[obj_num].num_vert;
  printf("* OBJECT: p_num (center) = %2d\n",object[obj_num].p_num[0]);
  /* the anion vertices */
  for (i=0; i<nv; i++) {

    Z = object[obj_num].Z[i];
    R = object[obj_num].R[i];
    p = object[obj_num].p_num[i];

    tmp = cart(&L,&object[obj_num].v[i]);
    ztoelm( c, Z );
    printf("atom %2s%d %15.9lf%15.9lf%15.9lf\n",
	   c,p,tmp.x,tmp.y,tmp.z);
    print_speZ_xbs( Z, p, R );

  }

  print_obj_bonds_xbs(obj_num);
  return;
}

/**************************************/
/*                                    */
/*  print Z specs for xbs             */
/*                                    */
/**************************************/
void print_speZ_xbs(int Z,int p,double R)
{
  char elm[2];
  struct vector rgb;

  ztoelm( elm, Z );
  rgb = ztorgb( Z );
  
  printf("spec %2s%d %15.9lf%15.9lf%15.9lf%15.9lf\n",elm,p,R,rgb.x,rgb.y,rgb.z);

  return;
}


/**************************************/
/*                                    */
/*  print species specs for xbs       */
/*                                    */
/**************************************/
void print_spe_xbs(char c,double R,double r,double g,double b)
{
  
  printf("spec %c %15.9lf%15.9lf%15.9lf%15.9lf\n",c,R,r,g,b);

  return;
}


/**************************************/
/*                                    */
/*    print obj bonds xbs             */
/*                                    */
/**************************************/
void print_obj_bonds_xbs(int n)
{
  int j,z0,zj,p0,pj;
  char e0[2],ej[2];
  double dist,small;

  z0 = object[n].Z[0];
  p0 = object[n].p_num[0];
  ztoelm(e0,z0);
  for(j=1; j<MAX_VERTS && object[n].vused[j]==1; j++){
    zj = object[n].Z[j];
    pj = object[n].p_num[j];
    ztoelm(ej,zj);
    dist = object[n].cvd[j];
    small = dist/100;
    /* Note: all bonds will be black here */
    printf("bonds %s%d     %s%d     %12.9lf%15.9lf%15.9lf    %s\n",
	   e0,p0,ej,pj,dist-small,dist+small,0.1,"Black");
  }

  return;
}


/**************************************/
/*                                    */
/*    print bond specs for xbs        */
/*                                    */
/**************************************/
void print_bnd_xbs(char c1,char c2,double m,double M,double R)
{
  /* Note: all bonds will be black here */
  printf("bonds %c %c %12.9lf%15.9lf%15.9lf    %s\n",c1,c2,m,M,R,"Black");

  return;
}

/**************************************/
/*                                    */
/*    print density of states         */
/*                                    */
/**************************************/
void print_dos(struct endos *DOS, int maxn, int L)
{
  int i;
  static int been_here=0;

  double loglogGE=0;

  pid_t processid;
  char pid_string[15];
  char dosfile[25]="DOS_",dos_string[4]=".dat";

  if ( !wl_prt_dos ) return;

  processid = getpid();
  if ( debug>1 && been_here==0 ) printf("* processid=%d\n",processid);
  been_here=1;
  snprintf(pid_string, 15, "%d", processid );
  strncat(dosfile, pid_string, 15);
  strncat(dosfile, dos_string, 4);
  dosfp = fopen(dosfile,"w");
  if ( dosfp==NULL ){
    printf("Can't open DOS output file.\n");
    exit_pack(0);
  }

  for(i=1;i<maxn;i++){
    if ( DOS[i].lngE == 1 ) continue; /* don't print empty states */
    if ( DOS[i].lngE < 1 ) loglogGE= -1;
    else loglogGE = log(DOS[i].lngE);
    if ( L == 0 ) fprintf(dosfp,"%d %30.20e\n",i,DOS[i].lngE);
    if ( L == 1 ) fprintf(dosfp,"%30.20e%30.20e%30.20e\n",
			  DOS[i].refen,
			  DOS[i].lngE,
			  loglogGE);
  }

  fclose(dosfp);
  return;
}

/**************************************/
/*                                    */
/*    print WL histogram              */
/*                                    */
/**************************************/
void print_hst(struct endos *HST, int maxn, int L, int hst_num)
{
  int i;
  static int been_here=0;

  double loglogHE=0;

  pid_t processid;
  char pid_string[15];
  char hstfile[25]="HST_",hst_underscore[1]="_",hst_string[4]=".dat";
  char Nhst[3]={0};

  if ( !wl_prt_hst ) return;

  processid = getpid();
  if ( debug>1 && been_here==0 ) printf("* processid=%d\n",processid);
  been_here=1;
  snprintf(pid_string, 15, "%d", processid );
  strncat(hstfile, pid_string, 15);
  strncat(hstfile, hst_underscore, 1);
  kr_itoa(hst_num,Nhst);
  strncat(hstfile, Nhst, 3 );
  strncat(hstfile, hst_string, 4);
  hstfp = fopen(hstfile,"w");
  if ( hstfp==NULL ){
    printf("Can't open HST output file.\n");
    exit_pack(0);
  }


  for(i=1;i<maxn;i++){
    if ( HST[i].lngE == 0 ) continue; /* don't print empty states */
    if ( HST[i].lngE < 1 ) loglogHE= -1;
    else loglogHE = log(HST[i].lngE);
    if ( L == 0 ) fprintf(hstfp,"%d %30.20e\n",i,HST[i].lngE);
    if ( L == 1 ) fprintf(hstfp,"%30.20e%30.20e%30.20e\n",
			  HST[i].refen,
			  HST[i].lngE,
			  loglogHE);
  }

  fclose(hstfp);
  return;
}

void print_moved_flags()
{
  int i,j;

  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0; j<MAX_VERTS && object[i].vused[j]==1; j++){
      if ( object[i].vmoved[j]!=0 ) printf("object[%2d %2d] is flagged.\n",i,j);
    }
  }

  return;
}

int get_ndummy( void )
{

  int i,j;
  int ndum;

  for(i=0; i<MAX_OBJS && object[i].used==1; i++){
    for(j=0;j<MAX_VERTS && object[i].vused[j]==1; j++){
      if ( tol(object[i].R[j], 0, 1e-6) ) ndum++; 
    }
  }

  return(ndum);
}
