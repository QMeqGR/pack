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

extern int debug;


int getz( char *ELMT )
{
  if ( 0==strcmp(ELMT,"H") )  return   1;
  if ( 0==strcmp(ELMT,"He") ) return   2;
  if ( 0==strcmp(ELMT,"Li") ) return   3;
  if ( 0==strcmp(ELMT,"Be") ) return   4;
  if ( 0==strcmp(ELMT,"B") )  return   5;
  if ( 0==strcmp(ELMT,"C") )  return   6;
  if ( 0==strcmp(ELMT,"N") )  return   7;
  if ( 0==strcmp(ELMT,"O") )  return   8;
  if ( 0==strcmp(ELMT,"F") )  return   9;
  if ( 0==strcmp(ELMT,"Ne") ) return  10;
  if ( 0==strcmp(ELMT,"Na") ) return  11;
  if ( 0==strcmp(ELMT,"Mg") ) return  12;
  if ( 0==strcmp(ELMT,"Al") ) return  13;
  if ( 0==strcmp(ELMT,"Si") ) return  14;
  if ( 0==strcmp(ELMT,"P") )  return  15;
  if ( 0==strcmp(ELMT,"S") )  return  16;
  if ( 0==strcmp(ELMT,"Cl") ) return  17;
  if ( 0==strcmp(ELMT,"Ar") ) return  18;
  if ( 0==strcmp(ELMT,"K") )  return  19;
  if ( 0==strcmp(ELMT,"Ca") ) return  20;
  if ( 0==strcmp(ELMT,"Sc") ) return  21;
  if ( 0==strcmp(ELMT,"Ti") ) return  22;
  if ( 0==strcmp(ELMT,"V") )  return  23;
  if ( 0==strcmp(ELMT,"Cr") ) return  24;
  if ( 0==strcmp(ELMT,"Mn") ) return  25;
  if ( 0==strcmp(ELMT,"Fe") ) return  26;
  if ( 0==strcmp(ELMT,"Co") ) return  27;
  if ( 0==strcmp(ELMT,"Ni") ) return  28;
  if ( 0==strcmp(ELMT,"Cu") ) return  29;
  if ( 0==strcmp(ELMT,"Zn") ) return  30;
  if ( 0==strcmp(ELMT,"Ga") ) return  31;
  if ( 0==strcmp(ELMT,"Ge") ) return  32;
  if ( 0==strcmp(ELMT,"As") ) return  33;
  if ( 0==strcmp(ELMT,"Se") ) return  34;
  if ( 0==strcmp(ELMT,"Br") ) return  35;
  if ( 0==strcmp(ELMT,"Kr") ) return  36;
  if ( 0==strcmp(ELMT,"Rb") ) return  37;
  if ( 0==strcmp(ELMT,"Sr") ) return  38;
  if ( 0==strcmp(ELMT,"Y") )  return  39;
  if ( 0==strcmp(ELMT,"Zr") ) return  40;
  if ( 0==strcmp(ELMT,"Nb") ) return  41;
  if ( 0==strcmp(ELMT,"Mo") ) return  42;
  if ( 0==strcmp(ELMT,"Tc") ) return  43;
  if ( 0==strcmp(ELMT,"Ru") ) return  44;
  if ( 0==strcmp(ELMT,"Rh") ) return  45;
  if ( 0==strcmp(ELMT,"Pd") ) return  46;
  if ( 0==strcmp(ELMT,"Ag") ) return  47;
  if ( 0==strcmp(ELMT,"Cd") ) return  48;
  if ( 0==strcmp(ELMT,"In") ) return  49;
  if ( 0==strcmp(ELMT,"Sn") ) return  50;
  if ( 0==strcmp(ELMT,"Sb") ) return  51;
  if ( 0==strcmp(ELMT,"Te") ) return  52;
  if ( 0==strcmp(ELMT,"I") )  return  53;
  if ( 0==strcmp(ELMT,"Xe") ) return  54;
  if ( 0==strcmp(ELMT,"Cs") ) return  55;
  if ( 0==strcmp(ELMT,"Ba") ) return  56;
  if ( 0==strcmp(ELMT,"La") ) return  57;
  if ( 0==strcmp(ELMT,"Ce") ) return  58;
  if ( 0==strcmp(ELMT,"Pr") ) return  59;
  if ( 0==strcmp(ELMT,"Nd") ) return  60;
  if ( 0==strcmp(ELMT,"Pm") ) return  61;
  if ( 0==strcmp(ELMT,"Sm") ) return  62;
  if ( 0==strcmp(ELMT,"Eu") ) return  63;
  if ( 0==strcmp(ELMT,"Gd") ) return  64;
  if ( 0==strcmp(ELMT,"Tb") ) return  65;
  if ( 0==strcmp(ELMT,"Dy") ) return  66;
  if ( 0==strcmp(ELMT,"Ho") ) return  67;
  if ( 0==strcmp(ELMT,"Er") ) return  68;
  if ( 0==strcmp(ELMT,"Tm") ) return  69;
  if ( 0==strcmp(ELMT,"Yb") ) return  70;
  if ( 0==strcmp(ELMT,"Lu") ) return  71;
  if ( 0==strcmp(ELMT,"Hf") ) return  72;
  if ( 0==strcmp(ELMT,"Ta") ) return  73;
  if ( 0==strcmp(ELMT,"W") )  return  74;
  if ( 0==strcmp(ELMT,"Re") ) return  75;
  if ( 0==strcmp(ELMT,"Os") ) return  76;
  if ( 0==strcmp(ELMT,"Ir") ) return  77;
  if ( 0==strcmp(ELMT,"Pt") ) return  78;
  if ( 0==strcmp(ELMT,"Au") ) return  79;
  if ( 0==strcmp(ELMT,"Hg") ) return  80;
  if ( 0==strcmp(ELMT,"Tl") ) return  81;
  if ( 0==strcmp(ELMT,"Pb") ) return  82;
  if ( 0==strcmp(ELMT,"Bi") ) return  83;
  if ( 0==strcmp(ELMT,"Po") ) return  84;
  if ( 0==strcmp(ELMT,"At") ) return  85;
  if ( 0==strcmp(ELMT,"Rn") ) return  86;
  if ( 0==strcmp(ELMT,"Fr") ) return  87;
  if ( 0==strcmp(ELMT,"Ra") ) return  88;
  if ( 0==strcmp(ELMT,"Ac") ) return  89;
  if ( 0==strcmp(ELMT,"Th") ) return  90;
  if ( 0==strcmp(ELMT,"Pa") ) return  91;
  if ( 0==strcmp(ELMT,"U") )  return  92;
  if ( 0==strcmp(ELMT,"Np") ) return  93;
  if ( 0==strcmp(ELMT,"Pu") ) return  94;
  if ( 0==strcmp(ELMT,"Am") ) return  95;
  if ( 0==strcmp(ELMT,"Cm") ) return  96;
  if ( 0==strcmp(ELMT,"Bk") ) return  97;
  if ( 0==strcmp(ELMT,"Cf") ) return  98;
  if ( 0==strcmp(ELMT,"Es") ) return  99;
  if ( 0==strcmp(ELMT,"Fm") ) return 100;
  if ( 0==strcmp(ELMT,"Md") ) return 101;
  if ( 0==strcmp(ELMT,"No") ) return 102;
  if ( 0==strcmp(ELMT,"Lr") ) return 103;
  if ( 0==strcmp(ELMT,"Rf") ) return 104;
  if ( 0==strcmp(ELMT,"Db") ) return 105;
  if ( 0==strcmp(ELMT,"Sg") ) return 106;
  if ( 0==strcmp(ELMT,"Bh") ) return 107;
  if ( 0==strcmp(ELMT,"Hs") ) return 108;
  if ( 0==strcmp(ELMT,"Mt") ) return 109;
  return(0);
}

double ztowt( int z )
{
  if ( z==1   ) return( 1.01  );
  if ( z==2   ) return( 4.00 );
  if ( z==3   ) return( 6.94 );
  if ( z==4   ) return( 9.01 );
  if ( z==5   ) return( 10.81  );
  if ( z==6   ) return( 12.01  );
  if ( z==7   ) return( 14.01  );
  if ( z==8   ) return( 16.00  );
  if ( z==9   ) return( 19.00  );
  if ( z==10  ) return( 20.18 );
  if ( z==11  ) return( 22.99 );
  if ( z==12  ) return( 24.30 );
  if ( z==13  ) return( 26.98 );
  if ( z==14  ) return( 28.09 );
  if ( z==15  ) return( 30.97  );
  if ( z==16  ) return( 32.06  );
  if ( z==17  ) return( 35.45 );
  if ( z==18  ) return( 39.95 );
  if ( z==19  ) return( 39.09  );
  if ( z==20  ) return( 40.08 );
  if ( z==21  ) return( 44.96 );
  if ( z==22  ) return( 47.90 );
  if ( z==23  ) return( 50.94  );
  if ( z==24  ) return( 52.00 );
  if ( z==25  ) return( 54.94 );
  if ( z==26  ) return( 55.85 );
  if ( z==27  ) return( 58.93 );
  if ( z==28  ) return( 58.71 );
  if ( z==29  ) return( 63.55 );
  if ( z==30  ) return( 65.38 );
  if ( z==31  ) return( 69.72 );
  if ( z==32  ) return( 72.59 );
  if ( z==33  ) return( 74.92 );
  if ( z==34  ) return( 78.96 );
  if ( z==35  ) return( 79.91 );
  if ( z==36  ) return( 83.80 );
  if ( z==37  ) return( 85.47 );
  if ( z==38  ) return( 87.62 );
  if ( z==39  ) return( 88.91  );
  if ( z==40  ) return( 91.22 );
  if ( z==41  ) return( 92.91 );
  if ( z==42  ) return( 95.94 );
  if ( z==43  ) return( 98.91 );
  if ( z==44  ) return( 101.07 );
  if ( z==45  ) return( 102.90 );
  if ( z==46  ) return( 106.40 );
  if ( z==47  ) return( 107.87 );
  if ( z==48  ) return( 112.40 );
  if ( z==49  ) return( 114.82 );
  if ( z==50  ) return( 118.69 );
  if ( z==51  ) return( 121.75 );
  if ( z==52  ) return( 127.60 );
  if ( z==53  ) return( 126.90  );
  if ( z==54  ) return( 131.30 );
  if ( z==55  ) return( 132.91 );
  if ( z==56  ) return( 137.34 );
  if ( z==57  ) return( 138.91 );
  if ( z==58  ) return( 140.12 );
  if ( z==59  ) return( 140.91 );
  if ( z==60  ) return( 144.24 );
  if ( z==61  ) return( 145.00 );
  if ( z==62  ) return( 150.35 );
  if ( z==63  ) return( 151.96 );
  if ( z==64  ) return( 157.25 );
  if ( z==65  ) return( 158.92 );
  if ( z==66  ) return( 162.50 );
  if ( z==67  ) return( 164.93 );
  if ( z==68  ) return( 167.26 );
  if ( z==69  ) return( 168.93 );
  if ( z==70  ) return( 173.04 );
  if ( z==71  ) return( 174.97 );
  if ( z==72  ) return( 178.49 );
  if ( z==73  ) return( 180.95 );
  if ( z==74  ) return( 183.85  );
  if ( z==75  ) return( 186.20 );
  if ( z==76  ) return( 190.20 );
  if ( z==77  ) return( 192.22 );
  if ( z==78  ) return( 195.09 );
  if ( z==79  ) return( 196.97 );
  if ( z==80  ) return( 200.59 );
  if ( z==81  ) return( 204.37 );
  if ( z==82  ) return( 207.20 );
  if ( z==83  ) return( 208.98 );
  if ( z==84  ) return( 210.00 );
  if ( z==85  ) return( 210.00 );
  if ( z==86  ) return( 222.00 );
  if ( z==87  ) return( 223.00 );
  if ( z==88  ) return( 226.02 );
  if ( z==89  ) return( 227.00 );
  if ( z==90  ) return( 232.04 );
  if ( z==91  ) return( 231.00 );
  if ( z==92  ) return( 238.03 );
  if ( z==93  ) return( 237.00 );
  if ( z==94  ) return( 244.00 );
  if ( z==95  ) return( 243.00 );

  return(0);
}

struct vector ztorgb( int z )
{
  double r,g,b,Z;
  struct vector rgb;

  Z = (double)z;

  if ( z == 1 ) {
    r = 0;
    g = 0.5;
    b = 0.9;
  } else if ( z%2==1 ) {
    r = 1;
    g = Z/120;
    b = 1-g;
  } else {
    r = Z/120;
    b = 1;
    g = 1-b;
  }
  
  rgb.x = r;
  rgb.y = g;
  rgb.z = b;
  
  return(rgb);

}


int wttoz( double wt ) {
  if ( tol(wt,1.01,0.1)==1 ) return( 1 );
  if ( tol(wt,4.00,0.1)==1 ) return( 2 );
  if ( tol(wt,6.94,0.1)==1 ) return( 3 );
  if ( tol(wt,9.01,0.1)==1 ) return( 4 );
  if ( tol(wt,10.81,0.1)==1 ) return( 5 );
  if ( tol(wt,12.01,0.1)==1 ) return( 6 );
  if ( tol(wt,14.01,0.1)==1 ) return( 7 );
  if ( tol(wt,16.00,0.1)==1 ) return( 8 );
  if ( tol(wt,19.00,0.1)==1 ) return( 9 );
  if ( tol(wt,20.18,0.1)==1 ) return( 10 );
  if ( tol(wt,22.99,0.1)==1 ) return( 11 );
  if ( tol(wt,24.30,0.1)==1 ) return( 12 );
  if ( tol(wt,26.98,0.1)==1 ) return( 13 );
  if ( tol(wt,28.09,0.1)==1 ) return( 14 );
  if ( tol(wt,30.97,0.1)==1 ) return( 15 );
  if ( tol(wt,32.06,0.1)==1 ) return( 16 );
  if ( tol(wt,35.45,0.1)==1 ) return( 17 );
  if ( tol(wt,39.95,0.1)==1 ) return( 18 );
  if ( tol(wt,39.09,0.1)==1 ) return( 19 );
  if ( tol(wt,40.08,0.1)==1 ) return( 20 );
  if ( tol(wt,44.96,0.1)==1 ) return( 21 );
  if ( tol(wt,47.90,0.1)==1 ) return( 22 );
  if ( tol(wt,50.94,0.1)==1 ) return( 23 );
  if ( tol(wt,52.00,0.1)==1 ) return( 24 );
  if ( tol(wt,54.94,0.1)==1 ) return( 25 );
  if ( tol(wt,55.85,0.1)==1 ) return( 26 );
  if ( tol(wt,58.93,0.1)==1 ) return( 27 );
  if ( tol(wt,58.71,0.1)==1 ) return( 28 );
  if ( tol(wt,63.55,0.1)==1 ) return( 29 );
  if ( tol(wt,65.38,0.1)==1 ) return( 30 );
  if ( tol(wt,69.72,0.1)==1 ) return( 31 );
  if ( tol(wt,72.59,0.1)==1 ) return( 32 );
  if ( tol(wt,74.92,0.1)==1 ) return( 33 );
  if ( tol(wt,78.96,0.1)==1 ) return( 34 );
  if ( tol(wt,79.91,0.1)==1 ) return( 35 );
  if ( tol(wt,83.80,0.1)==1 ) return( 36 );
  if ( tol(wt,85.47,0.1)==1 ) return( 37 );
  if ( tol(wt,87.62,0.1)==1 ) return( 38 );
  if ( tol(wt,88.91,0.1)==1 ) return( 39 );
  if ( tol(wt,91.22,0.1)==1 ) return( 40 );
  if ( tol(wt,92.91,0.1)==1 ) return( 41 );
  if ( tol(wt,95.94,0.1)==1 ) return( 42 );
  if ( tol(wt,98.91,0.1)==1 ) return( 43 );
  if ( tol(wt,101.07,0.1)==1 ) return( 44 );
  if ( tol(wt,102.90,0.1)==1 ) return( 45 );
  if ( tol(wt,106.40,0.1)==1 ) return( 46 );
  if ( tol(wt,107.87,0.1)==1 ) return( 47 );
  if ( tol(wt,112.40,0.1)==1 ) return( 48 );
  if ( tol(wt,114.82,0.1)==1 ) return( 49 );
  if ( tol(wt,118.69,0.1)==1 ) return( 50 );
  if ( tol(wt,121.75,0.1)==1 ) return( 51 );
  if ( tol(wt,127.60,0.1)==1 ) return( 52 );
  if ( tol(wt,126.90,0.1)==1 ) return( 53 );
  if ( tol(wt,131.30,0.1)==1 ) return( 54 );
  if ( tol(wt,132.91,0.1)==1 ) return( 55 );
  if ( tol(wt,137.34,0.1)==1 ) return( 56 );
  if ( tol(wt,138.91,0.1)==1 ) return( 57 );
  if ( tol(wt,140.12,0.1)==1 ) return( 58 );
  if ( tol(wt,140.91,0.1)==1 ) return( 59 );
  if ( tol(wt,144.24,0.1)==1 ) return( 60 );
  if ( tol(wt,145.00,0.1)==1 ) return( 61 );
  if ( tol(wt,150.35,0.1)==1 ) return( 62 );
  if ( tol(wt,151.96,0.1)==1 ) return( 63 );
  if ( tol(wt,157.25,0.1)==1 ) return( 64 );
  if ( tol(wt,158.92,0.1)==1 ) return( 65 );
  if ( tol(wt,162.50,0.1)==1 ) return( 66 );
  if ( tol(wt,164.93,0.1)==1 ) return( 67 );
  if ( tol(wt,167.26,0.1)==1 ) return( 68 );
  if ( tol(wt,168.93,0.1)==1 ) return( 69 );
  if ( tol(wt,173.04,0.1)==1 ) return( 70 );
  if ( tol(wt,174.97,0.1)==1 ) return( 71 );
  if ( tol(wt,178.49,0.1)==1 ) return( 72 );
  if ( tol(wt,180.95,0.1)==1 ) return( 73 );
  if ( tol(wt,183.85,0.1)==1 ) return( 74 );
  if ( tol(wt,186.20,0.1)==1 ) return( 75 );
  if ( tol(wt,190.20,0.1)==1 ) return( 76 );
  if ( tol(wt,192.22,0.1)==1 ) return( 77 );
  if ( tol(wt,195.09,0.1)==1 ) return( 78 );
  if ( tol(wt,196.97,0.1)==1 ) return( 79 );
  if ( tol(wt,200.59,0.1)==1 ) return( 80 );
  if ( tol(wt,204.37,0.1)==1 ) return( 81 );
  if ( tol(wt,207.20,0.1)==1 ) return( 82 );
  if ( tol(wt,208.98,0.1)==1 ) return( 83 );
  if ( tol(wt,210.00,0.1)==1 ) return( 84 );
  if ( tol(wt,210.00,0.1)==1 ) return( 85 );
  if ( tol(wt,222.00,0.1)==1 ) return( 86 );
  if ( tol(wt,223.00,0.1)==1 ) return( 87 );
  if ( tol(wt,226.02,0.1)==1 ) return( 88 );
  if ( tol(wt,227.00,0.1)==1 ) return( 89 );
  if ( tol(wt,232.04,0.1)==1 ) return( 90 );
  if ( tol(wt,231.00,0.1)==1 ) return( 91 );
  if ( tol(wt,238.03,0.1)==1 ) return( 92 );
  if ( tol(wt,237.00,0.1)==1 ) return( 93 );
  if ( tol(wt,244.00,0.1)==1 ) return( 94 );
  if ( tol(wt,243.00,0.1)==1 ) return( 95 );
  return(0);
}


/*************************************/
/*                                   */
/*      Z to element                 */
/*                                   */
/*************************************/
void ztoelm(char *elm, int z)
{

  if ( z==1   ) { strcpy(elm,"H"  );}
  if ( z==2   ) { strcpy(elm,"He" );}
  if ( z==3   ) { strcpy(elm,"Li" );}
  if ( z==4   ) { strcpy(elm,"Be" );}
  if ( z==5   ) { strcpy(elm,"B"  );}
  if ( z==6   ) { strcpy(elm,"C"  );}
  if ( z==7   ) { strcpy(elm,"N"  );}
  if ( z==8   ) { strcpy(elm,"O"  );}
  if ( z==9   ) { strcpy(elm,"F"  );}
  if ( z==10  ) { strcpy(elm,"Ne" );}
  if ( z==11  ) { strcpy(elm,"Na" );}
  if ( z==12  ) { strcpy(elm,"Mg" );}
  if ( z==13  ) { strcpy(elm,"Al" );}
  if ( z==14  ) { strcpy(elm,"Si" );}
  if ( z==15  ) { strcpy(elm,"P"  );}
  if ( z==16  ) { strcpy(elm,"S"  );}
  if ( z==17  ) { strcpy(elm,"Cl" );}
  if ( z==18  ) { strcpy(elm,"Ar" );}
  if ( z==19  ) { strcpy(elm,"K"  );}
  if ( z==20  ) { strcpy(elm,"Ca" );}
  if ( z==21  ) { strcpy(elm,"Sc" );}
  if ( z==22  ) { strcpy(elm,"Ti" );}
  if ( z==23  ) { strcpy(elm,"V"  );}
  if ( z==24  ) { strcpy(elm,"Cr" );}
  if ( z==25  ) { strcpy(elm,"Mn" );}
  if ( z==26  ) { strcpy(elm,"Fe" );}
  if ( z==27  ) { strcpy(elm,"Co" );}
  if ( z==28  ) { strcpy(elm,"Ni" );}
  if ( z==29  ) { strcpy(elm,"Cu" );}
  if ( z==30  ) { strcpy(elm,"Zn" );}
  if ( z==31  ) { strcpy(elm,"Ga" );}
  if ( z==32  ) { strcpy(elm,"Ge" );}
  if ( z==33  ) { strcpy(elm,"As" );}
  if ( z==34  ) { strcpy(elm,"Se" );}
  if ( z==35  ) { strcpy(elm,"Br" );}
  if ( z==36  ) { strcpy(elm,"Kr" );}
  if ( z==37  ) { strcpy(elm,"Rb" );}
  if ( z==38  ) { strcpy(elm,"Sr" );}
  if ( z==39  ) { strcpy(elm,"Y"  );}
  if ( z==40  ) { strcpy(elm,"Zr" );}
  if ( z==41  ) { strcpy(elm,"Nb" );}
  if ( z==42  ) { strcpy(elm,"Mo" );}
  if ( z==43  ) { strcpy(elm,"Tc" );}
  if ( z==44  ) { strcpy(elm,"Ru" );}
  if ( z==45  ) { strcpy(elm,"Rh" );}
  if ( z==46  ) { strcpy(elm,"Pd" );}
  if ( z==47  ) { strcpy(elm,"Ag" );}
  if ( z==48  ) { strcpy(elm,"Cd" );}
  if ( z==49  ) { strcpy(elm,"In" );}
  if ( z==50  ) { strcpy(elm,"Sn" );}
  if ( z==51  ) { strcpy(elm,"Sb" );}
  if ( z==52  ) { strcpy(elm,"Te" );}
  if ( z==53  ) { strcpy(elm,"I"  );}
  if ( z==54  ) { strcpy(elm,"Xe" );}
  if ( z==55  ) { strcpy(elm,"Cs" );}
  if ( z==56  ) { strcpy(elm,"Ba" );}
  if ( z==57  ) { strcpy(elm,"La" );}
  if ( z==58  ) { strcpy(elm,"Ce" );}
  if ( z==59  ) { strcpy(elm,"Pr" );}
  if ( z==60  ) { strcpy(elm,"Nd" );}
  if ( z==61  ) { strcpy(elm,"Pm" );}
  if ( z==61  ) { strcpy(elm,"Pm" );}
  if ( z==62  ) { strcpy(elm,"Sm" );}
  if ( z==63  ) { strcpy(elm,"Eu" );}
  if ( z==64  ) { strcpy(elm,"Gd" );}
  if ( z==65  ) { strcpy(elm,"Tb" );}
  if ( z==66  ) { strcpy(elm,"Dy" );}
  if ( z==67  ) { strcpy(elm,"Ho" );}
  if ( z==68  ) { strcpy(elm,"Er" );}
  if ( z==69  ) { strcpy(elm,"Tm" );}
  if ( z==70  ) { strcpy(elm,"Yb" );}
  if ( z==71  ) { strcpy(elm,"Lu" );}
  if ( z==72  ) { strcpy(elm,"Hf" );}
  if ( z==73  ) { strcpy(elm,"Ta" );}
  if ( z==74  ) { strcpy(elm,"W"  );}
  if ( z==75  ) { strcpy(elm,"Re" );}
  if ( z==76  ) { strcpy(elm,"Os" );}
  if ( z==77  ) { strcpy(elm,"Ir" );}
  if ( z==78  ) { strcpy(elm,"Pt" );}
  if ( z==79  ) { strcpy(elm,"Au" );}
  if ( z==80  ) { strcpy(elm,"Hg" );}
  if ( z==81  ) { strcpy(elm,"Tl" );}
  if ( z==82  ) { strcpy(elm,"Pb" );}
  if ( z==83  ) { strcpy(elm,"Bi" );}
  if ( z==84  ) { strcpy(elm,"Po" );}
  if ( z==85  ) { strcpy(elm,"At" );}
  if ( z==86  ) { strcpy(elm,"Rn" );}
  if ( z==87  ) { strcpy(elm,"Fr" );}
  if ( z==88  ) { strcpy(elm,"Ra" );}
  if ( z==89  ) { strcpy(elm,"Ac" );}
  if ( z==90  ) { strcpy(elm,"Th" );}
  if ( z==91  ) { strcpy(elm,"Pa" );}
  if ( z==92  ) { strcpy(elm,"U " );}
  if ( z==93  ) { strcpy(elm,"Np" );}
  if ( z==94  ) { strcpy(elm,"Pu" );}
  if ( z==95  ) { strcpy(elm,"Am" );}
  if ( z==96  ) { strcpy(elm,"Cm" );}
  if ( z==97  ) { strcpy(elm,"Bk" );}
  if ( z==98  ) { strcpy(elm,"Cf" );}
  if ( z==99  ) { strcpy(elm,"Es" );}
  if ( z==100 ) { strcpy(elm,"Fm" );}
  if ( z==101 ) { strcpy(elm,"Md" );}
  if ( z==102 ) { strcpy(elm,"No" );}
  if ( z==103 ) { strcpy(elm,"Lr" );}

  return(NULL);

}
