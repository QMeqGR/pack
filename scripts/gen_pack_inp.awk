#!/bin/awk -f

#
# Usage: gen_pack_inp.awk
#

#
# awk file to generate 'pack' input files
#


# version 3.1 Wed Feb  1 09:02:01 PST 2006
# -add dimer support.
#
# version 3.0 Fri Jan 20 12:27:39 PST 2006
# -clean up the mess.  Must edit this file.
#
# version 2.1 Thu Jan 19 15:37:40 PST 2006
# -added Rv_max depending on value of d
# and anion coordination; to prevent overlap
# of vertex atoms on the same anion.  See
# meeting book #2, p95, and alanate book #3, p135.
#
#
# version 2.0
#
# ehm SNL/CA 12 january 2006
#

BEGIN{

  DEBUG=0;

  CONVFMT="%.10g";
  OFMT="%.10g";

  count=0;
  DELTA=1e-3;       # a small number
  MAX_FILES=1000;   # don't generate more than X files.
  version = 4;      # pack major version

  cell= 10.0;       # starting cell length
  angles = 90.0;    # fixed to orth symmetry

##########################################
#    Formula unit details
##########################################

  fu_min = 1;       # minimum number of formula units
  fu_max = 2;

##########################################
#    Anion details
##########################################

  n_an = 3;         # number of anions
  an_coord = 0;     # 0=tetrahedral   1=octahedral  2=dimer
                    # these values are what I've been using
                    # from the beginning
  d = 1.23;         # center to vertex distance in cation
                    # B-H  = 1.23 (tet)
                    # Al-H = 1.65 (tet)

  Zcent = 5;       # atomic number of anion center
  Zvert = 1;

  Rv_min = 1.1;     # radius of anion vertex
  Rv_max = 1.3;      # if (Rv_max == -1) it will be calculated below!
  Rv_step = 0.05;

  q_an_c_min = 0.0; # anion center charge
  q_an_c_max = 0.2;
  q_an_c_step= 0.1;


##########################################
#    Cation details
##########################################

  n_type_ca = 2;   # number of types of cations 1, or 2

  n_ca_1 = 1;      # number of cations of type 1
  q_ca_1 = 2;      # charge on cations of type 1
  Rc1_min=  1.2;   # radius of cation 1
  Rc1_max=  1.4;
  Rc1_step= 0.05;
  Z_cation_1 = 20;

  n_ca_2 = 1;
  q_ca_2 = 1;
  Rc2_min=  0.9;
  Rc2_max=  1.3;
  Rc2_step= 0.1;
  Z_cation_2 = 19;



#########################################
# DON'T CHANGE ANYTHING BELOW THIS LINE!
#########################################



#####################################
# Set the appropriate Rv_max values
#####################################

  if ( an_coord==0 ) {
    printf("Generating inputs for tetrahedral anions.\n");
    num_vert=4;
  }
  if ( an_coord==1 ) {
    printf("Generating inputs for octahedral anions.\n");
    num_vert=6;
  }
  if ( an_coord==2 ) {
    printf("Generating inputs for dimer anions.\n");
    num_vert=1;
  }

  if (an_coord==0 && Rv_max == -1 ) {
    Rv_max = sqrt(2)*d/sqrt(3) - DELTA;
  }
  if (an_coord==1 && Rv_max == -1 ) {
    Rv_max = d/sqrt(2) - DELTA;
  }
  if (an_coord==2 && Rv_max == -1 ) {
    Rv_max = d - DELTA;
  }

  
  printf("Max vertex radius = %f\n",Rv_max);

  if (DEBUG){
    printf("fu_min = %d\n",fu_min);
    printf("fu_max = %d\n",fu_max);
    printf("q_an_c_min = %f\n",q_an_c_min);
    printf("q_an_c_max = %f\n",q_an_c_max);
    printf("Rv_min = %f\n",Rv_min);
    printf("Rv_max = %f\n",Rv_max);
    printf("Rc1_min = %f\n",Rc1_min);
    printf("Rc1_max = %f\n",Rc1_max);
  }

###################################################################
#     Main loops    
###################################################################
  
  for(fu=fu_min; fu<fu_max; fu++){
    for( q_an_c = q_an_c_min; q_an_c < q_an_c_max; q_an_c += q_an_c_step ){
      for(Rv=Rv_min; Rv<Rv_max; Rv+=Rv_step) {
	Rc = d - Rv - DELTA;
	
	for(Rc1=Rc1_min; Rc1<Rc1_max; Rc1+=Rc1_step) {

	  if ( DEBUG ) {
	    printf("fu=%d  q_an_c=%f  Rv=%f  Rc1=%f\n",fu,q_an_c,Rv,Rc1);
	  }
	  
	  if ( count > MAX_FILES ){
	    printf("Exceeded max file count.\n");
	    exit;
	  }
	  
	  Q_an   = -1.0 / n_an * ( n_ca_1 * q_ca_1 );
	  q_an_v =  1.0 / num_vert * ( -q_an_c + Q_an );
	  filename="pack_v" version "_" fu "fu_" q_an_c "qc_" Rc1 "Rc1_" Rv "Rv_" ".inp";
	  
	  if ( n_type_ca == 1 ) {
	    system("touch " filename);
	    count+=1;
	    
	    N_an = n_an * fu;
	    N_ca_1 = n_ca_1 * fu;
	    R_cation_1 = Rc1;
	    
	    
	    system("echo " cell " " cell " " cell " " angles " " angles " " angles   " >> " filename);
	    system("echo " an_coord " " N_an          " >> " filename);
	    system("echo " q_an_c " " q_an_v          " >> " filename);
	    system("echo " d " " Rc " " Rv            " >> " filename);
	    system("echo " Zcent " " Zvert            " >> " filename);
	    system("echo " n_type_ca                  " >> " filename);
	    system("echo " N_ca_1 " " q_ca_1          " >> " filename);
	    system("echo " R_cation_1 " " Z_cation_1  " >> " filename);
	    
	    continue;
	  }
	  
	  for(Rc2=Rc2_min; Rc2<Rc1; Rc2+=Rc2_step) {
	    
	    Q_an   = -1.0 / n_an * ( n_ca_1 * q_ca_1 + n_ca_2 * q_ca_2 );
	    q_an_v =  1.0 / num_vert * ( -q_an_c + Q_an );
	    filename="pack_v" version "_" fu "fu_" q_an_c "qc_" Rc1 "Rc1_" Rc2 "Rc2_" Rv "Rv_" ".inp";
	    
	    system("touch " filename);
	    count+=1;
	    
	    N_an = n_an * fu;
	    N_ca_1 = n_ca_1 * fu;
	    R_cation_1 = Rc1;
	    N_ca_2 = n_ca_2 * fu;
	    R_cation_2 = Rc2;
	    
	    system("echo " cell " " cell " " cell " " angles " " angles " " angles   " >> " filename);
	    system("echo " an_coord " " N_an          " >> " filename);
	    system("echo " q_an_c " " q_an_v          " >> " filename);
	    system("echo " d " " Rc " " Rv            " >> " filename);
	    system("echo " Zcent " " Zvert            " >> " filename);
	    system("echo " n_type_ca                  " >> " filename);
	    system("echo " N_ca_1 " " q_ca_1          " >> " filename);
	    system("echo " R_cation_1 " " Z_cation_1  " >> " filename);
	    system("echo " N_ca_2 " " q_ca_2          " >> " filename);
	    system("echo " R_cation_2 " " Z_cation_2  " >> " filename);
	  }
	}
      }
    }
  }
  
  
  printf("Done. Generated %d input files.\n",count);
  
}
