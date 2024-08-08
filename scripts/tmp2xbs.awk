#
# note the file format is hard coded here for
# line breaks at special locations. Need to
# code around this for easier use someday.
#

BEGIN{

  MAX=500;
  count=0;
  r[MAX]=0.0;
  Lt[9]=0.0;
  tiny = 0.001;

  rc1 = 1.0; gc1 = 0.0; bc1 = 0.0;
  rc2 = 0.5; gc2 = 0.5; bc2 = 0.0;

  rac = 0.0; gac = 0.0; bac = 1.0;
  rav = 0.0; gav = 1.0; bav = 0.0;


  # input variables

  a=0; b=0; c=0;
  alph=0; beta=0; gamm=0;
  
  num_vert=0; n_an=0;

  q_an_c=0; q_an_v=0;
  
  d=0; Rc=0; Rv=0;

  Zcent=0; Zvert=0;
  
  n_type_ca=0;
  
  n_ca_1=0; q_ca_1=0;
  R_cation_1=0; Z_cation_1=0;

  n_ca_2=0; q_ca_2=0;
  R_cation_2=0; Z_cation_2=0;

}

########################
#  functions
########################
function get_Lt(Lt,a,b,c,alph,beta,gamm){

  if ( debug ){
    printf("a b c = %20.10f%20.10f%20.10f\n",a,b,c);
    printf("angles = %20.10f%20.10f%20.10f\n",alph,beta,gamm);
  }
  
  ax = a;
  
  bx = b * cos(gamm);
  by = b * sin(gamm);

  cx = c * cos(beta);
  cy = ( 1.0 / ( b * sin(gamm) ) ) * ( c * b * cos(alph) - c * b * cos(beta) * cos(gamm) );
  cz = sqrt( c*c - cx*cx - cy*cy );

  if ( debug ){
    printf("L: not Lt!!\n");
    printf("%20.10f%20.10f%20.10f\n",a,0.0,0.0);
    printf("%20.10f%20.10f%20.10f\n",bx,by,0.0);
    printf("%20.10f%20.10f%20.10f\n",cx,cy,cz);
  }

  Lt[0] = ax;
  Lt[1] = bx;
  Lt[2] = cx;

  Lt[3] = ay;
  Lt[4] = by;
  Lt[5] = cy;

  Lt[6] = az;
  Lt[7] = bz;
  Lt[8] = cz;

  return;
}

function dist(at1,at2){
  d = sqrt( (r[at1+0]-r[at2+0])*(r[at1+0]-r[at2+0]) + (r[at1+1]-r[at2+1])*(r[at1+1]-r[at2+1]) + (r[at1+2]-r[at2+2])*(r[at1+2]-r[at2+2]) );
  return d;
}

function dist2(t,s){
  d = sqrt( (t[0]-s[0])*(t[0]-s[0]) + (t[1]-s[1])*(t[1]-s[1]) + (t[2]-s[2])*(t[2]-s[2]) );
  return d;
}

########################
# pattern action rules
########################

(NR==1){a=$1; b=$2; c=$3;}
(NR==2){alph=$1; beta=$2; gamm=$3;}
(NR==3){num_vert=$1; n_an=$2;}
(NR==4){q_an_c=$1; q_an_v=$2;}
(NR==5){d=$1; Rc=$2; Rv=$3; Zcent=$4; Zvert=$5;}
(NR==6){n_type_ca=$1;}
(NR==7){n_ca_1=$1; q_ca_1=$2;}
(NR==8){R_cation_1=$1; Z_cation_1=$2;}
(NR==9  && n_type_ca==2 ){n_ca_2=$1; q_ca_2=$2;}
(NR==10 && n_type_ca==2 ){R_cation_2=$1; Z_cation_2=$2;}
(n_type_ca==1 && NR>8){
  r[3*count+0] = $1;
  r[3*count+1] = $2;
  r[3*count+2] = $3;
  count++;
}
(n_type_ca==2 && NR>10){
  r[3*count+0] = $1;
  r[3*count+1] = $2;
  r[3*count+2] = $3;
  count++;
}

END{

  if ( debug ){
    printf("Input file:\n");
    printf("%15.10f%15.10f%15.10f\n",a,b,c);
    printf("%15.10f%15.10f%15.10f\n",alph,beta,gamm);
    printf("%15d%15d\n",num_vert,n_an);
    printf("%15.10f%15.10f\n",q_an_c,q_an_v);
    printf("%15.10f%15.10f%15.10f%5d%5d\n",d,Rc,Rv,Zcent,Zvert);
    printf("%15d\n",n_type_ca);
    printf("%15d%5d\n",n_ca_1,q_ca_1);
    printf("%15.10f%5d\n",R_cation_1,Z_cation_1);
    if (n_type_ca==2){
      printf("%15d%5d\n",n_ca_2,q_ca_2);
      printf("%15.10f%5d\n",R_cation_2,Z_cation_2);
    }
    
    for(i=0;i<count;i++){
      printf("%15.10f%15.10f%15.10f\n",r[3*i+0],r[3*i+1],r[3*i+2]);
    }
  }

  get_Lt(Lt,a,b,c,alph,beta,gamm);

  if (debug){
    printf("Lt: ");
    for(i=0;i<9;i++){
      printf("%8.3f",Lt[i]);
    }
    printf("\n");
  }

  for(i=0;i<count;i++){
    tx = r[3*i+0];
    ty = r[3*i+1];
    tz = r[3*i+2];
    r[3*i+0] = tx*Lt[0] + ty*Lt[1] + tz*Lt[2];
    r[3*i+1] = tx*Lt[3] + ty*Lt[4] + tz*Lt[5];
    r[3*i+2] = tx*Lt[6] + ty*Lt[7] + tz*Lt[8];
  }

# print cations first
  for(i=0; i<n_ca_1; i++){
    printf("atom%5s%15.10f%15.10f%15.10f\n","C1",r[3*i+0],r[3*i+1],r[3*i+2]);
  }
  if (n_ca_2>0){
    for(i=n_ca_1; i<(n_ca_1+n_ca_2); i++){
      printf("atom%5s%15.10f%15.10f%15.10f\n","C2",r[3*i+0],r[3*i+1],r[3*i+2]);
    }
  }
  
  k=3*n_ca_1+3*n_ca_2;
  for(i=0; i<n_an; i++){
    for(j=0;j<num_vert;j++){
      K=k+i*num_vert*3+3*j;
      if ( j==0 ) printf("atom%5s%15.10f%15.10f%15.10f\n","AC",r[K+0],r[K+1],r[K+2]);
      else printf("atom%5s%15.10f%15.10f%15.10f\n","AV",r[K+0],r[K+1],r[K+2]);
    }
  }

  printf("spec C1%10.6f%10.6f%10.6f%10.6f\n",R_cation_1,rc1,gc1,bc1);
  if ( n_type_ca == 2 ) printf("spec C2%10.6f%10.6f%10.6f%10.6f\n",R_cation_2,rc2,gc2,bc2);
  printf("spec AC%10.6f%10.6f%10.6f%10.6f\n",Rc,rac,gac,bac);
  printf("spec AV%10.6f%10.6f%10.6f%10.6f\n",Rv,rav,gav,bav);
  printf("inc 5.0\n");


  v0[0] = 0;
  v0[1] = 0;
  v0[2] = 0;

  v1[0] = ax;
  v1[1] = ay;
  v1[2] = az;

  v2[0] = ax+bx;
  v2[1] = ay+by;
  v2[2] = az+bz;

  v3[0] = bx;
  v3[1] = by;
  v3[2] = bz;

  v4[0] = cx;
  v4[1] = cy;
  v4[2] = cz;

  v5[0] = ax+cx;
  v5[1] = ay+cy;
  v5[2] = az+cz;

  v6[0] = ax+bx+cx;
  v6[1] = ay+by+cy;
  v6[2] = az+bz+cz;

  v7[0] = bx+cx;
  v7[1] = by+cy;
  v7[2] = bz+cz;

  printf("atom 0 %10.5f%10.5f%10.5f\n",0,0,0);
  printf("atom 1 %10.5f%10.5f%10.5f\n",ax,ay,az);
  printf("atom 2 %10.5f%10.5f%10.5f\n",ax+bx,ay+by,az+bz);
  printf("atom 3 %10.5f%10.5f%10.5f\n",bx,by,bz);
  printf("atom 4 %10.5f%10.5f%10.5f\n",cx,cy,cz);
  printf("atom 5 %10.5f%10.5f%10.5f\n",ax+cx,ay+cy,az+cz);
  printf("atom 6 %10.5f%10.5f%10.5f\n",ax+bx+cx,ay+by+cy,az+bz+cz);
  printf("atom 7 %10.5f%10.5f%10.5f\n",bx+cx,by+cy,bz+cz);

  printf("spec 0 0.05 0 0 0\n");
  printf("spec 1 0.05 0 0 0\n");
  printf("spec 2 0.05 0 0 0\n");
  printf("spec 3 0.05 0 0 0\n");
  printf("spec 4 0.05 0 0 0\n");
  printf("spec 5 0.05 0 0 0\n");
  printf("spec 6 0.05 0 0 0\n");
  printf("spec 7 0.05 0 0 0\n");

  d01 = dist2(v0,v1);
  d32 = dist2(v3,v2);
  d45 = dist2(v4,v5);
  d67 = dist2(v6,v7);
  d12 = dist2(v1,v2);
  d03 = dist2(v0,v3);
  d56 = dist2(v5,v6);
  d47 = dist2(v4,v7);
  d04 = dist2(v0,v4);
  d15 = dist2(v1,v5);
  d26 = dist2(v2,v6);
  d37 = dist2(v3,v7);
  printf("bonds 0 1 %10.6f%10.6f 0.01 Black\n",d01-tiny,d01+tiny);
  printf("bonds 3 2 %10.6f%10.6f 0.01 Black\n",d32-tiny,d32+tiny);
  printf("bonds 4 5 %10.6f%10.6f 0.01 Black\n",d45-tiny,d45+tiny);
  printf("bonds 6 7 %10.6f%10.6f 0.01 Black\n",d67-tiny,d67+tiny);
  printf("bonds 1 2 %10.6f%10.6f 0.01 Black\n",d12-tiny,d12+tiny);
  printf("bonds 0 3 %10.6f%10.6f 0.01 Black\n",d03-tiny,d03+tiny);
  printf("bonds 5 6 %10.6f%10.6f 0.01 Black\n",d56-tiny,d56+tiny);
  printf("bonds 4 7 %10.6f%10.6f 0.01 Black\n",d47-tiny,d47+tiny);
  printf("bonds 0 4 %10.6f%10.6f 0.01 Black\n",d04-tiny,d04+tiny);
  printf("bonds 1 5 %10.6f%10.6f 0.01 Black\n",d15-tiny,d15+tiny);
  printf("bonds 2 6 %10.6f%10.6f 0.01 Black\n",d26-tiny,d26+tiny);
  printf("bonds 3 7 %10.6f%10.6f 0.01 Black\n",d37-tiny,d37+tiny);

  printf("dup %10.5f%10.5f%10.5f\n",ax,ay,az);
  printf("dup %10.5f%10.5f%10.5f\n",bx,by,bz);
  printf("dup %10.5f%10.5f%10.5f\n",cx,cy,cz);

}
