#!/bin/awk -f

BEGIN{

  # input on command line

  # el=1.83;
  # fr=0.935;
  # see notebook series: E/T, book 4, page 119
  # if el=1, and fr=1 the distances will all be equal

  NUMAT=10;
  CONVFMT="%.15g";
  PI =  4.0*atan2(1.0,1.0);
  D2R = PI/180.0;
  R2D = 180.0/PI;
  
  theta_deg=45;
  theta_rad=45*(D2R);

  # x-y plane rotation matrix
  R11=  cos(theta_rad);
  R12=  sin(theta_rad);
  R21= -sin(theta_rad);
  R22=  cos(theta_rad);

  # declare the atom positions
  at[3*NUMAT];
  
  zcap= sqrt( el^2 * (fr^2-0.5) );
  zsft= sqrt(3)*el/4;
  ztop= zsft + zcap;
  zbot= -ztop;

  x= el/sqrt(2);
  at[0*3+0]=  x; at[0*3+1]=  0; at[0*3+2]=  zsft;
  at[1*3+0]= -x; at[1*3+1]=  0; at[1*3+2]=  zsft;
  at[2*3+0]=  0; at[2*3+1]=  x; at[2*3+2]=  zsft;
  at[3*3+0]=  0; at[3*3+1]= -x; at[3*3+2]=  zsft;

  # rotate the square in the plane
  for(i=4;i<9;i++){
    at[i*3+0]=  R11*at[(i-4)*3+0]+R12*at[(i-4)*3+1];
    at[i*3+1]=  R21*at[(i-4)*3+0]+R22*at[(i-4)*3+1];
    at[i*3+2]= -zsft;
  }

  at[8*3+0]=  0; at[8*3+1]=  0; at[8*3+2]=  ztop;
  at[9*3+0]=  0; at[9*3+1]=  0; at[9*3+2]=  zbot;

  for(i=0;i<NUMAT;i++){
    printf("%20.15f%20.15f%20.15f\n",
	   at[i*3+0],
	   at[i*3+1],
	   at[i*3+2]);
  }

  for(i=0;i<NUMAT;i++){
    printf("w[%2d].x =%20.15f;\n",i+1,at[i*3+0]);
    printf("w[%2d].y =%20.15f;\n",i+1,at[i*3+1]);
    printf("w[%2d].z =%20.15f;\n",i+1,at[i*3+2]);

  }

  exit;
}
