BEGIN{

#  debug=1;
  atom_line_count=0;
  cen_at_cnt=0;
  center_flag=0;

  # set on command line
  D=2.95;
  BHD=1.12;

# set in code
#  z[];
#  nz[];
#  at[];
#  atz[];
#  cen_at[];
#  hlist[];
#  vat;
#  bat[];
#  Dtol;


}

(NR==1 && $1~/Z:/){
  ntypat = NF-1;
  if ( debug>1 ) printf("ntypat = %d\n",ntypat);
  for(i=2;i<NF+1;i++){
    z[i-1] = $i;
    if ( debug>1 ) printf("Z line: z[%2d]= %2d\n",i-1,$i);
  }
}
(NR==2){ line2=$0; }
(NR==3){ line3=$0; }
(NR==4){ line4=$0; }
(NR==5){ line5=$0; }
(NR==6){
  for(i=1;i<NF+1;i++){
    if ( debug>1 ) printf("nat line: i=%2d, field=%2d, nz[%2d]=%2d\n",
			  i,$i,i,$i);
    nz[i] = $i;
    natoms += $i;

  }
}
(NR==7){
  line7=$0;
  if ( $1~/D/ || $1~/d/ ){
    printf("CONTCAR must be in cartesian coords!\n");
    exit;
  }
}
(NR>=8 && NF==3){

  at[3*atom_line_count+0]=$1;
  at[3*atom_line_count+1]=$2;
  at[3*atom_line_count+2]=$3;
  if ( debug ) printf("atom %3d: %s\n",atom_line_count,$0);

  tmp_at_count=0;
  for (k=1;k<ntypat+1;k++){
    tmp_at_count += nz[k];
    if ( debug>1 ) printf("setz: k=%d line_cnt=%2d tmp_at_count=%2d\
 z[%d]=%d\n",k,atom_line_count,tmp_at_count,k,z[k]);
    if ( atom_line_count == 0 ) { Z=z[1]; }
    if ( atom_line_count >= tmp_at_count ) {
      Z=z[k+1];
      if ( debug>1 ) printf("setz: setting Z=%d\n",Z);
    }
  }
  atz[atom_line_count]=Z;

  atom_line_count += 1;
}

END{

  Dtol= 0.02*D;
  if ( debug ) printf("D    = %f\n",D);
  if ( debug ) printf("Dtol = %f\n",Dtol);

  if ( debug ) printf("number of atoms= %d\n",natoms);
  for(i=0;i<natoms;i++){
    if ( debug ) printf("at[%3d]= %20.15f%20.15f%20.15f%5d\n",
			i,at[3*i+0],at[3*i+1],at[3*i+2],atz[i]);
  }

  for(i=0;i<natoms;i++){
    t1x=at[3*i+0];
    t1y=at[3*i+1];
    t1z=at[3*i+2];
    center_flag=0;
    for(j=0;j<natoms;j++){
      if (i==j) continue;
      t2x=at[3*j+0];
      t2y=at[3*j+1];
      t2z=at[3*j+2];
      dij= sqrt( (t1x-t2x)^2 + (t1y-t2y)^2 + (t1z-t2z)^2 );
      if ( debug>1 ) printf("d[%3d:%3d]= %f\n",i,j,dij);
      if ( dij > D-Dtol && dij < D+Dtol ){ center_flag++; }
    }
    if ( center_flag==12 ){ cen_at[cen_at_cnt++]=i; }
  }

  for(k=0;k<cen_at_cnt;k++){
    if ( debug ) printf("center atom %d = atom %d\n", k, cen_at[k]);
    t1x=at[3*cen_at[k]+0];
    t1y=at[3*cen_at[k]+1];
    t1z=at[3*cen_at[k]+2];
    for(j=0;j<natoms;j++){
      if (i==j) continue;
      t2x=at[3*j+0];
      t2y=at[3*j+1];
      t2z=at[3*j+2];
      dij= sqrt( (t1x-t2x)^2 + (t1y-t2y)^2 + (t1z-t2z)^2 );
      if ( debug>1 ) printf("d[%3d:%3d]= %f\n",i,j,dij);
      if ( dij > D-Dtol && dij < D+Dtol ){
	hlist[vat++]=j;
	if ( debug ) printf("hlist[%2d] = atom %d\n",vat-1,j);
	vecx = t1x - t2x;
	vecy = t1y - t2y;
	vecz = t1z - t2z;
	vecmag=sqrt( vecx^2 + vecy^2 + vecz^2 );
	nvecx = vecx/vecmag;
	nvecy = vecy/vecmag;
	nvecz = vecz/vecmag;
	bat[3*(vat-1)+0] = t2x + nvecx*BHD;
	bat[3*(vat-1)+1] = t2y + nvecy*BHD;
	bat[3*(vat-1)+2] = t2z + nvecz*BHD;
      }
    }
  }

  # print out the output POSCAR in cartesian coords
  if ( debug ) printf("\nBegin output\n\n");

  printf("Z: ");
  for(k=1;k<ntypat+1;k++){
    if ( z[k] != 5 ) printf("%d ",z[k]);
  }
  printf("5\n");

  printf("%s\n",line2);
  printf("%s\n",line3);
  printf("%s\n",line4);
  printf("%s\n",line5);

  for(k=1;k<ntypat+1;k++){
    if ( z[k] != 5 ) {
      printf("%d ",nz[k]);
    }
    if ( z[k] == 5 ) {
      printf("%d ",12*nz[k]);
    }
  }
  printf("\n");
  printf("%s\n",line7);
  
  for(k=1;k<ntypat+1;k++){
    if ( z[k] != 5 ) {
      for(j=0;j<natoms;j++){
	if ( atz[j]==z[k] ) printf("%20.15f%20.15f%20.15f\n",
				   at[3*j+0],at[3*j+1],at[3*j+2]);
      }
    }
  }
  # print B atoms last
  for(k=1;k<ntypat+1;k++){
    if ( z[k] == 5 ) {
      for(h=0;h<vat;h++){
	printf("%20.15f%20.15f%20.15f\n",
	       bat[3*h+0],bat[3*h+1],bat[3*h+2]);
      }
    }
  }




}
