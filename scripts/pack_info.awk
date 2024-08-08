BEGIN{

  #set formula units example: "-v fu=2" on command line.

  cvol = -1;
  ecc = -1;
  eng = -1;
  pf = -1;
  ar = -1;
  ortho = -1;
  an_chrg_c = -1e9;
  Z=1e10;
  max_trn_rej=1e-10;
  max_rot_rej=1e-10;
  max_lat_rej=1e-10;
  max_swp_rej=1e-10;
  min_trn_rej=1e10;
  min_rot_rej=1e10;
  min_lat_rej=1e10;
  min_swp_rej=1e10;
  version="0.0.0";
  metro_start_flag=0;
}

($0 ~/-- Starting Metropolis calculation --/){ metro_start_flag=1; }
($0 ~/Starting Metropolis calculation.../){ metro_start_flag=1; }
($3 == "(version"){ version = $4; }
($2 == "commandline"){ cline=$0; }
($2 == "ecc"){ ecc = $NF; }
($2 == "eng"){ eng = $NF; }
($2=="cell" && $3=="vol"){ cvol=$5; }
($2=="packing" && $3=="frac"){ pf=$5; }
($2=="aspect" && $3=="ratio"){ ar=$5; }
($2=="an_chrg_c"){ an_chrg_c=$4; }
($2=="cv/ortho"){ ortho=$4; }
($1=="*M" && metro_start_flag){
  trn=$11;
  rot=$12;
  lat=$13;
  swp=$14;

  if ( trn > max_trn_rej ) max_trn_rej=trn;
  if ( rot > max_rot_rej ) max_rot_rej=rot;
  if ( lat > max_lat_rej ) max_lat_rej=lat;
  if ( swp > max_swp_rej ) max_swp_rej=swp;

  if ( trn < min_trn_rej && trn>0.0 ) min_trn_rej=trn;
  if ( rot < min_rot_rej && rot>0.0 ) min_rot_rej=rot;
  if ( lat < min_lat_rej && lat>0.0 ) min_lat_rej=lat;
  if ( swp < min_swp_rej && swp>0.0 ) min_swp_rej=swp;
}

END{
  printf("%s\n",cline);

  Ncommand=split(cline,command,"_");
#  printf("Ncommand=%d\n",Ncommand);
  for( i=0; i<(Ncommand+1); i++){
# printf("command[%d] = %s\n",i,command[i]);
    if ( command[i] ~ /fu/ ) {
      
#     printf("found one!: %s\n",command[i]);
      zend=match(command[i],"fu");
#     printf("zend = %d\n",zend);
      Z = substr( command[i], 0 , zend-1 );
    }
  }



  printf("version= (%s\n",version);
  printf("Rejection rates: min ---> max\n");
  printf("trn: %6.2f ---> %6.2f\n",min_trn_rej,max_trn_rej);
  printf("rot: %6.2f ---> %6.2f\n",min_rot_rej,max_rot_rej);
  printf("lat: %6.2f ---> %6.2f\n",min_lat_rej,max_lat_rej);
  printf("swp: %6.2f ---> %6.2f\n",min_swp_rej,max_swp_rej);
  printf("an_chrg_c= %.5f\n",an_chrg_c);
  printf("ecc= %+.5e\n",ecc);
  printf("eng= %+.5e\n",eng);
  printf("fu = %d\n",fu);
  printf("e0pack= %.10f\n",ecc+eng);
  if ( fu > 0 ) printf("e0pack/fu= %.6f\n",(ecc+eng)/fu);
  printf("Z= %s\n",Z);
  if ( Z != 1e10 ) { printf("e_pack/fu= %.10f / %s = %f\n",ecc+eng,Z, (ecc+eng)/Z); }
  if ( Z == 1e10 ) { printf("e_pack/fu= %.10f /  =\n",ecc+eng); }
  printf("vol= %+.4f\n",cvol);
  printf("pf= %+.4f\n",pf);
  printf("ar= %+.4f\n",ar);
  printf("ortho= %+.4f\n",ortho);
}
