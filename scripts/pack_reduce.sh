#!/bin/bash

for file in *.bs; do
    cat $file | awk '($2=="Rcat_1"){Rc=$4}($2=="Rvert"){Rv=$4}($2=="an_chrg_c"){qc=$4}($2=="ecc"){ecc=$4} END{printf("%10.6f%10.6f%10.6f%10.6f\n",ecc,qc,Rc,Rv);}' >> pack_reduce.dat;
done
