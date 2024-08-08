#!/bin/sh

AWKFILE=/home/ehm/abinitio/pack/scripts/make_b12h12.awk

for file in POSCAR*; do

echo "processing file "$file
poscar_cnvrt -f $file -c > TEMPCAR
cat TEMPCAR | awk -f $AWKFILE > TEMPCAR_2
poscar_cnvrt -f TEMPCAR_2 -C > $file"_B"

done

rm -f TEMPCAR TEMPCAR_2

exit
