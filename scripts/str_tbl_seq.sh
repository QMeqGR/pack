#!/bin/sh

# for this run, total restart file length is 27
# top 10 lines are input file
# bottom 17 lines are positions

n=0

############################
# extract the input files and position data
for file in blah*.bs; do
    pack_run -o ${file%%.bs} -R &&
    head --lines=2 restart.dat > cel_$n ;
    head --lines=8 restart.dat | tail --lines=6 > inp_$n ;
    tail --lines=22 restart.dat > pos_$n ;
    n=$((n+1)) ;
done

rm -f restart.dat ;

#############################
# pair the files together
for file_inp in inp_*; do

    for file_pos in pos_*; do

	cat cel_${file_pos##*_} $file_inp > inp.dat ;
	cat cel_${file_pos##*_} $file_inp $file_pos > restart.dat ;

	outbs="out_"${file_inp##*_}"X"${file_pos##*_} ;
	echo "# Simplexing "$outbs

	pack -f inp.dat -X -E -V -q 5 > $outbs".bs" ;
	
    done

done

# rm -f inp_* pos_* cel_* restart.dat ;

exit
