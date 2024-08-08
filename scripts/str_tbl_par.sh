#!/bin/sh

# for this run, total restart file length is 27
# top 10 lines are input file
# bottom 17 lines are positions

n=0
node=48
directory=`pwd`
echo "running in directory "$directory

############################
# extract the input files and position data
echo "Extracting restart data and forming cel, inp and pos files"
for file in pack*.bs; do
    cat $file | awk '($1=="*R"){print $0}' | cut -d" " -f 2-300 > restart.dat &&
    head --lines=2 restart.dat > cel_$n ;
    head --lines=8 restart.dat | tail --lines=6 > inp_$n ;
    tail --lines=22 restart.dat > pos_$n ;
    n=$((n+1)) ;
done

rm -f restart.dat ;

inp_list=`ls inp_*`
pos_list=`ls pos_*`

# echo "inp list= "$inp_list
# echo "pos list= "$pos_list

#############################
# pair the files together
for file_inp in $inp_list; do

    for file_pos in $pos_list; do


        cat cel_${file_pos##*_} $file_inp > inp.dat ;
        cat cel_${file_pos##*_} $file_inp $file_pos > restart.dat ;

	outbs="out_"${file_inp##*_}"X"${file_pos##*_} ;
	echo "# Simplexing "$outbs

	pack_run -f $directory/inp.dat -N$node -o $outbs -O " -X -E -V -q 10000 " &
	cat restart.dat ;
	sleep 7

	echo "job sent to node "$node
	node=$((node-1)) ;
	if [ $node -eq 30 ]; then
	    node=48 ;
	    echo "one batch sent: will sleep for 30 minutes"
	    sleep 1800
	fi
	
    done

done

rm -f inp_* pos_* cel_* ;

exit
