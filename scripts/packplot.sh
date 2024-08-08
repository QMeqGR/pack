#!/bin/bash

count=0;
bdir=/home/users/ehmajzo/bin/
wdir=`pwd`
oname=xmgr

switch=$1

if [ -z $switch ]; then
    echo
    echo "You must enter a choice for output!"
    echo "r --- rm -f all xmgr.*.dat files"
    echo "a --- append all xmgr.*.dat"
    echo
    exit;
fi

nlist=`ls *.bs`
#echo $nlist

if [ 'r' = $switch ]; then
    rm -f xmgr.*.dat
fi

#################################
#      Metro Output             #
#################################
# header from version 2.2.4
# *  run      Temp       E_ion        E_ss   cel vol        pf        ar          Cv        Sig2            % Rejected Changes
# run number is column 2
# temp is column 3

# fetch information from one column in each file and plot
for num in 1 2 3 4 5 6 7 8 9 10 11 12; do
    case $num in
	1 ) lab="Temp" ;;
	2 ) lab="E_ion" ;;
	3 ) lab="E_rep" ;;
	4 ) lab="cel_vol" ;;
	5 ) lab="pf" ;;
	6 ) lab="ar" ;;
	7 ) lab="Cv" ;;
	8 ) lab="Esig2" ;;
	9 ) lab="trn" ;;
	10) lab="rot" ;;
	11) lab="lat" ;;
	12) lab="swp" ;;
    esac
    outfile=$oname.metro.$lab.dat
#    echo "@with g"$num >> $outfile
    echo "@     xaxis label \"Temperature\"" >> $outfile
    echo "@     yaxis label \"$lab \"" >> $outfile
    count=0
    for fil in $nlist; do
	col=$(($num+2))
	echo "@     s$count legend \"$fil \"" >> $outfile
	echo "@     s$count comment \"$fil\"" >> $outfile
#	echo "@target G"$num".S$count" >> $outfile
	echo "@type xy" >> $outfile
	cat $fil | awk -v A=$col '($1=="*M") { printf("%12.3e%20.5e\n",$3,$A); }'  >> $outfile
	echo "&" >> $outfile
	count=$(($count+1))
    done
done

#################################
#     Simplex Output            #
#################################
#*      nfunk           rtol        pf        ar             e_ion              e_ss             e_tot

# fetch information from one column in each file and plot
for num in 1 2 3 4 5 6; do
    case $num in
	1 ) lab="rtol" ;;
	2 ) lab="pf" ;;
	3 ) lab="ar" ;;
	4 ) lab="e_ion" ;;
	5 ) lab="e_rep" ;;
	6 ) lab="e_tot" ;;
    esac
    outfile=$oname.simplex.$lab.dat
#    echo "@with g0" >> $outfile
    echo "@     xaxis label \"Nfunk\"" >> $outfile
    echo "@     yaxis label \"$lab \"" >> $outfile
    count=0
    for fil in $nlist; do
	col=$(($num+2))
	echo "@     s$count legend \"$fil \"" >> $outfile
	echo "@     s$count comment \"$fil\"" >> $outfile
#	echo "@target G0.S$count" >> $outfile
	echo "@type xy" >> $outfile
	cat $fil | awk -v A=$col '($1=="*X") { printf("%12.3e%20.5e\n",$2,$A); }'  >> $outfile
	echo "&" >> $outfile
	count=$(($count+1))
    done
done


exit
