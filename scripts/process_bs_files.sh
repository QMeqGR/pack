#!/bin/sh

if [ "$1" = "" ]; then
    echo "Use: $ process_bs_files.sh r_tol"
    exit
fi

tolerance=$1;

# extract POSCARS
echo "Extracting POSCARS..."
for file in *.bs; do
    number=`cat $file | awk '($2~/commandline/){print $NF}'`
    echo "getting POSCAR "$number
    cat $file | awk '($1=="*P"){print $0}' | sed 's/*P //' > POSCAR-n$number
done

# get the energy data
echo "Getting energy data and symmetries..."
for file in *.bs; do
    number=`cat $file | awk '($2~/commandline/){print $NF}'`
    echo -e out-$number '\t' `cat $file | awk -f ~/awkfiles/pack_info.awk | grep e0` \
        "(tol=$tolerance)"\
        `symsearch -f POSCAR-n$number -r $tolerance | grep Space`
done | sort -n -k 3 | tee energies.dat

exit

