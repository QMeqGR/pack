#!/bin/sh

#
# This script should be called in the directory ~/pack
#
# Usage: script_run NaAlH4 v4831 10
#
# This will generate 10 runs in the directory NaAlH4, and the output
# string for each run will be out_v4831_n, where n is the instance number.
# 
# -assumes input file name is input.in
#

dir_name=$1;
outstring=$2;
instances=$3;


HOME=$HOME
PACK=$HOME/bin/pack


if [ $# -lt 3 ]; then
    echo "Number of command line arguments= "$#
    echo "Usage: script_run.sh DIRNAME OUTSTRING NUMBER"
    exit
fi

#########
# BEGIN #
#########

echo "Changing directory to "$dir_name;
echo "Using: outstring= "$outstring " and instances= "$instances
cd $dir_name;

if [ -f PACK_OPTIONS ]; then
    OPTIONS=`cat PACK_OPTIONS`;
else
    echo "WARN: PACK_OPTIONS file not found, using default in script setup."
    #OPTIONS="-r 1000 -t 10 -B 25";
    OPTIONS="-V -W -r 2 -t 2";
fi

directory=`pwd`;
input=input.in;

echo "BEGIN{for(i=0;i<$instances;i++){printf(\"%d \",i);}}" > tmp.awk;
i_list=`awk -f tmp.awk`

for n in $i_list; do

    script_n="$dir_name"_$n.sh;
    output=out_"$outstring"_$n.bs
    echo "#!/bin/sh" > $script_n;
    echo "#$ -o $directory" >> $script_n;
    echo "#$ -e $directory" >> $script_n;
    echo "#$ -cwd" >> $script_n;
    echo "cd $directory" >> $script_n;
    echo "$PACK -f $input  $OPTIONS -g $n > $output " >> $script_n;
    echo " " >> $script_n;

    qsub $script_n >> qsub.stdout;

done

rm -f tmp.awk



echo "Finished submitting jobs with qsub."
exit



