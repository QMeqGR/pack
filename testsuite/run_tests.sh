#!/bin/bash

PACK=../pack_exe

LIST=`ls *.inp`

for file in $LIST; do
    outfile=${file%%inp}bs
    echo "***************************************************"
    echo "testing file "$file
    echo "Without simplex"
    $PACK -f $file -V -W -r 5 -t 2 > $outfile;
    cat $outfile | grep dup | wc | awk '($1 != 3){printf("     Error on this file!\n");}';
    cat $outfile | grep dup | wc | awk '($1 == 3){printf("     Pass!\n");}';
    echo "With simplex"
    $PACK -f $file -V -q 5 -r 1 -t 2 > $outfile;
    cat $outfile | grep dup | wc | awk '($1 != 3){printf("     Error on this file!\n");}';
    cat $outfile | grep dup | wc | awk '($1 == 3){printf("     Pass!\n");}';
    echo "With DSM"
    $PACK -f $file -V -q 5 -r 1 -t 2 -B 1 > $outfile;
    cat $outfile | grep dup | wc | awk '($1 != 3){printf("     Error on this file!\n");}';
    cat $outfile | grep dup | wc | awk '($1 == 3){printf("     Pass!\n");}';
    echo
    echo
done

echo "finished tests."
# rm -f *.bs
rm -f PR*
exit;
