#!/bin/sh
# start bialk runs

# location of pack executable
pack=~/bin/pack

if [ $# -eq 0 ]; then
    echo
    echo "##########################"
    echo "#      pack_batch.sh     #"
    echo "##########################"
    echo
    echo " pack_batch version 1.0"
    echo
    echo "use: pack_batch.sh -o basename [-O \"pack options\"] [-FNnV]"
    echo
    echo "    -o --- output base name"
    echo "    -O --- put 'pack' options here in quotes"
    echo "    -V --- use pack-version number 2.x.x (ver 3.x.x is default)"
    echo "    -N --- specify the node to start runs on"
    echo "    -F --- specify the file count to start from"
    echo "    -m --- specify the number of files to generate"
    echo "    -n --- specify the number of processes per node"
    echo
    echo " example: "
    echo
    echo "$ pack_batch -o naalh4_0 -F 20 -m 18 -N 48 -n 3 -O \"-r 100 -t 50 -I 50 -G 0.5 -Y 2.1\""
    echo
    echo
    echo Eric Majzoub
    echo Sandia National Laboratories
    echo 11 october 2005
    echo
    exit
fi

# default values
opts=""
allo=-1
outn=tmp
post=-1
ecol=-1
clen=-1
node=-1
ver2=-1

# base file name
base=cabh
# formula units in file
fu=2
# processor count (goes down from here)
pc=48
# file count (goes up from here)
fc=20
fn=0

declare SWITCH
while getopts "o:O:N:VF:m:n:" SWITCH; do
    case $SWITCH in
    o) outn=$OPTARG ;;
    O) allo=$OPTARG ;;
    N) node=$OPTARG ;;
    V) ver2=1 ;;
    F) fc=$OPTARG ;;
    m) fn=$OPTARG ;;
    n) pc=$OPTARG ;;
    C) clen=1 ;;
    esac
done

# Check if we are on Mac OS X or Linux
onDarwin=`echo \`uname -a\` | grep "Darwin"`
onLinux=`echo \`uname -a\` | grep "Linux"`

for x in `echo $fc $fn | awk '{start=$1; for(i=0;i<$2;i++){printf("%d ",start++)}}'`
do
    echo pack_run -f pack.bialk_2fu.inp -N 48 -o bialk_$x -O "-r 100000 -t 25 -I 5 -J 0.02 -Y 1.05" &
done

