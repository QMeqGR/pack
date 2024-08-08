#! /bin/bash
#
# 'pack' data converter.
#
# Takes as input files:
#
#       1) pack output files
#
# Eric Majzoub
# Sandia National Laboratories
# 14 June 2005
#
# Last modified
#

name=$1

trunk=${name##/*/}
base=${name%/*}
extension=${name##*.}

# echo trunk = $trunk
# echo base = $base
# echo extension = $extension

namedat=${trunk%%bs}dat

switch=$2

if [ $# -eq 0 ]; then
    echo
    echo "####################"
    echo "#   packcnvrt.sh   #"
    echo "####################"
    echo
    echo "Version 1.0"
    echo 
    echo Eric Majzoub
    echo Sandia National Laboratories
    echo originally written 2005
    echo
    echo Converts:
    echo "1) pack files"
    echo
    echo Usage: packcnvrt.sh fname [r]
    echo
    echo " pack file conversion (expects .bs):"
    echo
    echo "    Abscissa: Temp vs."
    echo "        -r --- regular"
    echo 
    exit
fi

if [ -z $switch ]; then
    echo
    echo You must enter a choice for output!
    echo
    exit;
fi

if [ "bs" = $extension ]; then
    if [ 'r' = $switch ]; then
    cat $name | awk -f ~/awkfiles/pack.awk > $namedat
    fi
else
    echo
    echo "Extension does correspond to chosen output option."
    echo
fi
