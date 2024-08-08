#!/bin/sh

# automatically run a program with a set of different files
# ehm
# SNL/CA 21 Dec 2005

# version 2.0  16 jan 2006

directory=`pwd`
echo "Running in directory "$directory
echo

function get_node {
    ~/src/pack_machines g;
    node=`cat machines | awk '(NF>0){ printf("%d",substr( $0 , match($0,/node/)+4 ) ); }'`;
}

for input in `ls $directory | grep .inp`; do

    bsfile=${input%%.inp}.bs;
#    echo "bsfile = "$bsfile;
    if [ -e $bsfile ]; then
        echo "Found file "$bsfile
    fi

    # if there is not already a .bs file, then run pack on the .inp file
    if [ ! -e $bsfile ]; then

        get_node

        while [ -z $node ]; do
            echo "no nodes: sleeping for 30 minutes"
            sleep 1800  
            get_node
        done

        echo Running $input on node = $node
        pack_run -f $directory/$input -N$node -o ${input%%.inp} -O "-r 200 -t 4 -B 50" ;
        echo
        rm -f machines

    fi

done
