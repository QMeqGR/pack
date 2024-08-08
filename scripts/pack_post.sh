#!/bin/bash

dir=$1;
file=$2;

if [ -z "$dir" ] && [ -z "$file" ]; then
    echo "use: pack_post.sh dirname filename.bs"
    exit
fi

echo "using dir= "$dir;
echo "using file= "${file%%.bs};

mkdir $dir;
cd $dir;
cp ../$file .
pack_run -f ../input.in -o ${file%%.bs} -D ;
~/src/VASP_setup.sh -D ~/../rrstump/Dvasp/pawpot_PBE/ -k 4
cd ..

exit