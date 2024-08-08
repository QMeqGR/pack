#!/bin/sh

echo "rotate_anion.c" >  trouble_files.tmp;
echo "cent_of_mass.c" >> trouble_files.tmp;

make clean;

for file in `cat trouble_files.tmp`; do
    echo "Compiling "$file " without optimization"
    icc -c $file;
done

make icc=core2;
make install;

rm trouble_files.tmp
exit
