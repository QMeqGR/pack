#!/bin/sh

# E. Majzoub 07 Dec 2005
# v 1.1

# file defs
src=./src
main=$src/pack/main.c

#########################
#  BEGIN

cd $HOME

major=`cat $main | grep VERSION_MAJOR | grep define | cut -d " " -f 3`
minor=`cat $main | grep VERSION_MINOR | grep define | cut -d " " -f 3`
bugfx=`cat $main | grep VERSION_BUGFX | grep define | cut -d " " -f 3`
trivl=`cat $main | grep VERSION_TRIVL | grep define | cut -d " " -f 3`

gzipdir=pack-$major.$minor.$bugfx.$trivl

# Check if we are on Mac OS X
onDarwin=`echo \`uname -a\` | grep "Darwin"`
if [ -n "$onDarwin" ]; then
mach=darwin_ppc
fi

# Check if we are on Linux
onLinux=`echo \`uname -a\` | grep "Linux"`
if [ -n "$onLinux" ]; then
mach=linux_x86
fi

cd $src
if [ -d $gzipdir ]; then
    echo "packaging directory " $gzipdir " already exists!"
    echo "exiting."
    exit
fi
if [ ! -d $gzipdir ]; then
    mkdir $gzipdir
    cp -ra ./pack/* $gzipdir
    rm -f $gzipdir/*~
fi

# tar cvzf pack-$major.$minor.$bugfx-$mach.tgz ./pack/
tar cvzf pack-$major.$minor.$bugfx.$trivl.tgz $gzipdir
rm -rf $gzipdir

# mv pack-$major.$minor.$bugfx-$mach.tgz ~/
mv pack-$major.$minor.$bugfx.$trivl.tgz ~/

exit

