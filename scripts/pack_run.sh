#!/bin/bash

script_name="pack_run.sh"
script_version="2.9"
script_date="Tue Nov 28 09:56:13 PST 2006"

#script_version="2.8"
#script_date="Mon Aug 14 13:51:30 PDT 2006"

# location of pack executable
pack=~/bin/pack
pver=`$pack -h | grep version | awk '(NF==3){print $3}'`

NICE=19

# Changes:
#
# ver 2.8: added switch for xbsa and xmgrace so they don't start automatically
# ver 2.7: added calculation of pair dist function from poscar file
# ver 2.6: added extraction of restart file (26may2006)
# ver 2.5: added qposcar output file (see changes in pack.c 4.4.0.4)
# ver 2.4: nice value set to 19
# ver 2.3: updates movie file
#

if [ $# -eq 0 ]; then
    echo
    echo "########################"
    echo "#      "$script_name
    echo "########################"
    echo
    echo "version "$script_version
    echo $script_date
    echo
    echo "use: pack.sh  -f file.inp [-o outbase] [-O \"pack options\"] [-PSCVDQR]"
    echo
    echo "    -f --- input file name"
    echo "    -o --- output base name"
    echo "    -O --- put 'pack' options here in quotes"
    echo "    -V --- use pack-version number 3.x.x (ver 4.x.x is default)"
    echo "    -N --- specify the node to run on"
    echo "    -P --- post process the .bs file only"
    echo "    -S --- start xbsa and xmgrace to show xray and structure (must use with -P switch)"
    echo "    -C --- clean directory (delete all but the .bs file)"
    echo "    -D --- extract POSCAR file only"
    echo "    -R --- extract restart file only"
    echo "    -Q --- extract QPOSCAR file only, then run ewald"
    echo
    echo "This script will use pack version ("$pver
    echo
    echo " examples: "
    echo
    echo "$ pack_run -o naalh4_01 -f pack.tet.inp -O \"-r 100 -t 50 -I 50 -G 0.5 -Y 2.1\""
    echo "$ pack_run -o naalh4_01 -f pack.tet.inp -o jnk_01 -N 23 -O \"-r 100 -t 50 -I 50 -G 0.5 -Y 2.1\""
    echo "$ pack_run -o naalh4_01 -P"
    echo "$ pack_run -o naalh4_01 -C"
    echo "$ pack_run -o naalh4_01 -D"
    echo "$ pack_run -o naalh4_01 -Q"
    echo "$ pack_run -o naalh4_01 -R"
    echo " i.e. ... "
    echo "$ pack_run -o base_name -option"
    echo
    echo
    echo Eric Majzoub
    echo Sandia National Laboratories
    echo 22 March 2006
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
ver3=-1
posc=-1
qosc=-1
rsrt=-1
xbsa=0

declare SWITCH
while getopts "f:o:O:N:PCVDQRS" SWITCH; do
    case $SWITCH in
    f) name=$OPTARG ;;
#       outn=${OPTARG%%.inp} ;;
    o) outn=$OPTARG ;;
    O) allo=$OPTARG ;;
    N) node=$OPTARG ;;
    V) ver3=1 ;;
    P) post=1 ;;
    C) clen=1 ;;
    D) posc=1 ;;
    Q) qosc=1 ;;
    R) rsrt=1 ;;
    S) xbsa=1 ;;
    esac
done

# Check if we are on Mac OS X or Linux
onDarwin=`echo \`uname -a\` | grep "Darwin"`
onLinux=`echo \`uname -a\` | grep "Linux"`

# Check to see if the 80 column switch was
# given in the command line to pack
# must have a space before the minus sign, so
# that -e 1e-8 does not turn it on, hence the [" "]
eightcolflag=`echo "$allo" | grep [" "][\-][8]`
if [ ! -z "$eightcolflag" ]; then
    ecol=1
fi

trunk=${name##/*/}
base=${name%/*}
extension=${name##*.}
bs=$outn.bs
symin=$outn.sym.in
symout=$outn.sym.out
symlog=$outn.sym.log
fsymin=$outn.fsym.in
fsymout=$outn.fsym.out
fsymlog=$outn.fsym.log
xrin=$outn.xr.in
xrout=$outn.xr.dat
mvfile=$outn.mv
poscar=$outn.poscar
qposcar=$outn.qposcar
dxlout=$outn.str
restart=$outn.restart
rdf=$outn.rdf

hostn=`hostname`
hostbase=${hostn%%.*}

# echo trunk = $trunk
# echo base = $base
# echo extension = $extension
# echo bs = $bs
# echo symin = $symin
# echo symout = $symout
# echo xrin = $xrin
# echo xrout = $xrout

# echo $runs
# echo $trun
# echo $epsi
# echo $othr

#
# [test] expects integer expressions, for those not integer, quote them
# and use the string comparisons (!= is for strings, -ne is for integers).
#
if [ "$allo" != -1 ]; then
    opts=$opts" "$allo
fi

#################################
#   Set Environment Variables   #
#################################
if [ $hostbase = "master" -o $hostbase = "kelp" ] ; then
    newhome=`echo $HOME`
    export ISODATA=$newhome/bin/isotropy_2005/
    syms=$newhome/bin/isotropy_2005/symsearch
    findsym=$newhome/bin/isotropy_2005/findsym
    grace=`which xmgrace`
    nice=`which nice`
    powder=$newhome/bin/powder
else
    syms=/home/packages/binaries/isotropy_2005/symsearch
    findsym=/home/packages/binaries/isotropy_2005/findsym
#    grace=`locate xmgrace | grep 5.1.18 | grep -v '\.c' | grep -v '\.o'`
    grace=`which xmgrace`
    nice=`which nice`
    powder=/home/ehm/bin/powder
    ewald=/home/ehm/bin/ewald
    export ISODATA=/home/packages/binaries/isotropy_2005/
fi


###################################
#                                 #
#         Functions               #
#                                 #
###################################
function movie {

    if [ -f out.mv ]; then
	mv -f out.mv $mvfile;
    fi

}

function postprocess {

    echo "Parsing output bs file: "$bs ;
    cat $bs | awk '($1=="*F"){print $0}' | cut -d" " -f 2-300 > $fsymin ;
#    cat $bs | awk '($1=="*S"){print $0}' | cut -d" " -f 2-300 > $symin ;
    cat $bs | awk '($1=="*D"){print $0}' | cut -d" " -f 2-300 > $xrin ;
    cat $bs | awk '($1=="*P"){print $0}' | cut -d" " -f 2-300 > $poscar ;
    cat $bs | awk '($1=="*Q"){print $0}' | cut -d" " -f 2-300 > $qposcar ;
    cat $bs | awk '($1=="*R"){print $0}' | cut -d" " -f 2-300 > $restart ;
    
    cat $poscar | awk -f ~/awkfiles/poscar2vis.awk -v DXL=1 > $dxlout ;

    if [ -x `which powder` ]; then
	echo "########################################"
	echo "Generating powder diffraction pattern"
	echo "########################################"
	$powder -f $xrin -X -t 1.5405 -L -C > $xrout ;
    else
	echo "Powder executable not found"
    fi

    if [ -x `which contcar_pdf` ]; then
	echo "########################################"
	echo "Generating radial distribution function"
	echo "########################################"
	contcar_pdf -f $poscar > $rdf ;
    else
	echo "contcar_pdf executable not found"
    fi

    if [ -n "$onLinux" ]; then
	echo "########################################"
	echo "Running findsym"
	echo "########################################"
	$findsym < $fsymin > $fsymout ;
    elif [ -n "$onDarwin" ]; then
	echo "Findsym not available for OS X"
    else
	echo "Findsym executable not found"
    fi


    if [ $xbsa -eq 1 ]; then
    	if [ -x `which xbsa` ]; then
	    xbsa $bs &
	else
	    echo "xbsa executable not found"
	fi
	
	if [ -x $grace ]; then
	    $grace $xrout -free -geometry +592+0 &
	else
	    echo "grace executable not found"
	fi
    fi

    movie ;

    if [ -f symsearch.log ]; then
	mv -f symsearch.log $symlog
    fi
    
    if [ -f findsym.log ]; then
	mv -f findsym.log $fsymlog
    fi
}


##############################
#       BEGIN HERE           #
##############################

if [ $rsrt -eq 1 ]; then
    echo "# -----> Extracting restart.dat from file "$bs
    cat $bs | awk '($1=="*R"){print $0}' | cut -d" " -f 2-300 > restart.dat ;
    echo "# Done."
    exit
fi

if [ $posc -eq 1 ]; then
    echo "# ------> Extracting POSCAR from file "$bs
    cat $bs | awk '($1=="*P"){print $0}' | cut -d" " -f 2-300 > POSCAR ;
    echo "# Done."
    exit
fi

if [ $qosc -eq 1 ]; then
    echo "# ------> Extracting QPOSCAR from file "$bs
    cat $bs | awk '($1=="*Q"){print $0}' | cut -d" " -f 2-300 > QPOSCAR ;
    $ewald ;
    echo "# ------> Removing QPOSCAR file"
    echo "# Done."
    rm -f QPOSCAR ;
    exit
fi


if [ $clen -eq 1 ]; then
    rm -f $symin $symout $symlog $fsymin $fsymout $fsymlog $xrin $xrout $mvfile $poscar $restart $rdf $qposcar $dxlout;
    exit
fi

if [ $post -eq 1 ]; then
    postprocess ;
    exit
fi

if [ -f $mvfile ]; then
    rm -f $mvfile;
fi

if [ $ver3 == 1 ]; then
    pack=~/bin/pack_v3;
fi

# check if we are running on master, then farm out the work to a node
if [ $hostbase = "master" ]; then

    # where is the executable, and set up the directories
    # that we need to cd to

    cdir=`pwd`
    cdto=pack/${cdir##/*/}

    if [ $node -eq -1 ]; then
	echo
	echo "Looks like we are running on 'master'."
	echo "Node switch is not set, will use k_machines."
	echo "Using rsh... will farm out to first node found."
	echo

	echo "running k_machines..."
	k_machines 1 ;
	mach=(`cat machines`) ;
	echo
	echo "rsh will change directory to " $cdto " before running pack"
	rsh -n $mach "cd $cdto; $nice -n $NICE $pack -f $name $opts > ~/$cdto/$bs " &
	echo "Running pack...   job submitted to "$mach
	echo
	echo `date` $mach $bs>> jobs
	echo
    fi

    if [ $node -ne -1 ]; then
	echo
	echo "Looks like we are running on 'master'."
	echo "Node switch is set, will NOT use k_machines."
	echo
	echo "rsh will change directory to " $cdto " before running pack"
	rsh -n node$node "cd $cdto; $nice -n $NICE $pack -f $name $opts > ~/$cdto/$bs " &
	echo "Running pack...   job submitted to node"$node
	echo
	echo `date` node$node $bs>> jobs
	echo
    
    fi


    exit

else

    # if we are running locally, then nice the process and run it
    echo "Running pack locally (not on a node)..."
    if [ $ecol = -1 ]; then
	$nice -n $NICE $pack -f $name $opts > $bs & xterm -geometry 145x35 -e tail -f $bs ;
	movie ;
	xbsa $bs &
	exit
    fi

    if [ $ecol =  1 ]; then
	$nice -n $NICE $pack -f $name $opts > $bs & xterm -geometry 80x35 -e tail -f $bs ;
	movie ;
	xbsa $bs &
	exit
    fi

fi
