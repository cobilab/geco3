#!/bin/bash
###############################################################################
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
###############################################################################
RUN_GECO2=1;
RUN_GECO3=1;
###############################################################################
#
function RunGeCo2 {
   # 1 - params
   # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (/usr/bin/time -v ./GeCo2 $1 $2 ) &> ../../res/C_GECO2_LO_$2
   ls -la $2.co | awk '{ print $5;}' > ../../res/BC_GECO2_LO_$2
   #rm -f $2 $2.co;
fi
}

function RunGeCo3 {
    # 1 - params
    # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (/usr/bin/time -v ./GeCo3 $1 $2 ) &> ../../res/C_GECO3_LO_$2
   ls -la $2.co | awk '{ print $5;}' > ../../res/BC_GECO3_LO_$2
   #rm -f $2 $2.co;
fi
}
###############################################################################
# INSTALL
mkdir -p ds
mkdir -p progs
cd progs/
###############################################################################
# GET GECO2
#====================================================================
if [[ "$INSTALL_GECO2" -eq "1" ]]; then
   rm -fr geco2/
   git clone https://github.com/cobilab/geco2.git
   cd geco2/src/
   cp Makefile.linux Makefile
   sed -i 's/U32 garbage;//g' defs.h # fix for gcc 10
   sed -i 's/garbage                   =//g' gede2.c # fix for gcc 10
   make
   cp GeCo2 ../
   cp GeDe2 ../
   cd ../../
fi
###############################################################################
# GET GECO3
#====================================================================
if [[ "$INSTALL_GECO3" -eq "1" ]]; then
   rm -fr geco3/
   git clone https://github.com/cobilab/geco3.git
   cd geco3/src/
   make
   cp GeCo3 ../
   cp GeDe3 ../
   cd ../../
fi
##############################################################################
cd ..
###############################################################################
# DOWNLOAD
if [[ "$DOWNLOAD" -eq "1" ]]; then
    echo "Downloading ..."
    rm -rf ds/
    mkdir ds/
    cd ds
    wget https://tinyurl.com/DNAcorpus
    unzip DNAcorpus
    mv DNACorpus/* .
    rm -f DNAcorpus 
    rm -rf DNACorpus
    cd ..
fi
#
# RUN
#=========================================================================
if [[ "$RUN_GECO2" -eq "1" ]]; then
   echo "Running GeCo2 ...";
   mkdir -p res
   cd progs/geco2
   RunGeCo2 "-l 1 -lr 0.06 -hs 8" "BuEb"
   RunGeCo2 "-l 2 -lr 0.06 -hs 16" "AgPh"
   RunGeCo2 "-l 3 -lr 0.09 -hs 24" "YeMi"
   RunGeCo2 "-l 4 -lr 0.04 -hs 40" "HePy"
   RunGeCo2 "-l 5 -lr 0.04 -hs 16" "AeCa"
   RunGeCo2 "-l 5 -lr 0.04 -hs 40" "HaHi"
   RunGeCo2 "-l 6 -lr 0.03 -hs 40" "EsCo"
   RunGeCo2 "-l 7 -lr 0.03 -hs 40" "PlFa"
   RunGeCo2 "-l 8 -lr 0.03 -hs 40" "ScPo"
   RunGeCo2 "-l 9 -lr 0.05 -hs 64" "EnIn"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "DrMe"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "OrSa"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "DaRe"
   RunGeCo2 "-l 11 -lr 0.03 -hs 64" "GaGa"
   RunGeCo2 "-l 12 -lr 0.03 -hs 64" "HoSa"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p res
   cd progs/geco3
   RunGeCo3 "-l 1 -lr 0.06 -hs 8" "BuEb"
   RunGeCo3 "-l 2 -lr 0.06 -hs 16" "AgPh"
   RunGeCo3 "-l 3 -lr 0.09 -hs 24" "YeMi"
   RunGeCo3 "-l 4 -lr 0.04 -hs 40" "HePy"
   RunGeCo3 "-l 5 -lr 0.04 -hs 16" "AeCa"
   RunGeCo3 "-l 5 -lr 0.04 -hs 40" "HaHi"
   RunGeCo3 "-l 6 -lr 0.03 -hs 40" "EsCo"
   RunGeCo3 "-l 7 -lr 0.03 -hs 40" "PlFa"
   RunGeCo3 "-l 8 -lr 0.03 -hs 40" "ScPo"
   RunGeCo3 "-l 9 -lr 0.05 -hs 64" "EnIn"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "DrMe"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "OrSa"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "DaRe"
   RunGeCo3 "-l 11 -lr 0.03 -hs 64" "GaGa"
   RunGeCo3 "-l 12 -lr 0.03 -hs 64" "HoSa"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
