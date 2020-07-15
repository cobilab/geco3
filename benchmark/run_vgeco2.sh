#!/bin/bash
###############################################################################
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
###############################################################################
RUN_GECO2=0;
RUN_GECO3=0;
###############################################################################
#
function RunGeCo2 {
   # 1 - params
   # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (time ./GeCo2 $1 $2 ) &> ../../res/C_GECO2_$2
   ls -la $1.co | awk '{ print $5;}' > ../../res/BC_GECO2_$2
   rm -f $2 $2.co;
fi
}

function RunGeCo3 {
    # 1 - params
    # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (time ./GeCo3 $2 $1 ) &> ../../res/C_GECO3_$2
   ls -la $1.co | awk '{ print $5;}' > ../../res/BC_GECO3_$2
   rm -f $2 $2.co;
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
    cd ds
    wget https://tinyurl.com/DNAcorpus
    unzip DNAcorpus
    cd ..
fi
#
# RUN
#=========================================================================
if [[ "$RUN_GECO2" -eq "1" ]]; then
   echo "Running GeCo2 ...";
   mkdir -p res
   cd progs/geco2
   RunGeCo2 "-l 1 -lr 0.06 -hs 8" "DNACorpus/BuEb"
   RunGeCo2 "-l 2 -lr 0.06 -hs 16" "DNACorpus/AgPh"
   RunGeCo2 "-l 3 -lr 0.09 -hs 24" "DNACorpus/YeMi"
   RunGeCo2 "-l 4 -lr 0.04 -hs 40" "DNACorpus/HePy"
   RunGeCo2 "-l 5 -lr 0.04 -hs 16" "DNACorpus/AeCa"
   RunGeCo2 "-l 5 -lr 0.04 -hs 40" "DNACorpus/HaHi"
   RunGeCo2 "-l 6 -lr 0.03 -hs 40" "DNACorpus/EsCo"
   RunGeCo2 "-l 7 -lr 0.03 -hs 40" "DNACorpus/PlFa"
   RunGeCo2 "-l 8 -lr 0.03 -hs 40" "DNACorpus/ScPo"
   RunGeCo2 "-l 9 -lr 0.05 -hs 64" "DNACorpus/EnIn"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "DNACorpus/DrMe"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "DNACorpus/OrSa"
   RunGeCo2 "-l 10 -lr 0.03 -hs 64" "DNACorpus/DaRe"
   RunGeCo2 "-l 11 -lr 0.03 -hs 64" "DNACorpus/GaGa"
   RunGeCo2 "-l 12 -lr 0.03 -hs 64" "DNACorpus/HoSa"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p res
   cd progs/geco3
   RunGeCo3 "-l 1 -lr 0.06 -hs 8" "DNACorpus/BuEb"
   RunGeCo3 "-l 2 -lr 0.06 -hs 16" "DNACorpus/AgPh"
   RunGeCo3 "-l 3 -lr 0.09 -hs 24" "DNACorpus/YeMi"
   RunGeCo3 "-l 4 -lr 0.04 -hs 40" "DNACorpus/HePy"
   RunGeCo3 "-l 5 -lr 0.04 -hs 16" "DNACorpus/AeCa"
   RunGeCo3 "-l 5 -lr 0.04 -hs 40" "DNACorpus/HaHi"
   RunGeCo3 "-l 6 -lr 0.03 -hs 40" "DNACorpus/EsCo"
   RunGeCo3 "-l 7 -lr 0.03 -hs 40" "DNACorpus/PlFa"
   RunGeCo3 "-l 8 -lr 0.03 -hs 40" "DNACorpus/ScPo"
   RunGeCo3 "-l 9 -lr 0.05 -hs 64" "DNACorpus/EnIn"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "DNACorpus/DrMe"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "DNACorpus/OrSa"
   RunGeCo3 "-l 10 -lr 0.03 -hs 64" "DNACorpus/DaRe"
   RunGeCo3 "-l 11 -lr 0.03 -hs 64" "DNACorpus/GaGa"
   RunGeCo3 "-l 12 -lr 0.03 -hs 64" "DNACorpus/HoSa"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
