#!/bin/bash
###############################################################################
INSTALL_NAF=0;
INSTALL_XM=0;
INSTALL_JARVIS=0;
INSTALL_CMIX=0;
INSTALL_DEEPZIP=0;
INSTALL_ZPAQ=0;
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
###############################################################################
RUN_NAF=0;
RUN_XM=0;
RUN_JARVIS=0;
RUN_CMIX=0;
RUN_DEEPZIP=0;
RUN_ZPAQ=0;
RUN_GECO2=0;
RUN_GECO3=0;
#
function RunGeCo2 {
   # 1 - params
   # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (time ./GeCo2 $1 $2 ) &> ../../res/C_GECO2_l16_$2
   ls -la $2.co | awk '{ print $5;}' > ../../res/BC_GECO2_l16_$2
   rm -f $2 $2.co;
fi
}

function RunGeCo3 {
    # 1 - params
    # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (time ./GeCo3 $1 $2 ) &> ../../res/C_GECO3_l16_$2
   ls -la $2.co | awk '{ print $5;}' > ../../res/BC_GECO3_l16_$2
   rm -f $2 $2.co;
fi
}

function RunNAF {
    # 1 - seq
if test -f "../../ds/$1"; then
   echo ">" > $1.fasta
   cat ../../ds/$1 >> $1.fasta
   rm -f $1.fasta.naf
   (time ./ENNAF --temp-dir /tmp -22 $1.fasta -o $1.naf ) &> ../../res/C_NAF_$1
   ls -la $1.naf | awk '{ print $5;}' > ../../res/BC_NAF_$1
   rm -f $1.fasta $1.naf;
fi
}

function RunXM {
    # 1 - seq
if test -f "../../ds/$1"; then
   echo ">" > $1.fasta
   cat ../../ds/$1 >> $1.fasta

   rm -f $1.xm
   (time ./XM --real=$1.xm $1.fasta ) &> ../../res/C_XM_$1
   ls -la $1.xm | awk '{ print $5;}' > ../../res/BC_XM_$1
   rm -f $1.fasta $1.xm;
fi
}

function RunJARVIS {
    # 1 - params
    # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.jc
   (time ./JARVIS $1 $2 ) &> ../../res/C_JARVIS_$2
   ls -la $2.jc | awk '{ print $5;}' > ../../res/BC_JARVIS_$2
   rm -f $2 $2.jc;
fi
}

#
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
###############################################################################
# GET NAF
#====================================================================
if [[ "$INSTALL_NAF" -eq "1" ]]; then
   rm -fr naf/
   git clone --recurse-submodules https://github.com/KirillKryukov/naf.git
   cd naf && make
   cp ennaf/ennaf ENNAF # to avoid dir name collision
   cd ../
fi
###############################################################################
# GET XM
#====================================================================
if [[ "$INSTALL_XM" -eq "1" ]]; then
    rm -fr japsa/
    git clone https://github.com/mdcao/japsa.git
    cd japsa
    make install INSTALL_DIR=./xm MXMEM=7000m SERVER=true JLP=/usr/lib/jni
    cp xm/bin/jsa.xm.compress XM
    cd ../
fi
###############################################################################
# GET JARVIS
#====================================================================
if [[ "$INSTALL_JARVIS" -eq "1" ]]; then
    rm -fr jarvis/
    git clone https://github.com/cobilab/jarvis.git
    cd jarvis/src/
    sed -i 's/uint32_t garbage;//g' defs.h # fix for gcc 10
    make
    cp JARVIS ../
    cd ../../
fi
###############################################################################
# GET CMIX
#====================================================================
if [[ "$INSTALL_CMIX" -eq "1" ]]; then
    rm -fr cmix/
    git clone https://github.com/byronknoll/cmix.git
    cd cmix
    make
    cd ../
fi
###############################################################################
# GET DEEPZIP
#====================================================================
if [[ "$INSTALL_DEEPZIP" -eq "1" ]]; then
    rm -fr DeepZip/
    git clone --single-branch --branch noGPU https://github.com/mohit1997/DeepZip.git
    cd DeepZip
    python3.6 -m venv tf
    source tf/bin/activate
    bash install.sh
    deactivate
    cd ../
fi
###############################################################################
# GET ZPAQ
#====================================================================
if [[ "$INSTALL_ZPAQ" -eq "1" ]]; then
    rm -fr zpaq/
    git clone https://github.com/zpaq/zpaq.git
    cd zpaq
    g++ -O3 -march=native -Dunix zpaq.cpp libzpaq.cpp -pthread -o zpaq
    cd ../
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

    # ds4
    wget https://tinyurl.com/DNAcorpus
    unzip DNAcorpus
    mv DNACorpus/* .
    rm -f DNAcorpus 
    rm -rf DNACorpus

    #ds1
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/695/GCA_900067695.1_Pabies01/GCA_900067695.1_Pabies01_genomic.fna.gz
    
    cd ..
fi
#
# RUN
#=========================================================================
if [[ "$RUN_GECO2" -eq "1" ]]; then
    echo "Running GeCo2 ...";
    # l 16 equiv
    PARAM_DS4="-tm 1:1:0:0:0.9/0:0:0 -tm 2:1:0:0:0.78/0:0:0 -tm 3:1:0:0:0.9/0:0:0 -tm 4:1:1:0:0.78/0:0:0 -tm 5:10:1:1:0.90/0:0:0 -tm 6:1:1:0:0.85/0:0:0 -tm 7:1:1:0:0.85/0:0:0 -tm 8:1:1:0:0.91/0:0:0 -tm 9:10:0:0:0.9/0:0:0 -tm 10:10:0:0:0.9/0:0:0 -tm 11:10:0:0:0.9/0:0:0 -tm 12:20:1:1:0.94/0:0:0 -tm 13:10:1:0:0.95/0:0:0 -tm 14:50:1:1:0.95/0:0:0 -tm 16:200:1:10:0.95/1:50:0.95 -tm 17:100:1:20:0.9/3:10:0.9 -tm 20:500:1:30:0.95/2:20:0.95"
    
   mkdir -p res
   cd progs/geco2

   #ds4
   RunGeCo2 "-l 1 -lr 0.06 -hs 8" "BuEb"
   RunGeCo2 "-l 2 -lr 0.06 -hs 16" "AgPh"
   RunGeCo2 "-l 3 -lr 0.09 -hs 24" "YeMi"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "HePy"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "AeCa"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "HaHi"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "EsCo"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "PlFa"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "ScPo"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "EnIn"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "DrMe"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "OrSa"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "DaRe"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "GaGa"
   RunGeCo2 "$PARAM_DS4 -lr 0.03 -hs 64" "HoSa"

   #ds2 and ds3
   PARAM_DS23="-tm 3:1:1:1:0.8/0:0:0 -tm 6:1:1:1:0.85/0:0:0 -tm 9:1:1:1:0.85/0:0:0 -tm 12:10:0:1:0.85/0:0:0 -tm 15:200:1:10:0.85/2:1:0.85 -tm 17:200:1:10:0.85/2:1:0.85 -tm 20:500:1:40:0.85/5:20:0.85 -lr 0.03 -hs 64"

   RunGeCo2 PARAM_DS23 "HoSaY"
   RunGeCo2 PARAM_DS23 "Mitochondrion"

   RunGeCo2 PARAM_DS23 "VDB"
   RunGeCo2 PARAM_DS23 "Archaea"

   #ds1
   PARAM_DS1="-tm 3:1:1:1:0.70/0:0:0 -tm 8:1:1:1:0.85/0:0:0 -tm 13:10:0:1:0.85/0:0:0 -tm 19:500:1:40:0.85/5:20:0.85 -ls 0.03 -hs 64"
   RunGeCo2 PARAM_DS1 "GoGoC"
   RunGeCo2 PARAM_DS1 "PaTrC"
   RunGeCo2 PARAM_DS1 "HoSaC"
   RunGeCo2 PARAM_DS1 "PiAbC"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p res
   cd progs/geco3

   #ds4
   RunGeCo3 "-l 1 -lr 0.06 -hs 8" "BuEb"
   RunGeCo3 "-l 2 -lr 0.06 -hs 16" "AgPh"
   RunGeCo3 "-l 3 -lr 0.09 -hs 24" "YeMi"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "HePy"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "AeCa"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "HaHi"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "EsCo"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "PlFa"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "ScPo"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "EnIn"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "DrMe"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "OrSa"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "DaRe"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "GaGa"
   RunGeCo3 "-l 16 -lr 0.03 -hs 64" "HoSa"

   #ds2 and ds3
   PARAM_DS23="-tm 3:1:1:1:0.8/0:0:0 -tm 6:1:1:1:0.85/0:0:0 -tm 9:1:1:1:0.85/0:0:0 -tm 12:10:0:1:0.85/0:0:0 -tm 15:200:1:10:0.85/2:1:0.85 -tm 17:200:1:10:0.85/2:1:0.85 -tm 20:500:1:40:0.85/5:20:0.85 -lr 0.03 -hs 64"

   RunGeCo3 PARAM_DS23 "HoSaY"
   RunGeCo3 PARAM_DS23 "Mitochondrion"

   RunGeCo3 PARAM_DS23 "VDB"
   RunGeCo3 PARAM_DS23 "Archaea"

   #ds1
   PARAM_DS1="-tm 3:1:1:1:0.70/0:0:0 -tm 8:1:1:1:0.85/0:0:0 -tm 13:10:0:1:0.85/0:0:0 -tm 19:500:1:40:0.85/5:20:0.85 -ls 0.03 -hs 64"
   RunGeCo3 PARAM_DS1 "GoGoC"
   RunGeCo3 PARAM_DS1 "PaTrC"
   RunGeCo3 PARAM_DS1 "HoSaC"
   RunGeCo3 PARAM_DS1 "PiAbC"
   
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_NAF" -eq "1" ]]; then
   echo "Running NAF ...";
   mkdir -p res
   cd progs/naf

   #ds4
   RunNAF "BuEb"
   RunNAF  "AgPh"
   RunNAF  "YeMi"
   RunNAF  "HePy"
   RunNAF  "AeCa"
   RunNAF  "HaHi"
   RunNAF  "EsCo"
   RunNAF  "PlFa"
   RunNAF  "ScPo"
   RunNAF  "EnIn"
   RunNAF  "DrMe"
   RunNAF  "OrSa"
   RunNAF  "DaRe"
   RunNAF  "GaGa"
   RunNAF  "HoSa"

   #ds2 and ds3
   RunNAF  "HoSaY"
   RunNAF  "Mitochondrion"

   RunNAF  "VDB"
   RunNAF  "Archaea"

   #ds1
   RunNAF "GoGoC"
   RunNAF "PaTrC"
   RunNAF "HoSaC"
   RunNAF "PiAbC"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_XM" -eq "1" ]]; then
   echo "Running XM ...";
   mkdir -p res
   cd progs/japsa

   #ds4
   RunXM "BuEb"
   RunXM  "AgPh"
   RunXM  "YeMi"
   RunXM  "HePy"
   RunXM  "AeCa"
   RunXM  "HaHi"
   RunXM  "EsCo"
   RunXM  "PlFa"
   RunXM  "ScPo"
   RunXM  "EnIn"
   RunXM  "DrMe"
   RunXM  "OrSa"
   RunXM  "DaRe"
   RunXM  "GaGa"
   RunXM  "HoSa"

   #ds2 and ds3
   RunXM  "HoSaY"
   RunXM  "Mitochondrion"

   RunXM  "VDB"
   RunXM  "Archaea"

   #ds1
   RunXM "GoGoC"
   RunXM "PaTrC"
   RunXM "HoSaC"
   RunXM "PiAbC"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
