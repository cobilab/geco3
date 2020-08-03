#!/bin/bash
###############################################################################
INSTALL_NAF=0;
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
###############################################################################
RUN_NAF=0;
RUN_GECO2=1;
RUN_GECO3=1;
#
function RunGeCo2 {
   # 1 - params
   # 2 - seq
if test -f "../../ds/$2"; then
   cp ../../ds/$2 .
   rm -f $2.co
   (/usr/bin/time -v ./GeCo2 -v $1 $2 ) &> ../../res/C_GECO2_l16_$2
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
   (/usr/bin/time -v ./GeCo3 -v $1 $2 ) &> ../../res/C_GECO3_l16_$2
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
   (/usr/bin/time -v ./ENNAF --temp-dir /tmp -22 $1.fasta -o $1.naf ) &> ../../res/C_NAF_$1
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
   (/usr/bin/time -v ./XM --real=$1.xm $1.fasta ) &> ../../res/C_XM_$1
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
   (/usr/bin/time -v ./JARVIS $1 $2 ) &> ../../res/C_JARVIS_$2
   ls -la $2.jc | awk '{ print $5;}' > ../../res/BC_JARVIS_$2
   rm -f $2 $2.jc;
fi
}

function RunCMIX {
    # 1 - seq
if test -f "../../ds/$1"; then
   cp ../../ds/$1 .
   rm -f $1.cmix
   (/usr/bin/time -v ./cmix -c $1 $1.cmix ) &> ../../res/C_CMIX_$1
   ls -la $1.cmix | awk '{ print $5;}' > ../../res/BC_CMIX_$1
   rm -f $1 $1.cmix;
fi
}

function RunDEEPZIP {
    # 1 - seq
if test -f "../../ds/$1"; then
   source tf/bin/activate
   rm data/processed_files/*param.json data/processed_files/*npy
   rm -rf data/files_to_be_compressed/
   #rm -rf data/compressed/
   mkdir -p data/files_to_be_compressed/
   mkdir -p data/compressed/
   cp ../../ds/$1 data/files_to_be_compressed/$1
   cd data
   ./run_parser.sh
   cd ../src
   (/usr/bin/time -v ./run_experiments.sh biGRU ) &> ../../../res/C_DEEPZIP_$1
   cd ../data/compressed/
   ls -la $1/biGRU.compressed.combined | awk '{ print $5;}' > ../../../../res/BC_DEEPZIP_$1
   cd ../..
fi
}

function RunZPAQ {
    # 1 - seq
if test -f "../../ds/$1"; then
   cp ../../ds/$1 .
   rm -f $1.zpaq
   (/usr/bin/time -v ./zpaq a $1.zpaq $1  -m5 ) &> ../../res/C_ZPAQ_$1
   ls -la $1.zpaq | awk '{ print $5;}' > ../../res/BC_ZPAQ_$1
   rm -f $1 $1.zpaq;
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
    #params mentioned in paper
    sed -i 's/batch_size=128/batch_size=1024/g' src/trainer.py
    sed -i 's/num_epochs=20/num_epochs=3/g' src/trainer.py
    #don't decompress
    sed -i 's/\/usr\/bin\/time -v python decompressor.py/#\/usr\/bin\/time -v python decompressor.py/g' src/run_experiments.sh
    sed -i 's/   cmp $recon_file_name/    #cmp $recon_file_name/g' src/run_experiments.sh
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

    cd ..
fi
#
# RUN
#=========================================================================
###################################################################
#### Specific
if [[ "$RUN_NAF" -eq "1" ]]; then
   echo "Running NAF ...";
   mkdir -p res
   cd progs/naf

   RunNAF "virome"
   #
   cd ../../
   echo "Done!";
fi

###############################################################################
if [[ "$RUN_GECO2" -eq "1" ]]; then
   echo "Running GeCo2 ...";
   mkdir -p res
   cd progs/geco2
   PARAM_DS1="-tm 7:1:1:1:0.8/0:0:0 -tm 13:10:0:1:0.95/0:0:0 -tm 19:500:1:40:0.95/5:20:0.95 -ls 0.03 -hs 64"
   RunGeCo2 "$PARAM_DS1" "virome"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p res
   cd progs/geco3

   PARAM_DS1="-tm 7:1:1:1:0.8/0:0:0 -tm 13:10:0:1:0.95/0:0:0 -tm 19:500:1:40:0.95/5:20:0.95 -ls 0.03 -hs 64"
   RunGeCo3 "$PARAM_DS1" "virome"
   
   #
   cd ../../
   echo "Done!";
fi
