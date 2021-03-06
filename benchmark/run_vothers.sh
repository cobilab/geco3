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
RUN_CMIX=1;
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
## Gen purpose
if [[ "$RUN_CMIX" -eq "1" ]]; then
   echo "Running CMIX ...";
   mkdir -p res
   cd progs/cmix

   #ds4
   #RunCMIX  "BuEb"
   #RunCMIX  "AgPh"
   #RunCMIX  "YeMi"
   #RunCMIX  "HePy"
   #RunCMIX  "AeCa"
   #RunCMIX  "HaHi"
   #RunCMIX  "EsCo"
   #RunCMIX  "PlFa"
   #RunCMIX  "ScPo"
   #RunCMIX  "EnIn"
   RunCMIX  "HoSaY"
   RunCMIX  "DrMe"
   RunCMIX  "OrSa"
   RunCMIX  "DaRe"
   RunCMIX  "GaGa"
   RunCMIX  "HoSa"

   #ds2 and ds3
   #RunCMIX  "HoSaY"
   #RunCMIX  "Mitochondrion"

   #RunCMIX  "VDB"
   #RunCMIX  "Archaea"

   #ds1
   #RunCMIX "GoGoC"
   #RunCMIX "PaTrC"
   #RunCMIX "HoSaC"
   #RunCMIX "PiAbC"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_DEEPZIP" -eq "1" ]]; then
   echo "Running DEEPZIP ...";
   mkdir -p res
   cd progs/DeepZip

   #ds4 System reboots with larger files (>= AeCa)
   #RunDEEPZIP  "BuEb"
   #RunDEEPZIP  "AgPh"
   #RunDEEPZIP  "YeMi"
   #RunDEEPZIP  "HePy"
   #RunDEEPZIP  "AeCa"
   RunDEEPZIP  "HaHi"
   RunDEEPZIP  "EsCo"
   RunDEEPZIP  "PlFa"
   RunDEEPZIP  "ScPo"
   RunDEEPZIP  "EnIn"
   #RunDEEPZIP  "DrMe"
   #RunDEEPZIP  "OrSa"
   #RunDEEPZIP  "DaRe"
   #RunDEEPZIP  "GaGa"
   #RunDEEPZIP  "HoSa"

   #ds2 and ds3
   #RunDEEPZIP  "HoSaY"
   #RunDEEPZIP  "Mitochondrion"

   #RunDEEPZIP  "VDB"
   #RunDEEPZIP  "Archaea"

   #ds1
   #RunDEEPZIP "GoGoC"
   #RunDEEPZIP "PaTrC"
   #RunDEEPZIP "HoSaC"
   #RunDEEPZIP "PiAbC"
   #
   cd ../../
   echo "Done!";
fi


if [[ "$RUN_ZPAQ" -eq "1" ]]; then
   echo "Running ZPAQ ...";
   mkdir -p res
   cd progs/zpaq

   #ds4
   RunZPAQ  "BuEb"
   RunZPAQ  "AgPh"
   RunZPAQ  "YeMi"
   RunZPAQ  "HePy"
   RunZPAQ  "AeCa"
   RunZPAQ  "HaHi"
   RunZPAQ  "EsCo"
   RunZPAQ  "PlFa"
   RunZPAQ  "ScPo"
   RunZPAQ  "EnIn"
   RunZPAQ  "DrMe"
   RunZPAQ  "OrSa"
   RunZPAQ  "DaRe"
   RunZPAQ  "GaGa"
   RunZPAQ  "HoSa"

   #ds2 and ds3
   RunZPAQ  "HoSaY"
   RunZPAQ  "Mitochondrion"

   RunZPAQ  "VDB"
   RunZPAQ  "Archaea"

   #ds1
   RunZPAQ "GoGoC"
   RunZPAQ "PaTrC"
   RunZPAQ "HoSaC"
   RunZPAQ "PiAbC"
   #
   cd ../../
   echo "Done!";
fi
###################################################################
#### Specific
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

   #ds1 ERROR can't compress
   #RunXM "GoGoC"
   #RunXM "PaTrC"
   #RunXM "HoSaC"
   #RunXM "PiAbC"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_JARVIS" -eq "1" ]]; then
   echo "Running JARVIS ...";
   mkdir -p res
   cd progs/jarvis

   #ds4
   RunJARVIS "-l 1"  "BuEb"
   RunJARVIS "-l 2"  "AgPh"
   RunJARVIS "-l 2"  "YeMi"
   RunJARVIS "-l 3"  "HePy"
   RunJARVIS "-l 3"  "AeCa"
   RunJARVIS "-l 3"  "HaHi"
   RunJARVIS "-l 4"  "EsCo"
   RunJARVIS "-l 4"  "PlFa"
   RunJARVIS "-l 4"  "ScPo"
   RunJARVIS "-l 4"  "EnIn"
   RunJARVIS "-l 5"  "DrMe"
   RunJARVIS "-l 5"  "OrSa"
   RunJARVIS "-l 5"  "DaRe"
   RunJARVIS "-l 6"  "GaGa"
   RunJARVIS "-l 7"  "HoSa"

   #ds2 and ds3
   RunJARVIS "-l 7"  "HoSaY"
   RunJARVIS "-l 7"  "Mitochondrion"

   RunJARVIS "-l 7"  "VDB"
   RunJARVIS "-l 7"  "Archaea"

   #ds1 OOM with 32GB
   #RunJARVIS "-l 7" "GoGoC"
   #RunJARVIS "-l 7" "PaTrC"
   #RunJARVIS "-l 7" "HoSaC"
   #RunJARVIS "-l 7" "PiAbC"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
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

   RunGeCo2 "$PARAM_DS23" "HoSaY"
   RunGeCo2 "$PARAM_DS23" "Mitochondrion"

   RunGeCo2 "$PARAM_DS23" "VDB"
   RunGeCo2 "$PARAM_DS23" "Archaea"

   #ds1
   PARAM_DS1="-tm 3:1:1:1:0.70/0:0:0 -tm 8:1:1:1:0.85/0:0:0 -tm 13:10:0:1:0.85/0:0:0 -tm 19:500:1:40:0.85/5:20:0.85 -ls 0.03 -hs 64"
   RunGeCo2 "$PARAM_DS1" "GoGoC"
   RunGeCo2 "$PARAM_DS1" "PaTrC"
   RunGeCo2 "$PARAM_DS1" "HoSaC"
   RunGeCo2 "$PARAM_DS1" "PiAbC"
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

   RunGeCo3 "$PARAM_DS23" "HoSaY"
   RunGeCo3 "$PARAM_DS23" "Mitochondrion"

   RunGeCo3 "$PARAM_DS23" "VDB"
   RunGeCo3 "$PARAM_DS23" "Archaea"

   #ds1
   PARAM_DS1="-tm 3:1:1:1:0.70/0:0:0 -tm 8:1:1:1:0.85/0:0:0 -tm 13:10:0:1:0.85/0:0:0 -tm 19:500:1:40:0.85/5:20:0.85 -ls 0.03 -hs 64"
   RunGeCo3 "$PARAM_DS1" "GoGoC"
   RunGeCo3 "$PARAM_DS1" "PaTrC"
   RunGeCo3 "$PARAM_DS1" "HoSaC"
   RunGeCo3 "$PARAM_DS1" "PiAbC"
   
   #
   cd ../../
   echo "Done!";
fi
