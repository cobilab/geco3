#!/bin/bash
###############################################################################
INSTALL_GDC2=0;
INSTALL_IDOCOMP=0;
INSTALL_HRCM=0;
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
PARSE=1;
###############################################################################
RUN_GDC2=1;
RUN_IDOCOMP=1;
RUN_HRCM=1;
RUN_GECO2=1;
RUN_GECO3=1;
###############################################################################
#
function Parse {
if test -f "$1.fasta"; then
    cat $1.fasta | grep -v ">" | tr -d -c "ACGT" > XTMP
    cat HEADER XTMP > datasets/$1;
    cp datasets/$1 datasets-hrcm/$1.fasta
    fold -w80 datasets-hrcm/$1.fasta > XTMP
    cp XTMP datasets-hrcm/$1.fasta
    rm XTMP;
fi
}
#
function RunGDC2 {
   # 1 - TARGET
   # 2 - REFERENCE
if test -f "../../datasets/$1" && test -f "../../datasets/$2"; then
   cp ../../datasets/$1 .
   cp ../../datasets/$2 .
   rm -f xxx*
   (/usr/bin/time -v ./GDC2 c xxx $2 $1 ) &> ../../results/C_GDC_$1-$2
   ls -la xxx.gdc2_rc | awk '{ print $5;}' > ../../results/BC_GDC_$1-$2
   rm -f $2 $1 xxx*;
fi
}
#
function RunGeCo2 {
    PARAMR=" -rm 20:500:1:35:0.95/3:100:0.95 -rm 13:200:1:1:0.95/0:0:0 -rm 10:10:0:0:0.95/0:0:0 -lr 0.03 -hs 64 ";
    PARAMH=" -rm 20:500:1:35:0.95/3:100:0.95 -rm 13:200:1:1:0.95/0:0:0 -rm 10:10:0:0:0.95/0:0:0 -tm 4:1:0:1:0.9/0:0:0 -tm 17:100:1:10:0.95/2:20:0.95 -lr 0.03 -hs 64 ";
   # 1 - TARGET
   # 2 - REFERENCE
if test -f "../../datasets/$1" && test -f "../../datasets/$2"; then
   cp ../../datasets/$1 .
   cp ../../datasets/$2 .
   rm -f $1.co
   (/usr/bin/time -v ./GeCo2 $PARAMR -r $2 $1 ) &> ../../results/C_GECO2R_REFO_$1-$2
   ls -la $1.co | awk '{ print $5;}' > ../../results/BC_GECO2R_REFO_$1-$2

   (/usr/bin/time -v ./GeCo2 $PARAMH -r $2 $1 ) &> ../../results/C_GECO2H_REFO_$1-$2
   ls -la $1.co | awk '{ print $5;}' > ../../results/BC_GECO2H_REFO_$1-$2
   rm -f $2 $1 $1.co;
fi
}

function RunGeCo3 {
    PARAMR=" -rm 20:500:1:35:0.95/3:100:0.95 -rm 13:200:1:1:0.95/0:0:0 -rm 10:10:0:0:0.95/0:0:0 -lr 0.03 -hs 64 ";
    PARAMH=" -rm 20:500:1:35:0.95/3:100:0.95 -rm 13:200:1:1:0.95/0:0:0 -rm 10:10:0:0:0.95/0:0:0 -tm 4:1:0:1:0.9/0:0:0 -tm 17:100:1:10:0.95/2:20:0.95 -lr 0.03 -hs 64 ";
   # 1 - TARGET
   # 2 - REFERENCE
if test -f "../../datasets/$1" && test -f "../../datasets/$2"; then
   cp ../../datasets/$1 .
   cp ../../datasets/$2 .
   rm -f $1.co
   (/usr/bin/time -v ./GeCo3 $PARAMR -r $2 $1 ) &> ../../results/C_GECO3R_REFO_$1-$2
   ls -la $1.co | awk '{ print $5;}' > ../../results/BC_GECO3R_REFO_$1-$2

   (/usr/bin/time -v ./GeCo3 $PARAMH -r $2 $1 ) &> ../../results/C_GECO3H_REFO_$1-$2
   ls -la $1.co | awk '{ print $5;}' > ../../results/BC_GECO3H_REFO_$1-$2
   rm -f $2 $1 $1.co;
fi
}

function RunHRCM {
   # 1 - TARGET
   # 2 - REFERENCE
if test -f "../../datasets-hrcm/$1.fasta" && test -f "../../datasets-hrcm/$2.fasta"; then
   cp ../../datasets-hrcm/$1.fasta .
   cp ../../datasets-hrcm/$2.fasta .
   rm -f $1.7z
   (/usr/bin/time -v ./hrcm compress -r $2.fasta -t $1.fasta ) &> ../../results/C_HRCM_REF_$1-$2
   ls -la $1.7z | awk '{ print $5;}' > ../../results/BC_HRCM_REF_$1-$2
   rm -f $2.fasta $1.fasta $1.7z;
fi
}
#
function RunIDoComp {
   # 1 - TARGET
   # 2 - REFERENCE
if test -f "../../datasets/$1" && test -f "../../datasets/$2"; then
   cd sais-lite-2.4.1/
   rm -fr sa ref tar
   mkdir sa ref tar;
   mkdir tmp_oneline_ref;
   cp ../../../datasets/$1 tar/$1.fa
   cp ../../../datasets/$2 ref/$2.fa
   (./generateSA.sh ref sa ) &> TIME_SA
   TIMEOFSA=`cat TIME_SA | grep "..." | awk '{ print $5;}'`
   echo "ref/$2.fa tar/$1.fa sa/$2.sa" > f.txt;
   (/usr/bin/time -v ./iDoComp.run c f.txt OUT ) &> ../../../results/C_IDOCOMP_$1-$2
   cat ../../../results/C_IDOCOMP_$1-$2 | grep "Compressed Size:" | awk '{ print $3; }' > ../../../results/BC_IDOCOMP_$1-$2
   CTIME=`cat ../../../results/C_IDOCOMP_$1-$2 | grep "CPU T" | awk '{print $4;}'`
   echo "$TIMEOFSA+$CTIME" | bc -l > ../../../results/CT_IDOCOMP_$1-$2
   rm -fr sa ref tar tmp_oneline_ref
   cd ..
fi
}
#
###############################################################################
# INSTALL
mkdir -p datasets-hrcm
mkdir -p datasets
mkdir -p progs
cd progs/
###############################################################################
# GET iDoComp
#=================================================================
if [[ "$INSTALL_IDOCOMP" -eq "1" ]]; then
    rm -fr iDoComp/
    git clone https://github.com/mikelhernaez/iDoComp.git
    cd iDoComp/iDoComp_website_v1/sais-lite-2.4.1/source-code/
    gcc -o ../sa.run sa_generator.c sais.c -lm
    cd ../../simulations/source_code/
    gcc -o ../../../sais-lite-2.4.1/iDoComp.run idc_generate_mapping.c main.c stats.c arith.c \
	fasta_decompressor.c idc_load_chr.c os_stream.c fasta_compressor.c \
	sam_stream.c -lm
    cd ../../../../
fi
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
# GET GDC2
#=====================================================================
if [[ "$INSTALL_GDC2" -eq "1" ]]; then
    rm -fr GDC2/
    git clone https://github.com/refresh-bio/GDC2.git
    cd GDC2/gdc_2/Gdc2/
    # LIBRARIES ORDER ACCESS CREATE SOME PROBLES (WE ADD THEM TO THE END)
    printf '\nall: gdc2 \n\nCC      = g++\nCFLAGS  = -Wall -O3 -m64 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11\nCLINK   = -lm -O3 -m64 -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 -lz \n\n.cpp.o: \n\t$(CC) $(CFLAGS) -c $< -o $@ \n\ngdc2: c1stage.o c2stage.o fasta.o hasher.o main.o p1stage.o qsmodel.o queue.o rangecod.o timer.o \n\t$(CC) $(CLINK) -o gdc2 c1stage.o c2stage.o fasta.o hasher.o main.o p1stage.o qsmodel.o queue.o rangecod.o timer.o ../libs/libaelf64.a -lz -lpthread \n\nclean: \n\trm gdc2 \n\trm *.o \n' > Makefile;
    make clean ; make
    cp gdc2 ../../GDC2 # TO NOT OVERLAP FOLDER NAME
    cd ../../../
fi
###############################################################################
# GET HRCM
#====================================================================
if [[ "$INSTALL_HRCM" -eq "1" ]]; then
   rm -fr HRCM/
   git clone https://github.com/haicy/HRCM.git
   cd HRCM/
   chmod +x 7za
   make
   cd ../
fi
##############################################################################
cd ..
###############################################################################
# DOWNLOAD
if [[ "$DOWNLOAD" -eq "1" ]]; then
    echo "Downloading ..."
fi
# PARSE
#=======================================================================
if [[ "$PARSE" -eq "1" ]]; then
   echo "Parsing ..."
   echo ">X" > HEADER
   #
   Parse "HSK1_C1"
   Parse "HSK1_C2"
   Parse "HSK1_C3"
   Parse "HSK1_C4"
   Parse "HSK1_C5"
   Parse "HSK1_C6"
   Parse "HSK1_C7"
   Parse "HSK1_C8"
   Parse "HSK1_C9"
   Parse "HSK1_C10"
   Parse "HSK1_C11"
   Parse "HSK1_C12"
   Parse "HSK1_C13"
   Parse "HSK1_C14"
   Parse "HSK1_C15"
   Parse "HSK1_C16"
   Parse "HSK1_C17"
   Parse "HSK1_C18"
   Parse "HSK1_C19"
   Parse "HSK1_C20"
   Parse "HSK1_C21"
   Parse "HSK1_C22"
   Parse "HSK1_X"
   Parse "HSK1_Y"
   Parse "HSK1_M"
   #
   Parse "HSK2_C1"
   Parse "HSK2_C2"
   Parse "HSK2_C3"
   Parse "HSK2_C4"
   Parse "HSK2_C5"
   Parse "HSK2_C6"
   Parse "HSK2_C7"
   Parse "HSK2_C8"
   Parse "HSK2_C9"
   Parse "HSK2_C10"
   Parse "HSK2_C11"
   Parse "HSK2_C12"
   Parse "HSK2_C13"
   Parse "HSK2_C14"
   Parse "HSK2_C15"
   Parse "HSK2_C16"
   Parse "HSK2_C17"
   Parse "HSK2_C18"
   Parse "HSK2_C19"
   Parse "HSK2_C20"
   Parse "HSK2_C21"
   Parse "HSK2_C22"
   Parse "HSK2_X"
   Parse "HSK2_Y"
   Parse "HSK2_M"
   #
   echo "done!";
   fi
#
# RUN
#=========================================================================
if [[ "$RUN_HRCM" -eq "1" ]]; then
   echo "Running HRCM ...";
   mkdir -p results
   cd progs/HRCM
   # target $1, reference $2:
   RunHRCM "HSK1_C1" "HSK2_C1"
   #RunHRCM "HSK1_C2" "HSK2_C2"
   #RunHRCM "HSK1_C3" "HSK2_C3"
   #RunHRCM "HSK1_C4" "HSK2_C4"
   #RunHRCM "HSK1_C5" "HSK2_C5"
   #RunHRCM "HSK1_C6" "HSK2_C6"
   #RunHRCM "HSK1_C7" "HSK2_C7"
   #RunHRCM "HSK1_C8" "HSK2_C8"
   #RunHRCM "HSK1_C9" "HSK2_C9"
   #RunHRCM "HSK1_C10" "HSK2_C10"
   #RunHRCM "HSK1_C11" "HSK2_C11"
   #RunHRCM "HSK1_C12" "HSK2_C12"
   #RunHRCM "HSK1_C13" "HSK2_C13"
   #RunHRCM "HSK1_C14" "HSK2_C14"
   #RunHRCM "HSK1_C15" "HSK2_C15"
   #RunHRCM "HSK1_C16" "HSK2_C16"
   #RunHRCM "HSK1_C17" "HSK2_C17"
   #RunHRCM "HSK1_C18" "HSK2_C18"
   #RunHRCM "HSK1_C19" "HSK2_C19"
   #RunHRCM "HSK1_C20" "HSK2_C20"
   #RunHRCM "HSK1_C21" "HSK2_C21"
   #RunHRCM "HSK1_C22" "HSK2_C22"
   #RunHRCM "HSK1_X" "HSK2_X"
   #RunHRCM "HSK1_Y" "HSK2_Y"
   #RunHRCM "HSK1_M" "HSK2_M"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GDC2" -eq "1" ]]; then
   echo "Running GDC2 ...";
   mkdir -p results
   cd progs/GDC2
   # target $1, reference $2:
   RunGDC2 "HSK1_C1" "HSK2_C1"
   #RunGDC2 "HSK1_C2" "HSK2_C2"
   #RunGDC2 "HSK1_C3" "HSK2_C3"
   #RunGDC2 "HSK1_C4" "HSK2_C4"
   #RunGDC2 "HSK1_C5" "HSK2_C5"
   #RunGDC2 "HSK1_C6" "HSK2_C6"
   #RunGDC2 "HSK1_C7" "HSK2_C7"
   #RunGDC2 "HSK1_C8" "HSK2_C8"
   #RunGDC2 "HSK1_C9" "HSK2_C9"
   #RunGDC2 "HSK1_C10" "HSK2_C10"
   #RunGDC2 "HSK1_C11" "HSK2_C11"
   #RunGDC2 "HSK1_C12" "HSK2_C12"
   #RunGDC2 "HSK1_C13" "HSK2_C13"
   #RunGDC2 "HSK1_C14" "HSK2_C14"
   #RunGDC2 "HSK1_C15" "HSK2_C15"
   #RunGDC2 "HSK1_C16" "HSK2_C16"
   #RunGDC2 "HSK1_C17" "HSK2_C17"
   #RunGDC2 "HSK1_C18" "HSK2_C18"
   #RunGDC2 "HSK1_C19" "HSK2_C19"
   #RunGDC2 "HSK1_C20" "HSK2_C20"
   #RunGDC2 "HSK1_C21" "HSK2_C21"
   #RunGDC2 "HSK1_C22" "HSK2_C22"
   #RunGDC2 "HSK1_X" "HSK2_X"
   #RunGDC2 "HSK1_Y" "HSK2_Y"
   #RunGDC2 "HSK1_M" "HSK2_M"
   #
   cd ../../
   echo "Done!";
fi
#==============================================================================
if [[ "$RUN_IDOCOMP" -eq "1" ]]; then
   echo "Running iDoComp ...";
   mkdir -p results
   cd progs/iDoComp
   # target $1, reference $2:
   RunIDoComp "HSK1_C1" "HSK2_C1"
   #RunIDoComp "HSK1_C2" "HSK2_C2"
   #RunIDoComp "HSK1_C3" "HSK2_C3"
   #RunIDoComp "HSK1_C4" "HSK2_C4"
   #RunIDoComp "HSK1_C5" "HSK2_C5"
   #RunIDoComp "HSK1_C6" "HSK2_C6"
   #RunIDoComp "HSK1_C7" "HSK2_C7"
   #RunIDoComp "HSK1_C8" "HSK2_C8"
   #RunIDoComp "HSK1_C9" "HSK2_C9"
   #RunIDoComp "HSK1_C10" "HSK2_C10"
   #RunIDoComp "HSK1_C11" "HSK2_C11"
   #RunIDoComp "HSK1_C12" "HSK2_C12"
   #RunIDoComp "HSK1_C13" "HSK2_C13"
   #RunIDoComp "HSK1_C14" "HSK2_C14"
   #RunIDoComp "HSK1_C15" "HSK2_C15"
   #RunIDoComp "HSK1_C16" "HSK2_C16"
   #RunIDoComp "HSK1_C17" "HSK2_C17"
   #RunIDoComp "HSK1_C18" "HSK2_C18"
   #RunIDoComp "HSK1_C19" "HSK2_C19"
   #RunIDoComp "HSK1_C20" "HSK2_C20"
   #RunIDoComp "HSK1_C21" "HSK2_C21"
   #RunIDoComp "HSK1_C22" "HSK2_C22"
   #RunIDoComp "HSK1_X" "HSK2_X"
   #RunIDoComp "HSK1_Y" "HSK2_Y"
   #RunIDoComp "HSK1_M" "HSK2_M"
   #
   cd ../../
   echo "Done!";
fi
#==============================================================================
if [[ "$RUN_GECO2" -eq "1" ]]; then
   echo "Running GeCo2 ...";
   mkdir -p results
   cd progs/geco2
   # target $1, reference $2:
   RunGeCo2 "HSK1_C1" "HSK2_C1"
   RunGeCo2 "HSK1_C2" "HSK2_C2"
   RunGeCo2 "HSK1_C3" "HSK2_C3"
   RunGeCo2 "HSK1_C4" "HSK2_C4"
   RunGeCo2 "HSK1_C5" "HSK2_C5"
   RunGeCo2 "HSK1_C6" "HSK2_C6"
   RunGeCo2 "HSK1_C7" "HSK2_C7"
   RunGeCo2 "HSK1_C8" "HSK2_C8"
   RunGeCo2 "HSK1_C9" "HSK2_C9"
   RunGeCo2 "HSK1_C10" "HSK2_C10"
   RunGeCo2 "HSK1_C11" "HSK2_C11"
   RunGeCo2 "HSK1_C12" "HSK2_C12"
   RunGeCo2 "HSK1_C13" "HSK2_C13"
   RunGeCo2 "HSK1_C14" "HSK2_C14"
   RunGeCo2 "HSK1_C15" "HSK2_C15"
   RunGeCo2 "HSK1_C16" "HSK2_C16"
   RunGeCo2 "HSK1_C17" "HSK2_C17"
   RunGeCo2 "HSK1_C18" "HSK2_C18"
   RunGeCo2 "HSK1_C19" "HSK2_C19"
   RunGeCo2 "HSK1_C20" "HSK2_C20"
   RunGeCo2 "HSK1_C21" "HSK2_C21"
   RunGeCo2 "HSK1_C22" "HSK2_C22"
   RunGeCo2 "HSK1_X" "HSK2_X"
   RunGeCo2 "HSK1_Y" "HSK2_Y"
   RunGeCo2 "HSK1_M" "HSK2_M"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p results
   cd progs/geco3
   # target $1, reference $2:
   RunGeCo3 "HSK1_C1" "HSK2_C1"
   #RunGeCo3 "HSK1_C2" "HSK2_C2"
   #RunGeCo3 "HSK1_C3" "HSK2_C3"
   #RunGeCo3 "HSK1_C4" "HSK2_C4"
   #RunGeCo3 "HSK1_C5" "HSK2_C5"
   #RunGeCo3 "HSK1_C6" "HSK2_C6"
   #RunGeCo3 "HSK1_C7" "HSK2_C7"
   #RunGeCo3 "HSK1_C8" "HSK2_C8"
   #RunGeCo3 "HSK1_C9" "HSK2_C9"
   #RunGeCo3 "HSK1_C10" "HSK2_C10"
   #RunGeCo3 "HSK1_C11" "HSK2_C11"
   #RunGeCo3 "HSK1_C12" "HSK2_C12"
   #RunGeCo3 "HSK1_C13" "HSK2_C13"
   #RunGeCo3 "HSK1_C14" "HSK2_C14"
   #RunGeCo3 "HSK1_C15" "HSK2_C15"
   #RunGeCo3 "HSK1_C16" "HSK2_C16"
   #RunGeCo3 "HSK1_C17" "HSK2_C17"
   #RunGeCo3 "HSK1_C18" "HSK2_C18"
   #RunGeCo3 "HSK1_C19" "HSK2_C19"
   #RunGeCo3 "HSK1_C20" "HSK2_C20"
   #RunGeCo3 "HSK1_C21" "HSK2_C21"
   #RunGeCo3 "HSK1_C22" "HSK2_C22"
   #RunGeCo3 "HSK1_X" "HSK2_X"
   #RunGeCo3 "HSK1_Y" "HSK2_Y"
   #RunGeCo3 "HSK1_M" "HSK2_M"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
