#!/bin/bash
###############################################################################
INSTALL_GDC2=0;
INSTALL_IDOCOMP=0;
INSTALL_HRCM=0;
INSTALL_GECO2=0;
INSTALL_GECO3=0;
###############################################################################
DOWNLOAD=0;
PARSE=0;
###############################################################################
RUN_GDC2=0;
RUN_IDOCOMP=0;
RUN_HRCM=0;
RUN_GECO2=1;
RUN_GECO3=0;
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
   Parse "HS_C1"
   Parse "HS_C2"
   Parse "HS_C3"
   Parse "HS_C4"
   Parse "HS_C5"
   Parse "HS_C6"
   Parse "HS_C7"
   Parse "HS_C8"
   Parse "HS_C9"
   Parse "HS_C10"
   Parse "HS_C11"
   Parse "HS_C12"
   Parse "HS_C13"
   Parse "HS_C14"
   Parse "HS_C15"
   Parse "HS_C16"
   Parse "HS_C17"
   Parse "HS_C18"
   Parse "HS_C19"
   Parse "HS_C20"
   Parse "HS_C21"
   Parse "HS_C22"
   Parse "HS_X"
   Parse "HS_Y"
   Parse "HS_M"
   #
   Parse "PT_C1"
   Parse "PT_C2"
   Parse "PT_C3"
   Parse "PT_C4"
   Parse "PT_C5"
   Parse "PT_C6"
   Parse "PT_C7"
   Parse "PT_C8"
   Parse "PT_C9"
   Parse "PT_C10"
   Parse "PT_C11"
   Parse "PT_C12"
   Parse "PT_C13"
   Parse "PT_C14"
   Parse "PT_C15"
   Parse "PT_C16"
   Parse "PT_C17"
   Parse "PT_C18"
   Parse "PT_C19"
   Parse "PT_C20"
   Parse "PT_C21"
   Parse "PT_C22"
   Parse "PT_X"
   Parse "PT_Y"
   Parse "PT_M"
   #

   Parse "GG_C1"
   Parse "GG_C2"
   Parse "GG_C3"
   Parse "GG_C4"
   Parse "GG_C5"
   Parse "GG_C6"
   Parse "GG_C7"
   Parse "GG_C8"
   Parse "GG_C9"
   Parse "GG_C10"
   Parse "GG_C11"
   Parse "GG_C12"
   Parse "GG_C13"
   Parse "GG_C14"
   Parse "GG_C15"
   Parse "GG_C16"
   Parse "GG_C17"
   Parse "GG_C18"
   Parse "GG_C19"
   Parse "GG_C20"
   Parse "GG_C21"
   Parse "GG_C22"
   Parse "GG_X"
   Parse "GG_M"
   #
   Parse "PA_C1"
   Parse "PA_C2"
   Parse "PA_C3"
   Parse "PA_C4"
   Parse "PA_C5"
   Parse "PA_C6"
   Parse "PA_C7"
   Parse "PA_C8"
   Parse "PA_C9"
   Parse "PA_C10"
   Parse "PA_C11"
   Parse "PA_C12"
   Parse "PA_C13"
   Parse "PA_C14"
   Parse "PA_C15"
   Parse "PA_C16"
   Parse "PA_C17"
   Parse "PA_C18"
   Parse "PA_C19"
   Parse "PA_C20"
   Parse "PA_C21"
   Parse "PA_C22"
   Parse "PA_X"
   Parse "PA_M"
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
   RunHRCM "HS_C1" "GG_C1"
   RunHRCM "HS_C2" "GG_C2"
   RunHRCM "HS_C3" "GG_C3"
   RunHRCM "HS_C4" "GG_C4"
   RunHRCM "HS_C5" "GG_C5"
   RunHRCM "HS_C6" "GG_C6"
   RunHRCM "HS_C7" "GG_C7"
   RunHRCM "HS_C8" "GG_C8"
   RunHRCM "HS_C9" "GG_C9"
   RunHRCM "HS_C10" "GG_C10"
   RunHRCM "HS_C11" "GG_C11"
   RunHRCM "HS_C12" "GG_C12"
   RunHRCM "HS_C13" "GG_C13"
   RunHRCM "HS_C14" "GG_C14"
   RunHRCM "HS_C15" "GG_C15"
   RunHRCM "HS_C16" "GG_C16"
   RunHRCM "HS_C17" "GG_C17"
   RunHRCM "HS_C18" "GG_C18"
   RunHRCM "HS_C19" "GG_C19"
   RunHRCM "HS_C20" "GG_C20"
   RunHRCM "HS_C21" "GG_C21"
   RunHRCM "HS_C22" "GG_C22"
   RunHRCM "HS_X" "GG_X"
   RunHRCM "HS_M" "GG_M"
   #
   RunHRCM "PT_C1" "HS_C1"
   RunHRCM "PT_C2" "HS_C2"
   RunHRCM "PT_C3" "HS_C3"
   RunHRCM "PT_C4" "HS_C4"
   RunHRCM "PT_C5" "HS_C5"
   RunHRCM "PT_C6" "HS_C6"
   RunHRCM "PT_C7" "HS_C7"
   RunHRCM "PT_C8" "HS_C8"
   RunHRCM "PT_C9" "HS_C9"
   RunHRCM "PT_C10" "HS_C10"
   RunHRCM "PT_C11" "HS_C11"
   RunHRCM "PT_C12" "HS_C12"
   RunHRCM "PT_C13" "HS_C13"
   RunHRCM "PT_C14" "HS_C14"
   RunHRCM "PT_C15" "HS_C15"
   RunHRCM "PT_C16" "HS_C16"
   RunHRCM "PT_C17" "HS_C17"
   RunHRCM "PT_C18" "HS_C18"
   RunHRCM "PT_C19" "HS_C19"
   RunHRCM "PT_C20" "HS_C20"
   RunHRCM "PT_C21" "HS_C21"
   RunHRCM "PT_C22" "HS_C22"
   RunHRCM "PT_X" "HS_X"
   RunHRCM "PT_Y" "HS_Y"
   RunHRCM "PT_M" "HS_M"
   #
   RunHRCM "PA_C1" "HS_C1"
   RunHRCM "PA_C2" "HS_C2"
   RunHRCM "PA_C3" "HS_C3"
   RunHRCM "PA_C4" "HS_C4"
   RunHRCM "PA_C5" "HS_C5"
   RunHRCM "PA_C6" "HS_C6"
   RunHRCM "PA_C7" "HS_C7"
   RunHRCM "PA_C8" "HS_C8"
   RunHRCM "PA_C9" "HS_C9"
   RunHRCM "PA_C10" "HS_C10"
   RunHRCM "PA_C11" "HS_C11"
   RunHRCM "PA_C12" "HS_C12"
   RunHRCM "PA_C13" "HS_C13"
   RunHRCM "PA_C14" "HS_C14"
   RunHRCM "PA_C15" "HS_C15"
   RunHRCM "PA_C16" "HS_C16"
   RunHRCM "PA_C17" "HS_C17"
   RunHRCM "PA_C18" "HS_C18"
   RunHRCM "PA_C19" "HS_C19"
   RunHRCM "PA_C20" "HS_C20"
   RunHRCM "PA_C21" "HS_C21"
   RunHRCM "PA_C22" "HS_C22"
   RunHRCM "PA_X" "HS_X"
   RunHRCM "PA_M" "HS_M"
   #
   RunHRCM "GG_C1" "HS_C1"
   RunHRCM "GG_C2" "HS_C2"
   RunHRCM "GG_C3" "HS_C3"
   RunHRCM "GG_C4" "HS_C4"
   RunHRCM "GG_C5" "HS_C5"
   RunHRCM "GG_C6" "HS_C6"
   RunHRCM "GG_C7" "HS_C7"
   RunHRCM "GG_C8" "HS_C8"
   RunHRCM "GG_C9" "HS_C9"
   RunHRCM "GG_C10" "HS_C10"
   RunHRCM "GG_C11" "HS_C11"
   RunHRCM "GG_C12" "HS_C12"
   RunHRCM "GG_C13" "HS_C13"
   RunHRCM "GG_C14" "HS_C14"
   RunHRCM "GG_C15" "HS_C15"
   RunHRCM "GG_C16" "HS_C16"
   RunHRCM "GG_C17" "HS_C17"
   RunHRCM "GG_C18" "HS_C18"
   RunHRCM "GG_C19" "HS_C19"
   RunHRCM "GG_C20" "HS_C20"
   RunHRCM "GG_C21" "HS_C21"
   RunHRCM "GG_C22" "HS_C22"
   RunHRCM "GG_X" "HS_X"
   RunHRCM "GG_M" "HS_M"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GDC2" -eq "1" ]]; then
   echo "Running GDC2 ...";
   mkdir -p results
   cd progs/GDC2
   # target $1, reference $2:
   RunGDC2 "HS_C1" "GG_C1"
   RunGDC2 "HS_C2" "GG_C2"
   RunGDC2 "HS_C3" "GG_C3"
   RunGDC2 "HS_C4" "GG_C4"
   RunGDC2 "HS_C5" "GG_C5"
   RunGDC2 "HS_C6" "GG_C6"
   RunGDC2 "HS_C7" "GG_C7"
   RunGDC2 "HS_C8" "GG_C8"
   RunGDC2 "HS_C9" "GG_C9"
   RunGDC2 "HS_C10" "GG_C10"
   RunGDC2 "HS_C11" "GG_C11"
   RunGDC2 "HS_C12" "GG_C12"
   RunGDC2 "HS_C13" "GG_C13"
   RunGDC2 "HS_C14" "GG_C14"
   RunGDC2 "HS_C15" "GG_C15"
   RunGDC2 "HS_C16" "GG_C16"
   RunGDC2 "HS_C17" "GG_C17"
   RunGDC2 "HS_C18" "GG_C18"
   RunGDC2 "HS_C19" "GG_C19"
   RunGDC2 "HS_C20" "GG_C20"
   RunGDC2 "HS_C21" "GG_C21"
   RunGDC2 "HS_C22" "GG_C22"
   RunGDC2 "HS_X" "GG_X"
   RunGDC2 "HS_M" "GG_M"
   #
   RunGDC2 "PT_C1" "HS_C1"
   RunGDC2 "PT_C2" "HS_C2"
   RunGDC2 "PT_C3" "HS_C3"
   RunGDC2 "PT_C4" "HS_C4"
   RunGDC2 "PT_C5" "HS_C5"
   RunGDC2 "PT_C6" "HS_C6"
   RunGDC2 "PT_C7" "HS_C7"
   RunGDC2 "PT_C8" "HS_C8"
   RunGDC2 "PT_C9" "HS_C9"
   RunGDC2 "PT_C10" "HS_C10"
   RunGDC2 "PT_C11" "HS_C11"
   RunGDC2 "PT_C12" "HS_C12"
   RunGDC2 "PT_C13" "HS_C13"
   RunGDC2 "PT_C14" "HS_C14"
   RunGDC2 "PT_C15" "HS_C15"
   RunGDC2 "PT_C16" "HS_C16"
   RunGDC2 "PT_C17" "HS_C17"
   RunGDC2 "PT_C18" "HS_C18"
   RunGDC2 "PT_C19" "HS_C19"
   RunGDC2 "PT_C20" "HS_C20"
   RunGDC2 "PT_C21" "HS_C21"
   RunGDC2 "PT_C22" "HS_C22"
   RunGDC2 "PT_X" "HS_X"
   RunGDC2 "PT_Y" "HS_Y"
   RunGDC2 "PT_M" "HS_M"
   #
   RunGDC2 "PA_C1" "HS_C1"
   RunGDC2 "PA_C2" "HS_C2"
   RunGDC2 "PA_C3" "HS_C3"
   RunGDC2 "PA_C4" "HS_C4"
   RunGDC2 "PA_C5" "HS_C5"
   RunGDC2 "PA_C6" "HS_C6"
   RunGDC2 "PA_C7" "HS_C7"
   RunGDC2 "PA_C8" "HS_C8"
   RunGDC2 "PA_C9" "HS_C9"
   RunGDC2 "PA_C10" "HS_C10"
   RunGDC2 "PA_C11" "HS_C11"
   RunGDC2 "PA_C12" "HS_C12"
   RunGDC2 "PA_C13" "HS_C13"
   RunGDC2 "PA_C14" "HS_C14"
   RunGDC2 "PA_C15" "HS_C15"
   RunGDC2 "PA_C16" "HS_C16"
   RunGDC2 "PA_C17" "HS_C17"
   RunGDC2 "PA_C18" "HS_C18"
   RunGDC2 "PA_C19" "HS_C19"
   RunGDC2 "PA_C20" "HS_C20"
   RunGDC2 "PA_C21" "HS_C21"
   RunGDC2 "PA_C22" "HS_C22"
   RunGDC2 "PA_X" "HS_X"
   RunGDC2 "PA_M" "HS_M"
   #
   RunGDC2 "GG_C1" "HS_C1"
   RunGDC2 "GG_C2" "HS_C2"
   RunGDC2 "GG_C3" "HS_C3"
   RunGDC2 "GG_C4" "HS_C4"
   RunGDC2 "GG_C5" "HS_C5"
   RunGDC2 "GG_C6" "HS_C6"
   RunGDC2 "GG_C7" "HS_C7"
   RunGDC2 "GG_C8" "HS_C8"
   RunGDC2 "GG_C9" "HS_C9"
   RunGDC2 "GG_C10" "HS_C10"
   RunGDC2 "GG_C11" "HS_C11"
   RunGDC2 "GG_C12" "HS_C12"
   RunGDC2 "GG_C13" "HS_C13"
   RunGDC2 "GG_C14" "HS_C14"
   RunGDC2 "GG_C15" "HS_C15"
   RunGDC2 "GG_C16" "HS_C16"
   RunGDC2 "GG_C17" "HS_C17"
   RunGDC2 "GG_C18" "HS_C18"
   RunGDC2 "GG_C19" "HS_C19"
   RunGDC2 "GG_C20" "HS_C20"
   RunGDC2 "GG_C21" "HS_C21"
   RunGDC2 "GG_C22" "HS_C22"
   RunGDC2 "GG_X" "HS_X"
   RunGDC2 "GG_M" "HS_M"
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
   RunIDoComp "HS_C1" "GG_C1"
   RunIDoComp "HS_C2" "GG_C2"
   RunIDoComp "HS_C3" "GG_C3"
   RunIDoComp "HS_C4" "GG_C4"
   RunIDoComp "HS_C5" "GG_C5"
   RunIDoComp "HS_C6" "GG_C6"
   RunIDoComp "HS_C7" "GG_C7"
   RunIDoComp "HS_C8" "GG_C8"
   RunIDoComp "HS_C9" "GG_C9"
   RunIDoComp "HS_C10" "GG_C10"
   RunIDoComp "HS_C11" "GG_C11"
   RunIDoComp "HS_C12" "GG_C12"
   RunIDoComp "HS_C13" "GG_C13"
   RunIDoComp "HS_C14" "GG_C14"
   RunIDoComp "HS_C15" "GG_C15"
   RunIDoComp "HS_C16" "GG_C16"
   RunIDoComp "HS_C17" "GG_C17"
   RunIDoComp "HS_C18" "GG_C18"
   RunIDoComp "HS_C19" "GG_C19"
   RunIDoComp "HS_C20" "GG_C20"
   RunIDoComp "HS_C21" "GG_C21"
   RunIDoComp "HS_C22" "GG_C22"
   RunIDoComp "HS_X" "GG_X"
   RunIDoComp "HS_M" "GG_M"
   #
   RunIDoComp "PT_C1" "HS_C1"
   RunIDoComp "PT_C2" "HS_C2"
   RunIDoComp "PT_C3" "HS_C3"
   RunIDoComp "PT_C4" "HS_C4"
   RunIDoComp "PT_C5" "HS_C5"
   RunIDoComp "PT_C6" "HS_C6"
   RunIDoComp "PT_C7" "HS_C7"
   RunIDoComp "PT_C8" "HS_C8"
   RunIDoComp "PT_C9" "HS_C9"
   RunIDoComp "PT_C10" "HS_C10"
   RunIDoComp "PT_C11" "HS_C11"
   RunIDoComp "PT_C12" "HS_C12"
   RunIDoComp "PT_C13" "HS_C13"
   RunIDoComp "PT_C14" "HS_C14"
   RunIDoComp "PT_C15" "HS_C15"
   RunIDoComp "PT_C16" "HS_C16"
   RunIDoComp "PT_C17" "HS_C17"
   RunIDoComp "PT_C18" "HS_C18"
   RunIDoComp "PT_C19" "HS_C19"
   RunIDoComp "PT_C20" "HS_C20"
   RunIDoComp "PT_C21" "HS_C21"
   RunIDoComp "PT_C22" "HS_C22"
   RunIDoComp "PT_X" "HS_X"
   RunIDoComp "PT_Y" "HS_Y"
   RunIDoComp "PT_M" "HS_M"
   #
   RunIDoComp "PA_C1" "HS_C1"
   RunIDoComp "PA_C2" "HS_C2"
   RunIDoComp "PA_C3" "HS_C3"
   RunIDoComp "PA_C4" "HS_C4"
   RunIDoComp "PA_C5" "HS_C5"
   RunIDoComp "PA_C6" "HS_C6"
   RunIDoComp "PA_C7" "HS_C7"
   RunIDoComp "PA_C8" "HS_C8"
   RunIDoComp "PA_C9" "HS_C9"
   RunIDoComp "PA_C10" "HS_C10"
   RunIDoComp "PA_C11" "HS_C11"
   RunIDoComp "PA_C12" "HS_C12"
   RunIDoComp "PA_C13" "HS_C13"
   RunIDoComp "PA_C14" "HS_C14"
   RunIDoComp "PA_C15" "HS_C15"
   RunIDoComp "PA_C16" "HS_C16"
   RunIDoComp "PA_C17" "HS_C17"
   RunIDoComp "PA_C18" "HS_C18"
   RunIDoComp "PA_C19" "HS_C19"
   RunIDoComp "PA_C20" "HS_C20"
   RunIDoComp "PA_C21" "HS_C21"
   RunIDoComp "PA_C22" "HS_C22"
   RunIDoComp "PA_X" "HS_X"
   RunIDoComp "PA_M" "HS_M"
   #
   RunIDoComp "GG_C1" "HS_C1"
   RunIDoComp "GG_C2" "HS_C2"
   RunIDoComp "GG_C3" "HS_C3"
   RunIDoComp "GG_C4" "HS_C4"
   RunIDoComp "GG_C5" "HS_C5"
   RunIDoComp "GG_C6" "HS_C6"
   RunIDoComp "GG_C7" "HS_C7"
   RunIDoComp "GG_C8" "HS_C8"
   RunIDoComp "GG_C9" "HS_C9"
   RunIDoComp "GG_C10" "HS_C10"
   RunIDoComp "GG_C11" "HS_C11"
   RunIDoComp "GG_C12" "HS_C12"
   RunIDoComp "GG_C13" "HS_C13"
   RunIDoComp "GG_C14" "HS_C14"
   RunIDoComp "GG_C15" "HS_C15"
   RunIDoComp "GG_C16" "HS_C16"
   RunIDoComp "GG_C17" "HS_C17"
   RunIDoComp "GG_C18" "HS_C18"
   RunIDoComp "GG_C19" "HS_C19"
   RunIDoComp "GG_C20" "HS_C20"
   RunIDoComp "GG_C21" "HS_C21"
   RunIDoComp "GG_C22" "HS_C22"
   RunIDoComp "GG_X" "HS_X"
   RunIDoComp "GG_M" "HS_M"
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
   RunGeCo2 "HS_C1" "GG_C1"
   RunGeCo2 "HS_C2" "GG_C2"
   RunGeCo2 "HS_C3" "GG_C3"
   RunGeCo2 "HS_C4" "GG_C4"
   RunGeCo2 "HS_C5" "GG_C5"
   RunGeCo2 "HS_C6" "GG_C6"
   RunGeCo2 "HS_C7" "GG_C7"
   RunGeCo2 "HS_C8" "GG_C8"
   RunGeCo2 "HS_C9" "GG_C9"
   RunGeCo2 "HS_C10" "GG_C10"
   RunGeCo2 "HS_C11" "GG_C11"
   RunGeCo2 "HS_C12" "GG_C12"
   RunGeCo2 "HS_C13" "GG_C13"
   RunGeCo2 "HS_C14" "GG_C14"
   RunGeCo2 "HS_C15" "GG_C15"
   RunGeCo2 "HS_C16" "GG_C16"
   RunGeCo2 "HS_C17" "GG_C17"
   RunGeCo2 "HS_C18" "GG_C18"
   RunGeCo2 "HS_C19" "GG_C19"
   RunGeCo2 "HS_C20" "GG_C20"
   RunGeCo2 "HS_C21" "GG_C21"
   RunGeCo2 "HS_C22" "GG_C22"
   RunGeCo2 "HS_X" "GG_X"
   RunGeCo2 "HS_M" "GG_M"
   #
   RunGeCo2 "PT_C1" "HS_C1"
   RunGeCo2 "PT_C2" "HS_C2"
   RunGeCo2 "PT_C3" "HS_C3"
   RunGeCo2 "PT_C4" "HS_C4"
   RunGeCo2 "PT_C5" "HS_C5"
   RunGeCo2 "PT_C6" "HS_C6"
   RunGeCo2 "PT_C7" "HS_C7"
   RunGeCo2 "PT_C8" "HS_C8"
   RunGeCo2 "PT_C9" "HS_C9"
   RunGeCo2 "PT_C10" "HS_C10"
   RunGeCo2 "PT_C11" "HS_C11"
   RunGeCo2 "PT_C12" "HS_C12"
   RunGeCo2 "PT_C13" "HS_C13"
   RunGeCo2 "PT_C14" "HS_C14"
   RunGeCo2 "PT_C15" "HS_C15"
   RunGeCo2 "PT_C16" "HS_C16"
   RunGeCo2 "PT_C17" "HS_C17"
   RunGeCo2 "PT_C18" "HS_C18"
   RunGeCo2 "PT_C19" "HS_C19"
   RunGeCo2 "PT_C20" "HS_C20"
   RunGeCo2 "PT_C21" "HS_C21"
   RunGeCo2 "PT_C22" "HS_C22"
   RunGeCo2 "PT_X" "HS_X"
   RunGeCo2 "PT_Y" "HS_Y"
   RunGeCo2 "PT_M" "HS_M"
   #
   RunGeCo2 "PA_C1" "HS_C1"
   RunGeCo2 "PA_C2" "HS_C2"
   RunGeCo2 "PA_C3" "HS_C3"
   RunGeCo2 "PA_C4" "HS_C4"
   RunGeCo2 "PA_C5" "HS_C5"
   RunGeCo2 "PA_C6" "HS_C6"
   RunGeCo2 "PA_C7" "HS_C7"
   RunGeCo2 "PA_C8" "HS_C8"
   RunGeCo2 "PA_C9" "HS_C9"
   RunGeCo2 "PA_C10" "HS_C10"
   RunGeCo2 "PA_C11" "HS_C11"
   RunGeCo2 "PA_C12" "HS_C12"
   RunGeCo2 "PA_C13" "HS_C13"
   RunGeCo2 "PA_C14" "HS_C14"
   RunGeCo2 "PA_C15" "HS_C15"
   RunGeCo2 "PA_C16" "HS_C16"
   RunGeCo2 "PA_C17" "HS_C17"
   RunGeCo2 "PA_C18" "HS_C18"
   RunGeCo2 "PA_C19" "HS_C19"
   RunGeCo2 "PA_C20" "HS_C20"
   RunGeCo2 "PA_C21" "HS_C21"
   RunGeCo2 "PA_C22" "HS_C22"
   RunGeCo2 "PA_X" "HS_X"
   RunGeCo2 "PA_M" "HS_M"
   #
   RunGeCo2 "GG_C1" "HS_C1"
   RunGeCo2 "GG_C2" "HS_C2"
   RunGeCo2 "GG_C3" "HS_C3"
   RunGeCo2 "GG_C4" "HS_C4"
   RunGeCo2 "GG_C5" "HS_C5"
   RunGeCo2 "GG_C6" "HS_C6"
   RunGeCo2 "GG_C7" "HS_C7"
   RunGeCo2 "GG_C8" "HS_C8"
   RunGeCo2 "GG_C9" "HS_C9"
   RunGeCo2 "GG_C10" "HS_C10"
   RunGeCo2 "GG_C11" "HS_C11"
   RunGeCo2 "GG_C12" "HS_C12"
   RunGeCo2 "GG_C13" "HS_C13"
   RunGeCo2 "GG_C14" "HS_C14"
   RunGeCo2 "GG_C15" "HS_C15"
   RunGeCo2 "GG_C16" "HS_C16"
   RunGeCo2 "GG_C17" "HS_C17"
   RunGeCo2 "GG_C18" "HS_C18"
   RunGeCo2 "GG_C19" "HS_C19"
   RunGeCo2 "GG_C20" "HS_C20"
   RunGeCo2 "GG_C21" "HS_C21"
   RunGeCo2 "GG_C22" "HS_C22"
   RunGeCo2 "GG_X" "HS_X"
   RunGeCo2 "GG_M" "HS_M"
   #
   cd ../../
   echo "Done!";
fi

if [[ "$RUN_GECO3" -eq "1" ]]; then
   echo "Running GeCo3 ...";
   mkdir -p results
   cd progs/geco3
   # target $1, reference $2:
   RunGeCo3 "HS_C1" "GG_C1"
   RunGeCo3 "HS_C2" "GG_C2"
   RunGeCo3 "HS_C3" "GG_C3"
   RunGeCo3 "HS_C4" "GG_C4"
   RunGeCo3 "HS_C5" "GG_C5"
   RunGeCo3 "HS_C6" "GG_C6"
   RunGeCo3 "HS_C7" "GG_C7"
   RunGeCo3 "HS_C8" "GG_C8"
   RunGeCo3 "HS_C9" "GG_C9"
   RunGeCo3 "HS_C10" "GG_C10"
   RunGeCo3 "HS_C11" "GG_C11"
   RunGeCo3 "HS_C12" "GG_C12"
   RunGeCo3 "HS_C13" "GG_C13"
   RunGeCo3 "HS_C14" "GG_C14"
   RunGeCo3 "HS_C15" "GG_C15"
   RunGeCo3 "HS_C16" "GG_C16"
   RunGeCo3 "HS_C17" "GG_C17"
   RunGeCo3 "HS_C18" "GG_C18"
   RunGeCo3 "HS_C19" "GG_C19"
   RunGeCo3 "HS_C20" "GG_C20"
   RunGeCo3 "HS_C21" "GG_C21"
   RunGeCo3 "HS_C22" "GG_C22"
   RunGeCo3 "HS_X" "GG_X"
   RunGeCo3 "HS_M" "GG_M"
   #
   RunGeCo3 "PT_C1" "HS_C1"
   RunGeCo3 "PT_C2" "HS_C2"
   RunGeCo3 "PT_C3" "HS_C3"
   RunGeCo3 "PT_C4" "HS_C4"
   RunGeCo3 "PT_C5" "HS_C5"
   RunGeCo3 "PT_C6" "HS_C6"
   RunGeCo3 "PT_C7" "HS_C7"
   RunGeCo3 "PT_C8" "HS_C8"
   RunGeCo3 "PT_C9" "HS_C9"
   RunGeCo3 "PT_C10" "HS_C10"
   RunGeCo3 "PT_C11" "HS_C11"
   RunGeCo3 "PT_C12" "HS_C12"
   RunGeCo3 "PT_C13" "HS_C13"
   RunGeCo3 "PT_C14" "HS_C14"
   RunGeCo3 "PT_C15" "HS_C15"
   RunGeCo3 "PT_C16" "HS_C16"
   RunGeCo3 "PT_C17" "HS_C17"
   RunGeCo3 "PT_C18" "HS_C18"
   RunGeCo3 "PT_C19" "HS_C19"
   RunGeCo3 "PT_C20" "HS_C20"
   RunGeCo3 "PT_C21" "HS_C21"
   RunGeCo3 "PT_C22" "HS_C22"
   RunGeCo3 "PT_X" "HS_X"
   RunGeCo3 "PT_Y" "HS_Y"
   RunGeCo3 "PT_M" "HS_M"
   #
   RunGeCo3 "PA_C1" "HS_C1"
   RunGeCo3 "PA_C2" "HS_C2"
   RunGeCo3 "PA_C3" "HS_C3"
   RunGeCo3 "PA_C4" "HS_C4"
   RunGeCo3 "PA_C5" "HS_C5"
   RunGeCo3 "PA_C6" "HS_C6"
   RunGeCo3 "PA_C7" "HS_C7"
   RunGeCo3 "PA_C8" "HS_C8"
   RunGeCo3 "PA_C9" "HS_C9"
   RunGeCo3 "PA_C10" "HS_C10"
   RunGeCo3 "PA_C11" "HS_C11"
   RunGeCo3 "PA_C12" "HS_C12"
   RunGeCo3 "PA_C13" "HS_C13"
   RunGeCo3 "PA_C14" "HS_C14"
   RunGeCo3 "PA_C15" "HS_C15"
   RunGeCo3 "PA_C16" "HS_C16"
   RunGeCo3 "PA_C17" "HS_C17"
   RunGeCo3 "PA_C18" "HS_C18"
   RunGeCo3 "PA_C19" "HS_C19"
   RunGeCo3 "PA_C20" "HS_C20"
   RunGeCo3 "PA_C21" "HS_C21"
   RunGeCo3 "PA_C22" "HS_C22"
   RunGeCo3 "PA_X" "HS_X"
   RunGeCo3 "PA_M" "HS_M"
   #
   RunGeCo3 "GG_C1" "HS_C1"
   RunGeCo3 "GG_C2" "HS_C2"
   RunGeCo3 "GG_C3" "HS_C3"
   RunGeCo3 "GG_C4" "HS_C4"
   RunGeCo3 "GG_C5" "HS_C5"
   RunGeCo3 "GG_C6" "HS_C6"
   RunGeCo3 "GG_C7" "HS_C7"
   RunGeCo3 "GG_C8" "HS_C8"
   RunGeCo3 "GG_C9" "HS_C9"
   RunGeCo3 "GG_C10" "HS_C10"
   RunGeCo3 "GG_C11" "HS_C11"
   RunGeCo3 "GG_C12" "HS_C12"
   RunGeCo3 "GG_C13" "HS_C13"
   RunGeCo3 "GG_C14" "HS_C14"
   RunGeCo3 "GG_C15" "HS_C15"
   RunGeCo3 "GG_C16" "HS_C16"
   RunGeCo3 "GG_C17" "HS_C17"
   RunGeCo3 "GG_C18" "HS_C18"
   RunGeCo3 "GG_C19" "HS_C19"
   RunGeCo3 "GG_C20" "HS_C20"
   RunGeCo3 "GG_C21" "HS_C21"
   RunGeCo3 "GG_C22" "HS_C22"
   RunGeCo3 "GG_X" "HS_X"
   RunGeCo3 "GG_M" "HS_M"
   #
   cd ../../
   echo "Done!";
fi
###############################################################################
