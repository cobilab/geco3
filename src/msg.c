#include "msg.h"
#include <stdio.h>
#include <stdlib.h>

void ModelsExplanation(void){
  fprintf(stderr,
  "      -tm [NB_C]:[NB_D]:[NB_I]:[NB_H]:[NB_G]/[NB_S]:[NB_E]:[NB_A]       \n"
  "           Template of a target context model.                          \n"
  "           Parameters:                                                  \n"
  "           [NB_C]: (integer [1;20]) order size of the regular context   \n"
  "                   model. Higher values use more RAM but, usually, are  \n"
  "                   related to a better compression score.               \n"
  "           [NB_D]: (integer [1;5000]) denominator to build alpha, which \n"
  "                   is a parameter estimator. Alpha is given by 1/[NB_D].\n"
  "                   Higher values are usually used with higher [NB_C],   \n"
  "                   and related to confiant bets. When [NB_D] is one,    \n"
  "                   the probabilities assume a Laplacian distribution.   \n"
  "           [NB_I]: (integer {0,1,2}) number to define if a sub-program  \n"
  "                   which addresses the specific properties of DNA       \n"
  "                   sequences (Inverted repeats) is used or not. The     \n"
  "                   number 2 turns ON this sub-program without the       \n"
  "                   regular context model (only inverted repeats). The   \n"
  "                   number 1 turns ON the sub-program using at the same  \n"
  "                   time the regular context model. The number 0 does    \n"
  "                   not contemple its use (Inverted repeats OFF). The    \n"
  "                   use of this sub-program increases the necessary time \n"
  "                   to compress but it does not affect the RAM.          \n"
  "           [NB_H]: (integer [1;254]) size of the cache-hash for deeper  \n"
  "                   context models, namely for [NB_C] > 14. When the     \n"
  "                   [NB_C] <= 14 use, for example, 1 as a default. The   \n"
  "                   RAM is highly dependent of this value (higher value  \n"
  "                   stand for higher RAM).                               \n"
  "           [NB_G]: (real [0;1)) real number to define gamma. This value \n"
  "                   represents the decayment forgetting factor of the    \n"
  "                   regular context model in definition.                 \n"
  "           [NB_S]: (integer [0;20]) maximum number of editions allowed  \n"
  "                   to use a substitutional tolerant model with the same \n"
  "                   memory model of the regular context model with       \n"
  "                   order size equal to [NB_C]. The value 0 stands for   \n"
  "                   turning the tolerant context model off. When the     \n"
  "                   model is on, it pauses when the number of editions   \n"
  "                   is higher that [NB_C], while it is turned on when    \n"
  "                   a complete match of size [NB_C] is seen again. This  \n"
  "                   is probabilistic-algorithmic model very usefull to   \n"
  "                   handle the high substitutional nature of genomic     \n"
  "                   sequences. When [NB_S] > 0, the compressor used more \n"
  "                   processing time, but uses the same RAM and, usually, \n"
  "                   achieves a substantial higher compression ratio. The \n"
  "                   impact of this model is usually only noticed for     \n"
  "                   [NB_C] >= 14.                                        \n"
  "           [NB_E]: (integer [1;5000]) denominator to build alpha for    \n"
  "                   substitutional tolerant context model. It is         \n"
  "                   analogous to [NB_D], however to be only used in the  \n"
  "                   probabilistic model for computing the statistics of  \n"
  "                   the substitutional tolerant context model.           \n"
  "           [NB_A]: (real [0;1)) real number to define gamma. This value \n"
  "                   represents the decayment forgetting factor of the    \n"
  "                   substitutional tolerant context model in definition. \n"
  "                   Its definition and use is analogus to [NB_G].        \n"
  "                                                                        \n"
  "      ... (you may use several target models with custom parameters)    \n"
  "                                                                        \n"
  "      -rm [NB_C]:[NB_D]:[NB_I]:[NB_H]:[NB_G]/[NB_S]:[NB_E]:[NB_A]       \n"
  "           Template of a reference context model.                       \n"
  "           Use only when -r [FILE] is set (referential compression).    \n"
  "           Parameters: the same as in -tm.                              \n"
  "                                                                        \n"
  "      ... (you may use several reference models with custom parameters) \n"
  "                                                                        \n");
  }

void PrintMenuCompression(void){
  fprintf(stderr,
  "                                                                        \n"
  "                                                                        \n"
  "            ██████╗ ███████╗ ██████╗ ██████╗ ██████╗                    \n"
  "           ██╔════╝ ██╔════╝██╔════╝██╔═══██╗╚════██╗                   \n"
  "           ██║  ███╗█████╗  ██║     ██║   ██║ █████╔╝                   \n"
  "           ██║   ██║██╔══╝  ██║     ██║   ██║ ╚═══██╗                   \n"
  "           ╚██████╔╝███████╗╚██████╗╚██████╔╝██████╔╝                   \n"
  "            ╚═════╝ ╚══════╝ ╚═════╝ ╚═════╝ ╚═════╝                    \n"
  "                                                                        \n"
  "NAME                                                                    \n"
  "      GeCo3 v%u.%u,                                                     \n"
  "      efficient compression and analysis of genomic sequences.          \n"
  "                                                                        \n"
  "AUTHORS                                                                 \n"
  "      Diogo Pratas        pratas@ua.pt                                  \n"
  "      Milton Silva        teixeirasilva@ua.pt                           \n"
  "      Morteza Hosseini    seyedmorteza@ua.pt                            \n"
  "      Armando J. Pinho    ap@ua.pt                                      \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      ./GeCo3 [OPTION]... -r [FILE] [FILE]:[FILE]:[FILE]:[...]          \n"
  "                                                                        \n"
  "SAMPLE                                                                  \n"
  "      Run Compression         :  ./GeCo3 -v -l 3 sequence.txt           \n"
  "      Run Decompression       :  ./GeDe3 -v sequence.txt.co             \n"
  "      Run Information Profile :  ./GeCo3 -v -l 3 -e sequence.txt        \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      Compress and decompress genomic sequences for storage purposes.   \n"
  "      Measure an upper bound of the sequences entropy.                  \n"
  "      Compute information profiles of genomic sequences.                \n"
  "                                                                        \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu).                                     \n"
  "                                                                        \n"
  "      -V,  --version                                                    \n"
  "           Display program and version information.                     \n"
  "                                                                        \n"
  "      -F,  --force                                                      \n"
  "           force mode. Overwrites old files.                            \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information).                             \n"
  "                                                                        \n"
  "      -x,  --examples                                                   \n"
  "           show several running examples (parameter examples).          \n"
  "                                                                        \n"
  "      -s,  --show-levels                                                \n"
  "           show pre-computed compression levels (configured parameters).\n"
  "                                                                        \n",
  VERSION, RELEASE);

  #ifdef ESTIMATE
  fprintf(stderr,
  "      -e,  --estimate                                                   \n"
  "           it creates a file with the extension \".iae\" with the       \n"
  "           respective information content. If the file is FASTA or      \n"
  "           FASTQ it will only use the \"ACGT\" (genomic) sequence.      \n"
  "                                                                        \n");
  #endif

  fprintf(stderr,
  "      -l [NUMBER],  --level [NUMBER]                                    \n"
  "           Compression level (integer [1;16]).                          \n"
  "           Default level: %u.                                           \n"
  "           It defines compressibility in balance with computational     \n"
  "           resources (RAM & time). Use -s for levels perception.        \n"
  "                                                                        \n",
  DEFAULT_LEVEL);

  fprintf(stderr,
  "      -lr [NUMBER],  --learning-rate [NUMBER]                           \n"
  "           Learning rate (real).                                        \n"
  "           Default rate: %g.                                            \n"
  "           It defines learning rate the neural network uses.            \n"
  "                                                                        \n",
  DEFAULT_LEARNING_RATE);

  fprintf(stderr,
  "      -hs [NUMBER],  --hidden-size [NUMBER]                             \n"
  "           Hidden layer size (integer).                                 \n"
  "           Default size: %d.                                            \n"
  "           It defines number of hidden nodes for the neural network.    \n"
  "                                                                        \n",
  DEFAULT_HIDDEN_SIZE);

  ModelsExplanation();

  fprintf(stderr,
  "      -r [FILE], --reference [FILE]                                     \n"
  "           Reference sequence filename (\"-rm\" are trainned here).     \n"
  "           Example: -r file1.txt.                                       \n"
  "                                                                        \n"
  "      [FILE]                                                            \n"
  "           Input sequence filename (to compress) -- MANDATORY.          \n"
  "           File(s) to compress (last argument).                         \n"
  "           For more files use splitting \":\" characters.               \n"
  "           Example: file1.txt:file2.txt:file3.txt.                      \n"
  "                                                                        \n"
  "COPYRIGHT                                                               \n"
  "      Copyright (C) 2014-2020, IEETA, University of Aveiro.             \n"
  "      This is a Free software, under GPLv3. You may redistribute        \n"
  "      copies of it under the terms of the GNU - General Public          \n"
  "      License v3 <http://www.gnu.org/licenses/gpl.html>. There          \n"
  "      is NOT ANY WARRANTY, to the extent permitted by law.              \n"
  "                                                                        \n");
  }

void PrintMenuDecompression(void){
  fprintf(stderr,
  "                                                                        \n"
  "                                                                        \n"
  "             ██████╗ ███████╗██████╗ ███████╗██████╗                    \n"
  "            ██╔════╝ ██╔════╝██╔══██╗██╔════╝╚════██╗                   \n"
  "            ██║  ███╗█████╗  ██║  ██║█████╗   █████╔╝                   \n"
  "            ██║   ██║██╔══╝  ██║  ██║██╔══╝   ╚═══██╗                   \n"
  "            ╚██████╔╝███████╗██████╔╝███████╗██████╔╝                   \n"
  "             ╚═════╝ ╚══════╝╚═════╝ ╚══════╝╚═════╝                    \n"
  "                                                                        \n"
  "NAME                                                                    \n"
  "      GeDe3 v%u.%u,                                                     \n"
  "      decompress a genomic sequence compressed by GeCo3.                \n"
  "                                                                        \n"
  "AUTHORS                                                                 \n"
  "      Diogo Pratas        pratas@ua.pt                                  \n"
  "      Milton Silva        teixeirasilva@ua.pt                           \n"
  "      Morteza Hosseini    seyedmorteza@ua.pt                            \n"
  "      Armando J. Pinho    ap@ua.pt                                      \n"
  "                                                                        \n"
  "SYNOPSIS                                                                \n"
  "      ./GeDe3 [OPTION]... -r [FILE] [FILE]:[FILE]:[FILE]:[...]          \n"
  "                                                                        \n"
  "SAMPLE                                                                  \n"
  "      Run Decompression:  ./GeDe3 -v sequence.txt.co                    \n"
  "                                                                        \n"
  "DESCRIPTION                                                             \n"
  "      Decompress genomic sequences for compressed by GeCo3.             \n"
  "                                                                        \n"
  "      -h,  --help                                                       \n"
  "           usage guide (help menu).                                     \n"
  "                                                                        \n"
  "      -V,  --version                                                    \n"
  "           Display program and version information.                     \n"
  "                                                                        \n"
  "      -F,  --force                                                      \n"
  "           force mode. Overwrites old files.                            \n"
  "                                                                        \n"
  "      -v,  --verbose                                                    \n"
  "           verbose mode (more information).                             \n"
  "                                                                        \n"
  "      -r [FILE], --reference [FILE]                                     \n"
  "           Reference sequence filename (\"-rm\" are trainned here).     \n"
  "           Example: -r file1.txt.                                       \n"
  "                                                                        \n"
  "      [FILE]                                                            \n"
  "           Input compressed filename (to decompress) -- MANDATORY.      \n"
  "           File(s) to decompress (last argument).                       \n"
  "           For more files use splitting \":\" characters.               \n"
  "           Example: file1.txt:file2.txt:file3.txt.                      \n"
  "                                                                        \n"
  "COPYRIGHT                                                               \n"
  "      Copyright (C) 2014-2020, IEETA, University of Aveiro.             \n"
  "      This is a Free software, under GPLv3. You may redistribute        \n"
  "      copies of it under the terms of the GNU - General Public          \n"
  "      License v3 <http://www.gnu.org/licenses/gpl.html>. There          \n"
  "      is NOT ANY WARRANTY, to the extent permitted by law.              \n"
  "                                                                        \n",
  VERSION, RELEASE);
  }

void PrintVersion(void){
  fprintf(stderr,
  "                                                                        \n"
  "                          ==================                            \n"
  "                          |    GeCo3 %u.%u   |                          \n"
  "                          ==================                            \n"
  "                                                                        \n"
  "                  An efficient tool for compression                     \n"
  "                  and analysis of genomic sequences                     \n"
  "                                                                        \n"
  "             Copyright (C) 2014-2020 University of Aveiro.              \n"
  "                                                                        \n"
  "                This is a Free software, under GPLv3.                   \n"
  "                                                                        \n"
  "You may redistribute copies of it under the terms of the GNU - General  \n"
  "Public License v3 <http://www.gnu.org/licenses/gpl.html>. There is NOT  \n"
  "ANY WARRANTY, to the extent permitted by law. Developed and Written by  \n"
  "Diogo Pratas, Morteza Hosseini and Armando J. Pinho.\n\n", VERSION, RELEASE);
  }


void PrintExamples(void){
  fprintf(stderr,
  "                                                                       \n"
  "GeCo3 running examples:                                                \n"
  "                                                                       \n"
  "Considerations: the decompression is symmetric, therefore the same     \n"
  "resources, namely time and memory will be used as in the compression.  \n"
  "The memory used, after creating the models, will be constant, even in  \n"
  "deeper context models (cache-hash context model).                      \n"
  "                                                                       \n"
  "[A]=> Compressing sequences C(X) or C(X,Y):                            \n"
  "                                                                       \n"
  /*"1) Compression of a human genome (using 5.8 GB RAM memory):            \n"
  "   ./GeCo -tm 6:1:0:0/0 -tm 13:20:1:0/0 -tm 19:50:1:2/10 -c 35 -g 0.8 HS\n"
  "                                                                       \n"
  "5) Highly-redundant genomic sequence (full ACGT from fastq)            \n"
  "   ./GeCo -tm 4:1:0:0/0 -tm 11:1:0:0/0 -tm 14:20:0:0/0 -tm 20:100:0:1/10\n"
  "   -c 40 -g 0.8 SRR957627.fastq                                        \n"
  "   Decompression for A5: ./GeDe SRR957627.fastq.co                     \n"
  "   The decompressed file will be SRR957627.fastq.de                    \n"
  "                                                                       \n"
  "[B]=> Conditional (referential) exclusive compression C(X||Y):         \n"
  "                                                                       \n"
  "1) Compression of the gorilla (GG8) chromosome 8 given exclusively     \n"
  "   information from chimpanzee (PT8):                                  \n"
  "   ./GeCo -rm 4:1:0:0/0 -rm 20:1000:1:1/100 -c 20 -r PT8 GG8           \n"
  "   Decompression for B1: ./GeDe -r PT8 GG8.co                          \n"
  "   The decompressed file will be GG8.de                                \n"
  "                                                                       \n"
  "2) Compression of the same file (for identity studies):                \n"
  "   ./GeCo -rm 20:1000:0:0/0 -c 30 -r File1.txt File1.txt               \n"
  "   Decompression for B2: ./GeDe -r File1.txt File1.txt.co              \n"
  "   The decompressed file will be File1.txt.de                          \n"
  "                                                                       \n"
  "3) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
  "   chromsomes 5 given exclusively the gorilla (GG17) chromosome 17 as  \n"
  "   reference:                                                          \n"
  "   ./GeCo -rm 20:1000:1:1/100 -c 20 -r GG17 HS5:PT5:PA5                \n"
  "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
  "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"
  "                                                                       \n"
  "[C]=> Conditional compression C(X|Y) [use reference and target]:       \n"
  "                                                                       \n"
  "1) Compression of a human (HS5), chimpanzee (PT5) and orangutan (PA5)  \n"
  "   chromsomes 5 given the gorilla (GG17) chromosome 17 as reference:   \n"
  "   -rm 12:100:1:0/0 -rm 20:1000:1:1/100 -tm 4:1:0:0/0 -tm 14:20:1:1/10 \n"
  "   -c 20 -g 0.85 -r GG17 HS5:PT5:PA5                                   \n"
  "   Decompression for B3: ./GeDe -r GG17 HS5.co:PT5.co:PA5.co           \n"
  "   The decompressed files will be HS5.de, PT5.de and PA5.de            \n"*/
  "                                                                     \n");
  }
