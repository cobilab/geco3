#ifndef DEFS_H_INCLUDED
#define DEFS_H_INCLUDED

#include <stdint.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include <unistd.h>

#pragma pack(1)

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// RUNNING OPTIMIZATIONS : MEMORY / SPEED

#define ESTIMATE
#define PROGRESS
#define FSEARCHMODE

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// UNCOMMENT ONLY ONE:
//#define PREC32B // UNCOMMENT: CONTEXTS UP TO 28 (IT WILL USE HIGH MEMORY!)
#define PREC16B // UNCOMMENT: CONTEXTS UP TO 20 (IT WILL USE MEDIUM MEMORY!)
//#define PREC8B  // UNCOMMENT: CONTEXTS UP TO 16 (IT WILL USE LOW MEMORY!)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

typedef uint64_t ULL;
typedef uint64_t U64;
typedef uint32_t U32;
typedef uint16_t U16;
typedef uint8_t  U8;
typedef int64_t  I64;
typedef int32_t  I32;
typedef int16_t  I16;
typedef int8_t   I8;

typedef struct{
  U32      ctx;
  U32      den;
  U32      ir;
  U32      hashSize;
  double   gamma;
  U32      edits;
  U32      eDen;
  double   eGamma;
  U8       type;
  }
ModelPar;

typedef struct{
  U8       help;
  U8       verbose;
  U8       force;
  U8       estim;
  U8       level;
  ModelPar *model;
  char     *ref;
  char     **tar;
  U8       nTar;
  U64      checksum;
  U64      size;
  U32      watermark;
  U32      nModels;
  }
Parameters;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define RELEASE                0
#define VERSION                1

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define BITS_WATERMARK         32
#define BITS_CHECKSUM          46
#define BITS_SIZE              46
#define BITS_N_MODELS          16

#define BITS_TYPE              8
#define BITS_CTX               8
#define BITS_ALPHA_DEN         16
#define BITS_IR                8
#define BITS_GAMMA             32
#define BITS_HASH              8
#define BITS_EDITS             8
#define BITS_E_ALPHA_DEN       16
#define BITS_E_GAMMA           32

#define BITS_LR                32
#define BITS_HS                32

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define BUFFER_SIZE            262144
#define PROGRESS_MIN           1000000
#define DEF_VERSION            0
#define DEF_EXAMPLE            0
#define DEFAULT_HELP           0
#define DEFAULT_VERBOSE        0
#define DEFAULT_FORCE          0
#define DEFAULT_LEVEL          5
#define MAX_LEVEL              16
#define MIN_LEVEL              1
#define MAX_CTX                31
#define MIN_CTX                1
#define MAX_DEN                1000000
#define MIN_DEN                1
#define BGUARD                 32
#define DEFAULT_MAX_COUNT      ((1 << (sizeof(ACC) * 8)) - 1)
#define MX_PMODEL              65535
#define ALPHABET_SIZE          4
#define CHECKSUMGF             1073741824
#define WATERMARK              16112018
#define DEFAULT_GAMMA          0.90
#define MAX_HISTORYSIZE        1000000
#define REFERENCE              1
#define TARGET                 0
#define EXTRA_MOD_DEN          1
#define EXTRA_MOD_CTX          3
#define EXTRA_BIN_DEN          1
#define EXTRA_BIN_CTX          8
#define EXTRA_N_DEN            1
#define EXTRA_N_CTX            8
#define EXTRA_L_DEN            1
#define EXTRA_L_CTX            8
#define MAX_STR                2048

#define DEFAULT_LEARNING_RATE  0.03
#define DEFAULT_HIDDEN_SIZE    40
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#endif
