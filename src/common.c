#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include "defs.h"
#include "mem.h"
#include "common.h"
#include "msg.h"

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

float *logTable;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t FLog2(uint64_t i)
  {
  uint32_t n, m, k = 32, o = (i & (i - 1)) ? 1 : 0;
  static const uint64_t sizes[6] = 
  { 0x0000000000000002ull, 0x000000000000000Cull, 0x00000000000000F0ull, 
    0x000000000000FF00ull, 0x00000000FFFF0000ull, 0xFFFFFFFF00000000ull };

  for(n = 6 ; n-- ; )
    {
    o += (m = (i & *(sizes+n)) ? k : 0);
    i >>= m;
    k >>= 1;
    }

  return o;
  }


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//        Pow function from http://martin.ankerl.com/2007/10/04/
//        optimized-pow-approximation-for-java-and-c-c/
double Power (double base, double exponent) 
  {
  union 
    {
    double d;
    int x[2];
    } u = { base };
  u.x[1] = (int) (exponent * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void ShiftBuffer(uint8_t *buf, int size, uint8_t new){
  memmove(buf, buf + 1, size * sizeof(uint8_t));
  buf[size - 1] = new;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// REPLACE STRING
//
char *ReplaceSubStr(char *str, char *a, char *b){
  char *buf = (char *) Calloc(MAX_STR, sizeof(char));
  char *p;
  if(strlen(str) > MAX_STR){
    fprintf(stderr, "[x] Error: string too long!\n");
    exit(1);
    }
  if(!(p = strstr(str, a)))
    return str;
  strncpy(buf, str, p-str);
  buf[p-str] = '\0';
  sprintf(buf+(p-str), "%s%s", b, p+strlen(a));
  return buf;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t NBytesInFile(FILE *file){
  uint64_t size = 0;
  fseeko(file, 0, SEEK_END);
  size = ftello(file);
  rewind(file);
  return size;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t NDNASyminFile(FILE *file){
  uint8_t  buffer[BUFFER_SIZE];
  uint32_t k, idx;
  uint64_t nSymbols = 0;

  while((k = fread(buffer, 1, BUFFER_SIZE, file)))
    for(idx = 0 ; idx < k ; ++idx)
      switch(buffer[idx]){
        case 'A': ++nSymbols; break;
        case 'T': ++nSymbols; break;
        case 'C': ++nSymbols; break;
        case 'G': ++nSymbols; break;
        }

  rewind(file);
  return nSymbols;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t NDNASymInFasta(FILE *file){
  uint8_t  buffer[BUFFER_SIZE], sym = 0, header = 1;
  uint32_t k, idx;
  uint64_t nSymbols = 0;

  while((k = fread(buffer, 1, BUFFER_SIZE, file)))
    for(idx = 0 ; idx < k ; ++idx){
      sym = buffer[idx];
      if(sym == '>'){ header = 1; continue; }
      if(sym == '\n' && header == 1){ header = 0; continue; }
      if(sym == '\n') continue;
      if(sym == 'N' ) continue;
      if(header == 1) continue;
      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
        continue;
      ++nSymbols;
      }

  rewind(file);
  return nSymbols;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t NDNASymInFastq(FILE *file){
  uint8_t  buffer[BUFFER_SIZE], sym = 0, line = 0, dna = 0;
  uint32_t k, idx;
  uint64_t nSymbols = 0;

  while((k = fread(buffer, 1, BUFFER_SIZE, file)))
    for(idx = 0 ; idx < k ; ++idx){
      sym = buffer[idx];

      switch(line){
        case 0: if(sym == '\n'){ line = 1; dna = 1; } break;
        case 1: if(sym == '\n'){ line = 2; dna = 0; } break;
        case 2: if(sym == '\n'){ line = 3; dna = 0; } break;
        case 3: if(sym == '\n'){ line = 0; dna = 0; } break;
        }
      if(dna == 0 || sym == '\n') continue;
      if(dna == 1 && sym == 'N' ) continue;

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
        continue;
      ++nSymbols;
      }

  rewind(file);
  return nSymbols;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint64_t FopenBytesInFile(const char *fn)
  {
  uint64_t size = 0;
  FILE *file = Fopen(fn, "r");
  
  size = NBytesInFile(file);  
  fclose(file);

  return size;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FillLogTable(uint32_t nSym, uint32_t alphaDen, uint32_t maxCHigh)
  {
  uint32_t n, maxSize = nSym * maxCHigh * alphaDen;
 
  logTable = (float *) Malloc(maxSize * sizeof(float));
  for(n = 1 ; n != maxSize ; ++n)
    logTable[n] = FLog2(n);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double SearchLog(uint32_t idx)
  {
  return logTable[idx];
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t DNASymToNum(uint8_t symbol)
  {
  switch(symbol)
    {
    case 'A': return 0;
    case 'T': return 3;
    case 'C': return 1;
    case 'G': return 2;
    default : return 4;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t NumToDNASym(uint8_t symbol)
  {
  switch(symbol)
    {
    case 0: return 'A';
    case 3: return 'T';
    case 1: return 'C';
    case 2: return 'G';
    default: fprintf(stderr, "Error: unknown numerical symbols\n"); exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t GetCompSym(uint8_t symbol)
  {
  switch(symbol)
    {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default:  return symbol;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t GetCompNum(uint8_t symbol)
  {
  switch(symbol)
    {
    case 0: return 3;
    case 1: return 2;
    case 2: return 1;
    case 3: return 0;
    default:  fprintf(stderr, "symbol\n");
    return symbol;
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

FILE *Fopen(const char *path, const char *mode)
  {
  FILE *file = fopen(path, mode);

  if(file == NULL)
    {
    fprintf(stderr, "Error opening: %s (mode %s). Does the file exist?\n", 
    path, mode);
    exit(1);
    }

  return file;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

uint8_t *ReverseStr(uint8_t *str, uint32_t end)
  {
  uint32_t start = 0;

  while(start < end)
    {
    str[start] ^= str[end];
    str[end]   ^= str[start];
    str[start] ^= str[end];
    ++start;
    --end;
    }

  return str;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void SortString(char str[])
  {
  char     tmp;
  int32_t  i, j, length = strlen(str);

  for(i = 0 ; i != length-1 ; ++i)
    for (j = i+1 ; j != length ; ++j)
      if(str[i] > str[j])
        {
        tmp    = str[i];
        str[i] = str[j];
        str[j] = tmp;
        }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *CloneString(char *str)
  {
  char *clone;
  if(str == NULL)
    return NULL;
  strcpy((clone = (char*) Malloc((strlen(str) + 1) * sizeof(char))), str);
  return clone;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *concatenate(char *a, char *b)
  {
  char *base;

  base = (char *) Malloc(strlen(a) + strlen(b) + 1);
  strcpy(base, a);
  strcat(base, b);
  return base;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *RepString(const char *str, const char *old, const char *new)
  {
  size_t sLen = strlen(str) + 1;
  char *cout = 0, *p = 0, *tmp; 

  if(!(p = (cout = (char *) Malloc(sLen * sizeof(char)))))
    return 0;
  while(*str)
    if((*str & 0xc0) != 0x80 && !strncmp(str, old, strlen(old)))
      {
      p   -= (intptr_t) cout;
      tmp  = strcpy(p = (cout = (char *) Realloc(cout, sLen += strlen(new) - 
             strlen(old), strlen(new) - strlen(old))) + (intptr_t) p, new);
      p   += strlen(tmp);
      str += strlen(old);
      }
    else
      *p++ = *str++;
  *p = 0;
  return cout;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t ArgsNum(uint32_t d, char *a[], uint32_t n, char *s, char *s2, 
uint32_t l, uint32_t u)
  {
  uint32_t x;
  for( ; --n ; ) if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
    {
    if((x = atol(a[n+1])) < l || x > u)
      {
      fprintf(stderr, "[x] Invalid number! Interval: [%u;%u].\n", l, u);
      exit(EXIT_FAILURE);
      }
    return x;
    }
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double ArgsDouble(double d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      return atof(a[n+1]);
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t ArgsState(uint8_t d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      return d == 0 ? 1 : 0;
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *ArgsString(char *d, char *a[], uint32_t n, char *s, char *s2)
  {
  for( ; --n ; )
    if(!strcmp(s, a[n]) || !strcmp(s2, a[n]))
      return a[n+1];
  return d;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

ModelPar ArgsUniqModel(char *str, uint8_t type)
  {
  uint32_t  ctx, den, ir, hash, edits, eDen, models_exit = 0;
  double    gamma, eGamma;
  ModelPar  Mp;
 

  if(sscanf(str, "%u:%u:%u:%u:%lf/%u:%u:%lf", &ctx, &den, &ir, &hash, &gamma, 
    &edits, &eDen, &eGamma) == 8)
    {

    if(ctx > MAX_CTX || ctx < MIN_CTX)
      {
      fprintf(stderr, "ERROR: Invalid Context!\n");
      models_exit = 1;
      }

    if(den > MAX_DEN || den < MIN_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator!\n");
      models_exit = 1;
      }

    if(ir < 0 || ir > 2)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator!\n");
      models_exit = 1;
      }

    if(hash < 0 || hash > 255)
      {
      fprintf(stderr, "ERROR: Invalid cache-hash size!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma!\n");
      models_exit = 1;
      }

    if(edits < 0 || edits > 20)
      {
      fprintf(stderr, "ERROR: Invalid number of editions (substitutions)!\n");
      models_exit = 1;
      }

    if(eDen < 0 || eDen > MAX_DEN)
      {
      fprintf(stderr, "ERROR: Invalid Alpha denominator (substitutions)!\n");
      models_exit = 1;
      }

    if(gamma < 0 || gamma > 0.999999)
      {
      fprintf(stderr, "ERROR: Invalid gamma (substitutions)!\n");
      models_exit = 1;
      }

    if(models_exit == 1)
      {
      ModelsExplanation();
      fprintf(stderr, "\nReset models according to the above description.\n");
      exit(1);
      }

    Mp.ctx      = ctx;
    Mp.den      = den;
    Mp.ir       = ir;
    Mp.hashSize = hash;
    Mp.gamma    = ((int)( gamma * 65534)) / 65534.0;
    Mp.eGamma   = ((int)(eGamma * 65534)) / 65534.0;
    Mp.edits    = edits;
    Mp.eDen     = eDen;
    Mp.type     = type;
    return Mp;
    }
  else{
    fprintf(stderr, "Error: unknown scheme for model arguments!\n");
    ModelsExplanation();
    fprintf(stderr, "\nReset models according to the above description.\n");
    exit(1);
    }

  return Mp;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

char *ArgsFiles(char *arg[], uint32_t argc, char *str)
  {
  int32_t n = argc;

  for( ; --n ; )
    if(!strcmp(str, arg[n]))
      return CloneString(arg[n+1]);
  
  return concatenate(concatenate(arg[argc-2], arg[argc-1]), ".svg");
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void FAccessWPerm(char *fn)
  {
  if(access(fn, F_OK) != -1)
    {
    fprintf(stderr, "Error: file %s already exists!\n", fn);
    if(access(fn, W_OK) != -1)
      fprintf(stderr, "Note: file %s has write permission.\nTip: to force "
      "writing rerun with \"-f\" option.\nWarning: on forcing, the old (%s) "
      "file will be deleted permanently.\n", fn, fn);
    exit(1);
    }
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void TestReadFile(char *fn)
  {
  FILE *f = NULL;
  f = Fopen(fn, "r");
  fclose(f);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// GET STRING TO ARGV
//
int32_t StrToArgv(char *s, char ***v){
  int32_t n = 0, c = 2;
  char *d = strdup(s);   // STRDUP <=> CLONE STR
  while(d[n]) if(d[n++] == ' ') ++c;
  *v = (char **) Calloc(c, sizeof(char *));
  n = 0; (*v)[0] = d; c = 1;
  do if(d[n] == ' '){
       d[n] = '\0';
       (*v)[c++] = d+n+1;
       }
  while(d[++n]);
  return c;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint32_t ReadFNames(Parameters *P, char *arg)
  {
  uint32_t nFiles = 1, k = 0, argLen;
  
  argLen = strlen(arg);
  for( ; k != argLen ; ++k)
    if(arg[k] == ':')
      ++nFiles;
  P->tar = (char **) Malloc(nFiles * sizeof(char *));
  P->tar[0] = strtok(arg, ":");
  TestReadFile(P->tar[0]);
  for(k = 1 ; k != nFiles ; ++k)
    {
    P->tar[k] = strtok(NULL, ":");
    TestReadFile(P->tar[k]);
    }

  return nFiles;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void CalcProgress(uint64_t size, uint64_t i)
  {
  if(i % (size / 100) == 0 && size > PROGRESS_MIN)
    fprintf(stderr, "Progress:%3d %%\r", (uint8_t) (i / (size / 100)));
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

uint8_t CmpCheckSum(uint32_t cs, uint32_t checksum)
  {
  if(checksum != cs)
    { 
    fprintf(stderr, "Error: invalid reference file!\n"
    "Compression reference checksum ................. %u\n"
    "Decompression reference checksum ............... %u\n",
    cs, checksum);
    fprintf(stderr, "Tip: rerun with correct reference file\n");
    return 1;
    }
  return 0;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

void PrintArgs(Parameters *P)
  {
  uint32_t n;

  fprintf(stderr, "Force mode ......................... %s\n", P->force == 0 ? 
  "no" : "yes");

  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == 1)
      {
      fprintf(stderr, "Reference model %d:\n", n+1);
      fprintf(stderr, "  [+] Context order ................ %u\n", 
      P->model[n].ctx);
      fprintf(stderr, "  [+] Alpha denominator ............ %u\n", 
      P->model[n].den);
      switch(P->model[n].ir)
        {
        case 0:
        fprintf(stderr, "  [+] Inverted repeats ............. no (only regular)\n");
        break;
        case 1: 
        fprintf(stderr, "  [+] Inverted repeats ............. mix (regular and inverted)\n");
        break;
        case 2:
        fprintf(stderr, "  [+] Inverted repeats ............. inverted only\n");
        break;
        }
      fprintf(stderr, "  [+] Cache-hash size .............. %u\n",
      P->model[n].hashSize);
      fprintf(stderr, "  [+] Gamma ........................ %.3lf\n",
      P->model[n].gamma);
      fprintf(stderr, "  [+] Allowable substitutions ...... %u\n",
      P->model[n].edits);
      if(P->model[n].edits != 0){
        fprintf(stderr, "  [+] Substitutions alpha den ...... %u\n",
        P->model[n].eDen);
        fprintf(stderr, "  [+] Substitutions gamma .......... %.3lf\n",
        P->model[n].eGamma);
        }
      }

  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == 0)
      {
      fprintf(stderr, "Target model %d:\n", n+1);
      fprintf(stderr, "  [+] Context order ................ %u\n",
      P->model[n].ctx);
      fprintf(stderr, "  [+] Alpha denominator ............ %u\n",
      P->model[n].den);
      switch(P->model[n].ir)
        {
        case 0:
        fprintf(stderr, "  [+] Inverted repeats ............. no (only regular)\n");
        break;
        case 1:
        fprintf(stderr, "  [+] Inverted repeats ............. mix (regular and inverted)\n");
        break;
        case 2:
        fprintf(stderr, "  [+] Inverted repeats ............. inverted only\n");
        break;
        }
      fprintf(stderr, "  [+] Cache-hash size .............. %u\n",
      P->model[n].hashSize);
      fprintf(stderr, "  [+] Gamma ........................ %.3lf\n",
      P->model[n].gamma);
      fprintf(stderr, "  [+] Allowable substitutions ...... %u\n",
      P->model[n].edits);
      if(P->model[n].edits != 0)
        {
        fprintf(stderr, "  [+] Substitutions alpha den ...... %u\n",
        P->model[n].eDen);
        fprintf(stderr, "  [+] Substitutions gamma .......... %.3lf\n",
        P->model[n].eGamma);
        }
      }

  if(P->ref != NULL)
    fprintf(stderr, "Reference filename ................. %s\n", P->ref);
  fprintf(stderr, "Target files (%u):\n", P->nTar);
  for(n = 0 ; n < P->nTar ; ++n)
    fprintf(stderr, "  [+] Filename %-2u .................. %s\n", n + 1, 
    P->tar[n]);
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
