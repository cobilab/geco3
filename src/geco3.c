#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>
#include <time.h>
#include "mem.h"
#include "defs.h"
#include "msg.h"
#include "buffer.h"
#include "levels.h"
#include "common.h"
#include "pmodels.h"
#include "context.h"
#include "bitio.h"
#include "arith.h"
#include "arith_aux.h"
#include "mix.h"

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - C O M P R E S S O R - - - - - - - - - - - - - -

void Compress(Parameters *P, CModel **cModels, uint8_t id, INF *I, float lr, uint32_t hs){
  FILE        *Reader  = Fopen(P->tar[id], "r");
  char        *name    = concatenate(P->tar[id], ".co");
  FILE        *Writter = Fopen(name, "w");
  uint32_t    n, k, cModel, totModels, idxPos, j;
  uint64_t    compressed = 0, nSymbols = 0, nBases = 0;
  uint8_t     *readBUF, sym, irSym, *pos, type = 0,
              header = 1, line = 0, dna = 0;
  PModel      **pModel, *MX;
  FloatPModel *PT;
  CMWeight    *WM;
  CBUF        *symbBUF;

  #ifdef PROGRESS
  uint64_t    i = 0;
  #endif

  if(P->verbose)
    fprintf(stdout, "Analyzing data and creating models ...\n");

  #ifdef ESTIMATE
  FILE *IAE = NULL;
  char *IAEName = NULL;
  if(P->estim == 1){
    IAEName = concatenate(P->tar[id], ".iae");
    IAE = Fopen(IAEName, "w");
    }
  #endif

  sym = fgetc(Reader);
  switch(sym){
    case '>': type = 1; break;
    case '@': type = 2; break;
    default : type = 0;
    }
  rewind(Reader);

  switch(type){
    case 1:  nBases = NDNASymInFasta(Reader); break;
    case 2:  nBases = NDNASymInFastq(Reader); break;
    default: nBases = NDNASyminFile (Reader); break;
    }

  _bytes_output = 0;
  nSymbols      = NBytesInFile(Reader);

  // EXTRA MODELS DERIVED FROM EDITS
  totModels = P->nModels;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].edits != 0)
      totModels += 1;

  pModel        = (PModel  **) Calloc(totModels, sizeof(PModel *));
  for(n = 0 ; n < totModels ; ++n)
    pModel[n]   = CreatePModel(ALPHABET_SIZE);
  MX            = CreatePModel(ALPHABET_SIZE);
  PT            = CreateFloatPModel(ALPHABET_SIZE);
  WM            = CreateWeightModel(totModels);

  readBUF  = (uint8_t *) Calloc(BUFFER_SIZE, sizeof(uint8_t));
  symbBUF  = CreateCBuffer(BUFFER_SIZE, BGUARD);

  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == TARGET)
      cModels[n] = CreateCModel(TARGET, P->model[n].ctx, P->model[n].den,
                                P->model[n].ir, P->model[n].hashSize, P->model[n].gamma,
                                P->model[n].edits, P->model[n].eDen, P->model[n].eGamma);

  if(P->verbose)
    {
    fprintf(stdout, "Done!\nCompressing target sequence %d [bases:"
    "%"PRIu64"] ...\n", id + 1, nBases);
    }

  startoutputtingbits();
  start_encode();

  WriteNBits(WATERMARK,                               BITS_WATERMARK, Writter);
  WriteNBits(P->checksum,                              BITS_CHECKSUM, Writter);
  WriteNBits(nBases,                                       BITS_SIZE, Writter);
  WriteNBits(P->nModels,                               BITS_N_MODELS, Writter);
  for(n = 0 ; n < P->nModels ; ++n)
    {
    WriteNBits(P->model[n].type,                           BITS_TYPE, Writter);
    WriteNBits(cModels[n]->ctx,                             BITS_CTX, Writter);
    WriteNBits(cModels[n]->alphaDen,                  BITS_ALPHA_DEN, Writter);
    WriteNBits(cModels[n]->ir,                               BITS_IR, Writter);
    WriteNBits((int) (cModels[n]->gamma*65534),           BITS_GAMMA, Writter);
    WriteNBits(cModels[n]->hashSize,                       BITS_HASH, Writter);
    WriteNBits(cModels[n]->edits,                         BITS_EDITS, Writter);
    if(cModels[n]->edits > 0)
      {
      WriteNBits(cModels[n]->SUBS.eDen,             BITS_E_ALPHA_DEN, Writter);
      WriteNBits((int) (cModels[n]->SUBS.eGamma*65534), BITS_E_GAMMA, Writter);
      }
    }
  union {
    float    from;
    uint32_t to;
  } pun = { .from = lr };

  WriteNBits(pun.to,                                        BITS_LR, Writter);
  WriteNBits(hs,                                            BITS_HS, Writter);


  // GIVE SPECIFIC GAMMA:
  int pIdx = 0;
  for(n = 0 ; n < P->nModels ; ++n)
    {
    WM->gamma[pIdx++] = cModels[n]->gamma;
    if(P->model[n].edits != 0)
      {
      WM->gamma[pIdx++] = cModels[n]->SUBS.eGamma;
      }
    }

  I[id].header = _bytes_output;

  int nmodels = totModels + 1;
  float **probs = calloc(nmodels, sizeof(float *));
  for(n = 0; n < nmodels; ++n) {
    probs[n] = calloc(ALPHABET_SIZE, sizeof(float));
  }

  long **freqs = calloc(totModels, sizeof(long *));
  for (n = 0; n < totModels; ++n) {
    freqs[n] = calloc(ALPHABET_SIZE, sizeof(long));
  }
  long *sums = calloc(totModels, sizeof(long));

  mix_state_t *mxs = mix_init(nmodels, ALPHABET_SIZE, hs);

  while((k = fread(readBUF, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos)
      {

      #ifdef PROGRESS
      if(nSymbols > 100) {
	CalcProgress(nSymbols, ++i);
      }
      #endif

      sym = readBUF[idxPos];

      if(type == 1){  // IS A FAST[A] FILE
        if(sym == '>'){ header = 1; continue; }
        if(sym == '\n' && header == 1){ header = 0; continue; }
        if(sym == '\n') continue;
        if(header == 1) continue;
        }
      else if(type == 2){ // IS A FAST[Q] FILE
        switch(line){
          case 0: if(sym == '\n'){ line = 1; dna = 1; } break;
          case 1: if(sym == '\n'){ line = 2; dna = 0; } break;
          case 2: if(sym == '\n'){ line = 3; dna = 0; } break;
          case 3: if(sym == '\n'){ line = 0; dna = 0; } break;
          }
        if(dna == 0 || sym == '\n') continue;
        }

      // REMOVE SPECIAL SYMBOLS [WINDOWS TXT ISSUES]
      if(sym < 65 || sym > 122) continue;

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T'){
        #ifdef ESTIMATE
        if(P->estim != 0)
          fprintf(IAE, "0\n");
        #endif
        continue;
        }

      symbBUF->buf[symbBUF->idx] = sym = DNASymToNum(sym);
      memset((void *)PT->freqs, 0, ALPHABET_SIZE * sizeof(double));

      n = 0;
      pos = &symbBUF->buf[symbBUF->idx-1];
      for(cModel = 0 ; cModel < P->nModels ; ++cModel)
        {
        CModel *CM = cModels[cModel];
        GetPModelIdx(pos, CM);
        ComputePModel(CM, pModel[n], CM->pModelIdx, CM->alphaDen, freqs[n], &sums[n]);
        ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, 4);
        if(CM->edits != 0)
          {
          ++n;
          CM->SUBS.seq->buf[CM->SUBS.seq->idx] = sym;
          CM->SUBS.idx = GetPModelIdxCorr(CM->SUBS.seq->buf+
          CM->SUBS.seq->idx-1, CM, CM->SUBS.idx);
          ComputePModel(CM, pModel[n], CM->SUBS.idx, CM->SUBS.eDen, freqs[n], &sums[n]);
          ComputeWeightedFreqs(WM->weight[n], pModel[n], PT, 4);
          }
        ++n;
        }

      // Fill probabilities for all models

      for(n = 0 ; n < totModels ; ++n) {
        for(j = 0 ; j < ALPHABET_SIZE ; ++j) {
	  probs[n][j] = (float)freqs[n][j]/sums[n];
        }
      }

      for(j = 0 ; j < ALPHABET_SIZE ; ++j) {
	probs[totModels][j] = PT->freqs[j];
      }

      const float* y = mix(mxs, probs);

      for(n = 0 ; n < ALPHABET_SIZE ; ++n) {
        PT->freqs[n] = y[n];
      }

      ComputeMXProbs(PT, MX, 4);

      mix_update_state(mxs, probs, sym, lr);

      AESym(sym, (int *)(MX->freqs), (int) MX->sum, Writter);

      #ifdef ESTIMATE
      if(P->estim != 0)
        fprintf(IAE, "%.3g\n", PModelSymbolNats(MX, sym) / M_LN2);
      #endif

      CalcDecayment(WM, pModel, sym);

      for(n = 0 ; n < P->nModels ; ++n)
        {
        if(cModels[n]->ref == TARGET)
          {
          switch(cModels[n]->ir)
            {
            case 0:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            break;
            case 1:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            break;
            case 2:
            irSym = GetPModelIdxIR(symbBUF->buf+symbBUF->idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            break;
            default:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            break;
            }
          }
        }

      RenormalizeWeights(WM);

      n = 0;
      for(cModel = 0 ; cModel < P->nModels ; ++cModel)
        {
        if(cModels[cModel]->edits != 0)
          CorrectCModelSUBS(cModels[cModel], pModel[++n], sym);
        ++n;
        }


      UpdateCBuffer(symbBUF);
      ++compressed;
      }

  finish_encode(Writter);
  doneoutputtingbits(Writter);
  fclose(Writter);

  #ifdef ESTIMATE
  if(P->estim == 1)
    {
    fclose(IAE);
    Free(IAEName);
    }
  #endif

  RemovePModel(MX);
  Free(name);

  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
    else
      FreeCModel(cModels[n]);
  for(n = 0 ; n < totModels ; ++n){
    Free(pModel[n]->freqs);
    Free(pModel[n]);
    }
  Free(pModel);
  RemoveFPModel(PT);
  Free(readBUF);
  RemoveCBuffer(symbBUF);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stdout, "Done!                          \n");  // SPACES ARE VALID

  I[id].bytes = _bytes_output;
  I[id].size  = compressed;

  mix_free(mxs);
  free(sums);
  for (n = 0; n < totModels; ++n) {
    free(freqs[n]);
  }
  free(freqs);

  for(n = 0; n < nmodels; ++n) {
    free(probs[n]);
  }
  free(probs);

  RemoveWeightModel(WM);
}


//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - R E F E R E N C E - - - - - - - - - - - - -

CModel **LoadReference(Parameters *P)
  {
  FILE      *Reader = Fopen(P->ref, "r");
  uint32_t  n, k, idxPos;
  uint64_t  nBases = 0, y_bases = 0;
  int32_t   idx = 0;
  uint8_t   *readerBuffer, *symbolBuffer, sym, irSym = 0, type = 0, header = 1,
            line = 0, dna = 0;
  CModel    **cModels;
  #ifdef PROGRESS
  uint64_t  i = 0;
  #endif

  if(P->verbose == 1)
    fprintf(stdout, "Building reference model ...\n");

  readerBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + 1, sizeof(uint8_t));
  symbolBuffer  = (uint8_t *) Calloc(BUFFER_SIZE + BGUARD+1, sizeof(uint8_t));
  symbolBuffer += BGUARD;
  cModels       = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      cModels[n] = CreateCModel(REFERENCE, P->model[n].ctx, P->model[n].den,
      P->model[n].ir, P->model[n].hashSize, P->model[n].gamma,
      P->model[n].edits, P->model[n].eDen, P->model[n].eGamma);

  sym = fgetc(Reader);
  switch(sym){
    case '>': type = 1; break;
    case '@': type = 2; break;
    default : type = 0;
    }
  rewind(Reader);

  switch(type){
    case 1:  nBases = NDNASymInFasta(Reader); break;
    case 2:  nBases = NDNASymInFastq(Reader); break;
    default: nBases = NDNASyminFile (Reader); break;
    }

  y_bases = 0;
  P->checksum = 0;
  while((k = fread(readerBuffer, 1, BUFFER_SIZE, Reader)))
    for(idxPos = 0 ; idxPos < k ; ++idxPos){

      sym = readerBuffer[idxPos];
      if(type == 1){  // IS A FAST[A] FILE
        if(sym == '>'){ header = 1; continue; }
        if(sym == '\n' && header == 1){ header = 0; continue; }
        if(sym == '\n') continue;
        if(header == 1) continue;
        }
      else if(type == 2){ // IS A FAST[Q] FILE
        switch(line){
          case 0: if(sym == '\n'){ line = 1; dna = 1; } break;
          case 1: if(sym == '\n'){ line = 2; dna = 0; } break;
          case 2: if(sym == '\n'){ line = 3; dna = 0; } break;
          case 3: if(sym == '\n'){ line = 0; dna = 0; } break;
          }
        if(dna == 0 || sym == '\n') continue;
        }

      // FINAL FILTERING DNA CONTENT
      if(sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T')
        continue;

      symbolBuffer[idx] = sym = DNASymToNum(sym);
      P->checksum = (P->checksum + (uint8_t) sym);

      for(n = 0 ; n < P->nModels ; ++n)
        if(P->model[n].type == REFERENCE){

          GetPModelIdx(symbolBuffer+idx-1, cModels[n]);

	  // UPDATE ONLY IF IDX LARGER THAT CONTEXT
          switch(cModels[n]->ir)
            {
            case 0:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            break;
            case 1:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            break;
            case 2:
            irSym = GetPModelIdxIR(symbolBuffer+idx, cModels[n]);
            UpdateCModelCounter(cModels[n], irSym, cModels[n]->pModelIdxIR);
            break;
            default:
            UpdateCModelCounter(cModels[n], sym, cModels[n]->pModelIdx);
            break;
            }
          }
      ++y_bases;

      if(++idx == BUFFER_SIZE){
        memcpy(symbolBuffer - BGUARD, symbolBuffer + idx - BGUARD, BGUARD);
        idx = 0;
        }
      #ifdef PROGRESS
      if(nBases > 100) CalcProgress(nBases, ++i);
      #endif
      }

  P->checksum %= CHECKSUMGF;
  for(n = 0 ; n < P->nModels ; ++n)
    if(P->model[n].type == REFERENCE)
      ResetCModelIdx(cModels[n]);
  Free(readerBuffer);
  Free(symbolBuffer-BGUARD);
  fclose(Reader);

  if(P->verbose == 1)
    fprintf(stdout, "Done!                          \n");  // SPACES ARE VALID
  else
    fprintf(stdout, "                               \n");  // SPACES ARE VALID

  return cModels;
  }

//////////////////////////////////////////////////////////////////////////////
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - M A I N - - - - - - - - - - - - - - - - -
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

int32_t main(int argc, char *argv[]){
  char        **p = *&argv, **xargv, *xpl = NULL;
  CModel      **refModels;
  int32_t     n, xargc = 0;
  uint32_t    k, refNModels;
  uint64_t    totalBytes, headerBytes, totalSize;
  clock_t     stop = 0, start = clock();

  Parameters  *P;
  INF         *I;

  P = (Parameters *) Malloc(1 * sizeof(Parameters));
  if((P->help = ArgsState(DEFAULT_HELP, p, argc, "-h", "--help")) == 1
  || argc < 2){
    PrintMenuCompression();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_VERSION, p, argc, "-V", "--version")){
    PrintVersion();
    return EXIT_SUCCESS;
    }

  if(ArgsState(0, p, argc, "-s", "--show-levels")){
    PrintLevels();
    return EXIT_SUCCESS;
    }

  if(ArgsState(DEF_EXAMPLE, p, argc, "-x", "--examples")){
    PrintExamples();
    return EXIT_SUCCESS;
    }

  P->force    = ArgsState  (DEFAULT_FORCE,   p, argc, "-F", "--force");
  P->verbose  = ArgsState  (DEFAULT_VERBOSE, p, argc, "-v", "--verbose");
  P->estim    = ArgsState  (0,               p, argc, "-e", "--estimate");
  P->level    = ArgsNum    (0, p, argc, "-l", "--level", MIN_LEVEL, MAX_LEVEL);
  uint32_t hs = ArgsNum    (DEFAULT_HIDDEN_SIZE, p, argc, "-hs", "--hidden-size", 0, UINT32_MAX);
  float lr    = ArgsDouble (DEFAULT_LEARNING_RATE, p, argc, "-lr", "--learning-rate");

  printf("lr: %g, hs: %d\n", lr, hs);

  P->nModels  = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0 || strcmp(argv[n], "-tm") == 0)
      P->nModels += 1;

  if(P->nModels == 0 && P->level == 0)
    P->level = DEFAULT_LEVEL;

  if(P->level != 0){
    xpl = GetLevels(P->level);
    xargc = StrToArgv(xpl, &xargv);
    for(n = 1 ; n < xargc ; ++n) {
      if(strcmp(xargv[n], "-rm") == 0 || strcmp(xargv[n], "-tm") == 0)
        P->nModels += 1;
    }
  }

  if(P->nModels == 0){
    fprintf(stderr, "Error: at least you need to use a context model!\n");
    return 1;
    }

  P->model = (ModelPar *) Calloc(P->nModels, sizeof(ModelPar));

  k = 0;
  refNModels = 0;
  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-rm") == 0){
      P->model[k++] = ArgsUniqModel(argv[n+1], 1);
      ++refNModels;
      }
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-rm") == 0){
        P->model[k++] = ArgsUniqModel(xargv[n+1], 1);
        ++refNModels;
        }
    }

  for(n = 1 ; n < argc ; ++n)
    if(strcmp(argv[n], "-tm") == 0)
      P->model[k++] = ArgsUniqModel(argv[n+1], 0);
  if(P->level != 0){
    for(n = 1 ; n < xargc ; ++n)
      if(strcmp(xargv[n], "-tm") == 0)
        P->model[k++] = ArgsUniqModel(xargv[n+1], 0);
    }

  P->ref      = ArgsString (NULL, p, argc, "-r", "--reference");
  P->nTar     = ReadFNames (P, argv[argc-1]);
  P->checksum = 0;
  if(P->verbose)
    PrintArgs(P);

  if(refNModels == 0)
    refModels = (CModel **) Malloc(P->nModels * sizeof(CModel *));
  else{
    if(P->ref == NULL){
      fprintf(stderr, "Error: using reference model(s) in nonexistent "
      "reference sequence!\n");
      exit(1);
      }
    refModels = LoadReference(P);
    if(P->verbose)
      fprintf(stderr, "Checksum: %"PRIu64"\n", P->checksum);
    }

  I = (INF *) Calloc(P->nTar, sizeof(INF));

  totalSize   = 0;
  totalBytes  = 0;
  headerBytes = 0;
  for(n = 0 ; n < P->nTar ; ++n){
    Compress(P, refModels, n, I, lr, hs);
    totalSize   += I[n].size;
    totalBytes  += I[n].bytes;
    headerBytes += I[n].header;
    }

  if(P->nTar > 1)
    for(n = 0 ; n < P->nTar ; ++n){
      fprintf(stdout, "File %d compressed bytes: %"PRIu64" (", n+1, (uint64_t)
      I[n].bytes);
      PrintHRBytes(I[n].bytes);
      fprintf(stdout, ") , Normalized Dissimilarity Rate: %.6g\n",
      (8.0*I[n].bytes)/(2*I[n].size));
      }


  fprintf(stdout, "Total bytes: %"PRIu64" (", totalBytes);
  PrintHRBytes(totalBytes);
  fprintf(stdout, "), %.4g bpb, %.4g bps w/ no header, Normalized Dissimilarity"
  " Rate: %.6g\n", ((8.0*totalBytes)/totalSize), ((8.0*(totalBytes-headerBytes))
  /totalSize), (8.0*totalBytes)/(2.0*totalSize));
  stop = clock();
  fprintf(stdout, "Spent %g sec.\n", ((double)(stop-start))/CLOCKS_PER_SEC);

  Free(P->model);
  Free(P->tar);
  Free(P);

  if(refNModels == 0)
    Free(refModels);

  Free(I);

  return EXIT_SUCCESS;
  }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
