#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "common.h"
#include "mix.h"

mix_state_t* mix_init(uint32_t nmodels, uint32_t nsymbols, uint32_t hs) {
  mix_state_t *mxs = malloc(sizeof(mix_state_t));

  mxs->nmodels  = nmodels;
  mxs->nsymbols = nsymbols;
  int model_derived = 3; // hit, best, msym, bits
  int sequence_derived = 3; // last 7, last 15, last 50 symbols
  mxs->xs = (nmodels * (nsymbols + model_derived)) + (sequence_derived * nsymbols) + 1;
  mxs->x = calloc(mxs->xs, sizeof(double));
  mxs->y = calloc(4, sizeof(float));
  printf("xs: %d\n", mxs->xs);

  mxs->ann = genann_init(mxs->xs, 1, hs, nsymbols);//ann_init(xs, hs, nsymbols);

  // past symbols
  mxs->symlogs1 = 8; // empirically determined
  mxs->symlogs2 = 16; // empirically determined
  mxs->symlogs3 = 64; // empirically determined
  mxs->symlog = calloc(mxs->symlogs3, sizeof(uint8_t));

  
  // model performance
  mxs->hit = calloc(nmodels, sizeof(float));
  mxs->best = calloc(nmodels, sizeof(float));
  mxs->bits = calloc(nmodels, sizeof(float));

  mxs->nnbits = 0;

  float mean = 1.0 / mxs->nsymbols;
  mxs->smean = stretch(mean);
  mxs->lmean = fasterlog2(mean);

  return mxs;
}

float const* mix(mix_state_t* mxs, float **probs) {
  int i, j, k = 0;
  const float smean = mxs->smean;
  const uint32_t nmodels = mxs->nmodels;
  const uint32_t nsymbols = mxs->nsymbols;
  

  for(i = 0; i < nmodels; ++i) {
    for(j = 0; j < nsymbols; ++j) {
      mxs->x[k++] = stretch(probs[i][j]) - smean;
    }

    mxs->x[k++] = mxs->hit[i];
    mxs->x[k++] = mxs->best[i];
    mxs->x[k++] = mxs->bits[i];
  }

  mxs->x[k++] = mxs->nnbits;

  int sf1[nsymbols];
  int sf2[nsymbols];
  int sf3[nsymbols];
  for(i = 0; i < nsymbols; ++i) {
    sf1[i] = 0;
    sf2[i] = 0;
    sf3[i] = 0;
  }

  for(i = 0; i < mxs->symlogs1; ++i) {
    sf1[mxs->symlog[i]]++;
  }

  for(; i < mxs->symlogs2; ++i) {
    sf2[mxs->symlog[i]]++;
  }

  for(i = 0; i < mxs->symlogs3; ++i) {
    sf3[mxs->symlog[i]]++;
  }


  for(i = 0; i < nsymbols; ++i) {
    mxs->x[k++] = (((float)sf1[i] / mxs->symlogs1) - 0.5) * 2;
    mxs->x[k++] = (((float)sf2[i] / mxs->symlogs2) - 0.5) * 2;
    mxs->x[k++] = (((float)sf3[i] / mxs->symlogs3) - 0.5) * 2;
  }

  double const *prediction = genann_run(mxs->ann, mxs->x);
  
  for(i = 0; i < nsymbols; ++i) {
    mxs->y[i] = prediction[i];
  }

  return mxs->y;
}

void calc_aggregates(mix_state_t* mxs, float **probs, uint8_t sym, float psymbol) {
  const float lmean = mxs->lmean;
  const float a = 0.15;
  const float na = 1 - a;

  const float nnb = -fasterlog2(psymbol) + lmean;
  mxs->nnbits = (a * nnb) + (na * mxs->nnbits);

  const uint32_t nmodels = mxs->nmodels;
  const uint32_t nsymbols = mxs->nsymbols;

  // last N symbols (symbol log)
  int i, j;
  for(i = mxs->symlogs3 - 1; i > 0 ; --i) {
    mxs->symlog[i] = mxs->symlog[i-1];
  }
  mxs->symlog[0] = sym;

  const float alpha = 0.5;
  const float nalpha = 1 - alpha;

  float bestp = probs[0][sym];
  int ignorant[nmodels];
  int hit[nmodels];

  for(i = 0 ; i < nmodels; ++i) {
    const float psym = probs[i][sym];
    bestp = (bestp < psym) ? psym : bestp;

    const float p = probs[i][0];
    ignorant[i] = 1;
    hit[i] = (psym >= p);
    for(j = 1 ; j < nsymbols; ++j) {
      ignorant[i] &= (p == probs[i][j]);
      hit[i] &= (psym >= probs[i][j]);
    }
  }

  for(i = 0 ; i < nmodels; ++i) {
    const float psym = probs[i][sym];
    const float t = -fasterlog2(psym) + lmean;
    mxs->bits[i] = (alpha * t) + (nalpha * mxs->bits[i]);

    if(!ignorant[i]) {
      mxs->best[i] += (bestp == psym) ? 0.1 : -0.1;
      mxs->best[i] = fminf(1.0,  mxs->best[i]);
      mxs->best[i] = fmaxf(-1.0,  mxs->best[i]);

      mxs->hit[i] += hit[i] ? 0.1 : -0.1;
      mxs->hit[i] = fminf(1.0,  mxs->hit[i]);
      mxs->hit[i] = fmaxf(-1.0,  mxs->hit[i]);
    }
  }
}

void mix_update_state(mix_state_t* mxs, float **probs, uint8_t sym, float learning_rate, float psymbol) {
  calc_aggregates(mxs, probs, sym, psymbol);
  // Train NN
  double tdata[mxs->nsymbols];
  int i;
  for(i = 0 ; i < mxs->nsymbols; ++i) {
    tdata[i] = 0.0;
  }
  tdata[sym] = 1.0;

  genann_train(mxs->ann, mxs->x, tdata, learning_rate);
}

void mix_free(mix_state_t* mxs) {
  
  free(mxs->hit);
  free(mxs->best);
  free(mxs->bits);
  free(mxs->symlog);
  free(mxs);
}
