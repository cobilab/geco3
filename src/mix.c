#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "common.h"
#include "mix.h"
#include "fastonebigheader.h"

static inline float stretch(float p) {
  return fasterlog2(p / (1 - p));
}

mix_state_t* mix_init(uint32_t nmodels, uint32_t nsymbols, uint32_t hs) {
  mix_state_t *mxs = malloc(sizeof(mix_state_t));

  mxs->nmodels  = nmodels;
  mxs->nsymbols = nsymbols;
  int model_derived = 3; // hit, best, msym, bits
  int sequence_derived = 3; // last 7, last 15, last 50 symbols
  int xs = (nmodels * (nsymbols + model_derived)) + (sequence_derived * nsymbols) + 1;
  xs++; //bias neuron
  printf("xs: %d\n", xs);
  mxs->ann = ann_init(xs, hs, nsymbols);
  mxs->ann->x[xs - 1]  = 1.0;

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
  float *ret = mxs->ann->y;
  int i, j, k = 0;
  const float smean = mxs->smean;

  for(i = 0; i < mxs->nmodels; i++) {
    for(j = 0; j < mxs->nsymbols; j++) {
      mxs->ann->x[k++] = stretch(probs[i][j]) - smean;
    }

    mxs->ann->x[k++] = mxs->hit[i];
    mxs->ann->x[k++] = mxs->best[i];
    mxs->ann->x[k++] = mxs->bits[i];
  }

  mxs->ann->x[k++] = mxs->nnbits;

  int sf1[mxs->nsymbols];
  int sf2[mxs->nsymbols];
  int sf3[mxs->nsymbols];
  for(i = 0; i < mxs->nsymbols; ++i) {
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


  for(i = 0; i < mxs->nsymbols; ++i) {
    mxs->ann->x[k++] = (((float)sf1[i] / mxs->symlogs1) - 0.5) * 2;
    mxs->ann->x[k++] = (((float)sf2[i] / mxs->symlogs2) - 0.5) * 2;
    mxs->ann->x[k++] = (((float)sf3[i] / mxs->symlogs3) - 0.5) * 2;
  }

  ann_apply(mxs->ann);

  return ret;
}

void calc_aggregates(mix_state_t* mxs, float **probs, uint8_t sym) {
  int i, j;
  // Calculate aggregates
  // last N symbols (symbol log)
  for(i = mxs->symlogs3 - 1; i > 0 ; --i) {
    mxs->symlog[i] = mxs->symlog[i-1];
  }
  mxs->symlog[0] = sym;

  float bestp = probs[0][sym];
  int ignorant[mxs->nmodels];

  for(i = 0 ; i < mxs->nmodels; ++i) {
    if(bestp < probs[i][sym]) {
      bestp = probs[i][sym];
    }
    float p = probs[i][0];
    ignorant[i] = 1;
    for(j = 1 ; j < mxs->nsymbols; ++j) {
      if(p != probs[i][j]) {
	ignorant[i] = 0;
	break;
      }
    }
  }

  const float alpha = 0.5;
  const float nalpha = 1 - alpha;
  const float lmean = mxs->lmean;

  const float a = 0.15;
  const float na = 1 - a;

  float nnb = -fasterlog2(mxs->ann->y[sym]) + lmean;
  mxs->nnbits = (a * nnb) + (na * mxs->nnbits);

  for(i = 0 ; i < mxs->nmodels; ++i) {
    const float psym = probs[i][sym];
    const float t = -fasterlog2(psym) + lmean;
    mxs->bits[i] = (alpha * t) + (nalpha * mxs->bits[i]);

    if(!ignorant[i]) {
      float hit = 0.1;
      for(j = 0 ; j < mxs->nsymbols; ++j) {
	if(psym < probs[i][j]) {
	  hit = -0.1;
	  break;
	}
      }
      mxs->hit[i] += hit;
      mxs->hit[i] = fminf(1.0,  mxs->hit[i]);
      mxs->hit[i] = fmaxf(-1.0,  mxs->hit[i]);

      mxs->best[i] -= 0.1;
      mxs->best[i] += 0.2 * (bestp == psym);
      mxs->best[i] = fminf(1.0,  mxs->best[i]);
      mxs->best[i] = fmaxf(-1.0,  mxs->best[i]);

    }
  }
}

void mix_update_state(mix_state_t* mxs, float **probs, uint8_t sym, float learning_rate) {
  // Train NN
  float tdata[mxs->nsymbols];
  int i;
  for(i = 0 ; i < mxs->nsymbols; ++i) {
    tdata[i] = 0.0;
  }
  tdata[sym] = 1.0;

  ann_train(mxs->ann, tdata, learning_rate);
  calc_aggregates(mxs, probs, sym);
}

void mix_free(mix_state_t* mxs) {
  ann_free(mxs->ann);
  free(mxs->hit);
  free(mxs->best);
  free(mxs->bits);
  free(mxs->symlog);
  free(mxs);
}
