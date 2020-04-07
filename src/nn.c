#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nn.h"
#include "fastonebigheader.h"
#include "common.h"

ann_t* ann_init(uint64_t xs, uint64_t hs, uint64_t ys) {
  ann_t *ann = malloc(sizeof(ann_t));
  ann->xs = xs;
  ann->hs = hs;
  ann->ys = ys;
  ann->x = calloc(xs, sizeof(float));
  ann->h = calloc(hs, sizeof(float));
  ann->y = calloc(ys, sizeof(float));


  ann->wxh = calloc(xs * hs, sizeof(float));
  ann->why = calloc((hs+1) * ys, sizeof(float));

  float r;
  float *w1 = ann->wxh;
  float *w2 = ann->why;
  int i;
  for(i = 0; i < xs * hs; ++i) {
    r = (((float)rand() / RAND_MAX) * 2.0) - 1.0;
    *w1++ = r * sqrtf(6.0 / (xs + hs));
  }

  for(i = 0; i < (hs + 1) * ys; ++i) {
    r = (((float)rand() / RAND_MAX) * 2.0) - 1.0;
    *w2++ = r * sqrtf(6.0 / (ys + hs + 1));
  }

  return ann;
}

void ann_free(ann_t *ann) {
  free(ann->wxh);
  free(ann->why);
  free(ann->x);
  free(ann->h);
  free(ann->y);
  free(ann);
}

static inline float sig(float x) {
  if(x < -45) {
    return 0;
  } else {
    return fastsigmoid(x);
  }
}

void ann_apply(ann_t *ann) {
  int i, j;
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t ys = ann->ys;

  for(i = 0; i < hs; ++i) {
    ann->h[i] = 0;
  }

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    for(j = 0; j < hs; ++j) {
      ann->h[j] += *w1++ * ann->x[i];
    }
  }

  for(i = 0; i < hs; ++i) {
    ann->h[i] = sig(ann->h[i]);
  }

  for(i = 0; i < ys; ++i) {
    ann->y[i] = 0;
  }

  float *w2 = ann->why;
  for(i = 0; i < hs; ++i) {
    for(j = 0; j < ys; ++j) {
      ann->y[j] += *w2++ * ann->h[i];
    }
  }

  for(i = 0; i < ys; ++i) {
    ann->y[i] = sig(ann->y[i] + *w2++);
  }
}

void ann_train(ann_t *ann, float *t, float learning_rate) {
  int i, j;
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t ys = ann->ys;

  float d1[ys];
  for(i = 0; i < ys; ++i) {
    d1[i] = ann->y[i] * (1 - ann->y[i]) * (t[i] - ann->y[i]);
  }

  float *w2 = ann->why;
  float *w2t = ann->why;
  float d2[hs];
  for(i = 0; i < hs; ++i) {
    d2[i] = 0;
    for(j = 0; j < ys; ++j) {
      d2[i] += *w2t++ * d1[j];
      *w2++ += learning_rate * ann->h[i] * d1[j];
    }
    d2[i] = ann->h[i] * (1 - ann->h[i]) * d2[i];
  }

  for(i = 0; i < ys; ++i) {
    *w2++ += learning_rate * d1[i];
  }

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    for(j = 0; j < hs; ++j) {
      *w1++ += learning_rate * ann->x[i] * d2[j];
    }
  }
}
