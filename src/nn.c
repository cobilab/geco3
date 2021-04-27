#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nn.h"

ann_t* ann_init(uint64_t xs, uint64_t hs, uint64_t ys) {
  ann_t *ann = malloc(sizeof(ann_t));
  ann->xs = xs+1; // bias node
  ann->hs = hs+1; // bias node
  ann->ys = ys;

  ann->x = calloc(ann->xs, sizeof(float));
  ann->h = calloc(ann->hs, sizeof(float));
  ann->y = calloc(ann->ys, sizeof(float));

  ann->x[ann->xs - 1]  = 1.0;
  ann->h[ann->hs - 1]  = 1.0;
  
  ann->wxh = calloc(ann->xs * ann->hs, sizeof(float));
  ann->why = calloc(ann->hs * ann->ys, sizeof(float));

  float r;
  float *w1 = ann->wxh;
  float *w2 = ann->why;
  int i;
  for(i = 0; i < ann->xs * ann->hs; ++i) {
    r = (((float)rand() / RAND_MAX) * 2.0) - 1.0;
    *w1++ = r * sqrtf(6.0 / (ann->xs + ann->hs));
  }

  for(i = 0; i < ann->hs * ann->ys; ++i) {
    r = (((float)rand() / RAND_MAX) * 2.0) - 1.0;
    *w2++ = r * sqrtf(6.0 / (ann->ys + ann->hs));
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

void ann_apply(ann_t *ann) {
  int i, j;
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t xhs = ann->hs - 1;
  const uint32_t ys = ann->ys;

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    const float xi = ann->x[i];
    for(j = 0; j < xhs; ++j) {
      ann->h[j] += *w1++ * xi;
    }
  }

  for(i = 0; i < hs; ++i) {
    ann->h[i] = sig(ann->h[i]);
  }
  
  float *w2 = ann->why;
  for(i = 0; i < hs; ++i) {
    const float hi = ann->h[i];
    for(j = 0; j < ys; ++j) {
      ann->y[j] += *w2++ * hi;
    }
  }

  float maxy = ann->y[0];
  for(i = 1; i < ys; ++i) {
    maxy = fmaxf(ann->y[i], maxy);
  }

  float total = 0;
  for(i = 0; i < ys; ++i) {
    ann->y[i] = fasterexp(ann->y[i] - maxy);
    total += ann->y[i];
  }

  for(i = 0; i < ys; ++i) {
    ann->y[i] /= total;
  }
}

void ann_train(ann_t *ann, float *t, float learning_rate) {
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t xhs = ann->hs - 1;
  const uint32_t ys = ann->ys;

  int i, j;
  float d1[ys];
  for(i = 0; i < ys; ++i) {
    d1[i] = learning_rate * (t[i] - ann->y[i]);
    ann->y[i] = 0;
  }

  float *w2 = ann->why;
  float d2[hs];
  for(i = 0; i < hs; ++i) {
    d2[i] = 0;
    const float hi = ann->h[i];
    for(j = 0; j < ys; ++j) {
      d2[i] += *w2 * d1[j];
      *w2++ += hi * d1[j];
    }
    ann->h[i] = 0;
    d2[i] *= hi * (1 - hi);
  }

  ann->h[xhs] = 1.0;

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    const float xi = ann->x[i];
    for(j = 0; j < xhs; ++j) {
      *w1++ += xi * d2[j];
    }
  }
}
