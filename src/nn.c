#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nn.h"
#include "common.h"


static inline float Q_rsqrt( float number )
{	
	const float x2 = number * 0.5F;
	const float threehalfs = 1.5F;

	union {
		float f;
		unsigned long i;
	} conv  = { .f = number };
	conv.i  = 0x5f3759df - ( conv.i >> 1 );
	conv.f  *= ( threehalfs - ( x2 * conv.f * conv.f ) );
	return conv.f;
}

ann_t* ann_init(uint64_t xs, uint64_t hs, uint64_t ys) {
  ann_t *ann = malloc(sizeof(ann_t));
  ann->xs = xs;
  ann->hs = hs;
  ann->ys = ys;

  ann->x = calloc(xs, sizeof(float));
  ann->h = calloc(hs, sizeof(float));
  ann->y = calloc(ys, sizeof(float));

  ann->m = calloc(xs * (hs+1) * ys, sizeof(float));
  ann->v = calloc(xs * (hs+1) * ys, sizeof(float));

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

void ann_apply(ann_t *ann) {
  int i, j;
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t ys = ann->ys;

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    const float xi = ann->x[i];
    for(j = 0; j < hs; ++j) {
      ann->h[j] += *w1++ * xi;
    }
  }

  float *w2 = ann->why;
  for(i = 0; i < hs; ++i) {
    const float hi = sig(ann->h[i]);
    for(j = 0; j < ys; ++j) {
      ann->y[j] += *w2++ * hi;
    }
    ann->h[i] = hi;
  }

  for(i = 0; i < ys; ++i) {
    ann->y[i] = sig(ann->y[i] + *w2++);
  }
}

const float b1 = 0.9;
const float nb1 = 1- b1;
const float b2 = 0.999;
const float nb2 = 1 - b2;
const float eps = 1E-8;

const float bt1 = b1;
const float bt2 = b2;

const float ibt1 = 1 / (1 - bt1);
const float ibt2 = 1 / (1 - bt2);

void ann_train(ann_t *ann, float *t, float learning_rate) {

  float *m = ann->m;
  float *v = ann->v;
  int i, j;
  const uint32_t xs = ann->xs;
  const uint32_t hs = ann->hs;
  const uint32_t ys = ann->ys;

  float d1[ys];
  for(i = 0; i < ys; ++i) {
    d1[i] = ann->y[i] * (1 - ann->y[i]) * (t[i] - ann->y[i]);
  }

  float *w2 = ann->why;
  float d2[hs];
  for(i = 0; i < hs; ++i) {
    d2[i] = 0;
    const float hi = ann->h[i];
    
    for(j = 0; j < ys; ++j) {
      d2[i] += *w2 * d1[j];
      *m = (b1 * *m) + (nb1 * hi * d1[j]);
      *v = (b2 * *v) + (nb2 * hi * d1[j] * hi * d1[j]);

      float mhat = *m * ibt1;
      float vhat = *v * ibt2;
      
      *w2++ += learning_rate * mhat * Q_rsqrt(vhat);

      m++;
      v++;
    }
    ann->h[i] = 0;
    d2[i] *= hi * (1 - hi);
  }

  for(i = 0; i < ys; ++i) {
    *m = (b1 * *m) + (nb1 * d1[i]);
    *v = (b2 * *v) + (nb2 * d1[i] * d1[i]);

    float mhat = *m * ibt1;
    float vhat = *v * ibt2;
    
      
    *w2++ +=  learning_rate * mhat * Q_rsqrt(vhat);
    ann->y[i] = 0;
    m++;
    v++;
  }

  float *w1 = ann->wxh;
  for(i = 0; i < xs; ++i) {
    const float nxi1 = nb1 * ann->x[i];
    const float nxi2 = nb2 * ann->x[i] * ann->x[i];
    for(j = 0; j < hs; ++j) {
      *m = (b1 * *m) + (nxi1 * d2[j]);
      *v = (b2 * *v) + (nxi2 * d2[j] * d2[j]);
      
      float mhat = *m * ibt1;
      float vhat = *v * ibt2;
      *w1++ +=  learning_rate * mhat * Q_rsqrt(vhat);
      //printf("%f\n", (learning_rate) / (sqrt(vhat) + eps));
      m++;
      v++;
    }
  }

  //bt1 *= bt1;
  //bt2 *= bt2;
}
