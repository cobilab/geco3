#ifndef ANN_H_
# define ANN_H_

#include "defs.h"
#include "fastonebigheader.h"

typedef struct ann_t {
  uint64_t xs;
  uint64_t hs;
  uint64_t ys;
  float *x;
  float *h;
  float *y;
  float *wxh;
  float *why;
} ann_t;

ann_t* ann_init(uint64_t xs, uint64_t hs, uint64_t ys);
void ann_apply(ann_t *ann);
void ann_train(ann_t *ann, float *t, float learning_rate);
void ann_free(ann_t *ann);

static inline float sig(float x) {
  if(x < -45) {
    return 0;
  } else {
    return fastsigmoid(x);
  }
}

#endif
