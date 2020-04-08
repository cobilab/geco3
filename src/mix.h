#include "nn.h"

typedef struct mix_state_t {
  uint32_t nmodels;
  uint32_t nsymbols;

  ann_t *ann;

  float *hit;
  float *best;
  float *bits;
  float nnbits;
  uint8_t *symlog;

  float smean;
  float lmean;

  uint8_t symlogs1;
  uint8_t symlogs2;
  uint8_t symlogs3;
} mix_state_t;

mix_state_t* mix_init(uint32_t nmodels, uint32_t nsymbols, uint32_t hs);
float const* mix(mix_state_t* mxs, float **probs);
void mix_update_state(mix_state_t* mxs, float **probs, uint8_t sym, float learning_rate);
void mix_free(mix_state_t* mxs);

static inline float stretch(float p) {
  return fasterlog2(p / (1 - p));
}
