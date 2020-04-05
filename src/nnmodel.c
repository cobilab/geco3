#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "nnmodel.h"

nnmodel_t* nnmodel_init(uint32_t nsyms, uint32_t context, uint32_t hs) {
  nnmodel_t *nnm = malloc(sizeof(nnmodel_t));
  nnm->context = context;
  nnm->ema_err_short = 1;
  nnm->ema_err_long = 0;
  uint32_t xs = (nsyms*context) + 1;
  nnm->ann = ann_init(xs, hs, nsyms);
  nnm->ann->x[xs - 1] = 1.0;
  return nnm;
}

float const* nnmodel_predict(nnmodel_t *nnm) {
  float *ret = nnm->ann->y;
  ann_apply(nnm->ann);
  return ret;
}

void nnmodel_update(nnmodel_t *nnm, uint32_t sym, float lr) {
  ann_t *ann = nnm->ann;
  const uint32_t nsyms = ann->ys;
  float tdata[nsyms];

  int i;
  for(i = 0 ; i < nsyms; i++) {
    tdata[i] = 0.0;
  }
  tdata[sym] = 1.0;

  float err = 0;
  for(i = 0 ; i < nsyms; i++) {
    err += fabsf(tdata[i] - ann->y[i]);
  }

  const float gamma_short = 0.3;
  const float gamma_long = 0.05;

  nnm->ema_err_short = (gamma_short * err) + ((1 - gamma_short) * nnm->ema_err_short);
  nnm->ema_err_long = (gamma_long * err) + ((1 - gamma_long) * nnm->ema_err_long);


  //if(fabsf(nnm->ema_err_short - nnm->ema_err_long) > 0.0001) {
    ann_train(ann, tdata, lr);
    //}

  for(i = ann->xs - 2 ; i >= nsyms; i--) {
    ann->x[i] = ann->x[i - nsyms];
  }

  for(i = 0 ; i < nsyms; i++) {
    ann->x[i] = -1.0;
  }
  ann->x[sym] = 1.0;

}
