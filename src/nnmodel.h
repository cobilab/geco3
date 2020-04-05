#include "defs.h"
#include "nn.h"

typedef struct nnmodel_t {
  uint32_t context;
  float ema_err_short;
  float ema_err_long;
  ann_t *ann;
} nnmodel_t;

nnmodel_t* nnmodel_init(uint32_t nsyms, uint32_t context, uint32_t hs);
float const* nnmodel_predict(nnmodel_t *nnm);
void nnmodel_update(nnmodel_t *nnm, uint32_t sym, float lr);
