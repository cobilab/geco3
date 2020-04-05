/*------------------------------------------------------------------------------

Copyright 2010 IEETA / University of Aveiro, All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not be sold or incorporated into any commercial product. There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER. Use them at your own risk. If you do
happen to find a bug, or have modifications to suggest, please report
the same to Armando J. Pinho, ap@ua.pt. The copyright notice above
and this statement of conditions must remain an integral part of each
and every copy made of these files.

------------------------------------------------------------------------------*/

#ifndef ARITH_AUX_H_INCLUDED
#define ARITH_AUX_H_INCLUDED

#include "defs.h"

void GetInterval(int *low, int *high, int *count, int symbol);
int GetSymbol(int *low, int *high, int *count, int target, int nSymbols);
void WriteNBits(uint64_t bits, int nBits, FILE *oFp);
uint64_t ReadNBits(int nBits, FILE *iFp);
void AESym(int symbol, int *counters, int totalCount, FILE *oFp);
int ArithDecodeSymbol(int nSymbols, int *counters, int totalCount, FILE *iFp);

#endif /* ARITH_AUX_H_INCLUDED */

