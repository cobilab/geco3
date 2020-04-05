/******************************************************************************
File:      arith.h

Authors:   John Carpinelli   (johnfc@ecr.mu.oz.au)
           Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)
           Lang Stuiver      (langs@cs.mu.oz.au)

Purpose:   Data compression using revised arithmetic coding method.

Based on:  A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisted",
           Proc. IEEE Data Compression Conference, Snowbird, Utah, March 1995.

Copyright 1995 John Carpinelli and Wayne Salamonsen, All Rights Reserved.
Copyright 1996 Lang Stuiver, All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not sold or incorporated into any commercial product.  There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER.  Use them at your own risk.  If you do
happen to find a bug, or have modifications to suggest, please report
the same to Alistair Moffat, alistair@cs.mu.oz.au.  The copyright
notice above and this statement of conditions must remain an integral
part of each and every copy made of these files.

 ************************************************************************

Code adapted by Armando J. Pinho
ap@ua.pt
University of Aveiro, DETI/IEETA, 3810-193 Aveiro, Portugal
December 1999

******************************************************************************/

#ifndef CODER_H
#define CODER_H

/* ================= USER ADJUSTABLE PARAMETERS =================== */

	/* Default B_bits and F_bits */

#ifndef B_BITS
#define		B_BITS		32
#endif

#ifndef F_BITS
#define		F_BITS		27
#endif

/* Change these types for different precision calculations.  They may affect
 * the speed of the arithmetic operations (multiplcation, division, shift,
 * etc).
 * The way the stats module is implemented, the type of freq_value
 * must be able to accomodate f_bits+1 bits, instead of f_bits, to avoid
 * overflows.  Ie: For an f_bits of up to 31, type freq_value must be 32 bits.
 */
typedef unsigned long   code_value;	/* B_BITS of precision */
typedef unsigned long	freq_value;	/* F_BITS+1 of precision */
typedef unsigned long	div_value;	/* B_BITS-F_BITS of precision */


/* MAX_BITS_OUTSTANDING is a bound on bits_outstanding
 * If bits_outstanding ever gets above this number (extremely unlikely)
 * the program will abort with an error message.  (See arith.c for details).
 */
#define 	MAX_BITS_OUTSTANDING	((unsigned long)1<<31)


/* ================= END USER ADJUSTABLE PARAMETERS =================== */


/* Determine maximum bits allowed, based on size of the types used
 * to store them.  Also, that MAX_F_BITS <= MAX_B_BITS-2
 */

#define		MAX_B_BITS   (int)( sizeof(code_value) * 8)
#define		MAX_F_BITS   (int)((sizeof(freq_value)*8)-1 < MAX_B_BITS - 2\
				?  (sizeof(freq_value)*8)-1 : MAX_B_BITS - 2)

/* If varying bits, variables are B_bits, F_bits, Half and Quarter,
 *	otherwise #define them
 * These variables will be read (and possibly changed) by main.c and
 *  stats.c
 */

#define		B_bits		B_BITS
#define		F_bits  	F_BITS


extern char *coder_desc;


/* function prototypes */
void arithmetic_encode(freq_value l, freq_value h, freq_value t, FILE *s);
freq_value arithmetic_decode_target(freq_value t);
void arithmetic_decode(freq_value l, freq_value h, freq_value t, FILE *s);
void binary_arithmetic_encode(freq_value c0, freq_value c1, int bit, FILE *s);
int binary_arithmetic_decode(freq_value c0, freq_value c1, FILE *s);
void start_encode(void);
void finish_encode(FILE *s);
void start_decode(FILE *s);
void finish_decode(void);

#endif		/* ifndef arith.h */

