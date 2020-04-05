/******************************************************************************
File:      arith.c

Authors:   John Carpinelli   (johnfc@ecr.mu.oz.au)
           Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)
           Lang Stuiver      (langs@cs.mu.oz.au)
           Radford Neal      (radford@ai.toronto.edu)

Purpose:   Data compression using a revised arithmetic coding method.

Based on:  A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisted",
           Proc. IEEE Data Compression Conference, Snowbird, Utah, March 1995.

           Low-Precision Arithmetic Coding Implementation by Radford M. Neal

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

  $Log: arith.c,v $
  Revision 1.2  2008-01-23 15:59:09  ap
  No real changes...

  Revision 1.1  2005/07/01 09:22:27  ap

  New files.

  Revision 1.1  1996/08/07 01:34:11  langs
  Initial revision

 ************************************************************************

Code adapted by Armando J. Pinho
ap@ua.pt
University of Aveiro, DETI/IEETA, 3810-193 Aveiro, Portugal
December 1999

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "arith.h"
#include "bitio.h"

#define		Half		((code_value) 1 << (B_bits-1))
#define		Quarter		((code_value) 1 << (B_bits-2))

/* Separate input and output */
/* Input decoding state */
static code_value	in_R;				/* code range */
static code_value	in_D;				/* = V-L (V offset)*/
static div_value	in_r;				/* normalized range */

/* Output encoding state */
static code_value	out_L;				/* lower bound */
static code_value	out_R;				/* code range */
static unsigned long	out_bits_outstanding;		/* follow bit count */


/*
 * BIT_PLUS_FOLLOW(b, s)
 * responsible for outputting the bit passed to it and an opposite number of
 * bit equal to the value stored in bits_outstanding
 *
 */
#define ORIG_BIT_PLUS_FOLLOW(b, s)	\
do                                      \
{ 	  			        \
    OUTPUT_BIT((b), s);           	\
    while (out_bits_outstanding > 0)	\
    { 					\
	OUTPUT_BIT(!(b), s);      	\
	out_bits_outstanding--;    	\
    } 	                		\
} while (0)


#  define BIT_PLUS_FOLLOW(x, s)	ORIG_BIT_PLUS_FOLLOW(x, s)

/*
 * ENCODE_RENORMALISE
 * output code bits until the range has been expanded
 * to above QUARTER
 * With FRUGAL_BITS option, ignore first zero bit output
 * (a redundant zero will otherwise be emitted every time the encoder is
 * started)
 */
#define ENCODE_RENORMALISE(s)		\
do {					\
    while (out_R <= Quarter)		\
    {					\
        if (out_L >= Half)		\
    	{				\
    	    BIT_PLUS_FOLLOW(1, s);	\
    	    out_L -= Half;		\
    	}				\
    	else if (out_L+out_R <= Half)	\
    	{				\
    	    BIT_PLUS_FOLLOW(0, s);	\
    	}				\
    	else 				\
    	{				\
    	    out_bits_outstanding++;	\
    	    out_L -= Quarter;		\
    	}				\
    	out_L <<= 1;			\
    	out_R <<= 1;			\
    }					\
} while (0)


/*
 * DECODE_RENORMALISE
 * input code bits until range has been expanded to
 * more than QUARTER. Mimics encoder.
 * FRUGAL_BITS option also keeps track of bitstream input so it can work out
 * exactly how many disambiguating bits the encoder put out (1,2 or 3).
 */
#define DECODE_RENORMALISE(s)			\
do {						\
    while (in_R <= Quarter)			\
    {						\
    	in_R <<= 1;				\
    	ADD_NEXT_INPUT_BIT(in_D, 0, s);		\
    }						\
} while (0)


/*
 *
 * encode a symbol given its low, high and total frequencies
 *
 */
void arithmetic_encode(freq_value low, freq_value high, freq_value total,
  FILE *s)
{ 
  /* The following pseudocode is a concise (but slow due to arithmetic
   * calculations) description of what is calculated in this function.
   * Note that the division is done before the multiplication.  This is
   * one of the key points of this version of arithmetic coding.  This means
   * that much larger frequency values can be accomodated, but that the integer
   * ratio R/total (code range / frequency range) 
   * becomes a worse approximation as the total frequency count
   * is increased.  Therefore less than the whole code range will be
   * allocated to the frequency range.  The whole code range is used by
   * assigning the symbol at the end of the frequency range (where high==total)
   * this excess code range above and beyond its normal code range.  This
   * of course distorts the probabilities slightly, 
   * see the paper connected with this program for details and compression
   * results.
   *
   * Restricting the total frequency range means that the division and
   * multiplications can be done to low precision with shifts and adds.
   *
   * The following code describes this function:
   * (The actual code is only written less legibly to improve speed)
   *
   *    L += low*(R/total);			* Adjust low bound	*
   *
   *    if (high < total)			* 			*
   *    	R = (high-low) * (R/total);	* Restrict range	*
   *    else					* If symbol at end of	*
   *		R = R   -  low * (R/total);	*   range.  Restrict &	*
   *						*   allocate excess 	*
   * 						*   codelength to it.	*
   *    ENCODE_RENORMALISE;			* Expand code range	*
   *						*   and output bits	*
   *
   *    if (bits_outstanding > MAX_BITS_OUTSTANDING)
   *		abort();			* EXTREMELY improbable	*
   *						* (see comments below) *
   */

    code_value temp; 

   {
    div_value out_r;			
    out_r = out_R/total;		/* Calc range:freq ratio */
    temp = out_r*low;			/* Calc low increment */
    out_L += temp;			/* Increase L */
    if (high < total)
	out_R = out_r*(high-low);	/* Restrict R */
    else
	out_R -= temp;			/* If at end of freq range */
					/* Give symbol excess code range */
   }

    ENCODE_RENORMALISE(s);

    if (out_bits_outstanding > MAX_BITS_OUTSTANDING)
    {
/* 
 * For MAX_BITS_OUTSTANDING to be exceeded is extremely improbable, but
 * it is possible.  For this to occur the COMPRESSED file would need to
 * contain a sequence MAX_BITS_OUTSTANDING bits long (eg: 2^31 bits, or
 * 256 megabytes) of all 1 bits or all 0 bits.  This possibility is
 * extremely remote (especially with an adaptive model), and can only
 * occur if the compressed file is greater than MAX_BITS_OUTSTANDING
 * long.  Assuming the compressed file is effectively random,
 * the probability for any 256 megabyte section causing an overflow
 * would be 1 in 2^(2^31).  This is a number approximately 600 million
 * digits long (decimal).
 * 
 */
	
	fprintf(stderr,"Bits_outstanding limit reached - File too large\n");
	exit(1);
    }
}



/*
 * arithmetic decode_target()
 * decode the target value using the current total frequency
 * and the coder's state variables
 *
 * The following code describes this function:
 * The actual code is only written less legibly to improve speed,
 * including storing the ratio in_r = R/total for use by arithmetic_decode()
 *
 *  target = D / (R/total);	* D = V-L.  (Old terminology) 		*
 *				* D is the location within the range R	*
 *				* that the code value is located	*
 *				* Dividing by (R/total) translates it	*
 *				* to its correspoding frequency value	*
 *
 *				
 *  if (target < total) 	* The approximate calculations mean the *
 *	return target;		* encoder may have coded outside the	*
 *  else			* valid frequency range (in the excess	*
 *	return total-1;		* code range).  This indicates that the	*
 *				* symbol at the end of the frequency	*
 *				* range was coded.  Hence return end of *
 *				* range (total-1)			*
 *
 */
freq_value arithmetic_decode_target(freq_value total)
{
    freq_value target;
    
    in_r = in_R/total;
    target = (in_D)/in_r;

    return (target >= total ? total-1 : target);
}


/*
 * arithmetic_decode()
 *
 * decode the next input symbol
 *
 * The following code describes this function:
 * The actual code is only written less legibly to improve speed,
 * See the comments for arithmetic_encode() which is essentially the
 * same process.
 *
 * D -= low * (R/total);		* Adjust current code value 	   *
 *
 * if (high < total)
 *	R = (high-low) * (R/total);	* Adjust range 			   *
 * else
 *	R -= low * (R/total);		* End of range is a special case   *
 *
 * DECODE_RENORMALISE;			* Expand code range and input bits *
 *
 */

void arithmetic_decode(freq_value low, freq_value high, freq_value total,
  FILE *s)
{     
    code_value temp;

    /* assume r has been set by decode_target */
    temp = in_r*low;
    in_D -= temp;
    if (high < total)
	in_R = in_r*(high-low);
    else
	in_R -= temp;

    DECODE_RENORMALISE(s);

}


/*
 *
 * start the encoder
 *
 */
void start_encode(void)
{
    out_L = 0;				/* Set initial coding range to	*/
    out_R = Half;			/* [0,Half)			*/
    out_bits_outstanding = 0;
}


/*
 * finish_encode()
 * finish encoding by outputting follow bits and all the bits in L
 * (B_bits long) to make the last symbol unambiguous.  This also means the
 * decoder can read them harmlessly as it reads beyond end of what would
 * have been valid input.
 *
 */
void finish_encode(FILE *s)
{
  int nbits, i;
  code_value bits;

  nbits = B_bits;
  bits  = out_L;
  for (i = 1; i <= nbits; i++)        /* output the nbits integer bits */
        BIT_PLUS_FOLLOW(((bits >> (nbits-i)) & 1), s);
}

/* 
 * start_decode() start the decoder Fills the decode value in_D from the
 * bitstream.  If FRUGAL_BITS is defined only the first call reads in_D
 * from the bitstream.  Subsequent calls will assume the excess bits
 * that had been read but not used (sitting in in_V) are the start of
 * the next coding message, and it will put these into in_D.  (It also
 * initialises in_V to the start of the bitstream)
 * 
 * FRUGAL_BITS also means that the first bit (0) was not output, so only
 * take B_bits-1 from the input stream.  Since there are B_bits
 * "read-ahead" in the buffer, on second and subsequent calls
 * retrieve_excess_input_bits() is used and the last bit must be placed back
 * into the input stream.
 */
void start_decode(FILE *s)
{
 int i;
  in_D = 0;			/* Initial offset in range is 0 */
  in_R = Half;			/* Range = Half */

  for (i = 0; i<B_bits; i++)			/* Fill D */
	ADD_NEXT_INPUT_BIT(in_D, 0, s);

  if (in_D >= Half)
	{
	  fprintf(stderr,"Corrupt input file (start_decode())\n");
	  exit(1);
	}
}

/* 
 * finish_decode()
 *
 * Throw away B_bits in buffer by doing nothing
 * (encoder wrote these for us to consume.)
 * (They were mangled anyway as we only kept V-L, and cannot get back to V)
 */
void finish_decode(void)
{
	/* No action */
}

