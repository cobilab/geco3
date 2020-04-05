/******************************************************************************
File:      bitio.c

Authors:   John Carpinelli   (johnfc@ecr.mu.oz.au)          
           Wayne Salamonsen  (wbs@mundil.cs.mu.oz.au)       
           Lang Stuiver      (langs@cs.mu.oz.au)

Purpose:   Data compression using a revised arithmetic coding method.
Based on:  A. Moffat, R. Neal, I.H. Witten, "Arithmetic Coding Revisted",
           Proc. IEEE Data Compression Conference, Snowbird, Utah, March 1995.

           Low-Precision Arithmetic Coding Implementation by Radford M. Neal

Copyright 1996 Lang Stuiver, All Rights Reserved.

These programs are supplied free of charge for research purposes only,
and may not sold or incorporated into any commercial product.  There is
ABSOLUTELY NO WARRANTY of any sort, nor any undertaking that they are
fit for ANY PURPOSE WHATSOEVER.  Use them at your own risk.  If you do
happen to find a bug, or have modifications to suggest, please report
the same to Alistair Moffat, alistair@cs.mu.oz.au.  The copyright
notice above and this statement of conditions must remain an integral
part of each and every copy made of these files.

******************************************************************************

Code adapted by Armando J. Pinho
ap@ua.pt
University of Aveiro, DETI/IEETA, 3810-193 Aveiro, Portugal
December 1999

******************************************************************************/

#include <stdio.h>
#include "bitio.h"

/* 
 * The following variables are supposedly local, but actually global so they
 * can be referenced by macro
 */

unsigned int	_bytes_input = 0;
unsigned int	_bytes_output = 0;

int		_in_buffer;			/* I/O buffer */
unsigned char	_in_bit_ptr = 0;		/* bits left in buffer */
int		_in_garbage;			/* bytes read beyond eof */

int		_out_buffer;			/* I/O buffer */
int		_out_bits_to_go;		/* bits to fill buffer */

#ifndef FAST_BITIO
int		_bitio_tmp;			/* Used by some of the */
#endif						/* bitio.h macros */

/*
 *
 * initialize the bit output function
 *
 */
void startoutputtingbits(void)
{
    _out_buffer = 0;
    _out_bits_to_go = BYTE_SIZE;
}

/*
 *
 * start the bit input function
 *
 */
void startinputtingbits(void)
{
    _in_garbage = 0;	/* Number of bytes read past end of file */
    _in_bit_ptr = 0;	/* No valid bits yet in input buffer */
}

/*
 *
 * complete outputting bits
 *
 */
void doneoutputtingbits(FILE *s)
{
    if (_out_bits_to_go != BYTE_SIZE)
	OUTPUT_BYTE(_out_buffer << _out_bits_to_go, s);
    _out_bits_to_go = BYTE_SIZE;
}

/*
 *
 * complete inputting bits
 *
 */
void doneinputtingbits(void)
{
      _in_bit_ptr = 0;	      /* "Wipe" buffer (in case more input follows) */
}

