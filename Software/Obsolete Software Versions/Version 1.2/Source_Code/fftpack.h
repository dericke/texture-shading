/********************************************************************
 *
 * File: fftpack.h
 * Function: Fast discrete Fourier and cosine transforms and inverses
 *
 * Original author: Paul N. Swarztrauber
 * Last modification date: 1985 Apr (public domain)
 *
 * Modifications by: Monty <xiphmont@mit.edu>
 * Last modification date: 1996 Jul 01 (public domain)
 *
 * Modifications by: Leland Brown
 * Last modification date: 2011 Oct 08
 *
 * Copyright (c) 2011 Leland Brown.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************/

/*
 * These Fourier routines were originally based on the Fourier
 * routines of the same names from the NETLIB bihar and fftpack
 * Fortran libraries developed by Paul N. Swarztrauber at the
 * National Center for Atmospheric Research in Boulder, CO, USA.
 * They have been reimplemented in C and optimized in a few ways.
 */

#ifndef FFTPACK_H
#define FFTPACK_H

#include "compatibility.h"

#ifdef __cplusplus
extern "C" {
#endif

//#define FFTPACK_REAL float
#define FFTPACK_REAL double

//*******************************************************************************
//
//  costi initializes wsave and ifac, used in dcost.
//
//  Description:
//
//    The prime factorization of n together with a tabulation of the
//    trigonometric functions are computed and stored in ifac and wsave.
//
//  Parameters:
//
//    Input, int n, the length of the sequence to be transformed.  The
//    method is more efficient when n-1 is the product of small primes.
//
//    Output, REAL wsave[3*n], contains data, depending on n, and
//    required by the dcost algorithm.
//
//    Output, int ifac[].
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
//void costi(int n, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

//*******************************************************************************
//
//  cost computes the discrete Fourier cosine transform of an even sequence.
//
//  Description:
//
//    This routine is the unnormalized inverse of itself.  Two successive
//    calls will multiply the input sequence x by 2*(n-1).
//
//    The arrays wsave and ifac must be initialized by calling dcosti.
//
//    The transform is defined by:
//
//      x_out[i] = x_in[0] + (-1)^i * x_in[n-1] + sum ( 1 <= k <= n-2 )
//
//        2 * x_in[k] * cos ( k * i * pi / ( n - 1 ) )
//
//  Parameters:
//
//    Input, int n, the length of the sequence to be transformed.  The
//    method is more efficient when n-1 is the product of small primes.
//
//    Input/output, REAL x[n].
//    On input, the sequence to be transformed.
//    On output, the transformed sequence.
//
//    Input, REAL wsave[3*n].
//    The wsave array must be initialized by calling dcosti.  A different
//    array must be used for each different value of n.
//
//    Input, int ifac[].  The ifac array must be initialized by calling dcosti.
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
//void cost(int n, FFTPACK_REAL *RESTRICT x, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);


//*******************************************************************************
//
//  cosqi initializes wsave and ifac, used in dcosqf and dcosqb.
//
//  Description:
//
//    The prime factorization of n together with a tabulation of the
//    trigonometric functions are computed and stored in ifac and wsave.
//
//  Parameters:
//
//    Input, int n, the length of the array to be transformed.  The method
//    is more efficient when n is the product of small primes.
//
//    Output, REAL wsave[3*n], contains data, depending on n, and
//    required by the dcosqb and dcosqf algorithms.
//
//    Output, int ifac[].
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void cosqi(int n, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

//*******************************************************************************
//
//  cosqf computes the fast cosine transform of quarter wave data.
//
//  Description:
//
//    dcosqf computes the coefficients in a cosine series representation
//    with only odd wave numbers.
//
//    dcosqf is the unnormalized inverse of dcosqb since a call of dcosqf
//    followed by a call of dcosqb will multiply the input sequence x
//    by 4*n.
//
//    The arrays wsave and ifac must be initialized by calling dcosqi.
//
//    The transform is defined by:
//
//      x_out[i] = x_in[0] + sum ( 1 <= k <= n-1 )
//
//        2 * x_in[k] * cos ( ( 2 * i + 1 ) * k * pi / ( 2 * n ) )
//
//  Parameters:
//
//    Input, int n, the length of the array x.  The method is
//    more efficient when n is the product of small primes.
//
//    Input/output, REAL x[n].
//    On input, the data to be transformed.
//    On output, the transformed data.
//
//    Input, REAL wsave[3*n], contains data, depending on n, and
//    required by the algorithm.  The wsave array must be initialized by
//    calling dcosqi.  A different wsave array must be used for each different
//    value of n.
//
//    Input, int ifac[].  The ifac array must be initialized by calling dcosqi.
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void cosqf(int n, FFTPACK_REAL *RESTRICT x, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

//*******************************************************************************
//
//  cosqb computes the fast cosine transform of quarter wave data.
//
//  Description:
//
//    dcosqb computes a sequence from its representation in terms of a cosine
//    series with odd wave numbers.
//
//    The transform is defined by:
//
//      x_out[i] = sum ( 0 <= k <= n-1 ) 
//
//        4 * x_in[k] * cos ( ( 2 * k + 1 ) * i * pi / ( 2 * n ) )
//
//    dcosqb is the unnormalized inverse of dcosqf since a call of dcosqb
//    followed by a call of dcosqf will multiply the input sequence x by 4*n.
//
//    The arrays wsave and ifac must be initialized by calling dcosqi.
//
//  Parameters:
//
//    Input, int n, the length of the array x.  The method is
//    more efficient when n is the product of small primes.
//
//    Input/output, REAL x[n].
//    On input, the cosine series coefficients.
//    On output, the corresponding data vector.
//
//    Input, REAL wsave[3*n], contains data, depending on n, and
//    required by the algorithm.  The wsave array must be initialized by
//    calling dcosqi.  A different wsave array must be used for each different
//    value of n.
//
//    Input, int ifac[].  The ifac array must be initialized by calling dcosqi.
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void cosqb(int n, FFTPACK_REAL *RESTRICT x, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);


//*******************************************************************************
//
//  rffti initializes wsave and ifac, used in drfftf and drfftb.
//
//  Description:
//
//    The prime factorization of n together with a tabulation of the
//    trigonometric functions are computed and stored in ifac and wsave.
//
//  Parameters:
//
//    Input, int n, the length of the sequence to be transformed.
//
//    Output, REAL wsave[2*n], contains data, dependent on the value
//    of n, which is necessary for the drfftf and drfftb routines.
//
//    Output, int ifac[].
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void rffti(int n, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

//*******************************************************************************
//
//  rfftf computes the Fourier coefficients of a real periodic sequence.
//
//  Description:
//
//    This process is sometimes called Fourier analysis.
//
//    The transform is unnormalized.  A call to drfftf followed by a call
//    to drfftb will multiply the input sequence by n.
//
//    The transform is defined by:
//
//      r_out[0] = sum ( 0 <= i <= n-1 ) r_in[i]
//
//    For k = 1,...,(n-1)/2
//
//      r_out[2*k-1] = sum ( 0 <= i <= n-1 )
//
//        r_in[i] * cos ( k * i * 2 * pi / n )
//
//      r_out[2*k] = sum ( 0 <= i <= n-1 )
//
//        -r_in[i] * sin ( k * i * 2 * pi / n )
//
//    And, if n is even, then:
//
//      r_out[n-1] = sum ( 0 <= i <= n-1 ) (-1)^i * r_in[i]
//
//  Parameters:
//
//    Input, int n, the length of the array to be transformed.  The
//    method is more efficient when n is the product of small primes.
//
//    Input/output, REAL r[n].
//    On input, the sequence to be transformed.
//    On output, the transformed sequence.
//
//    Input, REAL wsave[2*n], a work array.  The wsave array
//    must be initialized by calling drffti.  A different wsave array must be
//    used for each different value of n.
//
//    Input, int ifac[].  The ifac array must be initialized by calling drffti.
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void rfftf(int n, FFTPACK_REAL *RESTRICT r, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

//*******************************************************************************
//
//  rfftb computes a real periodic sequence from its Fourier coefficients.
//
//  Description:
//
//    This process is sometimes called Fourier synthesis.
//
//    The transform is unnormalized.  A call to drfftf followed by a call to
//    drfftb will multiply the input sequence by n.
//
//    If n is even, the transform is defined by:
//
//      r_out[i] = r_in[0] + (-1)^i * r_in[n-1] + sum ( 1 <= k <= n/2-1 )
//
//        + 2 * r_in[2*k-1] * cos ( k * i * 2 * pi / n )
//
//        - 2 * r_in[2*k]   * sin ( k * i * 2 * pi / n )
//
//    If n is odd, the transform is defined by:
//
//      r_out[i] = r_in[0] + sum ( 1 <= k <= (n-1)/2 )
//
//        + 2 * r_in[2*k-1] * cos ( k * i * 2 * pi / n )
//
//        - 2 * r_in[2*k]   * sin ( k * i * 2 * pi / n )
//
//  Parameters:
//
//    Input, int n, the length of the array to be transformed.  The
//    method is more efficient when n is the product of small primes.
//
//    Input/output, REAL r[n].
//    On input, the sequence to be transformed.
//    On output, the transformed sequence.
//
//    Input, REAL wsave[2*n], a work array.  The wsave array must be
//    initialized by calling drffti.  A different wsave array must be used
//    for each different value of n.
//
//    Input, int ifac[].  The ifac array must be initialized by calling drffti.
//    ifac[0] = n, the number that was factored.
//    ifac[1] = nf, the number of factors.
//    ifac[2..1+nf], the factors.
//
//*******************************************************************************
void rfftb(int n, FFTPACK_REAL *RESTRICT r, FFTPACK_REAL *RESTRICT wsave, int *RESTRICT ifac);

#ifdef __cplusplus
}
#endif

#endif
