/*
 * dct_fftpack.c
 *
 * Created by Leland Brown on 2011 Feb 23.
 *
 * Copyright (c) 2011-2013 Leland Brown.
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
 */

//
// This file provides a particular implementation of the functions in dct.h,
// in this case using the DCT functions declared in myfftpack.h. Similar
// implementations are possible using other DCT or FFT libraries, if desired.
// These implementations should be functionally interchangeable as long as
// they conform to the interface in dct.h.
//

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include "dct.h"

#include "fftpack.h"

#include <stdlib.h>
#include <assert.h>

static const int max_factors = 30;

typedef struct {
    int     dct_type;   // 1 to 3 (DCT types I to III)
    int     nelems;     // length of inout_data0 and inout_data1 buffers
    double *inout_data0;// input/output buffer space
    double *inout_data1;// input/output buffer space
    double *wsave;      // workspace buffer
    int    *ifac;       // info on factorization of nelems
} Dct_Buffer;

int setup_dcts(
    // returns 0 on success, 1 if a memory allocation error occurred
    int       dct_type,     // 1, 2, or 3 (DCT types I, II, III)
    int       nelems,       // data length for each DCT
    void   *(*dct_buffer),  // internal buffer for use by perform_dcts()
    double *(in_data[2]),   // input  buffer locations for perform_dcts()
    double *(out_data[2])   // output buffer locations for perform_dcts()
    // on input, *dct_buffer, in_data[.], and out_data[.] must be NULL;
    // caller should NOT free the returned pointers - use cleanup_dcts() instead
)
// Specifies a DCT operation to be performed one or more times and
// allocates buffers to be used by perform_dcts().
// Note: returned in_data[.] and out_data[.] pointers may be the same.
{
    const int max_ifac = (int)( 1.8 * max_factors + 6.9 );
    
    Dct_Buffer *buf;
    double *data;

    assert( *dct_buffer == NULL );
    assert( in_data[0]  == NULL );
    assert( in_data[1]  == NULL );
    assert( out_data[0] == NULL );
    assert( out_data[1] == NULL );
    
    buf = (Dct_Buffer *)malloc( sizeof( Dct_Buffer ) );
    if (!buf) {
        return 1;
    }

    data = (double *)malloc
        ( 30 * nelems * sizeof( double ) + max_ifac * sizeof( int ) );
    if (!data) {
        free( buf );
        return 1;
    }
    

    buf->dct_type = dct_type;
    buf->nelems   = nelems;
    
    buf->inout_data0 =  data;
    buf->inout_data1 =  data + nelems;
    buf->wsave =        data + nelems * 2;
    buf->ifac = (int *)(data + nelems * 30);
    
    switch (dct_type) {
    //  case 1:
    //  //  costi( nelems, buf->wsave, buf->ifac );
    //      break;
        case 2: case 3:
            cosqi( nelems, buf->wsave, buf->ifac );
            break;
        default:
            assert( 0 );    // illegal or unsupported dct_type
    }
    
    assert( buf->ifac[1] <= max_factors );

    *dct_buffer = (void *)buf;
    in_data[0]  = buf->inout_data0;
    in_data[1]  = buf->inout_data1;
    out_data[0] = buf->inout_data0; // in-place transforms
    out_data[1] = buf->inout_data1; // in-place transforms
    return 0;
}

void perform_dcts(
    void   *dct_buffer,         // internal buffer from setup_dcts()
    double *const (in_data[2]), // input  data buffers (may be overwritten)
    double *const (out_data[2]) // output data buffers
)
// Performs two DCTs, each of size nelems (see setup_dcts()).
// All pointers must match those provided by setup_dcts().
// Note: input buffers may be overwritten, even if output buffers are different.
// Values in both data buffers must have similar magnitude to avoid roundoff error;
// for a single DCT, fill both buffers with the same values, or one with zeroes.
{
    Dct_Buffer *buf = (Dct_Buffer *)(dct_buffer);

    assert( in_data[0]  == buf->inout_data0 );
    assert( in_data[1]  == buf->inout_data1 );
    assert( out_data[0] == buf->inout_data0 );  // in-place transform
    assert( out_data[1] == buf->inout_data1 );  // in-place transform
    
    switch (buf->dct_type) {
    //  case 1:
    //  //  cost(  buf->nelems, buf->inout_data0, buf->wsave, buf->ifac );
    //  //  cost(  buf->nelems, buf->inout_data1, buf->wsave, buf->ifac );
    //      break;
        case 2:
            cosqb2(
                buf->nelems, buf->inout_data0, buf->inout_data1,
                buf->wsave, buf->ifac );
            break;
        case 3:
            cosqf2(
                buf->nelems, buf->inout_data0, buf->inout_data1,
                buf->wsave, buf->ifac );
            break;
        default:
            assert( 0 );    // illegal or unsupported dct_type
    }
}

void cleanup_dcts(
    void   *(*dct_buffer),  // internal buffer from setup_dcts()
    double *(in_data[2]),   // input  buffer locations from setup_dcts()
    double *(out_data[2])   // output buffer locations from setup_dcts()
    // on output, *dct_buffer, in_data[.], and out_data[.] will be NULL
)
// Frees memory allocated by setup_dcts().
{
    Dct_Buffer *buf = (Dct_Buffer *)(*dct_buffer);
    
    assert( in_data[0]  == buf->inout_data0 );
    assert( in_data[1]  == buf->inout_data1 );
    assert( out_data[0] == buf->inout_data0 );
    assert( out_data[1] == buf->inout_data1 );
    
    free( buf->inout_data0 );
    free( buf );
    
    in_data[0]  = NULL;
    in_data[1]  = NULL;
    out_data[1] = NULL;
    out_data[0] = NULL;
    *dct_buffer = NULL;
}
