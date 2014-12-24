/*
 * dct_fftpack.c
 *
 * Created by Leland Brown on 2011 Feb 23.
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
    int     nelems;     // length of in_data and aux_data arrays
    double *in_data;    // input array space
    double *aux_data;   // alternate output array space
    double *wsave;      // workspace array
    int    *ifac;       // info on factorization of nelems
} Dct_Info;

int setup_dct(
    // returns 0 on success, 1 if a memory allocation error occurred
    int       dct_type, // 1 to 3 (DCT types I to III)
    int       nelems,   // length for in_data array
    void   *(*dct_info),// internal data for DCT algorithm
    double *(*in_data)  // input array space for perform_dct()
    // on input, *dct_info and *in_data must be NULL
)
{
    Dct_Info *info;

    assert( *dct_info == NULL );
    assert( *in_data  == NULL );
    
    *dct_info = malloc( sizeof( Dct_Info ) );
    if (!*dct_info) {
        return 1;
    }

    *in_data = (double *)malloc
        ( 8 * nelems * sizeof( double ) + (max_factors+2) * sizeof( int ) );
    if (!*in_data) {
        free( *dct_info );
        *dct_info = NULL;
        return 1;
    }

    info = (Dct_Info *)(*dct_info);
    info->dct_type = dct_type;
    info->nelems   = nelems;
    info->in_data =     (*in_data);
    info->wsave =       (*in_data + 2 * nelems);
    info->ifac = (int *)(*in_data + 8 * nelems);
    switch (dct_type) {
        case 1:
            info->aux_data = (*in_data + 3 * nelems);
        //  costi( nelems, info->wsave, info->ifac );
            break;
        case 2: case 3:
            info->aux_data = (*in_data + 3 * nelems);
            cosqi( nelems, info->wsave, info->ifac );
            break;
        default:
            assert( 0 );    // illegal or unsupported dct_type
    }
    assert( info->ifac[1] <= max_factors );

    return 0;
}

double *perform_dct(
    // returns pointer to output array (may be same as in_data);
    // caller should NOT free returned pointer
    void   *dct_info,   // internal data for DCT algorithm
    double *in_data     // input array (may be overwritten)
)
{
    Dct_Info *info = (Dct_Info *)(dct_info);
    assert( in_data == info->in_data );
    switch (info->dct_type) {
        case 1:
        //  cost( info->nelems, info->in_data, info->wsave, info->ifac );
            break;
        case 2:
            cosqb( info->nelems, info->in_data, info->wsave, info->ifac );
            break;
        case 3:
            cosqf( info->nelems, info->in_data, info->wsave, info->ifac );
            break;
        default:
            assert( 0 );    // illegal or unsupported dct_type
    }
    return in_data;
}

void cleanup_dct(
    // frees memory allocated by setup_dct()
    void   *(*dct_info),// internal data for DCT algorithm
    double *(*in_data)  // input array space
)
{
    Dct_Info *info = (Dct_Info *)(*dct_info);
    assert( *in_data == info->in_data );
    free( *in_data );
    free( *dct_info );
    *in_data  = NULL;
    *dct_info = NULL;
}
