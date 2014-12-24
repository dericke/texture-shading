/*
 * dct.h
 *
 * Created by Leland Brown on 2010 Oct 31.
 *
 * Copyright (c) 2010-2011 Leland Brown.
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
// This file is intended to be a generic interface for DCT functions.
// The functions prototyped here may be implemented with appropriate
// calls to any DCT or FFT library. These implementations should be
// functionally interchangeable as long as they conform to this interface.
//

#ifndef DCT_H
#define DCT_H

#ifdef __cplusplus
extern "C" {
#endif

int setup_dct(
    // returns 0 on success, 1 if a memory allocation error occurred
    int       dct_type, // 1 to 3 (DCT types I to III)
    int       nelems,   // length for in_data array
    void   *(*dct_info),// internal data for DCT algorithm
    double *(*in_data)  // input array space for perform_dct()
    // on input, *dct_info and *in_data must be NULL
);

double *perform_dct(
    // returns pointer to output array (may be same as in_data);
    // caller should NOT free returned pointer
    void   *dct_info,   // internal data for DCT algorithm
    double *in_data     // input array (may be overwritten)
);

void cleanup_dct(
    // frees memory allocated by setup_dct()
    void   *(*dct_info),// internal data for DCT algorithm
    double *(*in_data)  // input array space
);

#ifdef __cplusplus
}
#endif

#endif
