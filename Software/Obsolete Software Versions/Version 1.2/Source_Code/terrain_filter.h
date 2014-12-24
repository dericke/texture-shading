/*
 * terrain_filter.h
 *
 * Created by Leland Brown on 2010 Oct 30.
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

#ifndef TERRAIN_FILTER_H
#define TERRAIN_FILTER_H

#ifdef __cplusplus
extern "C" {
#endif

enum Terrain_Reg {
    TERRAIN_REG_GRID = 1,   // nodes  - pixel centers on grid (adjacent regions overlap 1 pixel)
    TERRAIN_REG_CELL = 2    // pixels - pixel edges   on grid (adjacent regions don't overlap)
};

enum Terrain_Filter_Errors {
    TERRAIN_FILTER_SUCCESS       = 0,
    TERRAIN_FILTER_MALLOC_ERROR  = 1,   // memory allocation error occurred
    TERRAIN_FILTER_NULL_VALUES   = 2,   // input data contains NaN values
    TERRAIN_FILTER_INVALID_PARAM = 3,   // invalid data registration type
    TERRAIN_FILTER_CANCELED      = -1   // cancellation requested by progress callback function
};

struct Terrain_Progress_Callback {
    int (*callback)(void *state, int count, int total); // return nonzero to cancel operation
    void *state;                                        // optional state information
};

// Computes operator (-Laplacian)^(gain/2) applied to data array.
// Returns 0 on success, nonzero if an error occurred (see enum Terrain_Filter_Errors).
// Mean of data array is always (approximately) zero on output.
int terrain_filter(
    float *data,        // input/output: array of data to process (row-major order)
    double gain,        // input: "gain exponent" to be applied
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double xres,        // input: linear units (e.g., meters, not degrees) per column
    double yres,        // input: linear units (e.g., meters, not degrees) per row
//  enum Terrain_Reg registration,
    const struct Terrain_Progress_Callback
          *progress     // optional callback functor for status; NULL for none
);

// Determines X and Y scales at given latitude for geographic projection
void geographic_scale(
    double  latdeg, // input:  latitude in degrees
    double *xsize,  // output: meters per degree of longitude
    double *ysize   // output: meters per degree of latitude
);

// Determines graticule aspect ratio at given latitude
double geographic_aspect( double latdeg );

// Adjusts output of terrain_filter() for scale variation of Mercator-projected data
void fix_mercator(
    float *data,    // input/output: array of texture shading data (row-major order)
    double gain,    // input: "gain exponent" used to create texture shading
    int    nrows,   // input: number of rows    in data array
    int    ncols,   // input: number of columns in data array
    double lat1deg, // input: latitude at (center or) bottom edge of bottom pixels, degrees
    double lat2deg  // input: latitude at (center or) top    edge of top    pixels, degrees
//  enum Terrain_Reg registration
);

// Adjusts output of terrain_filter() for scale variation of polar stereographic projection
void fix_polar_stereographic(
    float *data,        // input/output: array of data to process (row-major order)
    double gain,        // input: "gain exponent" to be applied
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double center_res   // input: meters/pixel at pole (assumed to be center of array)
);

// Returns meters/pixel at pole (assumed to be center of array) for data in
// polar stereographic projection.
// Note: data is assumed to span no more than 90 degrees latitude - i.e., the center
// of the array must be in the same hemisphere (north or south) as the corners.
double polar_stereographic_center_res(
    int    nrows,           // input: number of rows    in data array
    int    ncols,           // input: number of columns in data array
    double corner_latdeg    // input: latitude at (center or) outer corner of corner pixels
);

// Adjusts the tone curve (contrast and saturation limits); assumes data mean is zero on input
void adjust_contrast(
    float *data,        // input/output: array of data to process (row-major order)
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double new_mean,    // input: adjust mean of data to this value
    double half_span,   // input: on return, num_sigmas standard deviations = half_span
    double num_sigmas   // input: # of std deviations to fit between new_mean +/- half_span
);

#ifdef __cplusplus
}
#endif

#endif
