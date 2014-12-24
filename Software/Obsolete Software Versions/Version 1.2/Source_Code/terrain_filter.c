/*
 * terrain_filter.c
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

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "terrain_filter.h"

#include "transpose_inplace.h"
#include "dct.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


static const double equatorial_radius = 6378137.0;      // WGS84 value in meters
static const double flattening = 1.0 / 298.257223563;   // WGS84 value
static const double ecc = 0.0818191908426215;           // first eccentricity
                     // = sqrt( (2.0 - flattening) * flattening )


// Utility function:

static int flt_isnan( float x )
{
    volatile float y = x;
    return y != y;
}


static double conformal_lat( double lat );

static double isometric_lat( double lat );

static double tan_lat_from_tan_conformal( double tan_conlat );

static double tan_lat_from_isometric( double isolat );

static double mercator_relscale_from_tan_lat( double tan_lat );


void adjust_contrast(
    float *data,        // input/output: array of data to process (row-major order)
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double new_mean,    // input: adjust mean of data to this value
    double half_span,   // input: on return, num_sigmas standard deviations = half_span
    double num_sigmas   // input: # of std deviations to fit between new_mean +/- half_span
)
// Adjusts the tone curve (contrast and saturation limits); assumes data mean is zero on input
{
    int i, j;
    float *ptr;

    double sum;
    double row_sum;
    double std_dev;
    double count;
    double factor;
    
    // Compute standard deviation of data, assuming mean = 0.0:

    count = (double)(ncols * nrows);

    sum = 0.0;
    for (i=0, ptr=data; i<nrows; ++i, ptr+=ncols) {
        row_sum = 0.0;
        for (j=0; j<ncols; ++j) {
            row_sum += ptr[j] * ptr[j];
        }
        sum += row_sum;
    }

    std_dev = sqrt( sum / count );
    
    // Transform data values so that num_sigmas * std_dev = half_span:
    
    factor = half_span / (num_sigmas * std_dev);

    for (i=0, ptr=data; i<nrows; ++i, ptr+=ncols) {
        for (j=0; j<ncols; ++j) {
            ptr[j] *= factor;

        //  // Saturate at +/-half_span:
        //
        //  if (ptr[j] > +half_span) {
        //      ptr[j] = +half_span;
        //  else if (ptr[j] < -half_span) {
        //      ptr[j] = -half_span;
        //  }
        
            ptr[j] += new_mean;
        }
    }
}


static double conformal_lat( double lat )
{
    // exact formula based on isometric latitude:
//  double isolat = isometric_lat( lat );
//  double temp = exp( fabs(isolat) );
//  double tan_conlat = copysign( temp - 1.0/temp, isolat ) * 0.5;
//  double conlat = atan( tan_conlat );

//  double conlat = lat;    // spherical earth approximation

    double geocentric_lat = atan( ((1.0-ecc)*(1.0+ecc)) * tan(lat) );
    double conlat = geocentric_lat; // better approximation

    return conlat;
}

static double isometric_lat( double lat )
{
    // exact formula, but difficult to invert:
//  double sin_lat = sin( lat );
//  double isolat =
//      ( log(1.0+    sin_lat) - log(1.0-    sin_lat) ) * 0.5 +
//      ( log(1.0-ecc*sin_lat) - log(1.0+ecc*sin_lat) ) * 0.5 * ecc

    double conlat = conformal_lat( lat );
    double sin_conlat = sin( conlat );
    double cos_conlat = cos( conlat );
    double isolat = copysign( log( cos_conlat / ( 1.0+fabs(sin_conlat) ) ), sin_conlat );
               // = ( log(1.0+sin_conlat) - log(1.0-sin_conlat) ) * 0.5
               // = atanh( tan(conlat*0.5) ) * 2.0
    
    return isolat;
}

static double tan_lat_from_tan_conformal( double tan_conlat )
{
//  double tan_lat = tan_conlat;    // spherical earth approximation

    double tan_geocentric_lat = tan_conlat; // better approximation
    double tan_lat = tan_geocentric_lat / ( (1.0-ecc) * (1.0+ecc) );

    return tan_lat;
}

static double tan_lat_from_isometric( double isolat )
{
    double temp = exp( fabs(isolat) );
    double tan_conlat = copysign( temp - 1.0/temp, isolat ) * 0.5;
    
    return tan_lat_from_tan_conformal( tan_conlat );
}

static double mercator_relscale_from_tan_lat( double tan_lat )
{
    double relscale = sqrt( 1.0 + ((1.0-ecc)*(1.0+ecc)) * tan_lat * tan_lat );
                 // = sqrt( (1.0-ecc*sin(lat)) * (1.0+ecc*sin(lat)) ) / cos(lat)
                 // = secant of parametric (reduced) latitude
    return relscale;
}

void geographic_scale(
    double  latdeg, // input:  latitude in degrees
    double *xsize,  // output: meters per degree of longitude
    double *ysize   // output: meters per degree of latitude
)
// Determines X and Y scales at given latitude for geographic projection
{
    double lat = (M_PI/180.0) * latdeg;

    double temp1 = ecc * sin( lat );
    double temp2 = (1.0 - temp1) * (1.0 + temp1);
    double temp3 = (1.0 - ecc)   * (1.0 + ecc);
    
    double normal_radius     = equatorial_radius / sqrt( temp2 );
    double meridional_radius = normal_radius * temp3 / temp2;

    *xsize = normal_radius * cos( lat );
    *ysize = meridional_radius;

    *xsize *= M_PI/180.0;
    *ysize *= M_PI/180.0;
}

double geographic_aspect( double latdeg )
// Determines graticule aspect ratio at given latitude
{
    double xsize, ysize;
    
    geographic_scale( latdeg, &xsize, &ysize );
    
    // Return graticule aspect ratio (width/height in linear units)
    return xsize / ysize;
}

void fix_mercator(
    float *data,    // input/output: array of data to process (row-major order)
    double gain,    // input: "gain exponent" to be applied
    int    nrows,   // input: number of rows    in data array
    int    ncols,   // input: number of columns in data array
    double lat1deg, // input: latitude at (center or) bottom edge of bottom pixels, degrees
    double lat2deg  // input: latitude at (center or) top    edge of top    pixels, degrees
//  enum Terrain_Reg registration
)
// Adjusts output of terrain_filter() for scale variation of Mercator-projected data
{
    enum Terrain_Reg registration = TERRAIN_REG_CELL;

    double ypix1;
    double ypix2;

    double pix2merc;
    double isolat0;
    double ypix0;
    
    int i, j;
    float *ptr;

    // convert latitudes from degrees to radians
    double lat1 = (M_PI/180.0) * lat1deg;
    double lat2 = (M_PI/180.0) * lat2deg;

    double isolat1 = isometric_lat( lat1 );
    double isolat2 = isometric_lat( lat2 );

    switch (registration) {
        case TERRAIN_REG_GRID:
            ypix1 = (double)(nrows - 1);    // center of bottom row of pixels
            ypix2 = 0.0;                    // center of top    row of pixels
            break;
        case TERRAIN_REG_CELL:
            ypix1 = (double)nrows - 0.5;    // bottom edge of bottom row of pixels
            ypix2 = -0.5;                   // top    edge of top    row of pixels
            break;
        default:
            // invalid data registration type
            ypix1 = 0.0;
            ypix2 = 0.0;
    }

    pix2merc = (isolat2 - isolat1) / (ypix2 - ypix1);
    isolat0  = (isolat1 + isolat2) / 2;
    ypix0    = (ypix1   + ypix2)   / 2;

    for (i=0, ptr=data; i<nrows; ++i, ptr+=ncols) {
        double ypix = (double)i;
        double isolat = isolat0 + (ypix - ypix0) * pix2merc;
        double tan_lat = tan_lat_from_isometric( isolat );
        double relscale = mercator_relscale_from_tan_lat( tan_lat );
        
        double zfactor = pow( relscale, gain );
        for (j=0; j<ncols; ++j) {
            ptr[j] *= zfactor;
        }
    }
}

void fix_polar_stereographic(
    float *data,        // input/output: array of data to process (row-major order)
    double gain,        // input: "gain exponent" to be applied
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double center_res   // input: meters/pixel at pole (assumed to be center of array)
)
// Adjusts output of terrain_filter() for scale variation of polar stereographic projection
{
    const double rfactor = sqrt( ((1+ecc)*(1-ecc)) * pow( (1+ecc)/(1-ecc), ecc ) ) * 0.5;

    double relscale;
    double zfactor;
    
    int i, j;
    float *ptr;

    double xres = center_res;
    double yres = center_res;

    // assume pole is in center of array
    double ipole = 0.5 * (double)(nrows-1);
    double jpole = 0.5 * (double)(ncols-1);

    double xfactor = xres / equatorial_radius;
    double yfactor = yres / equatorial_radius;

    for (i=0, ptr=data; i<nrows; ++i, ptr+=ncols) {
        double idiff = ((double)i - ipole) * yfactor;
        for (j=0; j<ncols; ++j) {
            double jdiff = ((double)j - jpole) * xfactor;
            double r = sqrt( idiff*idiff + jdiff*jdiff );
            if (r > 0.0) {
                double temp = r * rfactor;
                         // = exp( -isometric_lat ) for north polar projection
                         // = exp( +isometric_lat ) for south polar projection
                double tan_conlat = ( 1.0/temp - temp ) * 0.5;
                double tan_lat = tan_lat_from_tan_conformal( tan_conlat );
                relscale = r * mercator_relscale_from_tan_lat( tan_lat );
            } else {
                relscale = 1.0;
            }
            
            zfactor = pow( relscale, gain );
            ptr[j] *= zfactor;
        }
    }
}

double polar_stereographic_center_res(
    int    nrows,           // input: number of rows    in data array
    int    ncols,           // input: number of columns in data array
    double corner_latdeg    // input: latitude at (center or) outer corner of corner pixels
)
// Returns meters/pixel at pole (assumed to be center of array), if data is in
// polar stereographic projection.
// Note: data is assumed to span no more than 90 degrees latitude - i.e., the center
// of the array must be in the same hemisphere (north or south) as the corners.
{
    enum Terrain_Reg registration = TERRAIN_REG_CELL;

    const double rfactor = sqrt( ((1+ecc)*(1-ecc)) * pow( (1+ecc)/(1-ecc), ecc ) ) * 0.5;

    int idiff, jdiff;
    double conlat;
    double temp;
    double rcorner;
    double center_res;
    
    double corner_lat = (M_PI/180.0) * corner_latdeg;
    
    switch (registration) {
        case TERRAIN_REG_GRID:
            // assume pole is in center of array
            idiff = 0.5 * (double)(nrows-1);
            jdiff = 0.5 * (double)(ncols-1);
            break;
        case TERRAIN_REG_CELL:
            // assume pole is in center of array
            idiff = 0.5 * (double)nrows;
            jdiff = 0.5 * (double)ncols;
            break;
        default:
            // invalid data registration type
            return 0.0;
    }
    
    rcorner = sqrt( idiff*idiff + jdiff*jdiff );
    conlat = conformal_lat( corner_lat );
    temp = cos(conlat) / ( 1.0 + fabs( sin(conlat) ) );
             // = exp( -fabs(isometric_lat) )
    center_res = temp / rcorner * (equatorial_radius / rfactor);
    
    return center_res;
}


// 2-dimensional discrete cosine transform:

static int dct_2d_and_transpose( int nrows, int ncols, float *data, int dct_type );

static int dct_2d_and_transpose( int nrows, int ncols, float *data, int dct_type )
{
    // Perform 2-D DCT:
    
    int error;

    int i, j;
    float *ptr;
    
    void   *dct_info = NULL;
    double *in_data  = NULL;
    double *out_data;
    
    error = setup_dct( dct_type, ncols, &dct_info, &in_data );
    if (error) {
        return 1;
    }
    for (i=0, ptr=data; i<nrows; ++i, ptr+=ncols) {
        for (j=0; j<ncols; ++j) {
            in_data[j] = ptr[j];
        }
        out_data = perform_dct( dct_info, in_data );
        for (j=0; j<ncols; ++j) {
            ptr[j] = out_data[j];
        }
    }
    cleanup_dct( &dct_info, &in_data );

    error = transpose_inplace( data, nrows, ncols );
    if (error) {
        return 1;
    }

    error = setup_dct( dct_type, nrows, &dct_info, &in_data );
    if (error) {
        return 1;
    }
    for (i=0, ptr=data; i<ncols; ++i, ptr+=nrows) {
        for (j=0; j<nrows; ++j) {
            in_data[j] = ptr[j];
        }
        out_data = perform_dct( dct_info, in_data );
        for (j=0; j<nrows; ++j) {
            ptr[j] = out_data[j];
        }
    }
    cleanup_dct( &dct_info, &in_data );
    
    return 0;
}


// Thin-plate spline interpolator:

typedef struct {
    double a, b;
    double s;
    double temp1, temp2;
    double offset;
} Thin_Plate_Info;

static Thin_Plate_Info setup_thin_plate( double a, double b, double s );
static double thin_plate_spline( const Thin_Plate_Info  *info, double x, double y );

static Thin_Plate_Info setup_thin_plate( double a, double b, double s )
// Requires s > 1.0.
// Assumes a > 0.0 and b > 0.0.
{
    Thin_Plate_Info info;

    double temp3 = pow( 4.0, -s );
    
    info.a = a;
    info.b = b;
    info.s = s;
    
    info.temp1 = 4.0 * s * a * pow( a+0.25*b, -s-1.0 );
    info.temp2 = 4.0 * s * b * pow( b+0.25*a, -s-1.0 );

    info.offset  = pow( 0.25*a+2.25*b, -s ) + pow( 2.25*a+0.25*b, -s );
    info.offset += pow( 2.25*a+2.25*b, -s );
    info.offset += 0.0625 * temp3 * (info.temp1 + info.temp2);

    info.offset /= (0.25 - temp3);

    return info;
}

static double thin_plate_spline( const Thin_Plate_Info *info, double x, double y )
// Approximates sum(i=-inf..+inf, sum(j=-inf..+inf, (a*(x+i)^2 + b*(y+j)^2)^-s ) ).
// Requires 0.0 <= x < 1.0 and 0.0 <= y < 1.0.
// Assumes x != 0.0 OR y != 0.0.
// Assumes info->a > 0.0 and info->b > 0.0.
{
    double x1, y1, x2, y2;
    double result;
    
    double a = info->a;
    double b = info->b;
    double s = info->s;
    
//  x = x - floor(x);   // now 0.0 <= x < 1.0
//  y = y - floor(y);   // now 0.0 <= y < 1.0

    x1 = 1.0 - x;
    y1 = 1.0 - y;
    
    x2 = x - 0.5;
    y2 = y - 0.5;

    result  = pow( a*x1*x1 + b*y1*y1, -s ) + pow( a*x*x + b*y*y,   -s );
    result += pow( a*x1*x1 + b*y*y,   -s ) + pow( a*x*x + b*y1*y1, -s );
    result += (x2*x2) * info->temp1 + (y2*y2) * info->temp2;
    result += info->offset;
    
    return result;
}

// Higher-fidelity approximation to thin-plate spline:

typedef struct {
    Thin_Plate_Info info;
} Thin_Plate_Info2;

static Thin_Plate_Info2 setup_thin_plate2( double a, double b, double s );
static double thin_plate_spline2( const Thin_Plate_Info2 *info, double x, double y );

static Thin_Plate_Info2 setup_thin_plate2( double a, double b, double s )
{
    Thin_Plate_Info2 info;
    info.info = setup_thin_plate( 4.0*a, 4.0*b, s );
    return info;
}

static double thin_plate_spline2( const Thin_Plate_Info2 *info, double x, double y )
{
    double result =
        thin_plate_spline( &info->info, x*0.5,     y*0.5 ) +
        thin_plate_spline( &info->info, x*0.5+0.5, y*0.5+0.5 );
    result +=
        thin_plate_spline( &info->info, x*0.5,     y*0.5+0.5 ) +
        thin_plate_spline( &info->info, x*0.5+0.5, y*0.5 );
    return result;
}


// Fractional Laplacian operator:

static int apply_operator(
    float *data, double gain, int m, int n, double xscale, double yscale, enum Terrain_Reg registration );

static int apply_operator(
    float *data, double gain, int m, int n, double xscale, double yscale, enum Terrain_Reg registration )
{
    const int spline_type = 0;  // bicubic spline
//  const int spline_type = 1;  // thin-plate spline (fails for gain >= 2.0, unstable near 2.0)

    int m2, n2;
    int i, j;
    int j0;
    
    float *ptr;
    
    double *storage;
    double *nux, *xx, *splinex;
    double *nuy, *yy, *spliney;
    
    double power = gain * 0.5;

    double factor;
    double xfactor, yfactor;
    
    switch (registration) {
        case TERRAIN_REG_GRID:
            m2 = 2*m-2;
            n2 = 2*n-2;
            factor = 1.0 / (double)(n2*m2);     // DCT normalization factor
            break;
        case TERRAIN_REG_CELL:
            m2 = m+m;
            n2 = n+n;
            factor = 1.0 / (double)(n2*m2*4);   // DCT normalization factor
            break;
        default:
            return TERRAIN_FILTER_INVALID_PARAM;
    }
    
    factor *= pow( 2.0*M_PI, gain );    // constant factor for fractional Laplacian

    xfactor = 1.0 / (double)m2;
    yfactor = 1.0 / (double)n2;
    
    storage = (double *)malloc( sizeof( double ) * (m2+n2+2) * 3 );
    if (!storage) {
        return TERRAIN_FILTER_MALLOC_ERROR;
    }
    nux     = storage;          // size m2+1
    nuy     = (m2+1) + nux;     // size n2+1
    xx      = (n2+1) + nuy;     // size m2+1
    yy      = (m2+1) + xx;      // size n2+1
    splinex = (n2+1) + yy;      // size m2+1
    spliney = (m2+1) + splinex; // size n2+1

    for (i=0; i<=m2; ++i) {
        nux[i] = xfactor * (double)i;
    }

    for (j=0; j<=n2; ++j) {
        nuy[j] = yfactor * (double)j;
    }

    data[0] = 0.0;  // set "DC" component to zero
    
    if (spline_type == 1) {

        // Use approximation to thin-plate spline interpolator

        Thin_Plate_Info numer_info =
            setup_thin_plate( xscale * xscale, yscale * yscale, 2.0 - power );
        Thin_Plate_Info denom_info =
            setup_thin_plate( xscale * xscale, yscale * yscale, 2.0 );
        
        // skip i=0,j=0
        for (i=0, j0=1, ptr=data; i<m; ++i, j0=0, ptr+=n) {
            for (j=j0; j<n; ++j) {
                double numer = thin_plate_spline( &numer_info, nux[i], nuy[j] );
                double denom = thin_plate_spline( &denom_info, nux[i], nuy[j] );

                ptr[j] *= numer / denom;
                ptr[j] *= factor;   // normalization factor
            }
        }

    } else {
    
        // Use approximation to bicubic spline interpolator

        double xfreq, cosx, tempx;
        double yfreq, cosy, tempy;

        for (i=0; i<=m2; ++i) {
            xfreq = xscale * nux[i];
            xx[i] = xfreq * xfreq;
            cosx  = cos(M_PI * nux[i]);
            tempx = cosx * 0.5 + 0.5;
            splinex[i] = (tempx*tempx) * ( (cosx+2.0) / (cosx*cosx*2.0+1.0) );
        }

        for (j=0; j<=n2; ++j) {
            yfreq = yscale * nuy[j];
            yy[j] = yfreq * yfreq;
            cosy  = cos(M_PI * nuy[j]);
            tempy = cosy * 0.5 + 0.5;
            spliney[j] = (tempy*tempy) * ( (cosy+2.0) / (cosy*cosy*2.0+1.0) );
        }

        // skip i=0,j=0
        for (i=0, j0=1, ptr=data; i<m; ++i, j0=0, ptr+=n) {
            for (j=j0; j<n; ++j) {
                double mult1 = splinex[i]    * spliney[j];
                double mult2 = splinex[i]    * spliney[n2-j];
                double mult3 = splinex[m2-i] * spliney[j];
                double mult4 = splinex[m2-i] * spliney[n2-j];

                mult1 *= pow( xx[i]    + yy[j],    power );
                mult2 *= pow( xx[i]    + yy[n2-j], power );
                mult3 *= pow( xx[m2-i] + yy[j],    power );
                mult4 *= pow( xx[m2-i] + yy[n2-j], power );

                ptr[j] *= ( (mult1 + mult4) + (mult2 + mult3) );
                ptr[j] *= factor;   // normalization factor
            }
        }
        
    }   // spline_type

    free( storage );
    
    return TERRAIN_FILTER_SUCCESS;
}


// Main terrain_filter function:

int terrain_filter(
    float *data,        // input/output: array of data to process (row-major order)
    double gain,        // input: "gain exponent" to be applied
    int    nrows,       // input: number of rows    in data array
    int    ncols,       // input: number of columns in data array
    double xres,        // input: linear units (e.g., meters, not degrees) per column
    double yres,        // input: linear units (e.g., meters, not degrees) per row
    const struct Terrain_Progress_Callback
        *progress       // optional callback functor for status; NULL for none
//  enum Terrain_Reg registration
)
// Computes operator (-Laplacian)^(gain/2) applied to data array.
// Returns 0 on success, nonzero if an error occurred (see Terrain_Filter_Errors).
// Mean of data array is always (approximately) zero on output.
{
    enum Terrain_Reg registration = TERRAIN_REG_CELL;
    
    int error;

    int type_fwd, type_bwd;
    
    double xscale = 1.0 / xres;
    double yscale = 1.0 / yres;

    double dscale;

    xscale = fabs( xscale );
    yscale = fabs( yscale );

    // check for yscale > xscale, but allow a little slack (< 1/4 pixel relative error)
    dscale = yscale - xscale;
    if (dscale * ncols >= 0.25 * xscale || dscale * nrows >= 0.25 * yscale) {
        fprintf( stderr, "*** WARNING: " );
        fprintf( stderr, "Unusual pixel aspect ratio (> 1.0). Is this correct?\n" );
    }

    switch (registration) {
        case TERRAIN_REG_GRID:
            type_fwd = 1;
            type_bwd = 1;
            break;
        case TERRAIN_REG_CELL:
            type_fwd = 2;
            type_bwd = 3;
            break;
        default:
            return TERRAIN_FILTER_INVALID_PARAM;
    }

    // Perform 2-D DCT (and transpose matrix):

    if (progress && progress->callback( progress->state, 0, 3 )) {
        return TERRAIN_FILTER_CANCELED;
    }

    error = dct_2d_and_transpose( nrows, ncols, data, type_fwd );

    if (error) {
        return TERRAIN_FILTER_MALLOC_ERROR;
    }

    if (flt_isnan( data[0] )) {
        return TERRAIN_FILTER_NULL_VALUES;
    }

    if (progress && progress->callback( progress->state, 1, 3 )) {
        return TERRAIN_FILTER_CANCELED;
    }
    
    // Apply fractional Laplacian operator as Fourier multiplier

    error = apply_operator( data, gain, ncols, nrows, xscale, yscale, registration );

    if (error) {
        return error;
    }

    // Perform 2-D DCT (and transpose matrix):
    
    if (progress && progress->callback( progress->state, 2, 3 )) {
        return TERRAIN_FILTER_CANCELED;
    }
    
    error = dct_2d_and_transpose( ncols, nrows, data, type_bwd );

    if (error) {
        return TERRAIN_FILTER_MALLOC_ERROR;
    }

    if (progress && progress->callback( progress->state, 3, 3 )) {
        return TERRAIN_FILTER_CANCELED;
    }

    return TERRAIN_FILTER_SUCCESS;
}
