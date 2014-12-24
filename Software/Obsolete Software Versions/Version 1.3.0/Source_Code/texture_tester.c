/*
 * texture_tester.c
 *
 * Created by Leland Brown on 2013 Nov 17.
 *
 * Copyright (c) 2013 Leland Brown.
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

#include "terrain_filter.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

// CAUTION: This __DATE__ is only updated when THIS file is recompiled.
// If other source files are modified but this file is not touched,
// the version date may not be correct.
static const char version[] = "1.3, built " __DATE__;

static struct timeval tp1, tp2 = { 0, 0 };

static void report_time()
{
}

static const char *command_name;

static const char *get_command_name( const char *argv[] )
{
    const char *colon;
    const char *slash;
    const char *result;
    
    colon = strchr( argv[0], ':' );
    if (colon) {
        ++colon;
    } else {
        colon = argv[0];
    }
    slash = strrchr( colon, '/' );
    if (slash) {
        ++slash;
    } else {
        slash = colon;
    }
    result = strrchr( slash, '\\' );
    if (result) {
        ++result;
    } else {
        result = slash;
    }
    return result;
}

static void prefix_error()
{
    fprintf( stderr, "\n*** ERROR: " );
}

static void usage_exit( const char *message )
{
    if (message) {
        prefix_error();
        fprintf( stderr, "%s\n", message );
    }
    fprintf( stderr, "\n" );
    fprintf( stderr, "USAGE:    %s <nsize>\n", command_name );
    fprintf( stderr, "\n" );
    exit( EXIT_FAILURE );
}

static int print_progress( float portion, float steps_done, int total_steps, void *state )
{
    int *last_count = (int *)state;
    int  this_count = (int)steps_done;
    
    if (this_count > *last_count) {
        tp1 = tp2;
        gettimeofday(&tp2, NULL);
        float tdiff = (tp2.tv_sec - tp1.tv_sec) * 1000000.0 + (tp2.tv_usec - tp1.tv_usec);
        printf("Step time (microseconds) = %.0f\n", tdiff);
        printf( "Processing phase %d...\n", this_count + 1 );
        fflush( stdout );
        *last_count = this_count;
    }
    
    return 0;
}

#ifndef NOMAIN

int main( int argc, const char *argv[] )
{
    const int numargs = 2;  // including command name
    
    int last_count = -1;

    struct Terrain_Progress_Callback progress = { print_progress, &last_count };

    char *endptr;

    int nrows;
    int ncols;

    float *data;
    
    int error;

    printf( "\nTiming tester for terrain texture shading program - version %s\n", version );

    // Validate parameters:

//  command_name = "TEXTURE_TESTER";
    command_name = get_command_name( argv );

    if (argc == 1) {
        usage_exit( 0 );
    } else if (argc < numargs) {
        usage_exit( "Not enough command-line parameters." );
    } else if (argc > numargs) {
        usage_exit( "Too many command-line parameters." );
    }
        
    nrows = (int)strtol( argv[1], &endptr, 10 );
    ncols = nrows;
    if (endptr == argv[1] || *endptr != '\0') {
        usage_exit( "First parameter (detail) must be a number." );
    }
    
    data = (float *)malloc( (long)nrows * (long)ncols * sizeof( float ) );
    
    memset( data, 0, (long)nrows * (long)ncols * sizeof( float ) );

    // Process data:

    printf( "Processing %d column x %d row array...\n", ncols, nrows );
    fflush( stdout );

    gettimeofday(&tp2, NULL);
    error = terrain_filter(
        data, 1.0, nrows, ncols, 1.0, 1.0, TERRAIN_METERS, 0.0, &progress );

    if (error) {
        assert( error == TERRAIN_FILTER_MALLOC_ERROR );
        prefix_error();
        fprintf( stderr, "Memory allocation error occurred during processing of data.\n" );
        exit( EXIT_FAILURE );
    }
    
    free( data );
    
    printf( "DONE.\n" );

    return EXIT_SUCCESS;
}

#endif
