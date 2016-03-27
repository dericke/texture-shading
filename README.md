#User’s Guide for Texture Shading Program Version 1.3
(DRAFT)

*From User_Guide.pdf by Leland Brown*

## Purpose and Background Information

This software takes terrain elevation data and performs a transformation
I call “texture shading,” which consists primarily of a fractional
Laplacian operator. The output is then represented as a grayscale raster
image. Texture shading can function as an alternative to traditional
hillshading as a technique for producing shaded relief---or the two can
be used in combination to enhance standard hillshading. For a
description and some examples of the texture shading method, as well as
an overview of the mathematics behind it, see the annotated version of
the presentation I gave at NACIS 2010. A copy is available with this
software; the file is called "Texture_Shading_with_Notes.pdf".

## Compiling from Source

Using the GCC compiler, a recommended set of options for building the
texture shading programs is as follows:

```sh
gcc -O2 -funroll-loops -DNOMAIN -c *.c -lm
gcc -O2 -funroll-loops *.o texture.c -o texture -lm
gcc -O2 -funroll-loops *.o texture_image.c -o texture_image -lm
```

Use these commands on Mac/Linux/Unix. Or if using MinGW on Windows, add
".exe" to the output filenames “texture” and “texture_image”.The
programs should compile with other C compilers as well, with minimal or
no modifications, though this has not been tested.

## Source Files

The source code can be grouped into layers, with a middle layer that
implements the texture shading algorithm itself, a top layer that
provides the command-line interface and takes care of the file I/O, and
a bottom layer that does the FFT computations needed. The interface to
the FFT layer is designed to be generic enough to allow replacement of
the FFT implementation included here with any other general-purpose FFT
library if so desired. The included FFT implementation is also
interesting in its own right; it performs discrete cosine transforms
(DCT) on data of any size, and handles even prime sizes with reasonable
efficiency.

### Top Layer:
```
texture.c (main function for “texture” program)
texture_image.c (main function for “texture_image” program)
read_grid_files.h
read_grid_files.c
write_grid_files.h
write_grid_files.c
WriteGrayscaleTIFF.h
WriteGrayscaleTIFF.c
```
### Middle Layer:
```
terrain_filter.h (interface to texture shading algorithm)
terrain_filter.c
compatibility.h (helps deal with compiler dependencies)
dct.h (generic interface to FFT code in bottom layer)
```
### Bottom Layer:
```
compatibility.h (also used by middle layer)
dct_fftpack.c (adapts FFT library to dct.h interface)
fftpack.h (general-purpose DCT and real FFT library)
fftpack.ctranspose_inplace.h (general-purpose matrix transpose)
transpose_inplace.c
```
