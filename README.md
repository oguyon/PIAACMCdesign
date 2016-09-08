# PIAACMC type coronagraph design tool

## Downloading source code
Latest distribution is on [github](https://github.com/oguyon/PIAACMCdesign).
You can clone the repository, or download the latest .tar.gz distribution.

## Compilation
The source code follows the standard GNU build process:

    ./configure
    make
    make install

## Documentation, Getting Started
Please consult the documentation :  

 > src/PIAACMCsimul/doc/PIAACMCsimul.html


## Libraries
The following libraries are used:
- readline, for reading the command line input
- flex, for parsing the command line input
- bison, to interpret the command line input
- fftw, for performing Fourier Transforms
- gsl, for math functions and tools
- fitsio, for reading and writing FITS image files

If you use NVIDIA GPUs, install cuda and magma libraries, and add "--enable-cuda and --enable-magma" options to the configure command.


## Source Code Architecture 
Written in C.
The main is a command line interface (CLI). Source code is in CLIcore.c and CLIcore.h.
Key data structures (such as the image data structure) are declared in CLIcore.h.

