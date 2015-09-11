# PIAACMC type coronagraph design tool

## Downloading source code
Latest distribution is on [github](https://github.com/oguyon/PIAACMCdesign).
You can clone the repository, or download the latest .tar.gz distribution.

## Compilation
The source code follows the standard GNU build process:

    ./configure
    make
    make install

## Documentation 
Please consult the [online documentation]{http://oguyon.github.io/PIAACMCdesign/index.html} (generated by doxygen).


## Libraries
The following libraries are used:
- readline, for reading the command line input
- flex, for parsing the command line input
- bison, to interpret the command line input
- fftw, for performing Fourier Transforms
- gsl, for math functions and tools
- fitsio, for reading and writing FITS image files

## Source Code Architecture 
Written in C.
The main is a command line interface (CLI). Source code is in CLIcore.c and CLIcore.h.
Key data structures (such as the image data structure) are declared in CLIcore.h.

## Getting started
Example bash script to run PIAACMC design tool are provided in ./src/PIAACMCdesign/scripts/

To run design tool:
- copy scripts to working directory. A script is provided to do this:

    cd ~/src/PIAACMCdesign/src/PIAACMCsimul/scripts
    ./syncscripts <myworkdirectory>

- copy top level script run to your custom script name:

    cd <myworkdirectory>
    cp run <myrunscriptname>

- edit the top level script as needed

- execute top level script. Executing script with "help" as argument will give usage/example


