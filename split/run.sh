#!/bin/sh

rm *.fits xx
gcc -Wall -lm -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio split.c 
