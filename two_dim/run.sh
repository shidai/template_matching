#!/bin/sh

#gcc -lm -lgsl -lfftw3 simulate.c A7.c A9.c brent.c error.c fft.c find_peak.c find_peak_value.c get_toa.c preA7.c readfile.c -o get_toa.out
gcc -Wall -lm -lfftw3 -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio A7.c A9.c brent.c error.c fft.c find_peak.c find_peak_value.c get_toa.c preA7.c readfile.c main.c rms.c readfits.c simulatePPTA.c T2toolkit.c simulatePseudoBB.c tempo2pred.c cheby2d.c t1polyco.c -o get_toa.out
#gcc -Wall -lm -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio readfits.c 
