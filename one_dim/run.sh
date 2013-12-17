#!/bin/sh

gcc -Wall -lm -lgsl -lfftw3 simulate.c A7.c A9.c brent.c error.c fft.c find_peak.c find_peak_value.c get_toa.c preA7.c readfile.c -o get_toa.out
