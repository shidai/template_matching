#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

//extern double pi;

int readfile ( char *filename, int *ntxt, double *x, double *y );

int dft_profiles (int N, double *in, fftw_complex *out);

int simulate (int n, double SNR, double *s, double *p);

double find_peak_value (int n, double *s);

int preA7 (int *k, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double *s, double *p, int N);

//double A7 (double phase);
double A7 (double phase, double *a_s, double *a_p, double *p_s, double *p_p, int num);

//double A9 (double phase);
double A9 (double phase, double *a_s, double *a_p, double *p_s, double *p_p, int num);

//double zbrent(double (*func)(double), double x1, double x2, double tol);
double zbrent(double (*func)(double phase, double *a_s, double *a_p, double *p_s, double *p_p, int num), double x1, double x2, double tol, double *a_s, double *a_p, double *p_s, double *p_p, int num);

int find_peak (int n, double *s, int *position);

//int error (double phase, double b, double *errphase, double *errb);
int error (double phase, double b, double *errphase, double *errb, double *a_s, double *a_p, double *p_s, double *p_p, int num);
