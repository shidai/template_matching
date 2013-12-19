#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"

#define NP 512
#define PI 3.14159265359

long int stt_imjd ( char *name );
long int stt_smjd ( char *name );
double stt_offs ( char *name );

int get_nchan ( char *name );
int get_npol ( char *name );
int get_subint ( char *name );

int read_prof ( char *name, int subint, double *profile );
int print_t2pred ( char *name );
double read_freq ( char *name, int subint );
double read_offs ( char *name, int subint);

int readfile ( char *filename, int *ntxt, double *x, double *y );

int dft_profiles (int N, double *in, fftw_complex *out);

int simulate (int n, double SNR, double *s, double *p);

double find_peak_value (int n, double *s);

int preA7 (int *k, double amp_s[NP], double amp_p[NP], double phi_s[NP], double phi_p[NP], double *s, double *p, int nphase);

//double A7 (double phase);
double A7 (double phase, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num);

//double A9 (double phase);
double A9 (double phase, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num);

//double zbrent(double (*func)(double), double x1, double x2, double tol);
double zbrent(double (*func)(double phase, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num), double x1, double x2, double tol, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num);


int find_peak (int n, double *s, int *position);

//int error (double phase, double b, double *errphase, double *errb);
int error (double phase, double b, double *errphase, double *errb, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num);

int get_toa (double s[1024], double p[1024], double *phase, double *errphase);

