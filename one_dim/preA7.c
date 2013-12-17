// preparation for calculating A7 of Talyor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "get_toa.h"

double pi=3.1415926;

int preA7 (int *k, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double *s, double *p, int N)
{
	// k is the dimention of amp, N is the dimention of s
	int i;
	//printf ("%d\n", N);
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[N];  // initialize the system, don't know why....

	for (i=0;i<N;i++)
	{
		test[i]=s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	dft_profiles(N,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_s;
	fftw_complex *out_p;
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	dft_profiles(N,s,out_s);
	//printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	dft_profiles(N,p,out_p);

	int n=0;
	double r_s[N],im_s[N];
	double r_p[N],im_p[N];
	//double amp_s[N/2],phi_s[N/2];
	//double amp_p[N/2],phi_p[N/2];

	for (i=0;i<=N/2-1;i++)
	{
		r_s[i]=out_s[i+1][0];
		im_s[i]=out_s[i+1][1];
		r_p[i]=out_p[i+1][0];
		im_p[i]=out_p[i+1][1];
		//printf ("%lf %lf\n", r_p[i], im_p[i]);
		//printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		n++;
	}
	//printf ("%d\n", k);

	for (i=0;i<n;i++)
	{
		amp_s[i]=sqrt(r_s[i]*r_s[i]+im_s[i]*im_s[i]);
		amp_p[i]=sqrt(r_p[i]*r_p[i]+im_p[i]*im_p[i]);
		phi_s[i]=atan2(im_s[i],r_s[i]);
		phi_p[i]=atan2(im_p[i],r_p[i]);
		//printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		//printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		//printf ("%lf\n", amp_s[i]);
		//printf ("%lf\n", amp_p[i]);
	}
	
	(*k)=n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}

