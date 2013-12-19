// preparation for calculating A7 of Talyor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "get_toa.h"

double pi=3.1415926;

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn)
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////
	double test[nphase];  // initialize the system, don't know why....

	for (i=0;i<nphase;i++)
	{
		test[i]=s[i];
	}
	fftw_complex *out_t;
	out_t = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	dft_profiles(nphase,test,out_t);
	//////////////////////////////////////////////////////////////////////////////

    fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    amp_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    amp_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    phi_s[i][j]=atan2(im_s[j],r_s[j]);
		    phi_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	(*k)=n;

	fftw_free(out_s); 
	fftw_free(out_p); 
	fftw_free(out_t); 

	return 0;
}

