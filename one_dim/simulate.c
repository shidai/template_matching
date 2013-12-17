// simulate pulse profiles, adding white noise; return simulated profiles
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "get_toa.h"

int simulate (int n, double SNR, double *s, double *p)
{
	// simulate a profile with white noise
	///////////////////////////////////////////////////////////////////////
	// initialize gsl 
	
	int i;
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
  
	////////////////////////////////////////////////////////////////////////
	//  determine the amplitude of white noise according to SNR
	
	double scale;   // the scale multiply to white noise to get certain SNR
	double amp_noise, noise[n];
	double peak_s;

	for (i=0;i<n;i++)
	{
		noise[i]=gsl_ran_gaussian(r,1.0);
		amp_noise+=noise[i]*noise[i];
	}
	
	amp_noise=sqrt(amp_noise/n);

	peak_s=find_peak_value(n,s);   // find the peak flux of the std
    //printf ("peak of std: %g\n", peak_s);

	scale=peak_s/(SNR*amp_noise);
    //printf ("%g\n", scale);
	
	//////////////////////////////////////////////////////////////////////////
	//  add noise to std ==> p

	for (i=0;i<n;i++)
	{
		p[i]=(s[i]+scale*noise[i]);
	}

	/*
	double peak_p;
	peak_p=find_peak(n,p);  // find the peak flux of the profile
    //printf ("peak of profile: %g\n", peak_p);

	//  normalize the std and profile
	
	for (i=0;i<n;i++)
	{
		p[i]=p[i]/peak_p;
		s[i]=s[i]/peak_s;
		//printf ("%g %g\n", s[i], p[i]);
	}
	*/

	/*double rms=0.0;
	int m=0;

	for (i=300;i<700;i++)
	{
		rms+=p[i]*p[i];
		m++;
	}

	rms=sqrt(rms/m);
	//printf("rms is: %f\n", rms);*/

	gsl_rng_free (r);
  
	return 0;
}
