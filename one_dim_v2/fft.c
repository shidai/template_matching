// dft of profiles
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

int dft_profiles (int N, double *in, fftw_complex *out)
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

