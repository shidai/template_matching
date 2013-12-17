// calculate A9 of Talyor 1992, and get b  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

extern double *a_s,*a_p,*p_s,*p_p;
extern int num;

double A9 (double phase)
{
	double A9=0.0, sum=0.0;
	int i;

	for (i=0;i<num;i++)
	{
		A9+=a_s[i]*a_p[i]*cos(p_s[i]-p_p[i]+(i+1)*phase);
		sum+=a_s[i]*a_s[i];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
	}
	
	A9=A9/sum;

	return A9;
}
