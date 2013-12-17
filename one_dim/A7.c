// calculate A7 of Talyor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

extern double *a_s,*a_p,*p_s,*p_p;
extern int num;

double A7 (double phase)
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i;

	for (i=0;i<num;i++)
	{
		A7+=(i+1)*a_s[i]*a_p[i]*sin(p_s[i]-p_p[i]+(i+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
	}
	
	return A7;
}
