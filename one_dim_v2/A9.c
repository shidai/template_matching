// calculate A9 of Talyor 1992, and get b  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

//extern double *a_s,*a_p,*p_s,*p_p;
//extern int num;

double A9 (double phase, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num)
{
	double A9=0.0, sum=0.0;
	int j;

    for (j = 0; j < num; j++)
    {
	    A9+=a_s[j]*a_p[j]*cos(p_s[j]-p_p[j]+(j+1)*phase);
	    sum+=a_s[j]*a_s[j];
	    //printf ("%lf %lf\n", a_s[i], p_s[i]);
	}
	
	A9=A9/sum;

	return A9;
}

