// calculate A9 of Talyor 1992, and get b  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

//extern double *a_s,*a_p,*p_s,*p_p;
//extern int num;

double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    A9+=a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		    sum+=a_s[i][j]*a_s[i][j];
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}

double A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    A9+=(a_s[i][j]*a_p[i][j]*cos(p_s[i][j]-p_p[i][j]+(j+1)*phase))/(rms[i]*rms[i]);
		    sum+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}
