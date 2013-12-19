// calculate A7 of Talyor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

//extern double *a_s,*a_p,*p_s,*p_p;
//extern int num;

double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < num; j++)
	    {
			A7+=(j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}

double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < nchn; i++)
	{
		for (j = 0; j < num; j++)
	    {
			A7+=((j+1)*a_s[i][j]*a_p[i][j]*sin(p_s[i][j]-p_p[i][j]+(j+1)*phase))/(rms[i]*rms[i]);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}
