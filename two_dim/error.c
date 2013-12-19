// calculate the errors of phase, a and b according to Talyor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

//extern double *a_s,*a_p,*p_s,*p_p;
//extern int num;

//int error (double phase, double b, double a,  double *errphase, double *errb)
int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
{
	double rms,gk,s1,s2;
	int i,j,n;

	gk=0.0;
	n=0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=a_p[i][j]*a_p[i][j]+b*b*a_s[i][j]*a_s[i][j]-2.0*b*a_s[i][j]*a_p[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(j+1)*(j+1)*a_p[i][j]*a_s[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase);
		s2+=a_s[i][j]*a_s[i][j];
		n++;
		}
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*b*s1);
	(*errb)=rms/sqrt(2.0*s2);

	return 0;
}

int error_multi (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
{
	double s1,s2;
	int i,j,n;

	n=0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    s1+=((j+1)*(j+1)*a_p[i][j]*a_s[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase))/(rms[i]*rms[i]);
		    s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		n++;
		}
	}
	
	(*errphase)=1.0/sqrt(2.0*b*s1);
	(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}
