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
int error (double phase, double b, double *errphase, double *errb, double a_s[NP], double a_p[NP], double p_s[NP], double p_p[NP], int num)
{
	double rms,gk,s1,s2;
	int j,n;

	gk=0.0;
	n=0;

    for (j = 0; j < num; j++)
    {
	    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
	    gk+=a_p[j]*a_p[j]+b*b*a_s[j]*a_s[j]-2.0*b*a_s[j]*a_p[j]*cos(p_p[j]-p_s[j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(j+1)*(j+1)*a_p[j]*a_s[j]*cos(p_p[j]-p_s[j]-(j+1)*phase);
		s2+=a_s[j]*a_s[j];
		n++;
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*b*s1);
	(*errb)=rms/sqrt(2.0*s2);

	return 0;
}

