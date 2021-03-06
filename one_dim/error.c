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
int error (double phase, double b, double *errphase, double *errb, double *a_s, double *a_p, double *p_s, double *p_p, int num)
{
	double rms,gk,s1,s2;
	int i,n;

	gk=0.0;
	n=0;

	for (i=0;i<num;i++)
	{
		//gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(i+1)*(i+1)*a_p[i]*a_s[i]*cos(p_p[i]-p_s[i]-(i+1)*phase);
		s2+=a_s[i]*a_s[i];
		n++;
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*b*s1);
	(*errb)=rms/sqrt(2.0*s2);

	return 0;
}
