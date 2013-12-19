// calculate the rms of each subchannel  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

//extern double *a_s,*a_p,*p_s,*p_p;
//extern int num;

//int error (double phase, double b, double a,  double *errphase, double *errb)
int cal_rms (double phase, double b, double *rms, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
{
	double gk;
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
		n++;
		}
	}
	
	(*rms)=sqrt(gk/n);

	return 0;
}
