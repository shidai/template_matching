// find the peak flux of profiles
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "get_toa.h"

int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}
	peak=temp[n-1];

	for (i=0;i<n;i++)
	{
		if (fabs(peak-s[i])<1.0e-3)
		{
			(*position)=i;
		}
	}

	return 0;
}

/*int main (void)
{
	int i; 
	double s[1024];

	double x=0.0;
	for (i=0;i<1024;i++)
	{
		x+=1.0;
		s[i]=x;
	}

	printf ("%lf\n", find_peak(1024,s));

	return 0;
}*/
