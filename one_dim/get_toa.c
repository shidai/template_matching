// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "get_toa.h"

double *a_s,*a_p,*p_s,*p_p;
int num;

int main (int argc, char *argv[])
{
	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	double t[512],s[512];
	int n;

	readfile(argv[1],&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	double p[512];
	double SNR=atof(argv[2]);
	simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	int k;

	double amp_s[n/2],amp_p[n/2];  // elements for calculating A7
	double phi_s[n/2],phi_p[n/2];

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, n);
	
	a_s=amp_s;
	a_p=amp_p;
	p_s=phi_s;
	p_p=phi_p;
	num=k;    // k=n/2

	// initial guess of the phase
    int peak_s, peak_p;	

	find_peak(n,s,&peak_s);
	find_peak(n,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d=peak_p-peak_s;
	step=2.0*3.1415926/5120.0;

	if (d>=n/2)
	{
		ini_phase=2.0*3.1415926*(511-d)/512.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase)*A7(low_phase)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/512.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase)*A7(low_phase)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

    // calculate phase shift, a and b
    double phase,b;
    phase=zbrent(A7, low_phase, up_phase, 1.0e-16);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    b=A9(phase);
    //a=A4(b);

		
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
    double errphase, errb;	

	error(phase,b,&errphase,&errb);
	printf ("%.10lf %.10lf\n", ((phase/3.1415926)*5.75/2.0)*1.0e+3, ((errphase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);

	return 0;
}
