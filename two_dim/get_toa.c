// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "get_toa.h"

//double *a_s,*a_p,*p_s,*p_p;
//int num;

double get_toa (double s[1024], double p[1024], double psrfreq)
{
    int nphase=1024;
    int nchn=1;

	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	//double t[nphase*nchn],s[nphase*nchn];
	//int n;

	//readfile(name,&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	//double p[nphase*nchn];
	//double SNR=atof(argv[2]);
	//simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	double phi_s[nchn][NP],phi_p[nchn][NP];

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	
	// initial guess of the phase
    int peak_s, peak_p;	

	find_peak(nphase,s,&peak_s);
	find_peak(nphase,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d=peak_p-peak_s;
	step=2.0*3.1415926/10240.0;

	if (d>=nphase/2)
	{
		ini_phase=2.0*3.1415926*(1023-d)/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)*A7(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

    // calculate phase shift, a and b
    double phase,b;
    phase=zbrent(A7, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    b=A9(phase, amp_s, amp_p, phi_s, phi_p, k, nchn);
    //a=A4(b);

		
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
    double errphase, errb;	

	error(phase,b,&errphase,&errb, amp_s, amp_p, phi_s, phi_p, k,nchn);
	printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	// calculate the rms
	double rms;
	cal_rms(phase,b,&rms, amp_s, amp_p, phi_s, phi_p, k,nchn);

	return rms;
}

int get_toa_multi (double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq)
{
    int nphase=1024;

	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	//double t[nphase*nchn],s[nphase*nchn];
	//int n;

	//readfile(name,&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	//double p[nphase*nchn];
	//double SNR=atof(argv[2]);
	//simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	double phi_s[nchn][NP],phi_p[nchn][NP];

	preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	
	// initial guess of the phase
    int peak_s, peak_p;	

	find_peak(nphase,s,&peak_s);
	find_peak(nphase,p,&peak_p);

	int d;
	double step;
	double ini_phase,up_phase,low_phase;

	d=peak_p-peak_s;
	step=2.0*3.1415926/10240.0;

	if (d>=nphase/2)
	{
		ini_phase=2.0*3.1415926*(1023-d)/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*3.1415926*d/1024.0;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

    // calculate phase shift, a and b
    double phase,b;
    phase=zbrent_multi(A7_multi, low_phase, up_phase, 1.0e-16, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
    //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
    //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
    b=A9_multi(phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
    //a=A4(b);

		
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
    double errphase, errb;	

	error_multi(phase,b,&errphase,&errb, amp_s, amp_p, phi_s, phi_p, k, nchn, rms);
	printf ("multi-template\n");
	printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	return 0;
}
