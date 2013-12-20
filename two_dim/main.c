// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "get_toa.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
#include "simulatePseudoBB.h"


//double *a_s,*a_p,*p_s,*p_p;
//int num;

int main (int argc, char *argv[])
{
	// name of different extension
	char name_data[30]; 
	char name_predict[30]; 

	strcpy(name_data,argv[1]);
	strcpy(name_predict,argv[1]);

	char data[] = "[SUBINT]";
	char predict[] = "[T2PREDICT]";

	strcat(name_data, data);
	strcat(name_predict, predict);

	//puts(name_data);
	//puts(name_predict);
	
	////////////////////////////////////////////////
	long int imjd, smjd;
	double offs;
    int nphase=1024;
	int nchn;
	int nsub;
	int npol;
	
	imjd = stt_imjd(argv[1]);
	smjd = stt_smjd(argv[1]);
	offs = stt_offs(argv[1]);

    nchn = get_nchan(name_data);	
    npol = get_npol(name_data);	
    nsub = get_subint(name_data);	

	//printf ("%d\n", nchn);
	////////////////////////////////////////////////

	// read a std
	char std[30];
	strcpy(std,argv[2]);
	strcat(std, data);
	//puts(argv[1]);
	//puts(argv[2]);
	double s[nphase];
	//double tt[nphase];
	//int n;

	//readfile(argv[1],&n,tt,s);
	read_prof(std,1,s);

	/*
	int i;
	for (i = 0; i < nphase; i++)
	{
	    printf ("%lf\n", s[i]);
	}
	//puts(argv[1]);
	*/
	/////////////////////////////////////////////////////////////////////////////////
	FILE *fp;
	if ((fp = fopen("test.txt", "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	
	fprintf (fp, "FORMAT 1\n");

	////////////////////////////////////////////////////////////////////////////////
	double s_multi[nphase*nchn*npol];

	double p_multi[nchn*npol*nphase];
	double p_temp[nphase];
    //double SNR; 

	double rms[nchn];  // rms for each profile
	int h,i,j,z;
	double phase, e_phase;
	long double dt, e_dt;  
	long double t;     // TOA
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period, frequency, weight;
	double freq[nchn], wts[nchn];
	for (h = 1; h <= nsub; h++)
	{
	    //////////////////////////////////////////////////////////////////////////
	    // simulate data

		//SNR = 500.0 + 200.0*i;
	    //simulate(n,SNR,s,p_temp);

	    read_prof(name_data,h,p_multi);
	    //readfile(argv[2],&n,tt,p_multi);

		for (i = 0; i < nchn; i++)
		{
			for (j = 0; j < nphase; j++)
			{
				//printf ("%lf %lf\n", p_multi[j], s[j]);
				s_multi[i*nphase + j] = s[j];
				p_temp[j] = p_multi[i*nphase + j];
			}

			// calculate toa, rms for each profile
			rms[i] = get_toa(s, p_temp);
		}

		// do template matching, get the phase shift
		get_toa_multi(s_multi, p_multi, rms, nchn, &phase, &e_phase);

		////////////////////////////////////////////////////////////////////////////////////////

		// transform phase shift to TOAs
		// get the freq of the subint
		read_freq(name_data, h, freq, nchn);
		read_wts(name_data, h, wts, nchn);
		frequency = 0.0;
		weight = 0.0;
		for (z = 0; z < nchn; z++)
		{
			frequency += freq[z]*wts[z];
			weight += wts[z];
		}
		frequency = frequency/weight;
	    printf ("Frequency is %lf\n", frequency);

		// get the period
        print_t2pred(name_predict);   // output t2pred.dat

		T2Predictor_Init(&pred);  // prepare the predictor
		if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
	    {
			printf("Error: unable to read predictor\n");
			exit(1);
		}

		// get the offset of each subint
		offset = read_offs(name_data, h);

		// get the period at mjd0
		mjd0 = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) + (long double)(offset))/86400.0L;
		printf ("imjd is %ld \n", imjd);
		printf ("mjd0 is %.15Lf \n", mjd0);

		period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
	    printf ("Period is %.15lf\n", period);
	
		// transform phase shift to time shift
        //dt = (phase/PI)*period/2.0;
        //e_dt = (e_phase/PI)*period/2.0;
        dt = ((long double)(phase)/PI)*((long double)(period))/2.0L;
        e_dt = ((long double)(e_phase)/PI)*((long double)(period))/2.0L;
	    printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

		// calculate the TOA
        t = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
        //t = imjd;
		
	    printf ("offset is %lf\n", offset);
		fprintf (fp, "1713.pF  %lf  %.15Lf  %Lf  7\n", frequency, t, e_dt*1e+6);
	}

    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");


	return 0;
}
