// read txt file
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "get_toa.h"

int readfile ( char *filename, int *ntxt, double *x, double *y )
{
    FILE *fp;
	//double n[10];
	//double m[10];
	//char c[10000][100];
	//char d[10000][100];
	int i;
	if ((fp = fopen(filename, "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}
	i = 0;
	// printf ("%d\n", i);
	while (fscanf (fp, "%lf %lf", &x[i], &y[i]) == 2)
		{
			//printf ("%f %f\n", x[i], y[i]);
			i++;
			(*ntxt) = i;
		}
	//printf ("number of ogle %ld\n", *ntxt);
	if (fclose (fp) != 0)	
	    fprintf (stderr, "Error closing\n");
	
	return 0;
}
    	
   
