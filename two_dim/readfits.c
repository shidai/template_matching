// read PSRFITS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"
#include "get_toa.h"

long int stt_imjd ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	long int imjd;
    if ( fits_read_key(fptr, TLONG, (char *)"STT_IMJD", &imjd, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_IMJD: %ld\n", imjd);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return imjd;
}

long int stt_smjd ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	long int smjd;
    if ( fits_read_key(fptr, TLONG, (char *)"STT_SMJD", &smjd, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_SMJD: %ld\n", smjd);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return smjd;
}

double stt_offs ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    double stt_offs;
    if ( fits_read_key(fptr, TDOUBLE, (char *)"STT_OFFS", &stt_offs, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("STT_OFFS: %lf\n", stt_offs);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return stt_offs;
}

int get_nchan ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	int nchan;
    if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("number of nchan: %d\n", nchan);
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nchan;
}

int get_npol ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
	int npol;
    if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    printf ("number of npol: %d\n", npol);

	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return npol;
}

int get_subint ( char *name )
{  
//double *read_arrival_time( char *input, long *nrows )
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    printf ("number of subint: %ld\n", nrows);
    
	///////////////////////////////////////////////////////////////////////////

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return nrows;
}

int print_t2pred ( char *name )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "PREDICT", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = 1;
    int	anynull = 0;
    char nval[]="NULL";

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	// print predictor to t2pred.dat
	FILE *fp;
	char filename[]="t2pred.dat";
		    
	if ((fp = fopen(filename, "w+")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	int i;
    for (i = 1; i <= nrows; i++)                             // print the results
	{
		frow = i;

		fits_read_col(fptr, TSTRING, colnum, frow, felem, nelem, nval, line, &anynull, &status);           // read the column
	//fits_read_col_str(fptr,colnum,frow,felem,nelem,nval,line,&anynull,&status);

		//puts(line[0]);
        //printf("%s\n", line[0]);
		fprintf (fp, "%s\n", line[0]);
	}

    if (fclose (fp) != 0)
	    fprintf (stderr, "Error closing\n");


    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return 0;
}

double read_freq ( char *name, int subint )
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    if ( fits_get_colnum(fptr, CASEINSEN, "DAT_FREQ", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = 1;
    int null = 0;
    int	anynull = 0;

	double freq;

	//int i;
    //for (i = 1; i <= nrows; i++)                             // print the results
	//{
		frow = subint;

        fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &freq, &anynull, &status);           // read the column

		//puts(line[0]);
        //printf("%lf\n", freq);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return freq;
}

double read_offs ( char *name, int subint)
//int main (int argc, char *argv[] )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    if ( fits_get_colnum(fptr, CASEINSEN, "OFFS_SUB", &colnum, &status) )           // get the colnum number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    int frow;
    int	felem = 1;
    int nelem = 1;
    int null = 0;
    int	anynull = 0;

    double offs;

	//int i;
    //for (i = 1; i <= nrows; i++)                             // print the results
	//{
		frow = subint;

        fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, &offs, &anynull, &status);           // read the column

		//puts(line[0]);
        //printf("subint offset is %Lf\n", offs);
		//fprintf (fp, "%s\n", line[0]);
	//}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file \n" );
    }

    return offs;
}

//int main ( int argc, char *argv[] )
int read_prof ( char *name, int subint, double *profile )
{  
//double *read_arrival_time( char *input, long *nrows )
    //int subint = 1;
	//double profile[1024];
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

    if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    {
        printf( "error while getting the row number\n" );
    }
    //printf ("%ld\n", nrows);
    
    //if ( fits_get_colnum(fptr, CASEINSEN, "TSUBINT", &colnum, &status) )           // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

	//////////////////////////////////////////////////////////////////////////
	int npol;
    if ( fits_read_key(fptr, TINT, (char *)"NPOL", &npol, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", npol);

	int nchan;
    if ( fits_read_key(fptr, TINT, (char *)"NCHAN", &nchan, NULL, &status) )           // get the row number
    {
        printf( "error while getting the npol number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", nchan);
	///////////////////////////////////////////////////////////////////////////

	int nbin;
    int frow;
    int felem;
    int nelem;
    int null;
    int anynull;
    //double *profile;     // the array to store the profile   

	nbin = 1024;
    //profile = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value
    frow = subint;
    felem = 1;
    nelem = nbin*nchan*npol;
    //nelem = 1024;
    null = 0;
    anynull = 0;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status);           // read the column

	//int i;
    //for (i = 0; i < (nchan*npol*nbin); i++)                             // print the results
    //    printf("%d %lf \n", i, profile[i]);

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

    return 0;
}

