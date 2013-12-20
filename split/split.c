// read PSRFITS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fitsio.h"

int create (char *name_old, char *name_new );
int read_prof ( char *name, int nchan, double *p );
int write_prof ( char *name, double *p );
int print_prof ( char *name );

int main (int argc, char *argv[])
{
	// argv[1]: original psrfits file
	// argv[2]: new psrfits file
	// argv[3]: the profile of which frequency channel to be written in
	int nchan = atoi(argv[3]);

	// create new file for each profile in different channel
	create(argv[1], argv[2]);

	// read the profile of No. nchan channel
	double p[1024];

	read_prof(argv[1], nchan, p);

        /*
	int i;
	for (i = 0; i < 1024; i++)
	{
		printf ("%d %lf\n", i, p[i]);
	}
	*/

	// write the profile into the new file
	write_prof(argv[2], p);
	
	// print result
	print_prof(argv[2]);
	
	return 0;
}

int create (char *name_old, char *name_new )
// copy the original psrfits file to a new file, and delete the DATA colnum in SUBINT extension, and insert a new empty DATA colnum
{
	// open old file
    fitsfile *fptr_old;       // pointer to the FITS file, defined in fitsio.h 
    int status;

    status = 0;

    if ( fits_open_file(&fptr_old, name_old, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	//////////////////////////////////////////////////////////////////////////
    fitsfile *fptr_new;       // pointer to the FITS file, defined in fitsio.h 

	// create a new file
	if ( fits_create_file(&fptr_new, name_new, &status) )
    {
        printf( "error while creating file\n" );
    }

    fitsfile *fptr_e;       // pointer to the FITS file, defined in fitsio.h 

	// create an empty file
	if ( fits_create_file(&fptr_e, "!empty(psrheader)", &status) )
    {
        printf( "error while creating file\n" );
    }

	// copy old file to new file
	if ( fits_copy_file(fptr_old, fptr_new, 0, 1, 1, &status) )
    {
        printf( "error while copying file\n" );
    }

	///////////////////////////////////////////////////////////////////////////
	
	// move the SUBINT extension
	if ( fits_movabs_hdu(fptr_new, 6, NULL, &status) )
	{
        printf( "error while moving in file\n" );
	}

	// move the SUBINT extension
	if ( fits_movabs_hdu(fptr_e, 6, NULL, &status) )
	{
        printf( "error while moving in file\n" );
	}
	
	int colnum;
	if ( fits_get_colnum(fptr_new, CASEINSEN, "DATA", &colnum, &status) )
	{
	    printf( "error while getting the colnum number\n" );
	}

	//printf ("colnum number %d\n", colnum);

        /*
	// delete the DATA colnum
	if ( fits_delete_col(fptr_new, colnum, &status) )
	{
        printf( "error while deleting colnums in file\n" );
	}

	// insert a new empty DATA colnum
	char name[] = "DATA";
	char form[] = "1024I";
	if ( fits_insert_col(fptr_new, colnum, name, form, &status) )
	{
        printf( "error while inserting colnums in file\n" );
	}
	*/

	// copy the new empty DATA colnum
	if ( fits_copy_col(fptr_e, fptr_new, colnum, colnum, 0, &status) )
	{
        printf( "error while copying colnums in file\n" );
	}

	// close files
    if ( fits_close_file(fptr_old, &status) )
    {
        printf( " error while closing the file\n " );
    }

    if ( fits_close_file(fptr_new, &status) )
    {
        printf( " error while closing the file\n" );
    }

    if ( fits_close_file(fptr_e, &status) )
    {
        printf( " error while closing the file\n" );
    }

	return 0;
}

int read_prof ( char *name, int n, double *p )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    //long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	// move the SUBINT extension
	if ( fits_movabs_hdu(fptr, 6, NULL, &status) )
	{
        printf( "error while moving in file\n" );
	}
	
    //if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    //{
    //    printf( "error while getting the row number\n" );
    //}
    //printf ("%ld\n", nrows);
    
    // get the row number
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
    double *profile;     // the array to store the profile   

	nbin = 1024;
    profile = ( double *)malloc( (nchan*npol*nbin) * sizeof( double ) );               // allocate space for column value
    frow = 1;
    felem = 1;
    nelem = nbin*nchan*npol;
    //nelem = 1024;
    null = 0;
    anynull = 0;

    fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status);           // read the column

	int i;
    for (i = 0; i < nbin; i++)                             // print the results
	{
		p[i] = profile[nbin*(n - 1) + i];
        //printf("%d %lf \n", i, p[i]);
	}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

	//free(profile);

    return 0;
}

int write_prof ( char *name, double *p )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    //long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READWRITE, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	// move the SUBINT extension
	if ( fits_movabs_hdu(fptr, 6, NULL, &status) )
	{
        printf( "error while moving in file\n" );
	}
	
    // get the row number
    //if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    //{
    //    printf( "error while getting the row number\n" );
    //}
    //printf ("%ld\n", nrows);
    
    if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

    /*
	int i;
	for (i = 0; i < 1024; i++)
	{
		printf ("%d %lf\n", i, p[i]);
	}
	*/

	int nbin;
    int frow;
    int felem;
    int nelem;

	nbin = 1024;
    frow = 1;
    felem = 1;
    //nelem = 1;
    nelem = 1024;

    if ( fits_write_col(fptr, TDOUBLE, colnum, frow, felem, nelem, p, &status) )           // read the column
    {
        printf( " error while writing the file " );
    }

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file\n" );
    }

    return 0;
}

int print_prof ( char *name )
{  
    fitsfile *fptr;       // pointer to the FITS file, defined in fitsio.h 
    int status;
    int colnum;
    //long int nrows;

    status = 0;

    //if ( fits_open_file(&fptr, argv[1], READONLY, &status) )          // open the file
    if ( fits_open_file(&fptr, name, READONLY, &status) )          // open the file
    {
        printf( "error while openning file\n" );
    }

	// move the SUBINT extension
	if ( fits_movabs_hdu(fptr, 6, NULL, &status) )
	{
        printf( "error while moving in file\n" );
	}
	
    //if ( fits_get_num_rows(fptr, &nrows, &status) )           // get the row number
    //{
    //    printf( "error while getting the row number\n" );
    //}
    //printf ("%ld\n", nrows);
    
    // get the row number
    if ( fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status) )           // get the row number
    {
        printf( "error while getting the colnum number\n" );
		//fits_get_colnum(fptr, CASEINSEN, "DATA", &colnum, &status);
	}
    //printf ("%d\n", colnum);

	///////////////////////////////////////////////////////////////////////////

	int nbin;
    int frow;
    int felem;
    int nelem;
    int null;
    int anynull;
    double *profile;     // the array to store the profile   

	nbin = 1024;
    profile = ( double *)malloc( nbin * sizeof( double ) );               // allocate space for column value
    frow = 1;
    felem = 1;
    nelem = nbin;
    //nelem = 1024;
    null = 0;
    anynull = 0;

    if ( fits_read_col(fptr, TDOUBLE, colnum, frow, felem, nelem, &null, profile, &anynull, &status) )           // read the column
    {
        printf( " error while reading the file " );
    }

	int i;
    for (i = 0; i < nbin; i++)                             // print the results
	{
        printf("%d %lf \n", i, profile[i]);
	}

    if ( fits_close_file(fptr, &status) )
    {
        printf( " error while closing the file " );
    }

	free(profile);

    return 0;
}

