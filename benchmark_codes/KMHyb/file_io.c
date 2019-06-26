/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         file_io.c                                                 */
/*   Description:  This program reads point data from a file                 */
/*                 and write cluster output to files                         */
/*   Input file format:                                                      */
/*            First number represents Number of Objects (numObjs)            */
/*            Second number represents Number of Dimensions (numCoords)	     */
/*            After that each numCoords numbers represent a data point       */
/*                 coordinates) of each object                               */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                     			     */
/*   Modified by: Saiyedul Islam					     */
/*		  CSIS Department BITS Pilani				     */
/*		  saiyedul.islam@gmail.com		                     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* read(), close() */

#include "kmeans.h"

#define MAX_CHAR_PER_LINE 128


/*---< file_read() >---------------------------------------------------------*/
double** file_read(int   isBinaryFile,  /* flag: 0 or 1 */
                  char *filename,      /* input file name */
                  int  *numObjs,       /* no. data objects (local) */
                  int  *numCoords)     /* no. coordinates */
{
			//Input File Format ==>
				//First number represents Number of Objects (numObjs)	
				//Second number represents Number of Dimensions (numCoords)
				//After that each numCoords numbers represent a data point.

    double **objects;
    int     i, j, len;
    ssize_t numBytesRead;
//	printf("*** %d\n", isBinaryFile);
    if (isBinaryFile)
	{  /* input file is in raw binary format -------------*/
        int infile;
        if ((infile = open(filename, O_RDONLY, "0600")) == -1)
		{
            fprintf(stderr, "Error: no such file (%s)\n", filename);
            return NULL;
        }
        numBytesRead = read(infile, numObjs,    sizeof(int));
        assert(numBytesRead == sizeof(int));
        numBytesRead = read(infile, numCoords, sizeof(int));
        assert(numBytesRead == sizeof(int));
        if (_debug)
		{
            printf("File %s numObjs   = %d\n",filename,*numObjs);
            printf("File %s numCoords = %d\n",filename,*numCoords);
        }

        /* allocate space for objects[][] and read all objects */
        len = (*numObjs) * (*numCoords);
        objects    = (double**)malloc((*numObjs) * sizeof(double*));
        assert(objects != NULL);
        objects[0] = (double*) malloc(len * sizeof(double));
        assert(objects[0] != NULL);
        for (i=1; i<(*numObjs); i++)
            objects[i] = objects[i-1] + (*numCoords);

        numBytesRead = read(infile, objects[0], len*sizeof(double));
        assert(numBytesRead == len*sizeof(double));

        close(infile);
    }
	else 
	{		// input file is in following ASCII format
			
		FILE *infile;

		if ((infile = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr, "Error: no such file (%s)\n", filename);
			return NULL;
		}
		
		fscanf(infile, "%d", numObjs);
		fscanf(infile, "%d", numCoords);

		if (_debug) {
		    printf("File %s numObjs   = %d\n",filename,*numObjs);
		    printf("File %s numCoords = %d\n",filename,*numCoords);
		}

		/* allocate space for objects[][] and read all objects */
		len = (*numObjs) * (*numCoords);
		objects    = (double**)malloc((*numObjs) * sizeof(double*));
		assert(objects != NULL);
		objects[0] = (double*) malloc(len * sizeof(double));
		assert(objects[0] != NULL);
		for (i=1; i<(*numObjs); i++)
		    objects[i] = objects[i-1] + (*numCoords);

		i = 0;
		while(feof(infile)!= 1)
		{
			for (j=0; j<(*numCoords); j++)
		        {
				fscanf(infile, "%lf", &objects[i][j]);
				//printf("%lf ", objects[i][j]);
			}
			i++;
			//printf("\n");
		}
			
		fclose(infile);
		
		if((--i) != *numObjs)
		{
			fprintf(stderr, "Error: Broken file (%s) : Have lesser number of objects than claimed\n", filename);
			return NULL;
		}	
	}
    return objects;
}

/*---< file_write() >---------------------------------------------------------*/
int file_write(char      *filename,     /* input file name */
               int        numClusters,  /* no. clusters */
               int        numObjs,      /* no. data objects */
               int        numCoords,    /* no. coordinates (local) */
               double    **clusters,     /* [numClusters][numCoords] centers */
               int       *membership)   /* [numObjs] */
{
    FILE *fptr;
    int   i, j;
    char  outFileName[1024];

    /* output: the coordinates of the cluster centres ----------------------*/
    sprintf(outFileName, "%s.cluster_centres", filename);
    //printf("Writing coordinates of K=%d cluster centers to file \"%s\"\n", numClusters, outFileName);
    fptr = fopen(outFileName, "w");
    for (i=0; i<numClusters; i++) {
        fprintf(fptr, "%d ", i);
        for (j=0; j<numCoords; j++)
            fprintf(fptr, "%f ", clusters[i][j]);
        fprintf(fptr, "\n");
    }
    fclose(fptr);

    /* output: the closest cluster centre to each of the data points --------*/
    sprintf(outFileName, "%s.membership", filename);
    //printf("Writing membership of N=%d data objects to file \"%s\"\n", numObjs, outFileName);
    fptr = fopen(outFileName, "w");
    for (i=0; i<numObjs; i++)
        fprintf(fptr, "%d %d\n", i, membership[i]);
    fclose(fptr);

    return 1;
}
