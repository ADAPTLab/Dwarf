/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         mpi_main.c   (an MPI version)                             */
/*   Description:  This program shows an example on how to call a subroutine */
/*                 that implements a simple k-means clustering algorithm     */
/*                 based on Euclid distance.                                 */
/*   Input file format:                                                      */
/*                 ascii  file: each line contains 1 data object             */
/*                 binary file: first 4-byte integer is the number of data   */
/*                 objects and 2nd integer is the no. of features (or        */
/*                 coordinates) of each object                               */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department Northwestern University                         */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     /* strtok() */
#include <sys/types.h>  /* open() */
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>     /* getopt() */
#include <omp.h>
#include <mpi.h>
int      _debug;
#include "kmeans.h"

int     mpi_kmeans(double**, int, int, int, double, int*, double**, MPI_Comm, int);
double** mpi_read(int, char*, int*, int*, MPI_Comm);
int     mpi_write(int, char*, int, int, int, double**, int*, int, MPI_Comm);



/*---< usage() >------------------------------------------------------------*/
static void usage(char *argv0, double threshold) {
    char *help =
        "Usage: %s [switches] -i filename -n num_clusters\n"
        "       -i filename    : file containing data to be clustered\n"
        "       -b             : use MPI parallel I/O (default no)\n"
        "       -r             : output file in binary format (default no)\n"
        "       -n num_clusters: number of clusters (K must > 1)\n"
        "       -p nproc       : number of threads (default system allocated)\n"
        "       -t threshold   : threshold value (default %.4f)\n"
        "       -o             : output timing results (default no)\n"
        "       -d             : dimensions\n";
    fprintf(stderr, help, argv0, threshold);
}

/*---< main() >-------------------------------------------------------------*/
int main(int argc, char **argv) {
           int     opt;
    extern char   *optarg;
    extern int     optind;
           int     i, j;
           int     isInFileBinary,is_perform_atomic, isOutFileBinary;
           int     is_output_timing, is_print_usage;
int nthreads;
           int     numClusters, numCoords, numObjs, totalNumObjs;
           int    *membership;    /* [numObjs] */
           char   *filename;
           double **objects;       /* [numObjs][numCoords] data objects */
           double **clusters;      /* [numClusters][numCoords] cluster center */
           double   threshold;
           double  timing, io_timing, clustering_timing;

           int        rank, nproc, mpi_namelen;
           char       mpi_name[MPI_MAX_PROCESSOR_NAME];
           MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Get_processor_name(mpi_name,&mpi_namelen);

    /* some default values */
    _debug           = 0;
    threshold        = 0.001;
    numClusters      = 0;
    isInFileBinary   = 0;
    isOutFileBinary  = 0;
    is_output_timing = 1;
    is_print_usage   = 0;
    is_perform_atomic = 0;
    filename         = NULL;
	//nthreads=4;

    while ( (opt=getopt(argc,argv,"p:i:n:d:t:aborh"))!= EOF) {
        switch (opt) {
            case 'i': filename=optarg;
                      break;
            case 'b': isInFileBinary = 1;
                      break;
            case 'r': isOutFileBinary = 1;
                      break;

            case 't': threshold=atof(optarg);
                      break;
            case 'n': numClusters = atoi(optarg);
                      break;
            case 'p': nthreads = atoi(optarg);
                      break;
            case 'o': is_output_timing = 1;
                      break;
            case 'd': //_debug = 1;
			numCoords = atoi(optarg);
                      break;
            case 'h': is_print_usage = 1;
                      break;
            default: is_print_usage = 1;
                      break;
        }
    }

    if (filename == 0 || numClusters <= 1 || is_print_usage == 1) {
        if (rank == 0) usage(argv[0], threshold);
        MPI_Finalize();
        exit(1);
    }
    
    if(_debug)
	printf("Proc %d of %d running on %s\n", rank, nproc, mpi_name);

    int err = MPI_Barrier(MPI_COMM_WORLD);
	if (err != MPI_SUCCESS)
	{
        char errstr[MPI_MAX_ERROR_STRING];
        int  errlen;
        MPI_Error_string(err, errstr, &errlen);
        printf("Error at Barrier (%s)\n",errstr);
        MPI_Finalize();
        exit(1);
    }
    //io_timing = MPI_Wtime();

    /* read data points from file ------------------------------------------*/
    objects = mpi_read(isInFileBinary, filename, &numObjs, &numCoords,
                       MPI_COMM_WORLD);

	/*timing            = MPI_Wtime();
    io_timing         = timing - io_timing;
    clustering_timing = timing;*/

    if (_debug) { /* print the first 4 objects' coordinates */
        int num = (numObjs < 4) ? numObjs : 4;
        for (i=0; i<num; i++) {
            char strline[1024], strdouble[16];
            sprintf(strline,"%d: objects[%d]= ",rank,i);
            for (j=0; j<numCoords; j++) {
                sprintf(strdouble,"%10f",objects[i][j]);
                strcat(strline, strdouble);
            }
            strcat(strline, "\n");
            printf("%s",strline);
        }
    }

    /*timing            = MPI_Wtime();
    io_timing         = timing - io_timing;
    clustering_timing = timing;*/

    /* allocate a 2D space for clusters[] (coordinates of cluster centers)
       this array should be the same across all processes                  */
    clusters    = (double**) malloc(numClusters *             sizeof(double*));
    assert(clusters != NULL);
    clusters[0] = (double*)  malloc(numClusters * numCoords * sizeof(double));
    assert(clusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        clusters[i] = clusters[i-1] + numCoords;

    MPI_Allreduce(&numObjs, &totalNumObjs, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /* pick first numClusters elements in feature[] as initial cluster centers*/
    if (rank == 0) {
        for (i=0; i<numClusters; i++)
            for (j=0; j<numCoords; j++)
                clusters[i][j] = objects[i][j];
    }
    MPI_Bcast(clusters[0], numClusters*numCoords, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* membership: the cluster id for each data object */
    membership = (int*) malloc(numObjs * sizeof(int));
    assert(membership != NULL);
	if (is_output_timing) {
        timing            = wtime();
        io_timing         = timing - io_timing;
        clustering_timing = timing;
    }
    /* start the core computation -------------------------------------------*/
    mpi_kmeans(objects, numCoords, numObjs, numClusters, threshold, membership,
               clusters, MPI_COMM_WORLD, is_perform_atomic);

   // free(objects[0]);
   // free(objects);
	if (is_output_timing) {
        timing            = wtime();
        clustering_timing = timing - clustering_timing;
    }

    /* output: the coordinates of the cluster centres ----------------------*/
//    mpi_write(isOutFileBinary, filename, numClusters, numObjs, numCoords,
  //           clusters, membership, totalNumObjs, MPI_COMM_WORLD);
	
	//io_timing += MPI_Wtime() - timing;
    free(membership);
    free(clusters[0]);
    free(clusters);

    /*---- output performance numbers ---------------------------------------*/
    if (is_output_timing) {
        double max_io_timing, max_clustering_timing;

        io_timing += MPI_Wtime() - timing;

        /* get the max timing measured among all processes */
  /*      MPI_Reduce(&io_timing, &max_io_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&clustering_timing, &max_clustering_timing, 1, MPI_DOUBLE,
                   MPI_MAX, 0, MPI_COMM_WORLD);
*/
        if (rank == 0) {
            printf("\nPerforming **** Simple Kmeans  (MPI) ****\n");
            printf("Num of processes = %d\n", nproc);
            printf("Input file:        %s\n", filename);
            printf("numObjs          = %d\n", totalNumObjs);
            printf("numCoords        = %d\n", numCoords);
            printf("numClusters      = %d\n", numClusters);
            printf("threshold        = %.4f\n", threshold);

           // printf("I/O time           = %10.4f sec\n", io_timing);
            printf("Computation timing = %10.4f sec\n", clustering_timing);
        }
    }

    MPI_Finalize();
    return(0);
}

