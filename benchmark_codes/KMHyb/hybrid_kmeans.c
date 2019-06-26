/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   File:         seq_kmeans.c  (sequential version)                        */
/*   Description:  Implementation of simple k-means clustering algorithm     */
/*                 This program takes an array of N data objects, each with  */
/*                 M coordinates and performs a k-means clustering given a   */
/*                 user-provided value of the number of clusters (K). The    */
/*                 clustering results are saved in 2 arrays:                 */
/*                 1. a returned array of size [K][N] indicating the center  */
/*                    coordinates of K clusters                              */
/*                 2. membership[N] stores the cluster center ids, each      */
/*                    corresponding to the cluster a data object is assigned */
/*                                                                           */
/*   Author:  Wei-keng Liao                                                  */
/*            ECE Department, Northwestern University                        */
/*            email: wkliao@ece.northwestern.edu                             */
/*   Copyright, 2005, Wei-keng Liao                                          */
/*   See COPYRIGHT notice in top-level directory.                            */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include <stdlib.h>
#include<omp.h>
#include <math.h>
#include <mpi.h>
#include "kmeans.h"

#define MEAN_VARIATION			//remove this line to use "percentage change in membership" as the stopping criteria

/*----< euclid_dist_2() >----------------------------------------------------*/
/* square of Euclid distance between two multi-dimensional points            */
__inline static
double euclid_dist_2(int    numdims,  /* no. dimensions */
                    double *coord1,   /* [numdims] */
                    double *coord2)   /* [numdims] */
{
    int i;
    double ans=0.0;

    for (i=0; i<numdims; i++)
        ans += (coord1[i]-coord2[i]) * (coord1[i]-coord2[i]);

    return(sqrt(ans));
}

/*----< find_nearest_cluster() >---------------------------------------------*/
__inline static
int find_nearest_cluster(int     numClusters, /* no. clusters */
                         int     numCoords,   /* no. coordinates */
                         double  *object,      /* [numCoords] */
                         double **clusters)    /* [numClusters][numCoords] */
{
    int   index, i;
    double dist, min_dist;

    /* find the cluster id that has min distance to object */
    index    = 0;
    min_dist = euclid_dist_2(numCoords, object, clusters[0]);

    for (i=1; i<numClusters; i++) {
        dist = euclid_dist_2(numCoords, object, clusters[i]);
        /* no need square root */
        if (dist < min_dist) { /* find the min and its array index */
            min_dist = dist;
            index    = i;
        }
    }
    return(index);
}

/*----< mpi_kmeans() >-------------------------------------------------------*/
int mpi_kmeans(double    **objects,     /* in: [numObjs][numCoords] */
               int        numCoords,   /* no. coordinates */
               int        numObjs,     /* no. objects */
               int        numClusters, /* no. clusters */
               double      threshold,   /* % objects change membership */
               int       *membership,  /* out: [numObjs] */
               double    **clusters,    /* out: [numClusters][numCoords] */
               MPI_Comm   comm,        /* MPI communicator */
		int     is_perform_atomic)
{
    int      i, j, k, rank, index, loop=0, total_numObjs;
    int     *newClusterSize; /* [numClusters]: no. objects assigned in each
                                new cluster */
    int     *clusterSize;    /* [numClusters]: temp buffer for Allreduce */
    double    delta;          /* % of objects change their clusters */
    double    delta_tmp;
    double  **newClusters;    /* [numClusters][numCoords] */
	double meanVariation, temp;
	double t1, t2;
    extern int _debug;


    int      nthreads;             /* no. threads */
    int    **local_newClusterSize; /* [nthreads][numClusters] */
    double ***local_newClusters;    /* [nthreads][numClusters][numCoords] */

//int     is_perform_atomic=1;
    //nthreads = omp_get_max_threads();
    nthreads = 4;//omp_get_num_threads();
	//printf("threads = %d\n", nthreads);
   // if (_debug)
MPI_Comm_rank(comm, &rank);

    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    /* need to initialize newClusterSize and newClusters[0] to all 0 */
    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);
    clusterSize    = (int*) calloc(numClusters, sizeof(int));
    assert(clusterSize != NULL);

    newClusters    = (double**) malloc(numClusters *            sizeof(double*));
    assert(newClusters != NULL);
    newClusters[0] = (double*)  calloc(numClusters * numCoords, sizeof(double));
    assert(newClusters[0] != NULL);
    for (i=1; i<numClusters; i++)
        newClusters[i] = newClusters[i-1] + numCoords;

    MPI_Allreduce(&numObjs, &total_numObjs, 1, MPI_INT, MPI_SUM, comm);
    if (_debug) printf("%2d: numObjs=%d total_numObjs=%d numClusters=%d numCoords=%d\n",rank,numObjs,total_numObjs,numClusters,numCoords);
	if (!is_perform_atomic)
	{
        /* each thread calculates new centers using a private space,
           then thread 0 does an array reduction on them. This approach
           should be faster */
        local_newClusterSize    = (int**) malloc(nthreads * sizeof(int*));
        assert(local_newClusterSize != NULL);
        local_newClusterSize[0] = (int*)  calloc(nthreads*numClusters, sizeof(int));
        assert(local_newClusterSize[0] != NULL);
        for (i=1; i<nthreads; i++)
            local_newClusterSize[i] = local_newClusterSize[i-1]+numClusters;

		for(i=0; i< nthreads; i++)
			for (j=0; j < numClusters; j++)
		        local_newClusterSize[i][j] = 0;

        /* local_newClusters is a 3D array */
        local_newClusters    =(double***)malloc(nthreads * sizeof(double**));
        assert(local_newClusters != NULL);
        local_newClusters[0] =(double**) malloc(nthreads * numClusters * sizeof(double*));
        assert(local_newClusters[0] != NULL);
        for (i=1; i<nthreads; i++)
            local_newClusters[i] = local_newClusters[i-1] + numClusters;
        for (i=0; i<nthreads; i++)
        {
            for (j=0; j<numClusters; j++)
            {
                local_newClusters[i][j] = (double*)calloc(numCoords, sizeof(double));
                assert(local_newClusters[i][j] != NULL);
            }
        }
		for(i=0; i< nthreads; i++)
			for (j=0; j < numClusters; j++)
				for (k = 0; k < numCoords; k++)
			        local_newClusters[i][j][k] = 0.0;
    }   
double sum_time = 0;
    do {
        t1 = MPI_Wtime();
        delta = 0.0;
       if (is_perform_atomic)
		{
	            #pragma omp parallel for \
                    private(i,j,index) \
                    firstprivate(numObjs,numClusters,numCoords) \
                    shared(objects,clusters,membership,newClusters,newClusterSize) \
                    schedule(static) \
                    reduction(+:delta)
        for (i=0; i<numObjs; i++) {
            /* find the array index of nestest cluster center */
            index = find_nearest_cluster(numClusters, numCoords, objects[i],
                                         clusters);

            /* if membership changes, increase delta by 1 */
            if (membership[i] != index) delta += 1.0;

            /* assign the membership to object i */
            membership[i] = index;

             /* update new cluster centers : sum of objects located within */
		        #pragma omp atomic
			       	newClusterSize[index]++;

		        for (j=0; j<numCoords; j++)
				{
		        	#pragma omp atomic
					newClusters[index][j] += objects[i][j];
				}
			}
        }
 else
		{
			#pragma omp parallel shared(objects,clusters,membership,local_newClusters,local_newClusterSize)
            {
               	int tid = omp_get_thread_num();
               	#pragma omp for \
		        private(i,j,index) \
	            firstprivate(numObjs,numClusters,numCoords) \
	            schedule(static) \
	            reduction(+:delta)
	    	    for (i=0; i<numObjs; i++)
				{
		            /* find the array index of nestest cluster center */
		            index = find_nearest_cluster(numClusters, numCoords, objects[i], clusters);

		            /* if membership changes, increase delta by 1 */
		            if (membership[i] != index) delta += 1.0;

		            /* assign the membership to object i */
		            membership[i] = index;

		            /* update new cluster centers : sum of all objects located
		               within (average will be performed later) */
		            local_newClusterSize[tid][index]++;
		            for (j=0; j<numCoords; j++)
		                local_newClusters[tid][index][j] += objects[i][j];
				}
			} /* end of #pragma omp parallel */

		    /* let the main thread perform the array reduction */
		    for (i=0; i<numClusters; i++)
			{
		        for (j=0; j<nthreads; j++)
				{
		            newClusterSize[i] += local_newClusterSize[j][i];
		            local_newClusterSize[j][i] = 0.0;
		            for (k=0; k<numCoords; k++)
					{
		                newClusters[i][k] += local_newClusters[j][i][k];
		                local_newClusters[j][i][k] = 0.0;
		            }
		        }
		    }
        }

        /* sum all data objects in newClusters */

        MPI_Allreduce(MPI_IN_PLACE, newClusters[0], numClusters*numCoords,
                      MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(MPI_IN_PLACE, newClusterSize, numClusters, MPI_INT,
                      MPI_SUM, comm);

/*        MPI_Allreduce(newClusters[0], clusters[0], numClusters*numCoords,
                      MPI_DOUBLE, MPI_SUM, comm);
        MPI_Allreduce(newClusterSize, clusterSize, numClusters, MPI_INT,
                      MPI_SUM, comm);
*/
//	MPI_Barrier(comm);
#ifdef MEAN_VARIATION
        for (i=0; i<numClusters; i++)
        {
        	temp = 0.0f;
            for (j=0; j<numCoords; j++)
            {
               if (newClusterSize[i] > 0)
                    newClusters[i][j] = newClusters[i][j] / newClusterSize[i];
                temp += ( (clusters[i][j] - newClusters[i][j]) * (clusters[i][j] - newClusters[i][j]));
                clusters[i][j] = newClusters[i][j];
                newClusters[i][j] = 0.0;   //set back to 0
            }
            newClusterSize[i] = 0;   // set back to 0
	    meanVariation += sqrt(temp);
        }

        //delta = meanVariation / numObjs;
	delta = meanVariation / numClusters; 
	//if(rank == 0) printf("\nDelta ***: %lf\tK: %d\tD: %d", delta, numClusters, numCoords);
#else

        /* average the sum and replace old cluster centers with newClusters */
        for (i=0; i<numClusters; i++) {
            for (j=0; j<numCoords; j++) {
                if (clusterSize[i] > 1)
                    clusters[i][j] /= clusterSize[i];
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }
            
        MPI_Allreduce(&delta, &delta_tmp, 1, MPI_DOUBLE, MPI_SUM, comm);
        delta = delta_tmp / total_numObjs;
	if(rank == 0) printf("\nDelta: %lf", delta);
#endif
        meanVariation = 0.0;
		t2 = MPI_Wtime();
		if(rank == 0)
		printf("\niteration:%d\terror : %lf\tTime: %lf\n", loop+1, delta, t2-t1);
		sum_time += (t2 - t1);
		/*if (_debug) {
            double maxTime;
            curT = MPI_Wtime() - curT;
            MPI_Reduce(&curT, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
            if (rank == 0) printf("%2d: loop=%d time=%f sec\n",rank,loop,curT);
        }*/
    } while (delta > threshold && loop++ < 19);

    if (_debug && rank == 0) printf("%2d: delta=%f threshold=%f loop=%d\n",rank,delta,threshold,loop);

    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);
    free(clusterSize);

    return 1;
}

