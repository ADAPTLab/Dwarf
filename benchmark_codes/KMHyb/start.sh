#!/bin/sh

input=$1
data=../serialPublic

cluster=2000

#Kmeans error threshold
t=0.0001

#no of MPI processes to spawn
p=32
threads=4

#Profiler output will be appended into this file
output=results/kmeans.$input.2000.hybridPublic.$p.$threads

CC=mpic++


rr=mpirun


echo "======================================"
echo "Script for running Serial Kmeans Code"
echo "======================================"
echo
echo "\tCompiling Code..."
$CC -Ofast mpi_main.c hybrid_kmeans.c wtime.c file_io.c mpi_io.c -fopenmp -lrt -lm -o hybrid-kmeans 
if [ $? -eq 1 ] ;	then
	echo "\tError occured during compilation\nTerminating script\n"
	exit 1
else
	echo "Compilation done successfully."
fi

echo "\n#############################################################">>$output
echo "#############################################################">>$output
echo "Following code was profiled on : `date`">>$output
echo "\nInput File : $input\nMethod : Hybrid Kmeans">>$output
echo "#############################################################">>$output
echo "#############################################################\n">>$output

echo "--------------------------------------------------------------\n">>$output

echo "Performing Hybrid Kmeans clustring on $input"
#use following to run the code on cluster
$rr -np $p -hostfile ~/hostlist32 --map-by node -x OMP_NUM_THREADS=$threads -x GOMP_CPU_AFFINITY=0-7 ./hybrid-kmeans -n $cluster -i $data/$input -t $t | tee $output

#alternatively, use following to run the code locally
#$rr -np $p ./mpi-kmeans -n $cluster -i $data/$input -t $t | tee $output

if [ $? -ne 0 ] ;	then
	echo "\tError occured. Terminating script...\n"
	exit 1
fi
echo "================================================================"
