#!/bin/sh

# Assumption that dataset file in arff file format is available in this dirctory
CPP_OUT_DIR=cppsrc

# Change in this name will require change in C++ makefile in CPP_OUT_DIR
CPP_OUT_FILE=output

# Allowed levels: [0,5]
VERBOSITY_LEVEL=0

# Change "Model" parameter in the CONFIG_FILE to generate SERIAL (serial), DSITRIBUTED (mpi), HYBRID (hybrid) code
CONFIG_FILE=./dependencies/config.txt

if [ "$#" -ne 1 ]; then
  echo "Usage: sh dwarf.sh <dwarf file name>.dw" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi

echo "\tCompiling Dwarf Source Code using Dwarf Compiler"
java -jar dwarf.jar $1 $CPP_OUT_DIR/$CPP_OUT_FILE $VERBOSITY_LEVEL $CONFIG_FILE

if [ $? -eq 1 ] ;	then
	echo "\tError occured during compilation\nTerminating script\n"
	exit 1
else
	echo "\n\tDwarf Compilation Successful.\n"
fi

cd $CPP_OUT_DIR

echo "\tCompiling Dwarf Compiler generated C++ Code using C++ Compiler"
make clean all

if [ $? -eq 1 ] ;	then
	echo "\tError occured during compilation\nTerminating script\n"
	exit 1
else
	echo "\n\tC++ Compilation Successful.\n"
fi

# Uncomment one of the following
echo "\tStarting Sequential Execution on the local machine\n"
make localserial

#echo "\tStarting Distributed-memory Execution on the local machine using 4 local processes\n"
#make localpar p=4

#echo "\tStarting Distributed-memory Execution on the Cluster using 32 distributed processes\n"
#make mpircluster p=32

#echo "\tStarting Hybrid-memory Execution on the Cluster using 32 distributed processes with 4 threads in each\n"
#make hybridrcluster p=32 ARG=4