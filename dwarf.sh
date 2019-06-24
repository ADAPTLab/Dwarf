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

java -jar dwarf.jar $1 $CPP_OUT_DIR/$CPP_OUT_FILE $VERBOSITY_LEVEL $CONFIG_FILE