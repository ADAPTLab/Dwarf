# Dwarf
Dwarf is a Domain Specific Language for Representative-Based Clustering Algorithms.
Its Parallelizing Compiler takes the sequential Dwarf code and produces either Distributed-memory MPI C++ code or Hybrid-memory MPI+OpenMP C++ Code.

## Usage
```shell
sh dwarf.sh <dwarf file name>.dw
```
  It will produce output.cpp, outputslave.cpp, output.h, Point.cpp, and Point.h in cppsrc directory.
  

## Dependencies
The Dwarf Compiler requires Java 1.8.
It requires config.txt and weka.jar available in dependencies directory.

## Generated CPP Code
Check cppsrc directory for SERIAL, DISTRIBUTED, and HYBRID codes of K-means and EM algorithm generated from kmeans.dw and em.dw source codes.
