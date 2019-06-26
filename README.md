# Dwarf
Dwarf is a Domain Specific Language for Representative-Based Clustering Algorithms.
Its Parallelizing Compiler takes the sequential Dwarf code and produces either Distributed-memory MPI C++ code or Hybrid-memory MPI+OpenMP C++ Code.

## Usage
```shell
sh dwarf.sh <dwarf file name>.dw
```
  It will produce [output.cpp](cppsrc/output.cpp) , [outputslave.cpp](cppsrc/outputslave.cpp), [output.h](cppsrc/output.h), [Point.cpp](cppsrc/Point.cpp), and [Point.h](cppsrc/Point.h) in [cppsrc](cppsrc) directory.
  

## Dependencies
The Dwarf Compiler requires Java 1.8.
It requires [config.txt](dependencies/config.txt) and [weka.jar](dependencies/weka.jar) available in the [dependencies](dependencies) directory.

## Generated CPP Code
Compiler generates code in cppsrc directory. It already has a few codes.

| Platform\Algorithm  | K-means Algorithm | EM Algorithm | 
| ------------------- | ----------------- | ------------ |
| Serial C++ Code | [KDSeq](cppsrc/KDSeq)  | [EDSeq](cppsrc/EDSeq)  |
| Distributed-memory Parallel MPI C++ Code | [KDDis](cppsrc/KDDis)  | [EDDis](cppsrc/EDDis)  |
| Hybrid-memory Parallel MPI OpenMP C++ Code | [KDHyb](cppsrc/KDHyb)  | [EDHyb](cppsrc/EDHyb)  |

## Benchmark Codes
Comparison of various manual parallel implementations of Representative-based Clustering Algorithms with Dwarf Compiler generated automatically parallelized codes.
* K-means Algorithm
  * KMSeq (Manual C): [Prof. Wei-keng Liao's Repository](http://www.ece.northwestern.edu/~wkliao/Kmeans/index.html)
  * KMDis (Manual MPI C): [Prof. Wei-keng Liao's Repository](http://www.ece.northwestern.edu/~wkliao/Kmeans/index.html)
  * KMHyb (Manual MPI C): [KMHyb](benchmark_codes/KMHyb)
  * KSJav (Manual Spark Java): [AMP Camp Two - Big Data Bootcamp Strata 2013](http://ampcamp.berkeley.edu/exercises-strata-conf-2013/index.html)
  * KSScl (Manual Spark Scala): [AMP Camp Two - Big Data Bootcamp Strata 2013](http://ampcamp.berkeley.edu/exercises-strata-conf-2013/index.html)

* EM Algorithm
  * EMSeq (Manual C): [EMSeq](benchmark_codes/EMSeq)
  * EMDis (Manual MPI C++): [EMDis](benchmark_codes/EMDis)
  * EMHyb (Manual MPI OpenMP C++): [EMHyb](benchmark_codes/EMHyb)
  * ESJav (Manual Spark Java): [ESJav](benchmark_codes/ESJav)