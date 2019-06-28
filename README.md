# Dwarf
Dwarf is a Domain Specific Language for Representative-Based Clustering Algorithms.
Its Parallelizing Compiler takes the sequential Dwarf code and produces either Distributed-memory MPI C++ code or Hybrid-memory MPI+OpenMP C++ Code.

## Usage
```shell
sh dwarf.sh <dwarf file name>.dw
```
1. Compiles Dwarf Source code to generate C++ code as per [Dwarf Compiler Flags](dependencies/config.txt)
   * Master's C++ Code: [output.cpp](cppsrc/output.cpp)
   * Slave's C++ Code: [outputslave.cpp](cppsrc/outputslave.cpp)
   * C++ Header File: [output.h](cppsrc/output.h)
   * Point Type Implementation C++ Code: [Point.cpp](cppsrc/Point.cpp)
   * [Point.h](cppsrc/Point.h) in [cppsrc](cppsrc) directory.
2. Compiles C++ Code to generate executables
   * master
   * slave
3. Execute the executables using one of make targets
   * localserial : Sequential Execution on the local machine
   * localpar p=4 : Distributed-memory Execution by 4 processes on the local machine
   * localhybrid p=4 t=2 : Hybrid-memory Execution by 4 processes and 2 threads on the local machine
   * serial : Sequential Execution on the given *host* machine
   * mpircluster p=4 : Distributed-memory Execution by 4 processes on the *cluster* (assumes a hostlist for MPI)
   * hybridrcluster p=4 t=2 : Hybrid-memory Execution by 4 processes and 2 threads on the *cluster* (assumes a hostlist for MPI)


## Dependencies
* Java version 1.8
* C++ version C++11
* Weka Java API : [weka.jar](dependencies/weka.jar) available in the [dependencies](dependencies) directory.
* Dwarf Compiler Flag File : [config.txt](dependencies/config.txt) available in the [dependencies](dependencies) directory.

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

## Sample Test Codes
* [Basic Tests](dwarf_source_codes/test_cases/basic_tests)
* [Object Oriented Programs Tests](dwarf_source_codes/test_cases/oop_tests)
* [Loop Tests](dwarf_source_codes/test_cases/loop_tests)