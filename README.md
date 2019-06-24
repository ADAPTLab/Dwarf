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
