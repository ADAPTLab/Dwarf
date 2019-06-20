# Dwarf
Dwarf is a Domain Specific Language for Representative-Based Clustering Algorithms.

Its Parallelizing Compiler takes the sequential Dwarf code and produces either Distributed-memory MPI C++ code or Hybrid-memory MPI+OpenMP C++ Code.

Dependencies:
1. The Dwarf Compiler requires Java 1.8
2. Generated parallel C++ code requires Dwarf Core C++ Library to compile, requires C++11 compiler.
