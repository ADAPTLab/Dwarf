# Default values for Path and Name of generated executable
EXEC_NAME = master
EXEC_NAME_SLAVE = slave
host = localhost
CC_PATH=
ARFF_SRC=
# Location of Cereal library files
CEREAL_SRC = include

# C++ Compiler
CC = $(CC_PATH)mpic++

# C++ Compiler Flags
CC_FLAGS = -std=c++11 -w -Ofast -fstack-protector -I $(CEREAL_SRC) -fopenmp

SOURCES = $(wildcard arff_*.cpp) Point.cpp output.cpp
OBJECTS = $(SOURCES:.cpp=.o)

master: $(OBJECTS)
	$(CC) $(CC_FLAGS) $(OBJECTS) -o $(EXEC_PATH)$(EXEC_NAME)

SOURCES_SLAVE = Point.cpp outputslave.cpp
OBJECTS_SLAVE = $(SOURCES_SLAVE:.cpp=.o)

# Slave target
slave: $(OBJECTS_SLAVE)
	$(CC) $(CC_FLAGS) $(OBJECTS_SLAVE) -o $(EXEC_PATH)$(EXEC_NAME_SLAVE)

all: master slave

%.o: %.cpp
	@$(CC) -c $(CC_FLAGS) $< -o $@

clean:
	rm -f $(EXEC_NAME) $(OBJECTS) $(EXEC_NAME_SLAVE) $(OBJECTS_SLAVE)

# Targets for Executing code on a cluster
# Assumes that a hostlist32 file is available in the home directory with node names of all machines in the cluster

mpircluster:
	$(CC_PATH)mpirun -np 1 -host $(host) $(EXEC_PATH)$(EXEC_NAME) $(f) $(K) $(m) $(out) $(timefile) : -np $(p) --map-by node -hostfile ~/hostlist32 $(EXEC_PATH)$(EXEC_NAME_SLAVE)

serial:
	$(CC_PATH)mpirun -np 1 -host $(host) $(EXEC_PATH)$(EXEC_NAME) $(f) $(K) $(m) $(out) $(timefile)

hybridrcluster:
	$(CC_PATH)mpirun -np 1 --pernode -host $(host) $(EXEC_PATH)$(EXEC_NAME) $(f) $(K) $(m) $(out) $(timefile) : -np $(p) -x OMP_NUM_THREADS=$(t) -x GOMP_CPU_AFFINITY=0-7 -hostfile ~/hostlist32 $(EXEC_PATH)$(EXEC_NAME_SLAVE)


# Targets for local testing
localpar:
	$(CC_PATH)mpirun -oversubscribe -np 1 $(EXEC_PATH)$(EXEC_NAME) $(K) : -np $(p) $(EXEC_PATH)$(EXEC_NAME_SLAVE)

localserial:
	$(CC_PATH)mpirun -np 1 -host localhost $(EXEC_PATH)$(EXEC_NAME) $(K) $(m) $(out) $(timefile)

localhybrid:
	$(CC_PATH)mpirun -oversubscribe -np 1 $(EXEC_PATH)$(EXEC_NAME) $(K) : -np $(p) -x OMP_NUM_THREADS=$(t) -x GOMP_CPU_AFFINITY=0-7 $(EXEC_PATH)$(EXEC_NAME_SLAVE)

# targets to debug using xterm in parallel
debug:
	mpirun -np 1 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME)

# To debug using xterm in parallel
debugparallel:
	mpirun -np 1 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME) : -np 2 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME_SLAVE)

#for memory leak detection
memcheck:
	mpirun -np 1 $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --leak-check=full --log-file=MemoryProf $(EXEC_PATH)$(EXEC_NAME_SLAVE)
massif:
	mpirun -np 1 valgrind --tool=massif --stacks=yes $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --tool=massif --stacks=yes $(EXEC_PATH)$(EXEC_NAME_SLAVE)
callgrind:
	mpirun -np 1 valgrind --tool=callgrind $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --tool=callgrind $(EXEC_PATH)$(EXEC_NAME_SLAVE)
memcheckser:
	mpirun -np 1 valgrind --leak-check=full --log-file=mem_log.txt $(EXEC_PATH)$(EXEC_NAME)
callgrindserial:
	$(CC_PATH)mpirun -np 1 valgrind --tool=callgrind $(EXEC_PATH)$(EXEC_NAME) $(f) $(K) $(m) $(out) $(timefile)