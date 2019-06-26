# Default values for Path and Name of generated executable
EXEC_PATH = ./
HOSTFILE= ~/hostlist32

#CC_PATH=/share/apps/openmpi-1.10.0/bin/
CC_PATH=
# Location of C++ ARFF Parser Source Directory
#ARFF_SRC = arff
ARFF_SRC=
# Location of Cereal library files
CEREAL_SRC = include

# C++ Compiler
CC = $(CC_PATH)mpic++

# C++ Compiler Flags
CC_FLAGS = -std=c++11 -w -Ofast -fstack-protector -I $(CEREAL_SRC) -fopenmp
#CC_FLAGS = -std=c++11 -g -w -fstack-protector -I $(CEREAL_SRC) -fopenmp

SEXEC_NAME = soutput
SEXEC_NAME_SLAVE = soutputslave
MEXEC_NAME = moutput
MEXEC_NAME_SLAVE = moutputslave
HEXEC_NAME = houtput
HEXEC_NAME_SLAVE = houtputslave

# File names
SSOURCES = $(wildcard arff_*.cpp) Point.cpp serialoutput.cpp
SOBJECTS = $(SSOURCES:.cpp=.o)

MSOURCES = $(wildcard arff_*.cpp) Point.cpp mpioutput.cpp
MOBJECTS = $(MSOURCES:.cpp=.o)

# Slave File names
MSOURCES_SLAVE = Point.cpp mpioutputslave.cpp
MOBJECTS_SLAVE = $(MSOURCES_SLAVE:.cpp=.o)

HSOURCES = $(wildcard arff_*.cpp) Point.cpp hybridoutput.cpp
HOBJECTS = $(HSOURCES:.cpp=.o)

# Slave File names
HSOURCES_SLAVE = Point.cpp hybridoutputslave.cpp
HOBJECTS_SLAVE = $(HSOURCES_SLAVE:.cpp=.o)

# Main target
all: serialc mpic hybridc

# serial compilation head
serialc: $(SOBJECTS)
	@$(CC) $(CC_FLAGS) $(SOBJECTS) -o $(EXEC_PATH)$(SEXEC_NAME)

# mpi compilation head
mpic: mpimaster mpislave
mpimaster: $(MOBJECTS)
	$(CC) $(CC_FLAGS) $(MOBJECTS) -o $(EXEC_PATH)$(MEXEC_NAME)

mpislave: $(MOBJECTS_SLAVE)
	$(CC) $(CC_FLAGS) $(MOBJECTS_SLAVE) -o $(EXEC_PATH)$(MEXEC_NAME_SLAVE)

# hybrid compilation head
hybridc: hybridmaster hybridslave

hybridmaster: $(HOBJECTS)
	$(CC) $(CC_FLAGS) $(HOBJECTS) -o $(EXEC_PATH)$(HEXEC_NAME)

hybridslave: $(HOBJECTS_SLAVE)
	$(CC) $(CC_FLAGS) $(HOBJECTS_SLAVE) -o $(EXEC_PATH)$(HEXEC_NAME_SLAVE)

# To obtain object files
%.o: %.cpp
	@$(CC) $(CC_FLAGS) -c $< -o $@

# vampirtrace generated files
VTFILES = *.def *.events *.stats *.otf *.thumb *.z *.o

# To remove generated files
clean:
	rm -f $(SEXEC_NAME) $(OBJECTS) $(SEXEC_NAME_SLAVE) $(OBJECTS_SLAVE) $(VTFILES)  $(MEXEC_NAME) $(MEXEC_NAME_SLAVE) $(HEXEC_NAME) $(HEXEC_NAME_SLAVE)

# To run generated output executable (serial)
serialr:
#	$(CC_PATH)mpirun --mca opal_abort_print_stack 1 -np 1 -host compute-0-2 $(EXEC_PATH)$(SEXEC_NAME) $(f) $(K)
	$(CC_PATH)mpirun -np 1 -host compute-0-2 $(EXEC_PATH)$(SEXEC_NAME) $(f) $(K)

# To run generated output executable using mpi(parallel)
mpirlocal:
	$(CC_PATH)mpirun -np 1 $(EXEC_PATH)$(MEXEC_NAME) $(f) $(K) : -np 3 $(EXEC_PATH)$(MEXEC_NAME_SLAVE)

mpircluster:
	mpirun -np 1 -host localhost $(EXEC_PATH)$(MEXEC_NAME) $(f) $(K) : -np $(p) --map-by node -hostfile $(HOSTFILE) $(EXEC_PATH)$(MEXEC_NAME_SLAVE)

# To run generated output executable using mpi + openmp (hybrid)
hybridrlocal: export OMP_NUM_THREADS = $(ARG)
hybridrlocal:
	mpirun -np 1 $(EXEC_PATH)$(HEXEC_NAME) : -np 2 -x OMP_NUM_THREADS=$(ARG) -x GOMP_CPU_AFFINITY=0-7 $(EXEC_PATH)$(HEXEC_NAME_SLAVE)

hybridrcluster: export OMP_NUM_THREADS = $(ARG)
hybridrcluster:
	$(CC_PATH)mpirun -np 1 -host localhost $(EXEC_PATH)$(HEXEC_NAME) $(f) $(K) : -np $(p) --map-by node -hostfile $(HOSTFILE) -x OMP_NUM_THREADS=$(ARG) -x GOMP_CPU_AFFINITY=0-7 $(EXEC_PATH)$(HEXEC_NAME_SLAVE)


# To debug using xterm in parallel
debug:
	mpirun -np 1 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME)

# To debug using xterm in parallel
debugparallel:
	mpirun -np 1 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME) : -np 2 xterm -e gdb $(EXEC_PATH)$(EXEC_NAME_SLAVE)

#for memory leak detection
mem:
	mpirun -np 4 xterm -e valgrind --vgdb=yes --vgdb-error=0 $(EXEC_PATH)$(EXEC_NAME)

# temp
memcheck:
	mpirun -np 1 $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --leak-check=full --log-file=MemoryProf_house_1lac $(EXEC_PATH)$(EXEC_NAME_SLAVE)
massif:
	mpirun -np 1 valgrind --tool=massif --stacks=yes $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --tool=massif --stacks=yes $(EXEC_PATH)$(EXEC_NAME_SLAVE)
callgrind:
	mpirun -np 1 valgrind --tool=callgrind $(EXEC_PATH)$(EXEC_NAME) : -np 4 valgrind --tool=callgrind $(EXEC_PATH)$(EXEC_NAME_SLAVE)
memcheckser:
	mpirun -np 1 valgrind --leak-check=full --log-file=mem_log.txt $(EXEC_PATH)$(EXEC_NAME)

