#include "output.h"
#include <iostream>

List<Point> null;

void func() {
try {
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Status status;
}
catch(std::exception &e) {
	 std::cout << std::string(e.what()) +" At line no. 2, function "
	 + "func in file Basic_tests\/test-1.dw" << std::endl;
}
}
int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int signal = 0;
	while(signal != -1){
		MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(signal == -1) break;
		switch(signal){
			case -2 : null = null.getLocalList();
			break;
			default : cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}