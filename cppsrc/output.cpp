#include "output.h"
#include <iostream>
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
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	try{
	func();
	std::cout << "test successful\n" << std::endl;
		}
		catch(std::exception &e) {
			 std::cout << std::string(e.what()) +" At line no. 0, function"
	 + "" << std::endl;
}
		int xx_0_xx = -1;
		MPI_Bcast(&xx_0_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
