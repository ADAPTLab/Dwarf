#include "output.h"
#include <iostream>

List<Point> data;

KmeansCluster::KmeansCluster(Point p)
{
	rep = p;
	oldRep = p;
	count = 0;
}
double ConvergenceCriteria(List<KmeansCluster> CL)
{
	try
	{
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Status status;
		List<double> x;
		List<double> xx_0_xx;
		double y;

		for (auto &c : CL.Elements())
		{
			xx_0_xx.AddEle(distanceEuclidean(c.oldRep, c.rep));
		}
		x = xx_0_xx;
		y = *(double *)(Reduce(SUM, x.Elements(), true));
		return (y / CL.Size());
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 12, function " + "ConvergenceCriteria in file dwarf_source_codes\/kmeans.dw" << std::endl;
	}
}
int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int signal = 0;
	while (signal != -1)
	{
		MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (signal == -1)
			break;
		switch (signal)
		{
		case -2:
			data = data.getLocalList();
			break;
		default:
			cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}