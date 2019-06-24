#include "output.h"
#include <iostream>

List<Point> data;

EMCluster::EMCluster(Point p, int dim, double w)
{
	mu = p;
	oldmu = p;
	musum = p;
	sigma = p;
	wsum = w;
	pc = w;
	dim = sigma.InitializePoint(1.0);
	sigmasum = sigma;
}
double GMM(Point x, Point mu, Point sigma, double prob)
{
	try
	{
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Status status;
		double b;
		double prod;
		double a;
		List<int> xx_0_xx;
		double sum;
		double piv;
		int dim;
		double c;
		dim = x.Size();
		piv = 1.0;
		prod = 1.0;
		sum = 0.0;
		for (int xx_1_xx = 0; xx_1_xx <= dim - 1; xx_1_xx++)
		{
			xx_0_xx.AddEle(xx_1_xx);
		}
		for (auto &i : xx_0_xx.Elements())
		{
			sum = sum + ((x.GetEleAtIndex(i) - mu.GetEleAtIndex(i)) * (x.GetEleAtIndex(i) - mu.GetEleAtIndex(i)) / sigma.GetEleAtIndex(i));
			prod = prod * sigma.GetEleAtIndex(i);
			piv = piv * PI();
		}
		a = exponent(-sum / 2);
		b = sqrt(piv * (1 << dim) * prod);
		c = (a / b) * prob;
		return c;
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 21, function " + "GMM in file dwarf_source_codes\/em.dw" << std::endl;
	}
}
double convergence(List<EMCluster> clusList)
{
	try
	{
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Status status;
		List<double> x;
		List<double> xx_2_xx;

		for (auto &clus : clusList.Elements())
		{
			xx_2_xx.AddEle(distanceEuclideanSquared(clus.oldmu, clus.mu));
		}
		x = xx_2_xx;
		return *(double *)(Reduce(SUM, x.Elements(), true));
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 34, function " + "convergence in file dwarf_source_codes\/em.dw" << std::endl;
	}
}
double distanceEuclideanSquared(Point A, Point B)
{
	try
	{
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Status status;
		double t;
		t = distanceEuclidean(A, B);
		return t * t;
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 38, function " + "distanceEuclideanSquared in file dwarf_source_codes\/em.dw" << std::endl;
	}
}
Point pointSquare(Point A)
{
	try
	{
		int myrank;
		MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		int numprocs;
		MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
		MPI_Status status;
		List<int> xx_3_xx;
		int dim;
		Point B;
		B = A;
		dim = A.Size();
		for (int xx_4_xx = 0; xx_4_xx <= dim - 1; xx_4_xx++)
		{
			xx_3_xx.AddEle(xx_4_xx);
		}
		for (auto &i : xx_3_xx.Elements())
		{
			B.SetDoubleAt(i, A.GetEleAtIndex(i) * A.GetEleAtIndex(i));
		}
		return B;
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 42, function " + "pointSquare in file dwarf_source_codes\/em.dw" << std::endl;
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