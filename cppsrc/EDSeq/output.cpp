#include "output.h"
#include <iostream>
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
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	List<int> klist;
	List<int> xx_11_xx;
	List<Point> data;
	int dim;
	double EPS;
	int maxIter;
	List<Point> xx_5_xx;
	List<EMCluster> clusList;
	double error;
	int iter;
	int K;
	File inputFile;
	List<EMCluster> xx_8_xx;
	List<Point> seeds;
	try
	{
		K = 5;
		error = 0.2;
		EPS = 0.0001;
		maxIter = 10;
		inputFile.ReadDataset("datasets/3droad.arff", &data);
		dim = data.GetEleAtIndex(0).Size();

		List<int> xx_6_xx;
		for (int xx_7_xx = 0; xx_7_xx <= K - 1; xx_7_xx++)
		{
			xx_6_xx.AddEle(xx_7_xx);
		}
		for (auto &i : xx_6_xx.Elements())
		{
			xx_5_xx.AddEle(data.GetEleAtIndex(i));
		}
		seeds = xx_5_xx;

		List<int> xx_9_xx;
		for (int xx_10_xx = 0; xx_10_xx <= K - 1; xx_10_xx++)
		{
			xx_9_xx.AddEle(xx_10_xx);
		}
		for (auto &i : xx_9_xx.Elements())
		{
			xx_8_xx.AddEle(EMCluster(seeds.GetEleAtIndex(i), dim, 1.0 / K));
		}
		clusList = xx_8_xx;
		iter = 1;
		for (int xx_12_xx = 0; xx_12_xx <= K - 1; xx_12_xx++)
		{
			xx_11_xx.AddEle(xx_12_xx);
		}
		klist = xx_11_xx;

		while (error > EPS && iter <= maxIter)
		{
			for (auto &clus : clusList.Elements())
			{
				int b;
				int a;
				clus.oldmu = clus.mu;
				clus.wsum = 0.0;
				a = clus.musum.InitializePoint();
				b = clus.sigmasum.InitializePoint();
			}
			for (auto &xj : data.Elements())
			{
				List<double> numerators;
				List<double> xx_13_xx;
				double denominator;

				for (auto &clus : clusList.Elements())
				{
					xx_13_xx.AddEle(GMM(xj, clus.mu, clus.sigma, clus.pc));
				}
				numerators = xx_13_xx;
				denominator = *(double *)(Reduce(SUM, numerators.Elements(), true));
				for (auto &i : klist.Elements())
				{
					double gamma;
					Point temp;
					Point temp1;
					gamma = numerators.GetEleAtIndex(i) / denominator;
					clusList.GetEleAtIndex(i).wsum = clusList.GetEleAtIndex(i).wsum + gamma;
					temp = xj * gamma;
					clusList.GetEleAtIndex(i).musum = clusList.GetEleAtIndex(i).musum + temp;
					temp1 = pointSquare(xj) * gamma;
					clusList.GetEleAtIndex(i).sigmasum = clusList.GetEleAtIndex(i).sigmasum + temp1;
				}
				numerators.Free();
			}
			for (auto &clus : clusList.Elements())
			{
				clus.pc = clus.wsum / n;
				clus.mu = clus.musum / clus.wsum;
				clus.sigma = (clus.sigmasum / clus.wsum) - pointSquare(clus.mu);
			}
			error = convergence(clusList);
			iter = iter + 1;
		}
		for (auto &i : klist.Elements())
		{
			Point s;
			Point p;
			double c;
			p = clusList.GetEleAtIndex(i).mu;
			s = clusList.GetEleAtIndex(i).sigma;
			c = clusList.GetEleAtIndex(i).pc;
			std::cout << "Cluster centre: " << i << "\t" << p << "\t\tsigma: " << s << "\tpc: " << c << std::endl;
		}
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 0, function" + "" << std::endl;
	}
	int xx_14_xx = -1;
	MPI_Bcast(&xx_14_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
