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
		std::cout
			<< std::string(e.what()) + " At line no. 25, function " + "GMM in file //home\/master\/git\/dwarf\/CompilerDev\/Driver\/em9thOct17.dw.preprocessed"
			<< std::endl;
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

		if (clusList.Size() >= 5000)
		{
			int xx_41_xx = 0;
			MPI_Bcast(&xx_41_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			clusList.distributeList(numprocs);
			xx_2_xx.gatherList(numprocs);
		}
		else
		{
			for (auto &xx_5_xx : clusList.Elements())
			{
				auto xx_3_xx = distanceEuclideanSquared(xx_5_xx.oldmu,
														xx_5_xx.mu);
				xx_2_xx.AddEle(xx_3_xx);
			}
		}

		x = xx_2_xx;
		return *(double *)(Reduce(SUM, x.Elements(), true));
	}
	catch (std::exception &e)
	{
		std::cout
			<< std::string(e.what()) + " At line no. 41, function " + "convergence in file //home\/master\/git\/dwarf\/CompilerDev\/Driver\/em9thOct17.dw.preprocessed"
			<< std::endl;
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
		std::cout
			<< std::string(e.what()) + " At line no. 46, function " + "distanceEuclideanSquared in file //home\/master\/git\/dwarf\/CompilerDev\/Driver\/em9thOct17.dw.preprocessed"
			<< std::endl;
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
		List<int> xx_6_xx;
		int dim;
		Point B;
		B = A;
		dim = A.Size();
		for (int xx_7_xx = 0; xx_7_xx <= dim - 1; xx_7_xx++)
		{
			xx_6_xx.AddEle(xx_7_xx);
		}
		for (auto &i : xx_6_xx.Elements())
		{
			B.SetDoubleAt(i, A.GetEleAtIndex(i) * A.GetEleAtIndex(i));
		}
		return B;
	}
	catch (std::exception &e)
	{
		std::cout
			<< std::string(e.what()) + " At line no. 52, function " + "pointSquare in file //home\/master\/git\/dwarf\/CompilerDev\/Driver\/em9thOct17.dw.preprocessed"
			<< std::endl;
	}
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	List<double> numerators;
	List<double> wsuml;
	double denominator;
	List<List<double>> lll;
	List<int> xx_28_xx;
	File inputFile;
	List<int> xx_34_xx;
	int maxIter;
	double t4;
	int K;
	double t3;
	double t2;
	double t1;
	List<Point> munuml;
	List<Point> seeds;
	List<Point> ll;
	Point s;
	Point p;
	double EPS;
	int n;
	List<int> nlist;
	double error;
	int dim;
	int iter;
	List<int> klist;
	double c;
	List<EMCluster> clusList, _temp;
	int b;
	int a;
	List<Point> data;
	try
	{
		K = atoi(argv[2]);
		error = 100.0;
		EPS = 0.0001;
		maxIter = 5;
		t1 = MPI_Wtime();
		inputFile.ReadDataset(argv[1], &data);
		int xx_44_xx = -2;
		MPI_Bcast(&xx_44_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		data.distributeList(numprocs);
		MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		std::cout << "Time taken in File reading:\t" << t2 - t1 << std::endl;
		n = data.Size();
		dim = MaxEles;
		t2 = MPI_Wtime();
		List<Point> xx_8_xx;
		for (int xx_21_xx = 0; xx_21_xx <= K - 1; xx_21_xx++)
		{
			clusList.AddEle(EMCluster(data.GetEleAtIndex(xx_21_xx), dim, 1.0 / K));
		}

		t3 = MPI_Wtime();
		std::cout << "Time taken in Random Initial Points selection:\t"
				  << t3 - t2 << std::endl;
		iter = 1;

		t4 = MPI_Wtime();

		while (error > EPS && iter <= maxIter)
		{
			t2 = MPI_Wtime();

			for (auto &clus : clusList.Elements())
			{
				clus.oldmu = clus.mu;
				clus.wsum = 0.0;
				a = clus.musum.InitializePoint();
				b = clus.sigmasum.InitializePoint();
			}
			int xx_57_xx = 6;
			MPI_Bcast(&xx_57_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
			clusList.broadcastList(numprocs);

			for (int xx_61_xx = 1; xx_61_xx < numprocs; xx_61_xx++)
			{
				auto _temp = clusList.deserializeRecv(xx_61_xx);
				for (int xx_54_xx = 0; xx_54_xx < K; xx_54_xx++)
				{
					clusList[xx_54_xx].wsum += _temp[xx_54_xx].wsum;
					clusList[xx_54_xx].musum += _temp[xx_54_xx].musum;
					clusList[xx_54_xx].sigmasum += _temp[xx_54_xx].sigmasum;
				}
			}
			for (auto &clus : clusList.Elements())
			{
				clus.pc = clus.wsum / n;
				clus.mu = clus.musum / clus.wsum;
				clus.sigma = (clus.sigmasum / clus.wsum) - pointSquare(clus.mu);
			}

			error = convergence(clusList);
			t3 = MPI_Wtime();
			std::cout << "Iteration:" << iter << "\tError:" << error
					  << "\tTime taken:" << t3 - t2 << std::endl;
			iter = iter + 1;
		}
		t3 = MPI_Wtime();
		std::cout << "Time taken in all iterations:\t" << t3 - t4 << std::endl;
		std::cout << "Overall time taken:\t" << t3 - t1 << std::endl;

		for (int i = 0; i < K; i++)
		{
			p = clusList.GetEleAtIndex(i).mu;
			s = clusList.GetEleAtIndex(i).sigma;
			c = clusList.GetEleAtIndex(i).pc;
			std::cout << "Cluster centre: " << i << "\t" << p << "\t\tsigma: " << s << "\tpc: " << c << std::endl;
		}
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 0, function" + ""
				  << std::endl;
	}
	int xx_75_xx = -1;
	MPI_Bcast(&xx_75_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}