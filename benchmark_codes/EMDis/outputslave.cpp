#include "output.h"
#include <iostream>

List<Point> data;
int K;
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

void __ff__0()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<EMCluster> xx_4_xx;
	xx_4_xx = xx_4_xx.getLocalList();
	List<double> xx_42_xx;
	for (auto &xx_5_xx : xx_4_xx.Elements())
	{
		auto xx_3_xx = distanceEuclideanSquared(xx_5_xx.oldmu, xx_5_xx.mu);
		xx_42_xx.AddEle(xx_3_xx);
	}

	xx_42_xx.sendBack();
}

void __ff__1()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<Point> data;

	data = ::data;
	int xx_15_xx;
	int xx_16_xx;
	MPI_Recv(&xx_15_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_16_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_17_xx;
	xx_17_xx.makeRangeList(xx_15_xx, xx_16_xx);
	List<Point> xx_46_xx;
	for (auto &xx_11_xx : xx_17_xx.Elements())
	{
		auto xx_9_xx = data.GetEleAtIndex(xx_11_xx);
		xx_46_xx.AddEle(xx_9_xx);
	}

	xx_46_xx.sendBack();
}

void __ff__2()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<Point> seeds;

	seeds = seeds.getLocalList();
	int dim;
	MPI_Bcast(&dim, sizeof(dim), MPI_INT, 0, MPI_COMM_WORLD);
	int K;
	MPI_Bcast(&K, sizeof(K), MPI_INT, 0, MPI_COMM_WORLD);
	int xx_25_xx;
	int xx_26_xx;
	MPI_Recv(&xx_25_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_26_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_27_xx;
	xx_27_xx.makeRangeList(xx_25_xx, xx_26_xx);
	List<EMCluster> xx_49_xx;
	for (auto &xx_21_xx : xx_27_xx.Elements())
	{
		auto xx_19_xx = EMCluster(seeds.GetEleAtIndex(xx_21_xx), dim, 1.0 / K);
		xx_49_xx.AddEle(xx_19_xx);
	}

	xx_49_xx.sendBack();
}

void __ff__3()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	int xx_31_xx;
	int xx_32_xx;
	MPI_Recv(&xx_31_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_32_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_33_xx;
	xx_33_xx.makeRangeList(xx_31_xx, xx_32_xx);
	xx_33_xx.sendBack();
}

void __ff__4()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	int xx_37_xx;
	int xx_38_xx;
	MPI_Recv(&xx_37_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_38_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_39_xx;
	xx_39_xx.makeRangeList(xx_37_xx, xx_38_xx);
	xx_39_xx.sendBack();
}

void __ff__5()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	EMCluster clus;
	int a;
	int b;

	List<EMCluster> xx_56_xx;
	xx_56_xx = xx_56_xx.getLocalList();
	for (auto &clus : xx_56_xx.Elements())
	{
		clus.oldmu = clus.mu;
		clus.wsum = 0.0;
		a = clus.musum.InitializePoint();
		b = clus.sigmasum.InitializePoint();
	}

	if (myrank == (numprocs - 1))
	{
		MPI_Send(&a, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&b, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	xx_56_xx.sendBack();
}

void __ff__6()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	Point xj;
	List<double> numerators;
	List<EMCluster> clusList;
	List<int> klist;
	double denominator;

	clusList = clusList.getLocalList();

	for (auto &xj : data.Elements())
	{
		List<double> xx_40_xx;
		uint indexExp;
		denominator = 0.0;
		for (auto &clus : clusList.Elements())
		{
			auto _temp = GMM(xj, clus.mu, clus.sigma, clus.pc);
			numerators.AddEle(_temp);
			denominator += _temp;
		}
		double gamma;
		for (int i = 0; i < K; i++)
		{
			Point temp;
			Point temp1;
			uint indexExp;
			gamma = numerators.GetEleAtIndex(i) / denominator;
			clusList[i].wsum += gamma;

			temp = xj * gamma;
			clusList[i].musum += temp;

			temp1 = pointSquare(xj) * gamma;
			clusList[i].sigmasum += temp1;
		}
		numerators.Free();
	}

	clusList.sendBack();
}

void __ff__7()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	EMCluster clus;
	int n;

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	List<EMCluster> xx_72_xx;
	xx_72_xx = xx_72_xx.getLocalList();
	for (auto &clus : xx_72_xx.Elements())
	{
		clus.pc = clus.wsum / n;
		clus.mu = clus.musum / clus.wsum;
		clus.sigma = (clus.sigmasum / clus.wsum) - pointSquare(clus.mu);
	}

	if (myrank == (numprocs - 1))
	{
	}
	xx_72_xx.sendBack();
}

void __ff__8()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	Point p;
	Point s;
	double c;
	int i;
	List<EMCluster> clusList;

	clusList = clusList.getLocalList();

	List<int> xx_74_xx;
	xx_74_xx = xx_74_xx.getLocalList();
	for (auto &i : xx_74_xx.Elements())
	{
		p = clusList.GetEleAtIndex(i).mu;
		s = clusList.GetEleAtIndex(i).sigma;
		c = clusList.GetEleAtIndex(i).pc;
		std::cout << "Cluster centre: " << i << "\t" << p << "\t\tsigma: " << s
				  << "\tpc: " << c << std::endl;
	}

	xx_74_xx.sendBack();
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
			MPI_Bcast(&K, 1, MPI_INT, 0, MPI_COMM_WORLD);
			break;
		case 0:
			__ff__0();
			break;
		case 1:
			__ff__1();
			break;
		case 2:
			__ff__2();
			break;
		case 3:
			__ff__3();
			break;
		case 4:
			__ff__4();
			break;
		case 5:
			__ff__5();
			break;
		case 6:
			__ff__6();
			break;
		case 7:
			__ff__7();
			break;
		case 8:
			__ff__8();
			break;
		default:
			cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}