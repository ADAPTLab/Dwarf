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

		if (clusList.Size() >= 5000)
		{
			int xx_35_xx = 0;
			MPI_Bcast(&xx_35_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			clusList.distributeList(numprocs);
			xx_2_xx.gatherList(numprocs);
		}
		else
		{
			for (auto &xx_5_xx : clusList.Elements())
			{
				auto xx_3_xx = distanceEuclideanSquared(xx_5_xx.oldmu, xx_5_xx.mu);
				xx_2_xx.AddEle(xx_3_xx);
			}
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
	List<Point> data;
	List<int> xx_28_xx;
	int dim;
	double EPS;
	int maxIter;
	List<EMCluster> clusList;
	double error;
	int iter;
	int K;
	File inputFile;
	List<Point> seeds;
	try
	{
		K = 5;
		error = 0.2;
		EPS = 0.0001;
		maxIter = 10;
		inputFile.ReadDataset("datasets/3droad.arff", &data);
		int xx_39_xx = -2;
		MPI_Bcast(&xx_39_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		data.broadcastList(numprocs);
		dim = data.GetEleAtIndex(0).Size();
		List<Point> xx_8_xx;

		if ((K - 1) - (0) >= 5000)
		{
			int xx_40_xx = 1;
			MPI_Bcast(&xx_40_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			List<int> xx_12_xx;
			int xx_13_xx = 0;
			int xx_14_xx = K - 1;
			xx_12_xx.distributeRangeList(xx_13_xx, xx_14_xx, numprocs);

			xx_8_xx.gatherList(numprocs);
		}
		else
		{
			for (int xx_11_xx = 0; xx_11_xx <= K - 1; xx_11_xx++)
			{
				auto xx_9_xx = data.GetEleAtIndex(xx_11_xx);
				xx_8_xx.AddEle(xx_9_xx);
			}
		}

		seeds = xx_8_xx;
		List<EMCluster> xx_18_xx;

		if ((K - 1) - (0) >= 5000)
		{
			int xx_44_xx = 2;
			MPI_Bcast(&xx_44_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			List<int> xx_22_xx;
			seeds.broadcastList(numprocs);
			MPI_Bcast(&dim, sizeof(dim), MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&K, sizeof(K), MPI_INT, 0, MPI_COMM_WORLD);
			int xx_23_xx = 0;
			int xx_24_xx = K - 1;
			xx_22_xx.distributeRangeList(xx_23_xx, xx_24_xx, numprocs);

			xx_18_xx.gatherList(numprocs);
		}
		else
		{
			for (int xx_21_xx = 0; xx_21_xx <= K - 1; xx_21_xx++)
			{
				auto xx_19_xx = EMCluster(seeds.GetEleAtIndex(xx_21_xx), dim, 1.0 / K);
				xx_18_xx.AddEle(xx_19_xx);
			}
		}

		clusList = xx_18_xx;
		iter = 1;

		if ((K - 1) - (0) >= 5000)
		{
			int xx_48_xx = 3;
			MPI_Bcast(&xx_48_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
			int xx_29_xx = 0;
			int xx_30_xx = K - 1;
			xx_28_xx.distributeRangeList(xx_29_xx, xx_30_xx, numprocs);
			xx_28_xx.gatherList(numprocs);
		}
		else
		{
			for (int xx_49_xx = 0; xx_49_xx <= K - 1; xx_49_xx++)
			{
				xx_28_xx.AddEle(xx_49_xx);
			}
		}

		klist = xx_28_xx;

		while (error > EPS && iter <= maxIter)
		{

			if (clusList.Size() >= 5000)
			{
				int xx_50_xx = 4;
				MPI_Bcast(&xx_50_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				clusList.distributeList(numprocs);
				clusList.gatherList(numprocs);
			}
			else
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
			}

			if (data.Size() >= 5000)
			{
				int xx_53_xx = 5;
				MPI_Bcast(&xx_53_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				clusList.broadcastList(numprocs);
				klist.broadcastList(numprocs);

				HashMap<double> xx_54_xx = HashMap<double>();

				HashMap<Point> xx_55_xx = HashMap<Point>();

				HashMap<Point> xx_56_xx = HashMap<Point>();
				for (int xx_57_xx = 1; xx_57_xx < numprocs; xx_57_xx++)
				{
					xx_54_xx = xx_54_xx.deserializeRecv(xx_57_xx);
					for (auto iterator = xx_54_xx.htmap.begin(); iterator != xx_54_xx.htmap.end(); iterator++)
					{
						uint i = iterator->first;
						clusList[i].wsum += xx_54_xx.htmap[i];
					}
				}
				for (int xx_58_xx = 1; xx_58_xx < numprocs; xx_58_xx++)
				{
					xx_55_xx = xx_55_xx.deserializeRecv(xx_58_xx);
					for (auto iterator = xx_55_xx.htmap.begin(); iterator != xx_55_xx.htmap.end(); iterator++)
					{
						uint i = iterator->first;
						clusList[i].musum += xx_55_xx.htmap[i];
					}
				}
				for (int xx_59_xx = 1; xx_59_xx < numprocs; xx_59_xx++)
				{
					xx_56_xx = xx_56_xx.deserializeRecv(xx_59_xx);
					for (auto iterator = xx_56_xx.htmap.begin(); iterator != xx_56_xx.htmap.end(); iterator++)
					{
						uint i = iterator->first;
						clusList[i].sigmasum += xx_56_xx.htmap[i];
					}
				}
			}
			else
			{
				for (auto &xj : data.Elements())
				{
					List<double> numerators;
					List<double> xx_34_xx;
					double denominator;

					for (auto &clus : clusList.Elements())
					{
						xx_34_xx.AddEle(GMM(xj, clus.mu, clus.sigma, clus.pc));
					}
					numerators = xx_34_xx;
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
			}

			if (clusList.Size() >= 5000)
			{
				int xx_78_xx = 6;
				MPI_Bcast(&xx_78_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				clusList.distributeList(numprocs);
				clusList.gatherList(numprocs);
			}
			else
			{
				for (auto &clus : clusList.Elements())
				{
					clus.pc = clus.wsum / n;
					clus.mu = clus.musum / clus.wsum;
					clus.sigma = (clus.sigmasum / clus.wsum) - pointSquare(clus.mu);
				}
			}

			error = convergence(clusList);
			iter = iter + 1;
		}

		if (klist.Size() >= 5000)
		{
			int xx_81_xx = 7;
			MPI_Bcast(&xx_81_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
			clusList.broadcastList(numprocs);
			klist.distributeList(numprocs);
			klist.gatherList(numprocs);
		}
		else
		{
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
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 0, function" + "" << std::endl;
	}
	int xx_84_xx = -1;
	MPI_Bcast(&xx_84_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
