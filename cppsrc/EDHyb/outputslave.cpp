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

void __ff__0()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<EMCluster> xx_4_xx;
	xx_4_xx = xx_4_xx.getLocalList();
	List<double> xx_36_xx;
	xx_36_xx = List<double>(xx_4_xx.Size());
	int xx_37_xx = 0;
#pragma omp parallel for firstprivate(xx_37_xx) schedule(static)
	for (auto xx_38_xx = xx_4_xx.Elements().begin(); xx_38_xx < xx_4_xx.Elements().end(); xx_38_xx++)
	{
		auto &xx_5_xx = *xx_38_xx;
		auto xx_3_xx = distanceEuclideanSquared(xx_5_xx.oldmu, xx_5_xx.mu);
		xx_36_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_4_xx.Size()) + xx_37_xx, xx_3_xx);
		xx_37_xx++;
	}

	xx_36_xx.sendBack();
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
	List<Point> xx_41_xx;
	xx_41_xx = List<Point>(xx_17_xx.Size());
	int xx_42_xx = 0;
#pragma omp parallel for firstprivate(xx_42_xx) schedule(static)
	for (auto xx_43_xx = xx_17_xx.Elements().begin(); xx_43_xx < xx_17_xx.Elements().end(); xx_43_xx++)
	{
		auto &xx_11_xx = *xx_43_xx;
		auto xx_9_xx = data.GetEleAtIndex(xx_11_xx);
		xx_41_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_17_xx.Size()) + xx_42_xx, xx_9_xx);
		xx_42_xx++;
	}

	xx_41_xx.sendBack();
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
	List<EMCluster> xx_45_xx;
	xx_45_xx = List<EMCluster>(xx_27_xx.Size());
	int xx_46_xx = 0;
#pragma omp parallel for firstprivate(xx_46_xx) schedule(static)
	for (auto xx_47_xx = xx_27_xx.Elements().begin(); xx_47_xx < xx_27_xx.Elements().end(); xx_47_xx++)
	{
		auto &xx_21_xx = *xx_47_xx;
		auto xx_19_xx = EMCluster(seeds.GetEleAtIndex(xx_21_xx), dim, 1.0 / K);
		xx_45_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_27_xx.Size()) + xx_46_xx, xx_19_xx);
		xx_46_xx++;
	}

	xx_45_xx.sendBack();
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

	EMCluster clus;

	List<EMCluster> xx_51_xx;
	xx_51_xx = xx_51_xx.getLocalList();

#pragma omp parallel for firstprivate(clus)
	for (auto xx_52_xx = xx_51_xx.Elements().begin(); xx_52_xx < xx_51_xx.Elements().end(); xx_52_xx++)
	{
		auto &clus = *xx_52_xx;
		int b;
		int a;
		clus.oldmu = clus.mu;
		clus.wsum = 0.0;
		a = clus.musum.InitializePoint();
		b = clus.sigmasum.InitializePoint();
	}

	if (myrank == (numprocs - 1))
	{
	}
	xx_51_xx.sendBack();
}

void __ff__5()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	Point xj;
	List<EMCluster> clusList;
	List<int> klist;

	clusList = clusList.getLocalList();

	klist = klist.getLocalList();

	List<Point> xx_60_xx;
	int xx_61_xx = data.getDistributionCount(myrank, numprocs);
	int xx_62_xx = data.getLowerLimit(0, data.Size(), myrank, numprocs);
	for (int xx_63_xx = 0; xx_63_xx < xx_61_xx; xx_63_xx++)
	{
		xx_60_xx.AddEle(data[xx_62_xx + xx_63_xx]);
	}

	HashMap<double> xx_54_xx = HashMap<double>();

	HashMap<Point> xx_55_xx = HashMap<Point>();

	HashMap<Point> xx_56_xx = HashMap<Point>();

#pragma omp parallel firstprivate(xj)
	{
		HashMap<double> xx_64_xx = HashMap<double>();
		HashMap<Point> xx_65_xx = HashMap<Point>();
		HashMap<Point> xx_66_xx = HashMap<Point>();
#pragma omp for schedule(static)
		for (auto xx_67_xx = xx_60_xx.Elements().begin(); xx_67_xx < xx_60_xx.Elements().end(); xx_67_xx++)
		{
			auto &xj = *xx_67_xx;
			List<double> numerators;
			List<double> xx_34_xx;
			double denominator;
			uint indexExp;

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
				uint indexExp;
				gamma = numerators.GetEleAtIndex(i) / denominator;
				indexExp = i;
				auto xx_68_xx = xx_64_xx.htmap.find(indexExp);
				if (xx_68_xx != xx_64_xx.htmap.end())
				{
					xx_64_xx.htmap[indexExp] += gamma;
				}
				else
				{
					xx_64_xx.htmap[indexExp] = gamma;
				}
				temp = xj * gamma;
				indexExp = i;
				auto xx_69_xx = xx_65_xx.htmap.find(indexExp);
				if (xx_69_xx != xx_65_xx.htmap.end())
				{
					xx_65_xx.htmap[indexExp] += temp;
				}
				else
				{
					xx_65_xx.htmap[indexExp] = temp;
				}
				temp1 = pointSquare(xj) * gamma;
				indexExp = i;
				auto xx_70_xx = xx_66_xx.htmap.find(indexExp);
				if (xx_70_xx != xx_66_xx.htmap.end())
				{
					xx_66_xx.htmap[indexExp] += temp1;
				}
				else
				{
					xx_66_xx.htmap[indexExp] = temp1;
				}
			}
			numerators.Free();
		}
#pragma omp for ordered schedule(static, 1)
		for (int xx_71_xx = 0; xx_71_xx < omp_get_num_threads(); xx_71_xx++)
		{
#pragma omp ordered
			{
				for (auto xx_72_xx = xx_66_xx.htmap.begin(); xx_72_xx != xx_66_xx.htmap.end(); xx_72_xx++)
				{
					auto xx_73_xx = xx_56_xx.htmap.find(xx_72_xx->first);
					if (xx_73_xx != xx_56_xx.htmap.end())
					{
						xx_56_xx.htmap[xx_72_xx->first] = xx_56_xx.htmap[xx_72_xx->first] + xx_72_xx->second;
					}
					else
					{
						xx_56_xx.htmap[xx_72_xx->first] = xx_72_xx->second;
					}
				}
				for (auto xx_74_xx = xx_64_xx.htmap.begin(); xx_74_xx != xx_64_xx.htmap.end(); xx_74_xx++)
				{
					auto xx_75_xx = xx_54_xx.htmap.find(xx_74_xx->first);
					if (xx_75_xx != xx_54_xx.htmap.end())
					{
						xx_54_xx.htmap[xx_74_xx->first] = xx_54_xx.htmap[xx_74_xx->first] + xx_74_xx->second;
					}
					else
					{
						xx_54_xx.htmap[xx_74_xx->first] = xx_74_xx->second;
					}
				}
				for (auto xx_76_xx = xx_65_xx.htmap.begin(); xx_76_xx != xx_65_xx.htmap.end(); xx_76_xx++)
				{
					auto xx_77_xx = xx_55_xx.htmap.find(xx_76_xx->first);
					if (xx_77_xx != xx_55_xx.htmap.end())
					{
						xx_55_xx.htmap[xx_76_xx->first] = xx_55_xx.htmap[xx_76_xx->first] + xx_76_xx->second;
					}
					else
					{
						xx_55_xx.htmap[xx_76_xx->first] = xx_76_xx->second;
					}
				}
			}
		}
	}

	xx_54_xx.serializeSend(0);

	xx_55_xx.serializeSend(0);

	xx_56_xx.serializeSend(0);

	if (myrank == (numprocs - 1))
	{
	}
}

void __ff__6()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	EMCluster clus;

	List<EMCluster> xx_79_xx;
	xx_79_xx = xx_79_xx.getLocalList();

#pragma omp parallel for firstprivate(clus)
	for (auto xx_80_xx = xx_79_xx.Elements().begin(); xx_80_xx < xx_79_xx.Elements().end(); xx_80_xx++)
	{
		auto &clus = *xx_80_xx;
		clus.pc = clus.wsum / n;
		clus.mu = clus.musum / clus.wsum;
		clus.sigma = (clus.sigmasum / clus.wsum) - pointSquare(clus.mu);
	}

	if (myrank == (numprocs - 1))
	{
	}
	xx_79_xx.sendBack();
}

void __ff__7()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	int i;
	List<EMCluster> clusList;

	clusList = clusList.getLocalList();

	List<int> xx_82_xx;
	xx_82_xx = xx_82_xx.getLocalList();

#pragma omp parallel for firstprivate(i)
	for (auto xx_83_xx = xx_82_xx.Elements().begin(); xx_83_xx < xx_82_xx.Elements().end(); xx_83_xx++)
	{
		auto &i = *xx_83_xx;
		Point s;
		Point p;
		double c;
		p = clusList.GetEleAtIndex(i).mu;
		s = clusList.GetEleAtIndex(i).sigma;
		c = clusList.GetEleAtIndex(i).pc;
		std::cout << "Cluster centre: " << i << "\t" << p << "\t\tsigma: " << s << "\tpc: " << c << std::endl;
	}

	xx_82_xx.sendBack();
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
		default:
			cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}