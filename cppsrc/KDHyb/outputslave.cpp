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
		double y;
		List<double> xx_0_xx;

		if (CL.Size() >= 5000)
		{
			int xx_28_xx = 0;
			MPI_Bcast(&xx_28_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			CL.distributeList(numprocs);
			xx_0_xx.gatherList(numprocs);
		}
		else
		{
			for (auto &xx_3_xx : CL.Elements())
			{
				auto xx_1_xx = distanceEuclidean(xx_3_xx.oldRep, xx_3_xx.rep);
				xx_0_xx.AddEle(xx_1_xx);
			}
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

void __ff__0()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<KmeansCluster> xx_2_xx;
	xx_2_xx = xx_2_xx.getLocalList();
	List<double> xx_29_xx;
	xx_29_xx = List<double>(xx_2_xx.Size());
	int xx_30_xx = 0;
#pragma omp parallel for firstprivate(xx_30_xx) schedule(static)
	for (auto xx_31_xx = xx_2_xx.Elements().begin(); xx_31_xx < xx_2_xx.Elements().end(); xx_31_xx++)
	{
		auto &xx_3_xx = *xx_31_xx;
		auto xx_1_xx = distanceEuclidean(xx_3_xx.oldRep, xx_3_xx.rep);
		xx_29_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_2_xx.Size()) + xx_30_xx, xx_1_xx);
		xx_30_xx++;
	}

	xx_29_xx.sendBack();
}

void __ff__1()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<Point> data;

	data = ::data;
	int xx_11_xx;
	int xx_12_xx;
	MPI_Recv(&xx_11_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_12_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_13_xx;
	xx_13_xx.makeRangeList(xx_11_xx, xx_12_xx);
	List<Point> xx_34_xx;
	xx_34_xx = List<Point>(xx_13_xx.Size());
	int xx_35_xx = 0;
#pragma omp parallel for firstprivate(xx_35_xx) schedule(static)
	for (auto xx_36_xx = xx_13_xx.Elements().begin(); xx_36_xx < xx_13_xx.Elements().end(); xx_36_xx++)
	{
		auto &xx_7_xx = *xx_36_xx;
		auto xx_5_xx = data.GetEleAtIndex(xx_7_xx);
		xx_34_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_13_xx.Size()) + xx_35_xx, xx_5_xx);
		xx_35_xx++;
	}

	xx_34_xx.sendBack();
}

void __ff__2()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	List<Point> seeds;

	seeds = seeds.getLocalList();
	int xx_21_xx;
	int xx_22_xx;
	MPI_Recv(&xx_21_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&xx_22_xx, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	List<int> xx_23_xx;
	xx_23_xx.makeRangeList(xx_21_xx, xx_22_xx);
	List<KmeansCluster> xx_38_xx;
	xx_38_xx = List<KmeansCluster>(xx_23_xx.Size());
	int xx_39_xx = 0;
#pragma omp parallel for firstprivate(xx_39_xx) schedule(static)
	for (auto xx_40_xx = xx_23_xx.Elements().begin(); xx_40_xx < xx_23_xx.Elements().end(); xx_40_xx++)
	{
		auto &xx_17_xx = *xx_40_xx;
		auto xx_15_xx = KmeansCluster(seeds.GetEleAtIndex(xx_17_xx));
		xx_38_xx.SetEleAt(int_ceil(omp_get_thread_num(), omp_get_num_threads(), xx_23_xx.Size()) + xx_39_xx, xx_15_xx);
		xx_39_xx++;
	}

	xx_38_xx.sendBack();
}

void __ff__3()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	KmeansCluster clus;

	List<KmeansCluster> xx_42_xx;
	xx_42_xx = xx_42_xx.getLocalList();

#pragma omp parallel for firstprivate(clus)
	for (auto xx_43_xx = xx_42_xx.Elements().begin(); xx_43_xx < xx_42_xx.Elements().end(); xx_43_xx++)
	{
		auto &clus = *xx_43_xx;
		int a;
		clus.count = 0;
		clus.oldRep = clus.rep;
		a = clus.rep.InitializePoint();
	}

	if (myrank == (numprocs - 1))
	{
	}
	xx_42_xx.sendBack();
}

void __ff__4()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	Point i;
	List<KmeansCluster> CL;

	CL = CL.getLocalList();

	List<Point> xx_49_xx;
	int xx_50_xx = data.getDistributionCount(myrank, numprocs);
	int xx_51_xx = data.getLowerLimit(0, data.Size(), myrank, numprocs);
	for (int xx_52_xx = 0; xx_52_xx < xx_50_xx; xx_52_xx++)
	{
		xx_49_xx.AddEle(data[xx_51_xx + xx_52_xx]);
	}

	HashMap<Point> xx_45_xx = HashMap<Point>();

	HashMap<int> xx_46_xx = HashMap<int>();

#pragma omp parallel firstprivate(i)
	{
		HashMap<Point> xx_53_xx = HashMap<Point>();
		HashMap<int> xx_54_xx = HashMap<int>();
#pragma omp for schedule(static)
		for (auto xx_55_xx = xx_49_xx.Elements().begin(); xx_55_xx < xx_49_xx.Elements().end(); xx_55_xx++)
		{
			auto &i = *xx_55_xx;
			int xx_26_xx;
			int index;
			double xx_25_xx;
			double xx_24_xx;
			int xx_27_xx;
			uint indexExp;
			xx_24_xx = 0.0;
			xx_25_xx = DBL_MAX;
			xx_26_xx = -1;
			xx_27_xx = 0;
			for (auto &CC : CL.Elements())
			{
				xx_24_xx = distanceEuclidean(i, CC.oldRep);
				if (xx_24_xx < xx_25_xx)
				{
					KmeansCluster CC;
					xx_25_xx = xx_24_xx;
					xx_26_xx = xx_27_xx;
				}
				xx_27_xx = xx_27_xx + 1;
			}
			index = xx_26_xx;
			indexExp = index;
			auto xx_56_xx = xx_53_xx.htmap.find(indexExp);
			if (xx_56_xx != xx_53_xx.htmap.end())
			{
				xx_53_xx.htmap[indexExp] += i;
			}
			else
			{
				xx_53_xx.htmap[indexExp] = i;
			}
			indexExp = index;
			auto xx_57_xx = xx_54_xx.htmap.find(indexExp);
			if (xx_57_xx != xx_54_xx.htmap.end())
			{
				xx_54_xx.htmap[indexExp] += 1;
			}
			else
			{
				xx_54_xx.htmap[indexExp] = 1;
			}
		}
#pragma omp for ordered schedule(static, 1)
		for (int xx_58_xx = 0; xx_58_xx < omp_get_num_threads(); xx_58_xx++)
		{
#pragma omp ordered
			{
				for (auto xx_59_xx = xx_53_xx.htmap.begin(); xx_59_xx != xx_53_xx.htmap.end(); xx_59_xx++)
				{
					auto xx_60_xx = xx_45_xx.htmap.find(xx_59_xx->first);
					if (xx_60_xx != xx_45_xx.htmap.end())
					{
						xx_45_xx.htmap[xx_59_xx->first] = xx_45_xx.htmap[xx_59_xx->first] + xx_59_xx->second;
					}
					else
					{
						xx_45_xx.htmap[xx_59_xx->first] = xx_59_xx->second;
					}
				}
				for (auto xx_61_xx = xx_54_xx.htmap.begin(); xx_61_xx != xx_54_xx.htmap.end(); xx_61_xx++)
				{
					auto xx_62_xx = xx_46_xx.htmap.find(xx_61_xx->first);
					if (xx_62_xx != xx_46_xx.htmap.end())
					{
						xx_46_xx.htmap[xx_61_xx->first] = xx_46_xx.htmap[xx_61_xx->first] + xx_61_xx->second;
					}
					else
					{
						xx_46_xx.htmap[xx_61_xx->first] = xx_61_xx->second;
					}
				}
			}
		}
	}

	xx_45_xx.serializeSend(0);

	xx_46_xx.serializeSend(0);
}

void __ff__5()
{
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

	KmeansCluster clus;

	List<KmeansCluster> xx_64_xx;
	xx_64_xx = xx_64_xx.getLocalList();

#pragma omp parallel for firstprivate(clus)
	for (auto xx_65_xx = xx_64_xx.Elements().begin(); xx_65_xx < xx_64_xx.Elements().end(); xx_65_xx++)
	{
		auto &clus = *xx_65_xx;
		Point sum;
		sum = clus.rep;
		if (clus.count != 0)
		{
			clus.rep = sum / clus.count;
		}
	}

	xx_64_xx.sendBack();
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
		default:
			cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}