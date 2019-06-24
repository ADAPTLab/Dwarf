#include "output.h"
#include <iostream>
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

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	File ff;
	Point p;
	List<KmeansCluster> CL;
	List<Point> data;
	double threshold;
	double error;
	int iter;
	int K;
	List<Point> seeds;
	int xy;
	try
	{
		ff.ReadDataset("datasets/3droad.arff", &data);
		int xx_32_xx = -2;
		MPI_Bcast(&xx_32_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		data.broadcastList(numprocs);
		K = 100;
		List<Point> xx_4_xx;

		if ((K - 1) - (0) >= 5000)
		{
			int xx_33_xx = 1;
			MPI_Bcast(&xx_33_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			List<int> xx_8_xx;
			int xx_9_xx = 0;
			int xx_10_xx = K - 1;
			xx_8_xx.distributeRangeList(xx_9_xx, xx_10_xx, numprocs);

			xx_4_xx.gatherList(numprocs);
		}
		else
		{
			for (int xx_7_xx = 0; xx_7_xx <= K - 1; xx_7_xx++)
			{
				auto xx_5_xx = data.GetEleAtIndex(xx_7_xx);
				xx_4_xx.AddEle(xx_5_xx);
			}
		}

		seeds = xx_4_xx;
		List<KmeansCluster> xx_14_xx;

		if ((K - 1) - (0) >= 5000)
		{
			int xx_37_xx = 2;
			MPI_Bcast(&xx_37_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

			List<int> xx_18_xx;
			seeds.broadcastList(numprocs);
			int xx_19_xx = 0;
			int xx_20_xx = K - 1;
			xx_18_xx.distributeRangeList(xx_19_xx, xx_20_xx, numprocs);

			xx_14_xx.gatherList(numprocs);
		}
		else
		{
			for (int xx_17_xx = 0; xx_17_xx <= K - 1; xx_17_xx++)
			{
				auto xx_15_xx = KmeansCluster(seeds.GetEleAtIndex(xx_17_xx));
				xx_14_xx.AddEle(xx_15_xx);
			}
		}

		CL = xx_14_xx;
		threshold = 0.5;
		std::cout << "Number of points:" << data.Size() << "\nK :" << K << "\nError threshold :" << threshold << std::endl;
		std::cout << "\nseeds:" << seeds << std::endl;
		error = 100.0;
		iter = 0;

		while (error > threshold)
		{

			if (CL.Size() >= 5000)
			{
				int xx_41_xx = 3;
				MPI_Bcast(&xx_41_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				CL.distributeList(numprocs);
				CL.gatherList(numprocs);
			}
			else
			{
				for (auto &clus : CL.Elements())
				{
					int a;
					clus.count = 0;
					clus.oldRep = clus.rep;
					a = clus.rep.InitializePoint();
				}
			}

			if (data.Size() >= 5000)
			{
				int xx_44_xx = 4;
				MPI_Bcast(&xx_44_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				CL.broadcastList(numprocs);

				HashMap<Point> xx_45_xx = HashMap<Point>();

				HashMap<int> xx_46_xx = HashMap<int>();
				for (int xx_47_xx = 1; xx_47_xx < numprocs; xx_47_xx++)
				{
					xx_45_xx = xx_45_xx.deserializeRecv(xx_47_xx);
					for (auto iterator = xx_45_xx.htmap.begin(); iterator != xx_45_xx.htmap.end(); iterator++)
					{
						uint index = iterator->first;
						CL[index].rep += xx_45_xx.htmap[index];
					}
				}
				for (int xx_48_xx = 1; xx_48_xx < numprocs; xx_48_xx++)
				{
					xx_46_xx = xx_46_xx.deserializeRecv(xx_48_xx);
					for (auto iterator = xx_46_xx.htmap.begin(); iterator != xx_46_xx.htmap.end(); iterator++)
					{
						uint index = iterator->first;
						CL[index].count += xx_46_xx.htmap[index];
					}
				}
			}
			else
			{
				for (auto &i : data.Elements())
				{
					int xx_26_xx;
					int index;
					double xx_25_xx;
					double xx_24_xx;
					int xx_27_xx;
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
					CL.GetEleAtIndex(index).rep = CL.GetEleAtIndex(index).rep + i;
					CL.GetEleAtIndex(index).count = CL.GetEleAtIndex(index).count + 1;
				}
			}

			if (CL.Size() >= 5000)
			{
				int xx_63_xx = 5;
				MPI_Bcast(&xx_63_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
				CL.distributeList(numprocs);
				CL.gatherList(numprocs);
			}
			else
			{
				for (auto &clus : CL.Elements())
				{
					Point sum;
					sum = clus.rep;
					if (clus.count != 0)
					{
						clus.rep = sum / clus.count;
					}
				}
			}

			error = ConvergenceCriteria(CL);
		}
		std::cout << "\nClustering Results: " << std::endl;
		xy = 0;
		for (auto &c : CL.Elements())
		{
			xy = xy + 1;
			std::cout << "\nCluster: " << xy << ", Members: " << c.count << ", Centroid:" << c.rep << std::endl;
		}
	}
	catch (std::exception &e)
	{
		std::cout << std::string(e.what()) + " At line no. 0, function" + "" << std::endl;
	}
	int xx_66_xx = -1;
	MPI_Bcast(&xx_66_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
