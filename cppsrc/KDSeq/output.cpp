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
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;
	File ff;
	Point p;
	List<KmeansCluster> CL;
	List<Point> data;
	double threshold;
	List<Point> xx_1_xx;
	double error;
	int iter;
	int K;
	List<Point> seeds;
	List<KmeansCluster> xx_4_xx;
	int xy;
	try
	{
		ff.ReadDataset("datasets/3droad.arff", &data);
		K = 100;

		List<int> xx_2_xx;
		for (int xx_3_xx = 0; xx_3_xx <= K - 1; xx_3_xx++)
		{
			xx_2_xx.AddEle(xx_3_xx);
		}
		for (auto &i : xx_2_xx.Elements())
		{
			xx_1_xx.AddEle(data.GetEleAtIndex(i));
		}
		seeds = xx_1_xx;

		List<int> xx_5_xx;
		for (int xx_6_xx = 0; xx_6_xx <= K - 1; xx_6_xx++)
		{
			xx_5_xx.AddEle(xx_6_xx);
		}
		for (auto &i : xx_5_xx.Elements())
		{
			xx_4_xx.AddEle(KmeansCluster(seeds.GetEleAtIndex(i)));
		}
		CL = xx_4_xx;
		threshold = 0.5;
		std::cout << "Number of points:" << data.Size() << "\nK :" << K << "\nError threshold :" << threshold << std::endl;
		std::cout << "\nseeds:" << seeds << std::endl;
		error = 100.0;
		iter = 0;

		while (error > threshold)
		{
			for (auto &clus : CL.Elements())
			{
				int a;
				clus.count = 0;
				clus.oldRep = clus.rep;
				a = clus.rep.InitializePoint();
			}
			for (auto &i : data.Elements())
			{
				double xx_8_xx;
				int xx_10_xx;
				double xx_7_xx;
				int index;
				int xx_9_xx;
				xx_7_xx = 0.0;
				xx_8_xx = DBL_MAX;
				xx_9_xx = -1;
				xx_10_xx = 0;
				for (auto &CC : CL.Elements())
				{
					xx_7_xx = distanceEuclidean(i, CC.oldRep);
					if (xx_7_xx < xx_8_xx)
					{
						KmeansCluster CC;
						xx_8_xx = xx_7_xx;
						xx_9_xx = xx_10_xx;
					}
					xx_10_xx = xx_10_xx + 1;
				}
				index = xx_9_xx;
				CL.GetEleAtIndex(index).rep = CL.GetEleAtIndex(index).rep + i;
				CL.GetEleAtIndex(index).count = CL.GetEleAtIndex(index).count + 1;
			}
			for (auto &clus : CL.Elements())
			{
				Point sum;
				sum = clus.rep;
				if (clus.count != 0)
				{
					clus.rep = sum / clus.count;
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
	int xx_11_xx = -1;
	MPI_Bcast(&xx_11_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
