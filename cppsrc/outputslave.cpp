#include "output.h"
#include <iostream>

List<Point> data;

KmeansCluster::KmeansCluster(Point p) {
	rep = p;
	oldRep = p;
	count = 0;
}
double ConvergenceCriteria(List<KmeansCluster> CL) {
try {
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int numprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Status status;
	List<double> x;
	double y;
List<double> xx_0_xx;

	if (CL.Size() >= 5000) {
		int xx_28_xx = 0;
		MPI_Bcast(&xx_28_xx, 1, MPI_INT, 0, MPI_COMM_WORLD);

	CL.distributeList(numprocs);
	xx_0_xx.gatherList(numprocs);
	} else {
		for (auto &xx_3_xx : CL.Elements()) {
			auto xx_1_xx = distanceEuclidean(xx_3_xx.oldRep, xx_3_xx.rep);
			xx_0_xx.AddEle(xx_1_xx);
		}
	}


	x = xx_0_xx;
	y = *(double *)(Reduce(SUM, x.Elements(), true));
	return (y / CL.Size());
}
catch(std::exception &e) {
	 std::cout << std::string(e.what()) +" At line no. 12, function "
	 + "ConvergenceCriteria in file dwarf_source_codes\/kmeans.dw" << std::endl;
}
}

void __ff__0(){
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;

List<KmeansCluster> xx_2_xx;
	xx_2_xx = xx_2_xx.getLocalList();
	List< double > xx_29_xx;
	for(auto &xx_3_xx : xx_2_xx.Elements()){
		auto xx_1_xx = distanceEuclidean(xx_3_xx.oldRep, xx_3_xx.rep);
		xx_29_xx.AddEle(xx_1_xx);
	}

	xx_29_xx.sendBack();

}

void __ff__1(){
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
	List< Point > xx_33_xx;
	for(auto &xx_7_xx : xx_13_xx.Elements()){
		auto xx_5_xx = data.GetEleAtIndex(xx_7_xx);
		xx_33_xx.AddEle(xx_5_xx);
	}

	xx_33_xx.sendBack();

}

void __ff__2(){
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
	List< KmeansCluster > xx_36_xx;
	for(auto &xx_17_xx : xx_23_xx.Elements()){
		auto xx_15_xx = KmeansCluster(seeds.GetEleAtIndex(xx_17_xx));
		xx_36_xx.AddEle(xx_15_xx);
	}

	xx_36_xx.sendBack();

}

void __ff__3(){
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;



	KmeansCluster clus;




List<KmeansCluster> xx_39_xx;
xx_39_xx = xx_39_xx.getLocalList();
for( auto &clus : xx_39_xx.Elements()) {
		int a;
	clus.count = 0;
	clus.oldRep = clus.rep;
	a = clus.rep.InitializePoint();
	}
	 
if (myrank == (numprocs - 1)){
}
	xx_39_xx.sendBack();

}

void __ff__4(){
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;



	Point i;
		List<KmeansCluster> CL;


CL=CL.getLocalList();



List<Point> xx_45_xx;
int xx_46_xx = data.getDistributionCount(myrank, numprocs);
int xx_47_xx = data.getLowerLimit(0, data.Size(), myrank, numprocs);
for (int xx_48_xx = 0; xx_48_xx < xx_46_xx; xx_48_xx++) {
	xx_45_xx.AddEle(data[xx_47_xx + xx_48_xx]);
}

 HashMap<Point> xx_41_xx= HashMap<Point>();

 HashMap<int> xx_42_xx= HashMap<int>();
for( auto &i : xx_45_xx.Elements()) {
		int xx_26_xx;
	int index;
	double xx_25_xx;
	double xx_24_xx;
	int xx_27_xx;
uint indexExp ;
	xx_24_xx = 0.0;
	xx_25_xx = DBL_MAX;
	xx_26_xx = -1;
	xx_27_xx = 0;
for( auto &CC : CL.Elements()) {
		xx_24_xx = distanceEuclidean(i, CC.oldRep);
	if(xx_24_xx < xx_25_xx) {
	KmeansCluster CC;
	xx_25_xx = xx_24_xx;
	xx_26_xx = xx_27_xx;
	}
	xx_27_xx = xx_27_xx + 1;
	}
	index = xx_26_xx;
indexExp = index;
auto xx_49_xx = xx_41_xx.htmap.find(indexExp);
if(xx_49_xx != xx_41_xx.htmap.end()){
xx_41_xx.htmap[indexExp] += i;
 }
else{
xx_41_xx.htmap[indexExp] = i;
 }
indexExp = index;
auto xx_50_xx = xx_42_xx.htmap.find(indexExp);
if(xx_50_xx != xx_42_xx.htmap.end()){
xx_42_xx.htmap[indexExp] += 1;
 }
else{
xx_42_xx.htmap[indexExp] = 1;
 }
	}

 xx_41_xx.serializeSend(0);

 xx_42_xx.serializeSend(0);
	 

}

void __ff__5(){
	int numprocs, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Status status;



	KmeansCluster clus;




List<KmeansCluster> xx_52_xx;
xx_52_xx = xx_52_xx.getLocalList();
for( auto &clus : xx_52_xx.Elements()) {
		Point sum;
	sum = clus.rep;
	if(clus.count != 0) {
	clus.rep = sum / clus.count;
	}
	}
	 
	xx_52_xx.sendBack();

}
int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int signal = 0;
	while(signal != -1){
		MPI_Bcast(&signal, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if(signal == -1) break;
		switch(signal){
			case -2 : data = data.getLocalList();
			break;
			case 0 : __ff__0();
			break;
			case 1 : __ff__1();
			break;
			case 2 : __ff__2();
			break;
			case 3 : __ff__3();
			break;
			case 4 : __ff__4();
			break;
			case 5 : __ff__5();
			break;
			default : cout << "Unrecognized signal received: " << signal;
			break;
		}
	}
	MPI_Finalize();
	return 0;
}