#include <string>
#include <vector>
#include <cfloat>
#include <omp.h>
#include "Point.h"
#include "List.h"
#include "HashMap.h"
#include "src/graph"
#include "Library.h"
#include "mpi.h"

#include "RTree.h"
#include "kdtree++/kdtree_wrapper.hpp"
#include "clustering.h"

class EMCluster;

class EMCluster
{
public:
	Point mu;
	Point oldmu;
	Point musum;
	double pc;
	double wsum;
	Point sigma;
	Point sigmasum;
	EMCluster() {}
	EMCluster(Point p, int dim, double w);
	template <class Archive>
	void serialize(Archive &ar)
	{
		ar(mu, oldmu, musum, pc, wsum, sigma, sigmasum);
	}
};

double GMM(Point x, Point mu, Point sigma, double prob);
double convergence(List<EMCluster> clusList);
double distanceEuclideanSquared(Point A, Point B);
Point pointSquare(Point A);
