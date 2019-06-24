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

class KmeansCluster;

class KmeansCluster
{
public:
	Point rep;
	Point oldRep;
	int count;
	KmeansCluster() {}
	KmeansCluster(Point p);
	template <class Archive>
	void serialize(Archive &ar)
	{
		ar(rep, oldRep, count);
	}
};

double ConvergenceCriteria(List<KmeansCluster> CL);
