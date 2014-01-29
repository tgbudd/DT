#ifndef DISK_CIRCLE_PACKING_H
#define DISK_CIRCLE_PACKING_H

#include "boost/array.hpp"

#include "Triangulation.h"
#include "LinearAlgebra.h"
#include "Edge.h"
#include "Vertex.h"
#include "BabyUniverseDetector.h"

class DiskCirclePacking
{
public:
	DiskCirclePacking(Triangulation * const triangulation);
	~DiskCirclePacking() {}
	bool FindEmbedding(const std::list<const Edge*> & boundary, const Edge * centerEdge);
	void getCircles(std::vector<std::pair<Vertex*,std::pair<Vector2D,double> > > & circles);
private:
	bool FindDiskRadii(const std::list<const Edge*> & boundary);
	bool DiskLayout(const std::list<const Edge*> & boundary, const Edge * centerEdge);

	std::vector<Triangle *> disk_triangles_;
	std::vector<Vertex *> disk_vertices_;

	double Angle(Edge * edge);
	double AngleSum(int vertexId);

	const Triangulation * const triangulation_;
	BabyUniverseDetector babyuniversedetector_;

	std::vector<int> degree_;
	std::vector<double> radius_;
	std::vector<boost::array<double,3> > angles_;
	std::vector<double> euclidean_radius_;
	std::vector<Vector2D> coordinate_;
	std::vector<Vector2D> hyp_coordinate_;
	std::vector<bool> circle_fixed_;
	std::vector<int> circle_order_;

	int max_iterations_;
	double epsilon_, delta_;
};


#endif