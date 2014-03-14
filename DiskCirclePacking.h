#ifndef DISK_CIRCLE_PACKING_H
#define DISK_CIRCLE_PACKING_H

#include "boost/array.hpp"

#include <Eigen/Sparse>

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
	void getHyperbolicCoordinates(std::vector<std::pair<Vertex*,Vector2D> > & coor) const;
	void getBoundaryPositions(std::vector<double> & angles);
	double getCenterRadius() const;
	int getStepsUsed() const {
		return steps_;
	}
	std::vector<Triangle *> getDiskTriangles()
	{
		return disk_triangles_;
	}
private:
	bool FindDiskRadii(const std::list<const Edge*> & boundary);
	bool DiskLayout(const std::list<const Edge*> & boundary, const Edge * centerEdge);
	bool FindFlatEmbedding();
	void LayoutBoundary(std::vector<double> & radius, std::vector<Vector2D> & coordinate, int bLength); 
	void scaleRadii(std::vector<double> & radius, int bLength);

	std::vector<Triangle *> disk_triangles_;
	std::vector<Vertex *> disk_vertices_;

	std::vector<Eigen::Triplet<double> > lapl_rules_;
	std::vector<Eigen::Triplet<double> > boundary_rules_;
	std::vector<double> vertex_weight_;

	std::vector<boost::array<int,4> > edge_to_vert_;

	double Angle(Edge * edge);
	double AngleSum(int vertexId);
	double Angle(Edge * edge, const std::vector<int> & vertexPos, const std::vector<double> & radius);

	double boundaryAngle(std::vector<double> & radius, int bLength, double r);
	double boundaryAngleDerivative(std::vector<double> & radius, int bLength, double r);

	const Triangulation * const triangulation_;
	BabyUniverseDetector babyuniversedetector_;
	
	std::list<const Edge*> boundary_;
	const Edge * center_edge_;

	bool use_one_minus_radius_;
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
	bool reset_radius_;

	int steps_;
};


#endif