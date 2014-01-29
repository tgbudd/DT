#ifndef CIRCLE_PACKING_H
#define CIRCLE_PACKING_H

#include "boost/array.hpp"

#include "Embedding.h"
#include "CohomologyBasis.h"
#include "LinearAlgebra.h"
#include "BabyUniverseDetector.h"
#include "LaplacianMatrix.h"

class CirclePacking :
	public Embedding
{
public:
	CirclePacking(Triangulation * const triangulation, CohomologyBasis * const cohomologybasis);
	~CirclePacking() {}

	bool FindEdgeMeasure();
	void CalculateEdgeMeasure(std::vector<boost::array<double,3> > & measure);
	bool FindRadii();
	void RadiiToAngles();

	bool GetRadii(std::vector<double> & radii);
	bool FindDiskEmbedding(const std::list<const Edge*> & boundary, const Edge * centerEdge);
	void getCircles(std::vector<std::pair<Vector2D,double> > & circles);
private:
	bool FindDiskRadii(const std::list<const Edge*> & boundary);
	bool DiskLayout(const std::list<const Edge*> & boundary);
	void MobiusTransformation(const Edge * centerEdge);
	std::vector<Triangle *> disk_triangles_;
	std::vector<Vertex *> disk_vertices_;
	std::vector<boost::array<double,3> > disk_edge_measure_;

	std::pair<Vector2D,double> CircleMobius(Vector2D c, double radius, Vector2D z0, double angle);

//	void RadiiToAngles();
	double AngleSum(int vertexId);
	double Angle(Edge * edge);
	const Triangulation * const triangulation_;

	BabyUniverseDetector babyuniversedetector_;
	
	std::vector<int> degree_;
	std::vector<double> radius_;
	std::vector<boost::array<double,3> > angles_;
	std::vector<Vector2D> disk_coordinate_;


	int max_iterations_;
	double epsilon_, delta_;
};

#endif
