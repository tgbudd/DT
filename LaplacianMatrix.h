#ifndef LAPLACIAN_MATRIX_H
#define LAPLACIAN_MATRIX_H

#include <vector>
#include "boost/array.hpp"

#include "Triangulation.h"
#include "LinearAlgebra.h"
#include "Vertex.h"
#include "Edge.h"
#include "utilities.h"

class LaplacianMatrix : public linearalgebra::Matrix
{
public:
	LaplacianMatrix (const Triangulation * const triangulation);
	LaplacianMatrix (const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure);
	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
private:
	void Initialize(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure);
	std::vector<std::map<int,double> > laplacianRules_;
};

class DomainLaplacianMatrix : public linearalgebra::Matrix
{
public:
	DomainLaplacianMatrix (const Triangulation * const triangulation, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed);
	DomainLaplacianMatrix (const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed);
	void MultiplyVector(const std::vector<double> & from, std::vector<double> & to) const;
	bool FindHarmonic(std::vector<double> & x, double eps, int maxIterations = 100 );
private:
	void Initialize(const Triangulation * const triangulation, const std::vector<boost::array<double,3> > & edge_measure, const std::vector<const Vertex*> & vertices, const std::vector<const Vertex*> & fixed);
	void GetTarget(const std::vector<double> & x, std::vector<double> & target);
	std::vector<std::map<int,double> > laplacianRules_;
	std::vector<std::map<int,double> > targetRules_;
	const std::vector<const Vertex*> vertices_;
};


#endif