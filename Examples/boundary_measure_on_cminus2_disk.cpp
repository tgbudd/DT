#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "LinearAlgebra.h"
#include "BoundaryMeasure.h"

class CMinusTwoBoundaryMeasure : public BoundaryMeasure
{
public:
	CMinusTwoBoundaryMeasure(Triangulation * triangulation, CMinusTwoBuilder * builder);
	void DetermineBoundary();
private:
	CMinusTwoBuilder * builder_;
	Triangulation * triangulation_;
};

CMinusTwoBoundaryMeasure::CMinusTwoBoundaryMeasure(Triangulation * triangulation, CMinusTwoBuilder * builder)
	: BoundaryMeasure(triangulation),
	builder_(builder),
	triangulation_(triangulation)
{
}

void CMinusTwoBoundaryMeasure::DetermineBoundary()
{
	Vertex * v0 = triangulation_->getRandomVertex();
	Vertex * v1;
	do {
		v1 = triangulation_->getRandomVertex();
	} while( v0==v1 );
	std::list<const Edge *> path;
	builder_->getPath(v0,v1,path);
	boundary_ = path;
	for(std::list<const Edge *>::const_reverse_iterator it = path.rbegin();it!=path.rend();it++)
	{
		boundary_.push_back((*it)->getAdjacent());
	}
	//builder_->getDiskBoundary(boundary_);
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	int b = 0;//param.Read<int>("boundary length");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n,b);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	CMinusTwoBoundaryMeasure bmeasure(&triangulation,&builder);

	Simulation simulation( &triangulation, 0, secondsperoutput, output );
	simulation.AddObservable( &bmeasure, 1 );

	simulation.Run();

	return 0;
}