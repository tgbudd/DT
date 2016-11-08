#include <iostream>
#include <sstream>
#include <iomanip>

#include "Triangulation.h"
#include "Vertex.h"
#include "Simulation.h"
#include "CMinusTwoBuilder.h"
#include "utilities.h"
#include "Histogram.h"

class DegreeDistribution : public Observable 
{
public:
	DegreeDistribution(Triangulation & triangulation, int maxdegree);
	void Measure();
	std::string OutputData() const;
private:
	Triangulation & triangulation_;
	Histogram<int> degree_histogram_;
};

DegreeDistribution::DegreeDistribution(Triangulation & triangulation, int maxdegree) :
	triangulation_(triangulation), 
	degree_histogram_(0,maxdegree,maxdegree)
{
}

void DegreeDistribution::Measure()
{
	for(int i=0;i<triangulation_.NumberOfVertices();++i)
	{
		Vertex * vertex = triangulation_.getVertex(i);
		int degree = vertex->getDegree();
		degree_histogram_.Insert(degree);
	}
}

std::string DegreeDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << "degreedistribution -> ";
	degree_histogram_.PrintTo(stream);
	return stream.str();
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;
	
	int numtriangles = param.Read<int>("number of triangles");

	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	{
		// Use the c=-2 triangulation builder to produce a 
		// decent initial configuration
		CMinusTwoBuilder builder(&triangulation,0,numtriangles);
		triangulation.setDominantMatter( &builder );
		triangulation.DoSweep();
		triangulation.clearDominantMatter();
	}

	// The Simulation object takes care of producing thermalized samples of the triangulation
	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	DegreeDistribution degreedistribution(triangulation,100);
	
	// Make sure the observable degreedistribution is called for each sample (separated
	// by measurementSweeps sweeps of triangle flips).
	simulation.AddObservable( &degreedistribution, measurementSweeps );

	// directory for output
	simulation.SetDirectory("./");

	// Start simulation which will run indefinitely.
	simulation.Run();
	
	return 0;
}