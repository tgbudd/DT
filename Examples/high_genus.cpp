#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <set>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "SimplexAttribute.h"
#include "TriangulationProperties.h"
#include "Histogram.h"

class TriangulationSaver : public Observable {
public:
	TriangulationSaver(const Triangulation * triangulation, std::string filenameStart)
		: triangulation_(triangulation), filename_start_(filenameStart), index_(0)
	{}
	void Measure()
	{
		std::ostringstream filename;
		filename << filename_start_ << index_ << ".txt";
		index_++;
		std::ofstream file(filename.str());

		file << "{ \"genus\" -> " << triangulation_->CalculateGenus();
		file << ", \"tri\" -> {{";
		for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;++i)
		{
			file << (i>0?"},{":"");
			const Triangle * triangle = triangulation_->getTriangle(i);
			for(int j=0;j<3;++j)
			{
				int vertId = triangle->getEdge(j)->getOpposite()->getId()+1;
				file << (j>0?",":"") << vertId;
			}
		}
		file << "}}, \"adj\" -> {{";
		for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;++i)
		{
			file << (i>0?"},{":"");
			const Triangle * triangle = triangulation_->getTriangle(i);
			for(int j=0;j<3;++j)
			{
				int triId = triangle->getEdge(j)->getAdjacent()->getParent()->getId()+1;
				int edgeId = triangle->getEdge(j)->getAdjacent()->getId()+1;
				file << (j>0?",":"") << "{" << triId << "," << edgeId << "}";
			}
		}
		file << "}}}\n";
	}
	std::string OutputData()
	{
		return "{}";
	}
private:
	const Triangulation * triangulation_;
	std::string filename_start_;
	int index_;
};

void IncreaseGenus(Triangulation * triangulation, int genus)
{
	int currentGenus=0;
	while(currentGenus < genus)
	{
		Edge * edge1 = triangulation->getRandomEdge();
		Edge * edge2 = triangulation->getRandomEdge();
		if( edge1->getEnd() == edge1->getStart() ||
			edge2->getEnd() == edge2->getStart() ||
			edge1->getEnd() == edge2->getEnd() ||
			edge1->getStart() == edge2->getEnd() ||
			edge1->getEnd() == edge2->getStart() ||
			edge1->getStart() == edge2->getStart() )
		{
			continue;
		}

		Edge * edge1adj = edge1->getAdjacent();
		Edge * edge2adj = edge2->getAdjacent();

		edge1->bindAdjacent(edge2);
		edge1adj->bindAdjacent(edge2adj);

		triangulation->DetermineVertices();
		currentGenus++;

		BOOST_ASSERT( currentGenus == triangulation->CalculateGenus() );
	}
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;

	int n = param.Read<int>("triangles");
	int genus = param.Read<int>("genus");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = 1000000;
	bool output = param.UserInput();

	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter( &builder );
	triangulation.DoSweep();
	triangulation.clearDominantMatter();
	
	IncreaseGenus(&triangulation,genus);

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	std::ostringstream os;
	os << "./output/tri/genus" << genus << "-" << n << "-"; 
	TriangulationSaver saver( &triangulation, os.str() );
	
	simulation.AddObservable( &saver, measurementSweeps );
	simulation.SetDirectory("./output/");
	simulation.Run();
	return 0;
}