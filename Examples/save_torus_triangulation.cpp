#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"

int ProperFloor(double x)
{
	return static_cast<int>(std::floor(x));
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	std::string filename = param.Read<std::string>("file");

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,1,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding embedding( &triangulation, &cohom );
	embedding.SetAccuracy( 1.0e-8 );
	
	embedding.MakeUpToDate();

	std::pair<double,double> modulus = embedding.CalculateModuli();

	std::ofstream file(filename.c_str());

	file << triangulation.NumberOfTriangles() << "\n";
	file << std::fixed << modulus.first << " " << modulus.second << "\n";

	for(int i=0;i<triangulation.NumberOfTriangles();i++)
	{
		Triangle * t = triangulation.getTriangle(i);
		Vector2D center = embedding.getCoordinate(t->getEdge(0)->getOpposite());
		center = AddVectors2D(center,ScaleVector2D(embedding.getForm(i,2),1/3.0));
		center = AddVectors2D(center,ScaleVector2D(embedding.getForm(i,1),-1/3.0));

		Vector2D shift = {-std::floor(center[0]),-std::floor(center[1])};

		file << t->getEdge(0)->getOpposite()->getId() << " " << t->getEdge(1)->getOpposite()->getId() << " " << t->getEdge(2)->getOpposite()->getId() << " ";

		Vector2D x0 = AddVectors2D(embedding.getCoordinate(t->getEdge(0)->getOpposite()),shift);
		file << ProperFloor(x0[0]) << " " << ProperFloor(x0[1]) << " ";
		x0 = AddVectors2D(x0,embedding.getForm(i,2));
		file << ProperFloor(x0[0]) << " " << ProperFloor(x0[1]) << " ";
		x0 = AddVectors2D(x0,embedding.getForm(i,0));
		file << ProperFloor(x0[0]) << " " << ProperFloor(x0[1]) << "\n";
	}

	for(int i=0;i<triangulation.NumberOfVertices();i++)
	{
		Vector2D v = embedding.getCoordinate(i);
		file << v[0] << " " << v[1] << "\n";
	}
	return 0;
}