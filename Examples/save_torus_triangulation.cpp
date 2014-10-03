#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>

#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "ConnectivityRestrictor.h"

int ProperFloor(double x)
{
	return static_cast<int>(std::floor(x));
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n = param.Read<int>("N");
	int c = param.Read<int>("central charge (-2 or 0)");
	int sweeps=0;
	if( c==0 )
	{
		sweeps = param.Read<int>("sweeps");
	}
	std::string filename = param.Read<std::string>("file");

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,1,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	if( c==0 )
	{
		ConnectivityRestrictor conn(&triangulation,ConnectivityRestrictor::NO_DOUBLE_EDGES);
		triangulation.clearDominantMatter();
		triangulation.AddMatter(&conn);
		triangulation.DoSweep(sweeps);
	}

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);


	HarmonicEmbedding embedding( &triangulation, &cohom );
	embedding.SetAccuracy( 1.0e-8 );
	
	embedding.MakeUpToDate();

	std::pair<double,double> modulus = embedding.CalculateModuli();

	std::ofstream file(filename.c_str());

	file << "{\"numtriangles\" -> " << triangulation.NumberOfTriangles();
	file << ", \"modulus\" -> " << std::fixed << modulus.first << " + I*" << modulus.second;
	file << ", \"tri\" -> {";
	for(int i=0;i<triangulation.NumberOfTriangles();i++)
	{
		Triangle * t = triangulation.getTriangle(i);
		Vector2D center = embedding.getCoordinate(t->getEdge(0)->getOpposite());
		center = AddVectors2D(center,ScaleVector2D(embedding.getForm(i,2),1/3.0));
		center = AddVectors2D(center,ScaleVector2D(embedding.getForm(i,1),-1/3.0));

		Vector2D shift = {-std::floor(center[0]),-std::floor(center[1])};

		file << (i>0?",":"") << "{{" << t->getEdge(0)->getOpposite()->getId()+1 << "," << t->getEdge(1)->getOpposite()->getId()+1 << "," << t->getEdge(2)->getOpposite()->getId()+1 << "},";
		Vector2D x0 = AddVectors2D(embedding.getCoordinate(t->getEdge(0)->getOpposite()),shift);
		file << "{{" << x0[0] << "," << x0[1] << "},";
		x0 = AddVectors2D(x0,embedding.getForm(i,2));
		file << "{" << x0[0] << "," << x0[1] << "},";
		x0 = AddVectors2D(x0,embedding.getForm(i,0));
		file << "{" << x0[0] << "," << x0[1] << "}}}";
	}
	file << "}, \"coor\" -> {";
	for(int i=0;i<triangulation.NumberOfVertices();i++)
	{
		Vector2D v = embedding.getCoordinate(i);
		file << (i>0?",":"") << "{" << v[0] << "," << v[1] << "}";
	}
	file << "}}\n";
	return 0;
}