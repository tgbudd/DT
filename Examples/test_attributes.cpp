#include <iostream>
#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CMinusTwoBuilder.h"
#include "SimplexAttribute.h"


struct coor {
	int x,y;
};

int main(int argc, char* argv[])
{
	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,1,500);
	triangulation.setDominantMatter( &builder );
	triangulation.DoSweep();

	VertexAttribute<double> radius(&triangulation,1.0);
	const Vertex * v = triangulation.getVertex(12);
	radius[v] += 3.2;
	std::cout << radius[v] << "\n";

	TriangleAttribute<bool> visited(&triangulation,false);
	visited[triangulation.getTriangle(10)] = true;


	coor zerocoor = {0,0};
	EdgeAttribute<coor> coors(&triangulation,zerocoor);
	coors[triangulation.getRandomEdge()].x = 45;
	
	return 0;
}