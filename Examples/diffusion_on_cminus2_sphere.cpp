#include <iostream>
#include <sstream>

#include "Simulation.h"
#include "Triangulation.h"
#include "CMinusTwoBuilder.h"
#include "Diffusion.h"

int main(int argc, char* argv[])
{
	bool output = true;
	int n;
	int SecondsPerOutput;
	if( argc < 3 )
	{
		std::cout << "triangles = ";
		std::cin >> n;
		std::cout << "seconds per output = ";
		std::cin >> SecondsPerOutput;
	}else
	{
		std::istringstream is1(argv[1]);
		is1 >> n;
		std::istringstream is2(argv[2]);
		is2 >> SecondsPerOutput;
		output = false;
	}

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Diffusion diffusion( &triangulation );

	Simulation simulation( &triangulation, 0, SecondsPerOutput );
	simulation.AddObservable( &diffusion, 1 );
	simulation.Run();

	return 0;
}