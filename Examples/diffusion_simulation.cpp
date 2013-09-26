#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "Diffusion.h"

int main(int argc, char* argv[])
{
	bool output = true;
	int w,h;
	int thermalizationSweeps, MeasurementSweeps, SecondsPerOutput;
	if( argc < 6 )
	{
		std::cout << "width = ";
		std::cin >> w;
		std::cout << "height = ";
		std::cin >> h;
		std::cout << "thermalization sweeps = ";
		std::cin >> thermalizationSweeps;
		std::cout << "measurement sweeps = ";
		std::cin >> MeasurementSweeps;
		std::cout << "seconds per output = ";
		std::cin >> SecondsPerOutput;
	}else
	{
		std::istringstream is(argv[1]);
		is >> w;
		std::istringstream is2(argv[2]);
		is2 >> h;
		std::istringstream is3(argv[3]);
		is3 >> thermalizationSweeps;
		std::istringstream is4(argv[4]);
		is4 >> MeasurementSweeps;
		std::istringstream is5(argv[5]);
		is5 >> SecondsPerOutput;
		output = false;
	}

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);

	Diffusion diffusion( &triangulation );

	Simulation simulation( &triangulation, thermalizationSweeps, SecondsPerOutput );
	simulation.AddObservable( &diffusion, MeasurementSweeps );

	simulation.Run();

	return 0;
}