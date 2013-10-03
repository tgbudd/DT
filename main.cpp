#include <iostream>
#include <sstream>

#include "Simulation.h"
#include "Triangulation.h"
#include "CMinusTwoBuilder.h"
#include "Diffusion.h"
#include "DualScalarField.h"

class ParameterStream {
public:
	ParameterStream(int argc, char** argv) : argc_(argc), argv_(argv), current_(1) {}
	template<class T> T Read( std::string name )
	{
		T t;
		if( current_ < argc_ )
		{
			std::istringstream is(argv_[current_]);
			is >> t;
			std::cout << name << " = " << t << "\n";
		} else
		{
			std::cout << name << " = ";
			std::cin >> t;
		}
		current_++;
		return t;
	}
private: 
	int current_;
	int argc_;
	char** argv_;
};

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n =	param.Read<int>("triangles");
	double mass = param.Read<double>("mass");
	int thermalizationSweeps = param.Read<int>("thermalizationsweeps");
	int measurementSweeps = param.Read<int>("measurementsweeps");
	int SecondsPerOutput = param.Read<int>("seconds per output");

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();
	triangulation.clearDominantMatter();

	DualScalarField dualscalarfield( &triangulation );
	dualscalarfield.Initialize();
	triangulation.AddMatter( &dualscalarfield );

	Diffusion diffusion( &triangulation );

	Simulation simulation( &triangulation, thermalizationSweeps, SecondsPerOutput );
	simulation.AddObservable( &diffusion, measurementSweeps );
	simulation.Run();

	return 0;
}