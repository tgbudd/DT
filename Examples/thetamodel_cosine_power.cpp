#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "DualCohomologyBasis.h"
#include "ThetaModel.h"
#include "CirclePattern.h"
#include "ModuliObservable.h"
#include "AngleHistogram.h"
#include "ThetaHistogram.h"

int main(int argc, char* argv[])
{
	bool output = true;
	int w,h;
	int thermalizationSweeps, MeasurementSweeps, OutputSweeps;
	double cosinepower;
	if( argc < 7 )
	{
		std::cout << "width = ";
		std::cin >> w;
		std::cout << "height = ";
		std::cin >> h;
		std::cout << "power = ";
		std::cin >> cosinepower;
		std::cout << "thermalization sweeps = ";
		std::cin >> thermalizationSweeps;
		std::cout << "measurement sweeps = ";
		std::cin >> MeasurementSweeps;
		std::cout << "output sweeps = ";
		std::cin >> OutputSweeps;
	}else
	{
		std::istringstream is(argv[1]);
		is >> w;
		std::istringstream is2(argv[2]);
		is2 >> h;
		std::istringstream is3(argv[3]);
		is3 >> cosinepower;
		std::istringstream is4(argv[4]);
		is3 >> thermalizationSweeps;			// STUPID BUG!!!!!!!!!!!
		std::istringstream is5(argv[5]);
		is4 >> MeasurementSweeps;
		std::istringstream is6(argv[6]);
		is5 >> OutputSweeps;
		output = false;
	}

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	
	DualCohomologyBasis dualcohomologybasis( &triangulation );
	dualcohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &dualcohomologybasis );

	ThetaModel thetamodel( &triangulation, &dualcohomologybasis, 18000 );
	triangulation.setDominantMatter( &thetamodel );
	thetamodel.Initialize();
	thetamodel.setCosinePower( cosinepower );

	Simulation simulation( &triangulation, thermalizationSweeps, OutputSweeps );
	std::ostringstream config;
	config << "cosinepower -> " << std::fixed << cosinepower;
	simulation.AddConfigurationInfo(config.str());

	CohomologyBasis cohom( &triangulation, dualcohomologybasis );
	cohom.SetMakeUpToDateVia( &dualcohomologybasis );

	CirclePattern circlepattern( &triangulation, &cohom, &thetamodel );

	ModuliObservable moduli( &circlepattern );
	simulation.AddObservable( &moduli, MeasurementSweeps );
	AngleHistogram angle( &circlepattern );
	simulation.AddObservable( &angle, MeasurementSweeps );
	ThetaHistogram thetahistogram( &thetamodel );
	simulation.AddObservable( &thetahistogram, MeasurementSweeps );

	simulation.Run();

	return 0;
}