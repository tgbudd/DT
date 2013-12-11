#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CirclePacking.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "ConformalDistribution.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	ConnectivityRestrictor connect( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	triangulation.AddMatter( &connect );

	CirclePacking circlepacking( &triangulation, &cohom );

	ConformalDistribution conf( &circlepacking );

	simulation.AddObservable( &conf, measurementSweeps );

	simulation.Run();
	return 0;
}