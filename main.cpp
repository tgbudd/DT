#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "HarmonicDiffusion.h"
#include "CohomologyBasis.h"
#include "CMinusTwoBuilder.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n =	param.Read<int>("triangles");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,1,n);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding embedding( &triangulation, &cohom );

	HarmonicDiffusion diffusion( &triangulation, &embedding, &cohom );
	simulation.AddObservable( &diffusion, measurementSweeps );

	simulation.Run();
	return 0;
}