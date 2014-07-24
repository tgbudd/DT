#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "CMinusTwoBuilder.h"
#include "BallSizeDistribution.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;
	int matter = param.Read<int>("matter (0=none,1=minus2)");

	int thermalizationSweeps=0;
	int measurementSweeps=1;
	int n=2;
	if( matter == 0 )
	{
		int w =	param.Read<int>("width");
		int h = param.Read<int>("height");
		thermalizationSweeps = param.Read<int>("thermalization sweeps");
		measurementSweeps = param.Read<int>("measurement sweeps");
		triangulation.LoadRegularLattice(w,h);
	} else
	{
		n =	param.Read<int>("triangles");
	}
	int samples = param.Read<int>("samples");
	double factor = param.Read<double>("volume fraction factor (<1.0)");
	int fractions = param.Read<int>("max power of volume fraction");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	CMinusTwoBuilder builder(&triangulation,1,n);
	builder.setRemoveBabyUniverses(true);
	if( matter != 0 )
	{
		triangulation.setDominantMatter( &builder );
		triangulation.DoSweep();
	}
	ConnectivityRestrictor conn( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	if( matter == 0 )
	{
		triangulation.AddMatter( &conn );
	}

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding harm( &triangulation, &cohom );
	harm.SetAccuracy(1.0e-9);
	harm.setMaxIterations(20000);
	BallSizeDistribution ballsize( &harm, &triangulation, 18.0, 180, samples );

	ballsize.AddReferenceRadius(0.5);
	ballsize.AddReferenceRadius(0.25);
	ballsize.AddReferenceRadius(0.1);
	ballsize.AddReferenceRadius(0.05);
	double fraction = 1.0;
	for(int i=0;i<fractions;++i)
	{
		fraction *= factor; 
		ballsize.AddVolumeFraction(fraction);
	}
	ballsize.InitializeHistograms();
	
	simulation.AddObservable( &ballsize, measurementSweeps );
	simulation.SetDirectory("./output/");
	simulation.Run();
	return 0;
}