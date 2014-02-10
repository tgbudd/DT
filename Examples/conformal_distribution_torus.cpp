#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "ConformalDistribution.h"
#include "CMinusTwoBuilder.h"

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
	int method = param.Read<int>("method (0=SQUARE_ROOT_OF_AREA, 1=INSCRIBED_CIRCLE)");
	int ballsize = param.Read<int>("max ball size");
	double factor = param.Read<double>("ball size factor");
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
	ConformalDistribution conf( &harm, &triangulation );
	if( method == 0 )
	{
		harm.SetRadiusDefintion(HarmonicEmbedding::SQUARE_ROOT_OF_AREA);
		conf.setMethodString("SQUARE_ROOT_OF_AREA");
	}else
	{
		harm.SetRadiusDefintion(HarmonicEmbedding::INSCRIBED_CIRCLE);
		conf.setMethodString("INSCRIBED_CIRCLE");
	}
	//conf.setSaveEmbedding(true,"points.txt");
	conf.setMeasureEuclideanBallSize(true);
	conf.setRecordBallSizesExponentially(ballsize,factor);

	simulation.AddObservable( &conf, measurementSweeps );

	simulation.SetDirectory("./output/");
	//simulation.SetDirectory("C:\\Users\\Timothy\\Documents\\research\\projects\\fractaldimensions\\data\\conformaldist\\local\\");

	simulation.Run();
	return 0;
}