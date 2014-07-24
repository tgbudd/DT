#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "MetricGraphObservable.h"
#include "CMinusTwoBuilder.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;
	int matter = param.Read<int>("matter (0=none,1=minus2)");
	int thermalizationSweeps=0;
	int measurementSweeps=1;
	int n = param.Read<int>("triangles");
	if( matter == 0 )
	{
		thermalizationSweeps = param.Read<int>("thermalization sweeps");
		measurementSweeps = param.Read<int>("measurement sweeps");
	}
	int samples = param.Read<int>("samples");
	double maxlength = param.Read<int>("max length");
	int lengthbins = param.Read<int>("length bins");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter( &builder );
	triangulation.DoSweep();
	if( matter == 0 )
	{
		triangulation.clearDominantMatter();
	}

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	MetricGraphObservable metric( &triangulation, maxlength, lengthbins, samples );
	
	simulation.AddObservable( &metric, measurementSweeps );
	simulation.SetDirectory("./output/");
	simulation.Run();
	return 0;
}