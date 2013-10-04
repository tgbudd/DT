#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "SpanningTree.h"
#include "DualScalarField.h"
#include "utilities.h"
#include "Diffusion.h"

#include "CMinusTwoBuilder.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int method = param.Read<int>(" --- Models ---\n 0: spanning tree\n 1: c=-2\n 2: pure gravity\n 3: spanning tree + 1 scalar\n 4: spanning tree + 2 scalars\n 5: spanning tree + 3 scalars\n 6: pure gravity + 1 scalar\n\nmodel");
	double masssquared = 0.0;
	if( method == 3 || method == 4 || method == 5 || method == 6)
	{
		masssquared = param.Read<double>("mass^2");
	}
	int n =	param.Read<int>("triangles");
	int thermalizationSweeps = 0;
	int measurementSweeps = 1;
	if( method != 1 )
	{
		thermalizationSweeps = param.Read<int>("thermalizationsweeps");
		measurementSweeps = param.Read<int>("measurementsweeps");
	}
	int secondsPerOutput = param.Read<int>("seconds per output");

	bool output = param.UserInput();

	Triangulation triangulation;
	CMinusTwoBuilder builder(&triangulation,0,n);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	DominantMatter * dom;
	if( method == 0 || method == 3 || method == 4 || method == 5 )
	{	
		SpanningTree * spanningtree = new SpanningTree( &triangulation );
		spanningtree->Initialize();
		dom = spanningtree;
		triangulation.setDominantMatter( dom );
	}
	if( method == 2 || method == 6 )
	{
		triangulation.clearDominantMatter();
	}

	DualScalarField * dualscalarfield;
	if( method == 3 || method == 4 || method == 5 || method == 6 )
	{
		dualscalarfield = new DualScalarField( &triangulation, masssquared );
		dualscalarfield->Initialize();
		triangulation.AddMatter( dualscalarfield );
	}
	DualScalarField * dualscalarfield2;
	if( method == 4 || method == 5 )
	{
		dualscalarfield2 = new DualScalarField( &triangulation, masssquared );
		dualscalarfield2->Initialize();
		triangulation.AddMatter( dualscalarfield2 );
	}
	DualScalarField * dualscalarfield3;
	if( method == 5 )
	{
		dualscalarfield3 = new DualScalarField( &triangulation, masssquared );
		dualscalarfield3->Initialize();
		triangulation.AddMatter( dualscalarfield3 );
	}


	Simulation simulation( &triangulation, thermalizationSweeps, secondsPerOutput, output );

	Diffusion diffusion( &triangulation );
	simulation.AddObservable( &diffusion, measurementSweeps );

	simulation.Run();

	return 0;
}