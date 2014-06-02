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
#include "utilities.h"
#include "HyperbolicStructure.h"
//#include "LaplacianDeterminant.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	//double maxtheta = param.Read<double>("maxtheta (x PI)");
	double maxlength = param.Read<double>("maxlength");
//	double pow = param.Read<double>("cosine power");
//	double c = param.Read<double>("central charge");
	//int seed = param.Read<int>("seed");
	bool output = param.UserInput();


	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	
	//////
	//triangulation.SeedRandom(seed);
	//////
	
	DualCohomologyBasis dualcohomologybasis( &triangulation );
	dualcohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &dualcohomologybasis );

	ThetaModel thetamodel( &triangulation, &dualcohomologybasis, 18000 );
	triangulation.setDominantMatter( &thetamodel );
	thetamodel.Initialize();

/*	if( std::fabs(pow) > 0.001 )
	{
		thetamodel.setCosinePower(pow);
	}
	*/
/*	LaplacianDeterminant lapl(&triangulation, c );
	if( std::fabs(c) > 0.0001 )
	{
		triangulation.AddMatter(&lapl);
	}
*/

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation, dualcohomologybasis );
	cohom.SetMakeUpToDateViaReinitialization(true);

	CirclePattern circlepattern( &triangulation, &cohom, &thetamodel );

	//ModuliObservable moduli( &circlepattern );
	//simulation.AddObservable( &moduli, measurementSweeps );
	AngleHistogram angle( &circlepattern, &triangulation );
	simulation.AddObservable( &angle, measurementSweeps );
	//ThetaHistogram theta( &thetamodel );
	//simulation.AddObservable( &theta, measurementSweeps );
	
	HyperbolicStructure hyp(&triangulation,&thetamodel,&circlepattern,&dualcohomologybasis);
	//hyp.setMaxTheta(maxtheta*PI);
	hyp.setMaxLength(maxlength);
	simulation.AddObservable( &hyp, measurementSweeps );

	simulation.Run();

	return 0;
}