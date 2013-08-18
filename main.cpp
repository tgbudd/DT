#include <iostream>

#include "triangulation.h"
#include "DualCohomologyBasis.h"
#include "ThetaModel.h"

int main()
{
	Triangulation triangulation;
	DualCohomologyBasis dualcohomologybasis(&triangulation);
	ThetaModel thetamodel(&triangulation,&dualcohomologybasis);
	
	triangulation.AddDecoration( &dualcohomologybasis );  
	triangulation.setDominantMatter( &thetamodel );

	triangulation.LoadRegularLattice(5,6);
	dualcohomologybasis.Initialize(5,6);
	thetamodel.Initialize();

	while(true)
	{
		triangulation.DoSweep(2);
		std::cout << "state: " << thetamodel.PrintState() << "\n";
	}

	return 0;
}