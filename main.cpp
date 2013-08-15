#include <iostream>

#include "triangulation.h"
//#include "PottsModel.h"
#include "CohomologyBasis.h"
#include "ThetaModel.h"

int main()
{
	Triangulation triangulation;
	CohomologyBasis cohomologybasis(&triangulation);
	ThetaModel thetamodel(&triangulation,&cohomologybasis);
	
	//PottsModel * potts = new PottsModel(&triangulation,2);
	//potts->setInteraction(3.0);

	triangulation.AddDecoration( &cohomologybasis );  
	triangulation.AddMatter( &thetamodel );

	triangulation.LoadRegularLattice(5,6);

	cohomologybasis.Initialize(5,6);
	thetamodel.Initialize();

	while(true)
	{
		triangulation.DoSweep(10);
		std::cout << thetamodel.PrintState() << "\n";
	}

	return 0;
}