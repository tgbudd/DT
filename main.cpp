#include <iostream>

#include "triangulation.h"
#include "PottsModel.h"

int main()
{
	Triangulation triangulation;
	PottsModel * potts = new PottsModel(&triangulation,2);
	potts->setInteraction(3.0);

	triangulation.AddMatter( potts );  // Add an Ising model

	triangulation.LoadRegularLattice(5,6);
	potts->Initialize();
	
	while(true)
	{
		triangulation.DoSweep(1000);
		std::cout << potts->PrintState() << "\n";
	}

	return 0;
}