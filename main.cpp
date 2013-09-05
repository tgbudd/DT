#include <iostream>

#include "triangulation.h"
#include "CohomologyBasis.h"
#include "DualCohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "BitmapDrawer.h"
#include "ShortestLoop.h"

int main()
{
	/*Triangulation triangulation;
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
	}*/

	Triangulation triangulation;
	CohomologyBasis cohomologybasis( &triangulation );
	triangulation.AddDecoration( &cohomologybasis );  

	triangulation.LoadRegularLattice(7,10);
	cohomologybasis.Initialize(7,10);

	triangulation.DoSweep(20);

	BOOST_ASSERT(cohomologybasis.CheckClosedness());
	DualCohomologyBasis dualcohom(cohomologybasis);
	BOOST_ASSERT(dualcohom.CheckClosedness());
	CohomologyBasis cohomologybasis2( &triangulation, dualcohom );
	BOOST_ASSERT(cohomologybasis2.CheckClosedness());

	ShortestLoop shortestloop( &triangulation, &cohomologybasis);

	shortestloop.FindGenerators();
	std::vector<std::list<Edge*> > generators = shortestloop.getGenerators();
	std::vector<IntForm2D> integrals = shortestloop.getGeneratorIntegrals();

	std::vector<IntForm2D> integrals2,integrals3;
	for(int i=0;i<2;i++)
	{
		integrals2.push_back(cohomologybasis.Integrate(generators[i]));
		integrals3.push_back(cohomologybasis2.Integrate(generators[i]));
	}

	BOOST_ASSERT(integrals2[0][0] == integrals[0][0]);
	BOOST_ASSERT(integrals2[0][1] == integrals[0][1]);
	BOOST_ASSERT(integrals2[1][0] == integrals[1][0]);
	BOOST_ASSERT(integrals2[1][1] == integrals[1][1]);
	BOOST_ASSERT(integrals3[0][0] == integrals[0][0]);
	BOOST_ASSERT(integrals3[0][1] == integrals[0][1]);
	BOOST_ASSERT(integrals3[1][0] == integrals[1][0]);
	BOOST_ASSERT(integrals3[1][1] == integrals[1][1]);

	cohomologybasis2.Simplify();
	std::vector<IntForm2D> integrals4;
	for(int i=0;i<2;i++)
	{
		integrals4.push_back(cohomologybasis2.Integrate(generators[i]));
	}
	int det = integrals[0][0] * integrals[1][1] - integrals[1][0] * integrals[0][1];
	BOOST_ASSERT(integrals4[0][0] == integrals[0][0]);
	BOOST_ASSERT(integrals4[0][1] == integrals[0][1]);
	BOOST_ASSERT(integrals4[1][0] == integrals[1][0]);
	BOOST_ASSERT(integrals4[1][1] == integrals[1][1]);
	BOOST_ASSERT(cohomologybasis2.CheckClosedness());

	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );

	if( harmonicembedding.FindEmbedding() )
	{
		std::cout << "Embedding found\n";
	} else
	{
		std::cout << "Embedding not found\n";
	}

	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	BitmapDrawer bitmap(1024,768,4);
	bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.4);
	bitmap.setPenWidth(5);
	bitmap.setPenColor(180,180,180);
	tridrawer.Draw(bitmap);

	ShortestLoopDrawer loopdrawer(&shortestloop,&harmonicembedding);
	shortestloop.FindGenerators();
	shortestloop.FindShortestLoop();
	bitmap.setPenWidth(7);
	bitmap.setPenColor(0,140,10);
	loopdrawer.DrawGenerators(bitmap);
	bitmap.setPenColor(10,10,110);
	loopdrawer.Draw(bitmap);

	bitmap.SaveImage("test.bmp");

	return 0;
}