#include <iostream>

#include "triangulation.h"
#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "BitmapDrawer.h"
#include "ShortestLoop.h"

int main()
{
	Triangulation triangulation;
	CohomologyBasis cohomologybasis( &triangulation );
	triangulation.AddDecoration( &cohomologybasis );  

	triangulation.LoadRegularLattice(175,175);
	cohomologybasis.Initialize(175,175);

	for(int i=0;i<10;i++)
	{
		triangulation.DoSweep(1);
		cohomologybasis.Simplify(false);
		BOOST_ASSERT( cohomologybasis.CheckClosedness() );
		std::cout << i << "\n";
	}
	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );

	if( harmonicembedding.FindEmbedding() )
	{
		std::cout << "Embedding found\n";
	} else
	{
		std::cout << "Embedding not found\n";
	}

	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	FundamentalDomainDrawer funddrawer;

	BitmapDrawer bitmap(1000,1000,4);
	
	bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.5);
	bitmap.setPenWidth(8);
	bitmap.setPenColor(200,20,20);
	funddrawer.Draw(bitmap);
	bitmap.setPenWidth(5);
	bitmap.setPenColor(30,30,100);
	tridrawer.Draw(bitmap);


	/*ShortestLoop shortestloop(&triangulation,&cohomologybasis);
	ShortestLoopDrawer loopdrawer(&shortestloop,&harmonicembedding);
	shortestloop.FindGenerators();
	shortestloop.FindShortestLoop();
	bitmap.setPenWidth(7);
	bitmap.setPenColor(0,140,10);
	loopdrawer.DrawGenerators(bitmap);
	bitmap.setPenColor(10,10,110);
	loopdrawer.Draw(bitmap);*/

	bitmap.SaveImage("test.bmp");

	return 0;
}