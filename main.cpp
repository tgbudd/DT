#include <iostream>
#include <sstream>
#include <iomanip>

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

	triangulation.LoadRegularLattice(30,40);
	cohomologybasis.Initialize(30,40);

	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );

	
	// add harmonic embedding as decoration to keep track of the embedding after flip moves
	triangulation.AddDecoration( &harmonicembedding );

	/*HotEdgesDrawer hotedgesdrawer( &harmonicembedding );
	triangulation.AddDecoration( &hotedgesdrawer );
	hotedgesdrawer.setNumberOfEdges(120);
	hotedgesdrawer.setColor(0.9,0.6,0.6);
	*/

	TextDrawer textdrawer("lucida");
	textdrawer.setColor(200,0,0);
	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	FundamentalDomainDrawer funddrawer;

	BitmapDrawer bitmap(1280,720,2);

	harmonicembedding.FindEmbedding();

	for(int i=0;i<14;i++)
	{
		//harmonicembedding.setMaxIterations(20);
		if( i % 20 == 1 )
		{
			cohomologybasis.Simplify(true);
		}
		harmonicembedding.FindEmbedding();

		bitmap.Clear();

		bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.55);
		bitmap.setPenWidth(25);
		bitmap.setPenColor(170,170,170);
		funddrawer.Draw(bitmap);
		//hotedgesdrawer.Draw(bitmap);
		bitmap.setPenWidth(5);
		bitmap.setPenColor(0,0,70);
		tridrawer.Draw(bitmap);

		textdrawer.DrawText("Pure gravity\n2400 triangles\n2000 moves/s",bitmap,20,20,0.9);

		std::ostringstream os;
		os << "D:\\temp\\dt\\output\\test-" << std::setw( 4 ) << std::setfill( '0' ) << i << ".bmp";
		bitmap.SaveImage(os.str());
		std::cout << os.str() << " saved\n";

		int SuccesfulMoves = 0;
		while( SuccesfulMoves < 100 )
		{
			if( triangulation.TryFlipMove() )
				SuccesfulMoves++;
		}

	}


	



	/*ShortestLoop shortestloop(&triangulation,&cohomologybasis);
	ShortestLoopDrawer loopdrawer(&shortestloop,&harmonicembedding);
	shortestloop.FindGenerators();
	shortestloop.FindShortestLoop();
	bitmap.setPenWidth(7);
	bitmap.setPenColor(0,140,10);
	loopdrawer.DrawGenerators(bitmap);
	bitmap.setPenColor(10,10,110);
	loopdrawer.Draw(bitmap);*/


	return 0;
}