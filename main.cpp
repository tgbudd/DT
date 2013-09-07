#include <iostream>
#include <sstream>
#include <iomanip>

#include "triangulation.h"
#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "BitmapDrawer.h"
#include "ShortestLoop.h"
#include "ConnectivityRestrictor.h"


int main()
{

	Triangulation triangulation;
	CohomologyBasis cohomologybasis( &triangulation );
	triangulation.AddDecoration( &cohomologybasis );  
	ConnectivityRestrictor connect( &triangulation, &cohomologybasis, ConnectivityRestrictor::NO_CONTRACTIBLE_DOUBLE_EDGES );
	triangulation.AddMatter( &connect );

	triangulation.LoadRegularLattice(8,6);
	cohomologybasis.Initialize(30,25);

	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );

	
	// add harmonic embedding as decoration to keep track of the embedding after flip moves
	triangulation.AddDecoration( &harmonicembedding );

	ShortestLoop shortestloop( &triangulation, &cohomologybasis );

	TextDrawer textdrawer("lucida",2);
	textdrawer.setColor(160,0,30);
	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	FundamentalDomainDrawer funddrawer;

	BitmapDrawer bitmap(1280,720,4);

	harmonicembedding.FindEmbedding();

	int framespersecond = 10;
	int movesperframe = 1;

	for(int i=0;i<1000;i++)
	{
		if( i % 20 == 1 )
		{
			cohomologybasis.Simplify(true);
		}
		if( i % framespersecond == 0 )
		{
			movesperframe = std::min( 100, (int)(1.1*movesperframe + 1.0) );
		}

		harmonicembedding.FindEmbedding();

		bitmap.Clear();

		bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.55);
		bitmap.setPenWidth(30);
		bitmap.setPenColor(170,170,170);
		funddrawer.Draw(bitmap);
		bitmap.setPenWidth(8);
		bitmap.setPenColor(0,0,70);
		tridrawer.Draw(bitmap);

		shortestloop.FindGenerators();
		shortestloop.FindShortestLoop();
		std::ostringstream text;
		text << "Pure gravity\n";
		text << "Triangles: " << triangulation.NumberOfTriangles() << "\n";
		text << "Moves/s: " << movesperframe * framespersecond << "\n";
		text << "Shortest loop: " << shortestloop.getShortestLoop().size();
		textdrawer.DrawText(text.str(),bitmap,30,30,0.93);

		std::ostringstream os;
		os << "D:\\temp\\dt\\output\\test-" << std::setw( 4 ) << std::setfill( '0' ) << i << ".bmp";
		bitmap.SaveImage(os.str());
		std::cout << os.str() << " saved\n";

		int SuccesfulMoves = 0;
		while( SuccesfulMoves < movesperframe )
		{
			if( triangulation.TryFlipMove() )
				SuccesfulMoves++;
		}


	}

	return 0;
}