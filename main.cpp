#include <iostream>
#include <sstream>
#include <iomanip>

#include "triangulation.h"
#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "BitmapDrawer.h"
#include "ShortestLoop.h"
#include "ConnectivityRestrictor.h"
#include "LaplacianDeterminant.h"

int main()
{

	Triangulation triangulation;
	triangulation.LoadRegularLattice(10,9);
	
	CohomologyBasis cohomologybasis( &triangulation );
	cohomologybasis.Initialize(10,9);
	triangulation.AddDecoration( &cohomologybasis );  

	ConnectivityRestrictor connect( &triangulation, &cohomologybasis, ConnectivityRestrictor::NO_CONTRACTIBLE_DOUBLE_EDGES );
	triangulation.AddMatter( &connect );

	double centralcharge = 0.0;
	LaplacianDeterminant laplaciandeterminant( &triangulation, centralcharge );
	triangulation.AddMatter( &laplaciandeterminant );


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
	int movesperframe = 20;

	for(int i=0;i<1000;i++)
	{
		if( i % 20 == 1 )
		{
			cohomologybasis.Simplify(true);
		}
		if( i % framespersecond == 0 )
		{
			//movesperframe = std::min( 100, (int)(1.1*movesperframe + 1.0) );
		}
		if( i == 40 )
		{
			centralcharge = -20.0;
			laplaciandeterminant.setCentralCharge( centralcharge );
		}
		if( i == 100 )
		{
			centralcharge = -40.0;
			laplaciandeterminant.setCentralCharge( centralcharge );
		}
		if( i == 160 )
		{
			centralcharge = -80.0;
			laplaciandeterminant.setCentralCharge( centralcharge );
		}

		harmonicembedding.FindEmbedding();

		bitmap.Clear();

		bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.4);
		bitmap.setPenWidth(30);
		bitmap.setPenColor(170,170,170);
		//funddrawer.Draw(bitmap);
		bitmap.setPenWidth(8);
		bitmap.setPenColor(0,0,70);
		tridrawer.Draw(bitmap);

		shortestloop.FindGenerators();
		shortestloop.FindShortestLoop();
		std::ostringstream text;
		text << "Central charge: " << centralcharge << "\n";
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