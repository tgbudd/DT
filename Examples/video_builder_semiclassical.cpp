/* 
Source file used to make the video dt-semiclassical.avi on youtube 
(http://www.youtube.com/watch?v=TB9jpCRLTx0). The bitmaps were converted 
to an uncompressed avi using EasyBMPtoAVI with command:

EasyBMPtoAVI -start output/test-0000.bmp -end output/test-0699.bmp -framerate 10 -output out1.avi

This avi-file was then encoded in h264-codec using avidemux (including a vertical flip filter, because
for some reason they end up upside down.
*/

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

	for(int i=0;i<20;i++)
	{
		triangulation.DoSweep(10);
		cohomologybasis.Simplify(false);
	}

	int framespersecond = 10;
	int movesperframe = 10;

	for(int i=0;i<1400;i++)
	{
		if( i % 20 == 0 )
		{
			cohomologybasis.Simplify(true);
		}
		if( (i+1)%80 == 0 )
		{
			if( centralcharge > -0.1 )
				centralcharge = -1.0;
			else if( centralcharge > - 500.0 )
				centralcharge = (int)(2.0*centralcharge);
			laplaciandeterminant.setCentralCharge( centralcharge );

			if( centralcharge < -30.0 )
			{
				movesperframe--;
				if( movesperframe < 5 )
					movesperframe = 5;
			}
		}

		harmonicembedding.FindEmbedding();

		bitmap.Clear();

		bitmap.SetPeriodicDomain(harmonicembedding.CalculateModuli(),0.3);
		tridrawer.DrawShading(bitmap);
		bitmap.setPenWidth(4);
		bitmap.setPenColor(0,0,0);
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