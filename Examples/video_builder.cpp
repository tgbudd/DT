/* 
Source file used to make the video dt-1500t-acc.avi on youtube 
(http://www.youtube.com/watch?v=c3NdgSIe030). The bitmaps were converted 
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


int main()
{

	Triangulation triangulation;
	CohomologyBasis cohomologybasis( &triangulation );
	triangulation.AddDecoration( &cohomologybasis );  

	triangulation.LoadRegularLattice(30,25);
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