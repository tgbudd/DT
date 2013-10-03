#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "SpanningTree.h"
#include "DualScalarField.h"
#include "utilities.h"
#include "Diffusion.h"

#include "CMinusTwoBuilder.h"

#include "CohomologyBasis.h"
#include "HarmonicEmbedding.h"
#include "BitmapDrawer.h"


class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void Save();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;

	std::string info_;
	std::string prefix_;
	int counter_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  bitmap_(1280,720,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  counter_(0)
{
	textdrawer_.setColor(160,0,30);
}

void Snapshot::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		bitmap_.Clear();

		bitmap_.SetPeriodicDomain(embedding_->CalculateModuli(),0.24);
		bitmap_.setPenWidth(4);
		bitmap_.setPenColor(0,0,0);
		tridrawer_.Draw(bitmap_);
		std::ostringstream text;
		text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
		text << info_;
		textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

		Save();
	}
}

void Snapshot::Save()
{
	std::ostringstream os;
	os << prefix_ << std::setw( 4 ) << std::setfill( '0' ) << counter_++ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int method = param.Read<int>("spanning tree (0), c=-2 (1), c=0 (2), c=-1 (3)");
	int h =	param.Read<int>("height");
	int w = param.Read<int>("width");
	int thermalizationSweeps = param.Read<int>("thermalizationsweeps");
	int measurementSweeps = param.Read<int>("measurementsweeps");
	int secondsPerOutput = param.Read<int>("seconds per output");

	Triangulation triangulation;
	triangulation.SeedRandom( 123 );
	triangulation.LoadRegularLattice(w,h);

	DominantMatter * dom;
	if( method == 0 || method == 3 )
	{	
		SpanningTree * spanningtree = new SpanningTree( &triangulation );
		spanningtree->Initialize();
		dom = spanningtree;
	} else if( method == 1 )
	{
		dom = new CMinusTwoBuilder(&triangulation,0,w*h*2);
	}
	if( method != 2 )
	{
		triangulation.setDominantMatter( dom );
	}
	triangulation.DoSweep();

	DualScalarField * dualscalarfield;
	if( method == 3 )
	{
		dualscalarfield = new DualScalarField( &triangulation );
		dualscalarfield->Initialize();
		triangulation.AddMatter( dualscalarfield );
	}
	/*CohomologyBasis cohomologybasis( &triangulation );
	cohomologybasis.Initialize(w,h);
	triangulation.AddDecoration( &cohomologybasis );  

	//BOOST_ASSERT( cohomologybasis.CheckClosedness() );

	HarmonicEmbedding harmonicembedding( &triangulation, &cohomologybasis );
	
	// add harmonic embedding as decoration to keep track of the embedding after flip moves
	triangulation.AddDecoration( &harmonicembedding );

	DualScalarField dualscalarfield( &triangulation );
	dualscalarfield.Initialize();
	triangulation.AddMatter( &dualscalarfield );
	

	TextDrawer textdrawer("lucida",2);
	textdrawer.setColor(160,0,30);
	TriangulationDrawer tridrawer( &triangulation, &harmonicembedding );
	FundamentalDomainDrawer funddrawer;

	BitmapDrawer bitmap(1280,720,4);

	BOOST_ASSERT( harmonicembedding.FindEmbedding() );
	*/

	Simulation simulation( &triangulation, thermalizationSweeps, secondsPerOutput );

	Diffusion diffusion( &triangulation );
	simulation.AddObservable( &diffusion, measurementSweeps );

/*	Snapshot snapshot( &triangulation, &harmonicembedding );
	std::ostringstream prefix;
	prefix << "./output/snapshot-";
	snapshot.SetPrefix( prefix.str() );
	std::ostringstream info;
	info << "Spanning tree";
	snapshot.SetInfo( info.str() );

	snapshot.Measure();

	simulation.AddObservable( &snapshot, measurementSweeps );
*/

	simulation.Run();

	return 0;
}