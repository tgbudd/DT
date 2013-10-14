#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "Embedding.h"
#include "HarmonicEmbedding.h"
#include "CohomologyBasis.h"
#include "DualCohomologyBasis.h"
#include "SpanningTree.h"
#include "DualScalarField.h"
#include "utilities.h"
#include "Diffusion.h"
#include "CMinusTwoBuilder.h"
#include "BitmapDrawer.h"



class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, SpanningTree * spanningtree, DualScalarField * scalar );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void Save();
	void SetShade();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;
	SpanningTreeDrawer spanningtreedrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;
	SpanningTree * spanningtree_;
	DualScalarField * dualscalar_;

	std::string info_;
	std::string prefix_;
	int counter_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, SpanningTree * spanningtree, DualScalarField * scalar ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  spanningtree_(spanningtree),
	  dualscalar_(scalar),
	  bitmap_(1280,720,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  spanningtreedrawer_( triangulation, spanningtree, embedding),
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

		SetShade();
		tridrawer_.DrawShading(bitmap_);

		//bitmap_.setPenWidth(9);
		//bitmap_.setPenColor(220,120,120);
		//funddrawer_.Draw(bitmap_);
		bitmap_.setPenWidth(4);
		bitmap_.setPenColor(210,210,210);
		tridrawer_.Draw(bitmap_);
		bitmap_.setPenWidth(10);
		bitmap_.setPenColor(10,10,10);
		spanningtreedrawer_.Draw(bitmap_);
		std::ostringstream text;
		text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
		text << info_;
		textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

		Save();
	}
}

void Snapshot::SetShade()
{
	double average=0.0, sqaverage=0.0;
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		double field = dualscalar_->getField(i);
		average += field;
		sqaverage += field*field;
	}
	average /= triangulation_->NumberOfTriangles();
	sqaverage = std::sqrt( sqaverage/triangulation_->NumberOfTriangles() );
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		tridrawer_.SetTriangleShade(i,std::max(0.0,std::min(1.0,0.5+0.1*(dualscalar_->getField(i)-average)/sqaverage)));
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

	int w = param.Read<int>("width");
	int h = param.Read<int>("height");
	int sweepsize = param.Read<int>("sweep size");
	int therm = param.Read<int>("Thermalization sweeps");
	int snap = param.Read<int>("Snapshot sweeps");

	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.SetCustomSweepSize( sweepsize );
	triangulation.LoadRegularLattice(w,h);

	CohomologyBasis cohom( &triangulation );
	cohom.Initialize();
	cohom.SetStayInSameClass( true );
	//cohom.SetMakeUpToDateViaReinitialization(true);
	triangulation.AddDecoration( &cohom );

	SpanningTree spanningtree( &triangulation );
	spanningtree.InitializeWithLoopErasedRandomWalk();
	triangulation.setDominantMatter( &spanningtree );
	spanningtree.SetCustomSweepSize( sweepsize );

	DualScalarField dualscalar( &triangulation );
	dualscalar.Initialize();
	triangulation.AddMatter( &dualscalar );
	dualscalar.SetCustomSweepSize( 20*sweepsize );

	HarmonicEmbedding harm(  &triangulation, &cohom );
	harm.FindEmbedding();
	harm.setWorkWithHarmonicForms( true );
	triangulation.AddDecoration( &harm );

	Snapshot snapshot( &triangulation, &harm, &spanningtree, &dualscalar );
	std::ostringstream prefix;
	prefix << "./output/snapshot-";
	snapshot.SetPrefix( prefix.str() );
	std::ostringstream info;
	info << "Spanning tree & scalar";
	snapshot.SetInfo( info.str() );

	snapshot.Measure();
	
	Simulation simulation( &triangulation, therm, 10000000, output );
	simulation.AddObservable( &snapshot, snap );

	simulation.Run();

	return 0;
}