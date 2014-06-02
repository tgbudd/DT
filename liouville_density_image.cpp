#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "CohomologyBasis.h"
#include "CMinusTwoBuilder.h"
#include "BitmapDrawer.h"
#include "TriangulationProperties.h" 
#include "ShortestLoop.h"

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, ShortestLoop * shortestloop = NULL, double radius = 10.0 );
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
	ShortestLoop * shortestloop_;

	Triangulation * triangulation_;
	Embedding * embedding_;

	std::string info_;
	std::string prefix_;
	int counter_;

	double radius_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, ShortestLoop * shortestloop, double radius ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  bitmap_(1000,600,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  counter_(0),
	  radius_(radius),
	  shortestloop_(shortestloop)
{
	textdrawer_.setColor(160,0,30);
}

void Snapshot::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		std::pair<double,double> moduli = embedding_->CalculateModuli();

		if( moduli.first*moduli.first + (moduli.second-1.0)*(moduli.second-1.0) > radius_*radius_ )
		{
			std::cout << "to far from \tau = i\n";
			return;
		}

		ColorScheme::Scheme scheme = ColorScheme::RAINBOW;

		bitmap_.Clear();

		bitmap_.SetPeriodicDomain(moduli,0.45,0.5,0.5,0.5,0.5);

		tridrawer_.DrawShading(bitmap_,scheme);
		bitmap_.setPenWidth(3);
		bitmap_.setPenColor(20,20,20);
		tridrawer_.Draw(bitmap_);

		if( shortestloop_ != NULL )
		{
			shortestloop_->FindGenerators();
			shortestloop_->FindShortestLoop();
			ShortestLoopDrawer loopdrawer(shortestloop_,embedding_);
			bitmap_.setPenColor(250,250,250);
			bitmap_.setPenWidth(18);
			loopdrawer.Draw(bitmap_);
		}

		FundamentalDomainDrawer fund;
		bitmap_.setPenWidth(14);
		bitmap_.setPenColor(200,200,200);
		fund.Draw(bitmap_);

		Save();
	}
}

void Snapshot::Save()
{
	std::ostringstream os;
	os << prefix_ << std::setw( 4 ) << std::setfill( '0' ) << counter_ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
	counter_++;
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n =	param.Read<int>("triangles");
	int seed = param.Read<int>("seed");
	bool drawloop = (param.Read<int>("draw loop (1=yes)")==1);
	double radius = param.Read<double>("tau radius");

	int thermalizationSweeps = 0;
	int	measurementSweeps = 1; 
	int secondsperoutput = 999999; 
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.SeedRandom(seed);
	CMinusTwoBuilder builder(&triangulation,1,n);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding embedding( &triangulation, &cohom );
	embedding.SetAccuracy( 1.0e-8 );

	ShortestLoop * loop = NULL;
	if( drawloop )
	{
		loop = new ShortestLoop(&triangulation,&cohom);
	}

	Snapshot snapshot( &triangulation, &embedding, loop, radius );
	std::ostringstream prefix;
	prefix << "D:/temp/output/snapshot-" << simulation.GetIdentifier() << "-";
	snapshot.SetPrefix( prefix.str() );
	std::ostringstream info;
	info << "c = -2";
	snapshot.SetInfo( info.str() );
	
	simulation.AddObservable( &snapshot, measurementSweeps );
	simulation.Run();

	return 0;
}
