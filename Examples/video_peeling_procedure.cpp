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
#include "PeelingProcedure.h"

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, CohomologyBasis * cohomologybasis );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void DrawShading();
	void Save();
	void RunCommand();
	double Area(Triangle * triangle);

	DiscreteColorScheme scheme_;

	Triangle * SelectTriangle();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;
	PeelingProcedure peelingprocedure_;

	std::string info_;
	std::string prefix_;
	int counter_;
	int batch_counter_;
	std::string lastfile_;
	std::string firstfile_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, CohomologyBasis * cohomologybasis ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  bitmap_(1920,1080,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  counter_(0),
	  batch_counter_(0),
	  peelingprocedure_(triangulation,cohomologybasis),
	  scheme_(DiscreteColorScheme::Scheme::COLORS8)
{
	textdrawer_.setColor(160,0,30);
}

void Snapshot::DrawShading()
{
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		if( peelingprocedure_.in_mother_universe_[i] )
		{
			Edge * edge = triangulation_->getTriangle(i)->getEdge(0);
			int dist = std::min(peelingprocedure_.distance_[edge->getOpposite()->getId()],std::min(peelingprocedure_.distance_[edge->getNext()->getOpposite()->getId()],peelingprocedure_.distance_[edge->getPrevious()->getOpposite()->getId()]));
			tridrawer_.DrawTriangleShading( bitmap_, i, scheme_.getColor( dist ) );
		}
	}
}

void Snapshot::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		Triangle * triangle = SelectTriangle();
		std::pair<double,double> moduli = embedding_->CalculateModuli();
		scheme_ = DiscreteColorScheme(static_cast<DiscreteColorScheme::Scheme>(triangulation_->RandomInteger(0,5)));
		
		peelingprocedure_.InitializePeeling( triangle );
		Vector2D center = embedding_->getCoordinate(peelingprocedure_.start_vertex_);

		double totalarea = 0.0;
		double lastarea = -1.0;
		int steps = 0;
		do {
			if( totalarea - lastarea > 0.0015 )
			{
				bitmap_.Clear();

				bitmap_.SetPeriodicDomain(moduli,0.85,0.5,0.5,center[0],center[1]);

				DrawShading();
				bitmap_.setPenWidth(2);
				bitmap_.setPenColor(50,50,50);
				tridrawer_.Draw(bitmap_);
				std::ostringstream text;
				text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
				text << "Volume: " << peelingprocedure_.volume_within_frontier_ << "\n";
				text << "Frontier: " << peelingprocedure_.frontier_size_ << "\n";
				text << "Steps: " << steps;
				textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

				Save();
				lastarea = totalarea;
			}
			totalarea += Area(peelingprocedure_.last_triangle_);
			steps++;
		} while( !peelingprocedure_.DoPeelingStepOnTorus() );

		RunCommand();
		batch_counter_++;
		counter_=0;
	}
}
double Snapshot::Area(Triangle * triangle)
{
	Vector2D v0 = embedding_->getForm(triangle->getId(),0), v1 = NegateVector2D(embedding_->getForm(triangle->getId(),2));
	return std::fabs(0.5*(v0[0]*v1[1] - v0[1]*v1[0]));
}

Triangle * Snapshot::SelectTriangle() 
{
	Triangle * triangle;
	double minarea = 1.0;
	std::vector<double> areas;
	areas.reserve(triangulation_->NumberOfTriangles());
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Vector2D v0 = embedding_->getForm(i,0), v1 = NegateVector2D(embedding_->getForm(i,2));
		double area = std::fabs(0.5*(v0[0]*v1[1] - v0[1]*v1[0]));
		areas.push_back(area);
	}
	std::vector<double> areas2(areas);
	std::sort(areas2.begin(),areas2.end());
	double goal = areas2[static_cast<int>(0.6*areas2.size())];
	for(int i=0,endi=areas.size();i<endi;i++)
	{
		if( areas[i] == goal )
		{
			triangle = triangulation_->getTriangle(i);
		}
	}
	return triangle;
}

void Snapshot::Save()
{
	std::ostringstream os;
	os << prefix_ << batch_counter_ << "-" << std::setw( 4 ) << std::setfill( '0' ) << counter_ << ".bmp";
	lastfile_ = os.str();
	if( firstfile_.size() < 5 )
	{
		firstfile_ = lastfile_;
	}
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
	counter_++;
}

void Snapshot::RunCommand()
{
	std::ostringstream os;
	std::string prefix = lastfile_.substr(lastfile_.find_last_of("/\\")+1);
	prefix = prefix.substr(0,prefix.find_last_of("-"));
	os << "D:\\temp\\bmptoavi.bat " << firstfile_ << " " << lastfile_ << " " << prefix;
	std::cout << "run: " << os.str() << "\n";
	std::system(os.str().c_str());
	firstfile_.clear();
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int n =	param.Read<int>("triangles");
	int seed = param.Read<int>("seed");
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

	Snapshot snapshot( &triangulation, &embedding, &cohom );
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
