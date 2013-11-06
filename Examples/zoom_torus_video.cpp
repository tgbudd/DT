#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "HarmonicDiffusion.h"
#include "CohomologyBasis.h"
#include "CMinusTwoBuilder.h"
#include "BitmapDrawer.h"
#include "TriangulationProperties.h" 

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void SetShading(Vertex * vertex);
	void Save();
	void RunCommand();

	Vertex * SelectVertex();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;

	std::string info_;
	std::string prefix_;
	int counter_;
	int batch_counter_;
	std::string lastfile_;
	std::string firstfile_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  bitmap_(1920,1080,4), 
	  textdrawer_("lucida",2), 
	  tridrawer_( triangulation, embedding ),
	  counter_(0),
	  batch_counter_(0)
{
	textdrawer_.setColor(160,0,30);
}

void Snapshot::SetShading(Vertex * vertex)
{
	std::vector<int> distance;
	properties::VertexDistanceList(triangulation_,vertex,distance);
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Edge * edge = triangulation_->getTriangle(i)->getEdge(0);
		int dist = std::min(distance[edge->getOpposite()->getId()],std::min(distance[edge->getNext()->getOpposite()->getId()],distance[edge->getPrevious()->getOpposite()->getId()]));
		tridrawer_.SetTriangleColorIndex(i,dist);
	}
}

void Snapshot::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		Vertex * vertex = SelectVertex();
		Vector2D center = embedding_->getCoordinate(vertex);
		std::pair<double,double> moduli = embedding_->CalculateModuli();
		SetShading(vertex);
		DiscreteColorScheme::Scheme scheme = static_cast<DiscreteColorScheme::Scheme>(triangulation_->RandomInteger(0,5));
		for(int i=0;i<200;i++)
		{
			bitmap_.Clear();

			bitmap_.SetPeriodicDomain(moduli,0.1*std::exp(0.06*i),0.5,0.5,center[0],center[1]);

			tridrawer_.DrawShading(bitmap_,scheme);
			bitmap_.setPenWidth(3);
			bitmap_.setPenColor(0,0,0);
			tridrawer_.Draw(bitmap_);
			std::ostringstream text;
			text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
			text << info_;
			textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

			Save();
		}
		RunCommand();
		batch_counter_++;
	}
}

Vertex * Snapshot::SelectVertex()
{
	Vertex * vertex;
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
			vertex = triangulation_->getTriangle(i)->getEdge(triangulation_->RandomInteger(0,2))->getOpposite();
		}
	}
	return vertex;
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

	Snapshot snapshot( &triangulation, &embedding );
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
