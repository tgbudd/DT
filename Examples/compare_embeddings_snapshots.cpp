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
#include "CirclePacking.h"

class Snapshot : public Observable {
public:
	Snapshot( Triangulation * triangulation, Embedding * embedding, CirclePacking * circlepacking );
	void Measure();
	std::string OutputData() { return ""; }
	void SetInfo(std::string info) { info_=info; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }
private:
	void SetShading(Vertex * vertex);
	void Save(bool isCirclePacking);
	void RunCommand();

	Vertex * SelectVertex();

	BitmapDrawer bitmap_;
	FundamentalDomainDrawer funddrawer_;
	TextDrawer textdrawer_;
	TriangulationDrawer tridrawer_;

	Triangulation * triangulation_;
	Embedding * embedding_;
	CirclePacking * circlepacking_;

	std::vector<int> distance_;

	std::string info_;
	std::string prefix_;
	int counter_;
	int batch_counter_;
	std::string lastfile_;
	std::string firstfile_;
};

Snapshot::Snapshot(Triangulation * triangulation, Embedding * embedding, CirclePacking * circlepacking ) 
	: triangulation_(triangulation), 
	  embedding_(embedding), 
	  circlepacking_(circlepacking),
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
	properties::VertexDistanceList(triangulation_,vertex,distance_);
	for(int i=0;i<triangulation_->NumberOfTriangles();i++)
	{
		Edge * edge = triangulation_->getTriangle(i)->getEdge(0);
		int dist = std::min(distance_[edge->getOpposite()->getId()],std::min(distance_[edge->getNext()->getOpposite()->getId()],distance_[edge->getPrevious()->getOpposite()->getId()]));
		tridrawer_.SetTriangleColorIndex(i,dist);
	}
}

void Snapshot::Measure()
{
	if( (embedding_->IsUpToDate() || embedding_->MakeUpToDate()) &&
		(circlepacking_->IsUpToDate() || circlepacking_->MakeUpToDate() ) )
	{
		Vertex * vertex = SelectVertex();
		Vector2D center = embedding_->getCoordinate(vertex);
		std::pair<double,double> moduli = embedding_->CalculateModuli();
		std::vector<double> radii;
		circlepacking_->GetRadii(radii);

		SetShading(vertex);
		//DiscreteColorScheme::Scheme scheme = static_cast<DiscreteColorScheme::Scheme>(triangulation_->RandomInteger(0,5));
		DiscreteColorScheme::Scheme scheme = DiscreteColorScheme::COLORS2;
		DiscreteColorScheme circlescheme(scheme);
		//for(int i=0;i<200;i++)
		for(int i=10;i<30;i+=5)
		{
			for(int picnum=0;picnum<2;++picnum)
			{
				bool circlepack = picnum==1;
				bitmap_.Clear();

				bitmap_.SetPeriodicDomain(moduli,0.1*std::exp(0.06*i),0.5,0.5,center[0],center[1]);

				if( !circlepack )
				{
					tridrawer_.DrawShading(bitmap_,scheme);
					bitmap_.setPenWidth(1);
					//bitmap_.setPenWidth(3);
					bitmap_.setPenColor(0,0,0);
					tridrawer_.Draw(bitmap_);
				} else
				{

					bitmap_.Clear(0);
					for(int j=0,endj=radii.size();j<endj;j++)
					{
						boost::array<unsigned char,3> c = circlescheme.getColor(distance_[j]);
						bitmap_.setPenColor(c[0],c[1],c[2]);
						Edge * edge = triangulation_->getVertex(j)->getParent()->getPrevious();
						double fraction = radii[j]/(radii[j]+radii[edge->getEnd()->getId()]);
						Vector2D center = circlepacking_->getCoordinate(j);
						Vector2D through = AddScaledVectors2D(1.0,center,fraction,circlepacking_->getForm(edge));
						bitmap_.domainDisk(center,through);
					}
				}
				std::ostringstream text;
				text << "Triangles: " << triangulation_->NumberOfTriangles() << "\n";
				text << info_;
				//textdrawer_.DrawText(text.str(),bitmap_,30,30,0.93);

				Save(circlepack);
			}
		}
		//RunCommand();
		batch_counter_++;
	}
}

Vertex * Snapshot::SelectVertex()
{
	Vertex * vertex;
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

void Snapshot::Save(bool isCirclePacking)
{
	std::ostringstream os;
	os << prefix_ << batch_counter_ << (isCirclePacking? "-c-":"-h-") << std::setw( 4 ) << std::setfill( '0' ) << counter_/2 << ".bmp";
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
	std::string path = param.Read<std::string>("path+prefix (\"0\"=\"/tmp/snapshots/snapshot-\")");
	if( path == "0" )
		path = "/tmp/snapshots/snapshot-";
	int thermalizationSweeps = 0;
	int	measurementSweeps = 1; 
	int secondsperoutput = 999999; 
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.SeedRandom(seed);
	CMinusTwoBuilder builder(&triangulation,1,n);
	builder.setRemoveBabyUniverses(true);
	triangulation.setDominantMatter(&builder);
	triangulation.DoSweep();

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding embedding( &triangulation, &cohom );
	embedding.SetAccuracy( 1.0e-8 );

	CirclePacking circle( &triangulation, &cohom );
	circle.SetAccuracy(1.0e-8);
	circle.setMaxIterations(50000);

	Snapshot snapshot( &triangulation, &embedding, &circle );
	std::ostringstream prefix;
	prefix << path << simulation.GetIdentifier() << "-";
	snapshot.SetPrefix( prefix.str() );
	std::ostringstream info;
	info << "c = -2";
	snapshot.SetInfo( info.str() );
	
	simulation.AddObservable( &snapshot, measurementSweeps );
	simulation.Run();

	return 0;
}
