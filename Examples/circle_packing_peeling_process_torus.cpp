#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "DiskCirclePacking.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "PeelingProcedure.h"
#include "BitmapDrawer.h"

class PeelingMeasurement : public Observable {
public:
	PeelingMeasurement( Triangulation * triangulation, CohomologyBasis * cohomologybasis );
	void Measure();
	std::string OutputData() { return ""; }
	void SetPrefix(std::string prefix) { prefix_=prefix; }

private:
	Triangle * SelectTriangle();
	void SaveBitmap();

	Triangulation * triangulation_;
	DiskCirclePacking diskcirclepacking_;
	PeelingProcedure peelingprocedure_;
	CohomologyBasis * cohomologybasis_;

	BitmapDrawer bitmap_;
	int fullwidth_, fullheight_;

	std::string prefix_;
	int counter_;
};

PeelingMeasurement::PeelingMeasurement(Triangulation * triangulation, CohomologyBasis * cohomologybasis ) 
	: triangulation_(triangulation), 
	  diskcirclepacking_(triangulation),
	  peelingprocedure_(triangulation,cohomologybasis),
	  cohomologybasis_(cohomologybasis),
	  fullwidth_(1000),
	  fullheight_(1000),
	  bitmap_(500,500,2),
	  counter_(0)
{
}

void PeelingMeasurement::Measure()
{
	if( cohomologybasis_->IsUpToDate() || cohomologybasis_->MakeUpToDate() )
	{

		Triangle * triangle = SelectTriangle();
		
		peelingprocedure_.InitializePeeling( triangle );
		DiscreteColorScheme scheme(DiscreteColorScheme::COLORS8);
		int a=0;
		do {
			if( a%2 == 0 )
			{
				if( peelingprocedure_.volume_within_frontier_ > 50 && diskcirclepacking_.FindEmbedding(peelingprocedure_.frontier_,triangle->getEdge(0)) )
				{
					std::vector<std::pair<Vertex*,std::pair<Vector2D,double> > > circles;
					diskcirclepacking_.getCircles(circles);

					bitmap_.Clear();
					bitmap_.setPenColor(200,200,200);
					bitmap_.Disk(std::pair<int,int>(fullwidth_/2,fullheight_/2),static_cast<int>(0.48*fullheight_));
					bitmap_.setPenColor(20,20,200);
					for(int i=0,endi=circles.size();i<endi;i++)
					{
						boost::array<unsigned char,3> c = scheme.getColor(peelingprocedure_.distance_[circles[i].first->getId()]);
						bitmap_.setPenColor(c[0],c[1],c[2]);				
						std::pair<int,int> x;
						x.first = static_cast<int>(fullwidth_/2 + fullheight_ * 0.48 * circles[i].second.first[0] );
						x.second = static_cast<int>(fullheight_/2 + fullheight_ * 0.48 * circles[i].second.first[1] );
						bitmap_.Disk(x, static_cast<int>(fullheight_ * 0.48 * circles[i].second.second) );
					}
					SaveBitmap();
					int hy=0;
				}
			}
			a++;
		} while( !peelingprocedure_.DoPeelingStepOnTorus() );
	}
}

Triangle * PeelingMeasurement::SelectTriangle() 
{
	return triangulation_->getRandomTriangle();
}

void PeelingMeasurement::SaveBitmap()
{
	std::ostringstream os;
	os << prefix_ << "-" << std::setw( 4 ) << std::setfill( '0' ) << counter_ << ".bmp";
	bitmap_.SaveImage(os.str());
	std::cout << os.str() << " saved\n";
	counter_++;
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	//int n = param.Read<int>("N = ");
	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	bool output = param.UserInput();

	Triangulation triangulation;
	//triangulation.LoadSphericalBySubdivision(n);
	triangulation.LoadRegularLattice(w,h);

	Simulation simulation( &triangulation, thermalizationSweeps, 10000, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	ConnectivityRestrictor connect( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	triangulation.AddMatter( &connect );

	PeelingProcedure peelingprocedure( &triangulation, & cohom );
	PeelingMeasurement peelingmeasurement( &triangulation, & cohom );

	std::ostringstream prefix;
	prefix << "D:/temp/output/circ-" << simulation.GetIdentifier() << "-";
	peelingmeasurement.SetPrefix( prefix.str() );
	
	simulation.AddObservable( &peelingmeasurement, measurementSweeps );

	simulation.Run();
	return 0;
}