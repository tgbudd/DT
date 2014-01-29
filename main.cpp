#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CirclePacking.h"
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
	CirclePacking circlepacking_;
	PeelingProcedure peelingprocedure_;
	CohomologyBasis * cohomologybasis_;

	BitmapDrawer bitmap_;
	int fullwidth_, fullheight_;

	std::string prefix_;
	int counter_;
};

PeelingMeasurement::PeelingMeasurement(Triangulation * triangulation, CohomologyBasis * cohomologybasis ) 
	: triangulation_(triangulation), 
	  circlepacking_(triangulation,cohomologybasis),
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

		do {
			if( peelingprocedure_.volume_within_frontier_ > 20 && circlepacking_.FindDiskEmbedding(peelingprocedure_.frontier_,triangle->getEdge(0)) )
			{
				std::vector<std::pair<Vector2D,double> > circles;
				circlepacking_.getCircles(circles);

				bitmap_.Clear();
				bitmap_.setPenColor(20,20,200);
				for(int i=0,endi=circles.size();i<endi;i++)
				{
					std::pair<int,int> x;
					x.first = static_cast<int>(fullwidth_/2 + fullheight_ * ( 0.5 + 0.48 * circles[i].first[0] ));
					x.second = static_cast<int>(fullheight_/2 + fullheight_ * ( 0.5 + 0.48 * circles[i].first[1] ));
					bitmap_.Disk(x, static_cast<int>(fullheight_ * 0.48 * circles[i].second) );
				}
				SaveBitmap();
			}
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

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	ConnectivityRestrictor connect( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	triangulation.AddMatter( &connect );

	PeelingProcedure peelingprocedure( &triangulation, &cohom );
	PeelingMeasurement peelingmeasurement( &triangulation, & cohom );

	std::ostringstream prefix;
	prefix << "D:/temp/output/circ-" << simulation.GetIdentifier() << "-";
	peelingmeasurement.SetPrefix( prefix.str() );
	
	simulation.AddObservable( &peelingmeasurement, measurementSweeps );

	simulation.Run();
	return 0;
}