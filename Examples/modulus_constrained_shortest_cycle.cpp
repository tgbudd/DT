#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "HarmonicEmbedding.h"
#include "CohomologyBasis.h"
#include "CMinusTwoBuilder.h"
#include "ShortestLoop.h"

class ConstrainedLoopLength : public Observable 
{
public:
	ConstrainedLoopLength(Embedding * embedding, ShortestLoop * shortestloop, std::pair<double,double> modulus, double maxradius, int bins );
	void Measure();
	std::string OutputData() const;
private:
	Embedding * embedding_;
	ShortestLoop * shortestloop_;
	std::pair<double,double> modulus_;
	double maxradius_;

	std::vector<std::vector<int> > histogram_;
	int attempts_;
	int bins_;
};

ConstrainedLoopLength::ConstrainedLoopLength(Embedding * embedding, ShortestLoop * shortestloop, std::pair<double,double> modulus, double maxradius, int bins)
	: embedding_(embedding),
	 shortestloop_(shortestloop),
	 modulus_(modulus),
	 maxradius_(maxradius),
	 bins_(bins),
	 attempts_(0)
{
	histogram_.resize(bins);
}

void ConstrainedLoopLength::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		std::pair<double,double> tau = embedding_->CalculateModuli();
		double r = std::sqrt(  (tau.first - modulus_.first)*(tau.first - modulus_.first) + (tau.second - modulus_.second)*(tau.second - modulus_.second) );
		if( r < maxradius_ )
		{
			shortestloop_->FindGenerators();
			shortestloop_->FindShortestLoop();
			int length = shortestloop_->getShortestLoop().size();
			BOOST_ASSERT( length >= 1 );

			for(int i=static_cast<int>(bins_*r/maxradius_);i<bins_;i++)
			{
				if( length > static_cast<int>(histogram_[i].size()) )
				{
					histogram_[i].resize(length,0);
				}
				histogram_[i][length-1]++;
			}
		}
		attempts_++;
	}
}

std::string ConstrainedLoopLength::OutputData() const
{
	std::ostringstream stream;
	stream << "constrainedlooplength -> {tau -> " << modulus_.first << " + " << modulus_.second << "*I";
	stream << ", maxradius -> " << maxradius_ << ", bins -> " << bins_; 
	stream << ", measurementattempts -> " << attempts_ << ", histogram -> ";
	PrintToStream2D(stream,histogram_.begin(),histogram_.end());
	stream << "}";
	return stream.str();
}

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);
	
	Triangulation triangulation;
	int matter = param.Read<int>("matter (0=none,1=minus2)");

	int thermalizationSweeps=0;
	int measurementSweeps=1;
	int n=2;
	if( matter == 0 )
	{
		int w =	param.Read<int>("width");
		int h = param.Read<int>("height");
		thermalizationSweeps = param.Read<int>("thermalization sweeps");
		measurementSweeps = param.Read<int>("measurement sweeps");
		triangulation.LoadRegularLattice(w,h);
	} else
	{
		n =	param.Read<int>("triangles");
	}
	double tau1 = param.Read<double>("tau_1");
	double tau2 = param.Read<double>("tau_2");
	std::pair<double,double> tau(tau1,tau2);
	double maxradius = param.Read<double>("maxradius");
	int bins = param.Read<int>("bins");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	CMinusTwoBuilder builder(&triangulation,1,n);
	if( matter != 0 )
	{
		triangulation.setDominantMatter( &builder );
		triangulation.DoSweep();
	}

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding harm( &triangulation, &cohom );

	ShortestLoop shortestloop( &triangulation, &cohom );
	ConstrainedLoopLength loopobservable( &harm, &shortestloop, tau, maxradius, bins );

	simulation.AddObservable( &loopobservable, measurementSweeps );

	simulation.SetDirectory("./output/");
	//simulation.SetDirectory("C:\\Users\\Timothy\\Documents\\research\\projects\\fractaldimensions\\data\\conformaldist\\local\\");

	simulation.Run();
	return 0;
}