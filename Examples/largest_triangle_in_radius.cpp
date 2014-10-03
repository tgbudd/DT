#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "Histogram.h"
#include "CirclePacking.h"

class LargestTriangle : public Observable {
public:
	LargestTriangle(Triangulation * triangulation, CirclePacking * embedding, int samples, double maxInitialRadius)
		: triangulation_(triangulation), embedding_(embedding), samples_(samples),
		max_dist_(0.3),max_radius_ratio_(200),radius_bins_(200), max_initial_radius_(maxInitialRadius),
		min_initial_radius_(1.0e-6)
	{
		triangle_sizes_.resize(radius_bins_,Histogram<double>(0.0,static_cast<double>(max_radius_ratio_),radius_bins_));
	}
	void Measure();
	std::string OutputData() const;
private:
	void LookAround(const Vertex * v);

	CirclePacking * embedding_;
	Triangulation * triangulation_;
	std::pair<double,double> modulus_;
	int samples_;
	int max_radius_ratio_;
	double max_dist_;
	
	int radius_bins_;
	double max_initial_radius_;
	double min_initial_radius_;
	std::vector<Histogram<double> > triangle_sizes_;
	std::vector<double> radii_;
};

void LargestTriangle::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		modulus_ = embedding_->CalculateModuli();
		embedding_->GetAbsoluteRadii(radii_);

		int numsamples = 0;
		while(numsamples < samples_ )
		{
			const Vertex * vertex = triangulation_->getRandomVertex();
			double initialRadius = radii_[vertex->getId()];
			if( initialRadius < max_initial_radius_ &&
				initialRadius > min_initial_radius_ )
			{
				LookAround(vertex);
				numsamples++;
			}
		}
	}
}

double NormInDomain(Vector2D v, std::pair<double,double> modulus)
{
	v[0] = properfmod(v[0]+0.5,1.0)-0.5;
	v[1] = properfmod(v[1]+0.5,1.0)-0.5;

	double minsquared = 100.0;
	for(int x=-1;x<=1;++x)
	{
		for(int y=-1;y<=1;++y)
		{
			double normsq = NormSquaredTransformedByModulus(MakeVector2D(v[0]+x,v[1]+y),modulus);
			if( normsq < minsquared )
			{
				minsquared = normsq;
			}
		}
	}
	return std::sqrt(minsquared);
}

void LargestTriangle::LookAround(const Vertex * v)
{
	// Find 
	
	Vector2D x0 = embedding_->getCoordinate(v);
	double initialRadius = radii_[v->getId()];
	std::vector<double> maxsize(radius_bins_,0.0);
	for(int i=0,endi=triangulation_->NumberOfVertices();i<endi;++i)
	{
		Vector2D coor = embedding_->getCoordinate(i);
		double dist = NormInDomain(SubtractVectors2D(coor,x0),modulus_);
		if( dist < max_dist_ )
		{
			for(int bin = static_cast<int>(radius_bins_*dist/initialRadius/max_radius_ratio_),lastbin=std::min(radius_bins_,static_cast<int>(radius_bins_*max_dist_/initialRadius/max_radius_ratio_));bin<lastbin;++bin)
			{
				if( radii_[i] > maxsize[bin] )
				{
					maxsize[bin] = radii_[i];
				}
			}
		}
	}
	for(int i=0;i<radius_bins_;++i)
	{
		triangle_sizes_[i].Insert(maxsize[i]/initialRadius);
	}
}

std::string LargestTriangle::OutputData() const
{
	std::ostringstream os;
	os << "largesttriangle -> { samples -> " << std::fixed << samples_;
	os << ", maxinitialradius -> " << max_initial_radius_;
	os << ", mininitialradius -> " << min_initial_radius_;
	os << ", maxdist -> " << max_dist_;
	os << ", maxradiusratio -> " << max_radius_ratio_;
	os << ", radiusbins -> " << radius_bins_;
	os << ", trianglesizes -> {";
	for(int i=0;i<radius_bins_;++i)
	{
		os << (i>0?",":"");
		triangle_sizes_[i].PrintTo(os);
	}
	os << "}}";
	return os.str();
}


int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int w =	param.Read<int>("width");
	int h = param.Read<int>("height");
	double maxinitialradius = param.Read<double>("max initial radius");
	int samples = param.Read<int>("samples");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	std::string path = param.Read<std::string>("path (\"0\"=\"../../output/\")");
	if( path == "0" )
		path = "../../output/";
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);
	ConnectivityRestrictor conn( &triangulation, ConnectivityRestrictor::NO_DOUBLE_EDGES );
	triangulation.AddMatter( &conn );
	
	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	CirclePacking circle( &triangulation, &cohom );
	circle.SetAccuracy(1.0e-8);
	circle.setMaxIterations(50000);

	LargestTriangle largest( &triangulation, &circle, samples, maxinitialradius );

	simulation.AddObservable( &largest, measurementSweeps );
	simulation.SetDirectory(path);
	simulation.Run();
	return 0;
}
