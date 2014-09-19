#include <iostream>
#include <sstream>
#include <iomanip>

#include "Simulation.h"
#include "Triangulation.h"
#include "utilities.h"
#include "CohomologyBasis.h"
#include "ConnectivityRestrictor.h"
#include "Histogram.h"
#include "HarmonicEmbedding.h"

class LargestTriangle : public Observable {
public:
	LargestTriangle(Triangulation * triangulation, Embedding * embedding, int samples)
		: triangulation_(triangulation), embedding_(embedding), samples_(samples),
		maxradius_(0.3), radius_bins_(150)
	{
		triangle_sizes_.resize(radius_bins_,Histogram<double>(0.0,maxradius_,radius_bins_));
	}
	void Measure();
	std::string OutputData() const;
private:
	void LookAround(Vertex * v);

	Embedding * embedding_;
	Triangulation * triangulation_;
	std::pair<double,double> modulus_;
	int samples_;
	double maxradius_;
	int radius_bins_;
	std::vector<Histogram<double> > triangle_sizes_;
	std::vector<double> radii_;
};

void LargestTriangle::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		modulus_ = embedding_->CalculateModuli();
		embedding_->GetRadii(radii_);

		for(int i=0;i<samples_;++i)
		{
			LookAround(triangulation_->getRandomVertex());
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

void LargestTriangle::LookAround(Vertex * v)
{
	Vector2D x0 = embedding_->getCoordinate(v);
	std::vector<double> maxsize(radius_bins_,0.0);
	for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;++i)
	{
		const Triangle * triangle = triangulation_->getTriangle(i);
		Vector2D coor = embedding_->GetCentroid(triangle);
		double dist = NormInDomain(SubtractVectors2D(coor,x0),modulus_);
		for(int bin = static_cast<int>(dist/maxradius_*radius_bins_);bin<radius_bins_;++bin)
		{
			if( radii_[i] > maxsize[bin] )
			{
				maxsize[bin] = radii_[i];
			}
		}
	}
	for(int i=0;i<radius_bins_;++i)
	{
		triangle_sizes_[i].Insert(maxsize[i]);
	}
}

std::string LargestTriangle::OutputData() const
{
	std::ostringstream os;
	os << "{ samples -> " << samples_;
	os << ", maxradius -> " << maxradius_;
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
	int method = param.Read<int>("method (0=SQUARE_ROOT_OF_AREA, 1=INSCRIBED_CIRCLE)");
	int samples = param.Read<int>("samples");
	int thermalizationSweeps = param.Read<int>("thermalization sweeps");
	int	measurementSweeps = param.Read<int>("measurement sweeps");
	int secondsperoutput = param.Read<int>("seconds per output");
	bool output = param.UserInput();

	Triangulation triangulation;
	triangulation.LoadRegularLattice(w,h);

	Simulation simulation( &triangulation, thermalizationSweeps, secondsperoutput, output );

	CohomologyBasis cohom( &triangulation );
	cohom.SetMakeUpToDateViaReinitialization(true);

	HarmonicEmbedding harm( &triangulation, &cohom );
	harm.SetAccuracy(1.0e-9);
	harm.setMaxIterations(20000);
	if( method == 0 )
	{
		harm.SetRadiusDefintion(HarmonicEmbedding::SQUARE_ROOT_OF_AREA);
		simulation.AddConfigurationInfo("method -> \"SQUARE_ROOT_OF_AREA\"");
	}else
	{
		harm.SetRadiusDefintion(HarmonicEmbedding::INSCRIBED_CIRCLE);
		simulation.AddConfigurationInfo("method -> \"INSCRIBED_CIRCLE\"");
	}

	LargestTriangle largest( &triangulation, &harm, samples );

	simulation.AddObservable( &largest, measurementSweeps );
	simulation.SetDirectory("./output/");
	simulation.Run();
	return 0;
}