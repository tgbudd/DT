#include <cmath>
#include <vector>

#include "utilities.h"
#include "ModuliObservable.h"

void ToFundamentalDomain( std::pair< double, double > & x)
{
	if( x.second < 0.0 )
		x.second = -x.second;
		
	x.first = fmod(200.0+x.first+0.5,1.0)-0.5;
	double absv = x.first * x.first + x.second * x.second;

	while( absv < 1.0 )
	{
		x.first = fmod(200.0 -x.first / absv + 0.5,1.0) -0.5;
		x.second = x.second / absv;
		absv = x.first * x.first + x.second * x.second;
	}
}

ModuliObservable::ModuliObservable(Embedding * const embedding) : embedding_(embedding)
{
	tau1_min_ = -0.5;
	tau1_max_ = 0.5;
	tau1_bins_ = 10;
	tau2_min_ =  0.85;
	tau2_max_ = 10.00;
	tau2_bins_ = (1000 - 85)/5;
	measurements_ = 0;

	histogram_.resize(tau1_bins_,std::vector<int>(tau2_bins_,0));
}

ModuliObservable::~ModuliObservable(void)
{
}

void ModuliObservable::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		std::pair<double,double> tau = embedding_->CalculateModuli();
		ToFundamentalDomain(tau);

		int tau1bin = (int)(tau1_bins_*(tau.first - tau1_min_)/(tau1_max_ - tau1_min_));
		int tau2bin = (int)(tau2_bins_*(tau.second - tau2_min_)/(tau2_max_ - tau2_min_));

		if( tau1bin >= 0 && tau1bin < tau1_bins_ && tau2bin >=0 && tau2bin < tau2_bins_ )
		{
			histogram_[tau1bin][tau2bin]++;
		}else
		{
			int hy=0;
		}
		measurements_++;
	}
}

std::string ModuliObservable::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "moduliobservable -> {tau1min -> " << tau1_min_;
	stream << ", tau1max -> " << tau1_max_;
	stream << ", tau1bins -> " << tau1_bins_;
	stream << ", tau2min -> " << tau2_min_;
	stream << ", tau2max -> " << tau2_max_;
	stream << ", tau2bins -> " << tau2_bins_;
	stream << ", measurements -> " <<  measurements_;
	stream << ", histogram -> ";
	PrintToStream2D(stream,histogram_.begin(),histogram_.end());
	stream << "}";
	return stream.str();
}