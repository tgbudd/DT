#include "ConformalDistribution.h"

ConformalDistribution::ConformalDistribution(Embedding * const embedding) :
	embedding_(embedding),
	measurements_(0),
	minlograd_(-40.0),
	maxlograd_(0.0),
	bins_(1600)
{
	histogram_.resize(bins_,0);
}

ConformalDistribution::~ConformalDistribution()
{
}

void ConformalDistribution::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		std::vector<double> radius;
		
		if( embedding_->GetRadii(radius) )
		{
			for(int i=0,endi=radius.size();i<endi;i++)
			{
				if( radius[i] > 0.0 )
				{
					double logradius = std::log(radius[i]);
					if( logradius >= minlograd_ && logradius < maxlograd_ )
					{
						histogram_[static_cast<int>((logradius - minlograd_)/(maxlograd_-minlograd_)*bins_)]++;
					}
				}
			}
			measurements_++;
		}
	}
}

std::string ConformalDistribution::OutputData() const
{
	std::ostringstream stream;
	stream << std::fixed << "conformaldistribution -> {minlograd -> " << minlograd_;
	stream << ", maxlograd -> " << maxlograd_;
	stream << ", bins -> " << bins_;
	stream << ", measurements -> " <<  measurements_;
	stream << ", histogram -> ";
	PrintToStream(stream,histogram_.begin(),histogram_.end());
	stream << "}";
	return stream.str();
}