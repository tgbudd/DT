#include "ThetaHistogram.h"
#include "utilities.h"

ThetaHistogram::ThetaHistogram(const ThetaModel * const thetamodel, int bins, double min, double max)
	: thetamodel_(thetamodel), bins_(bins), min_(min), max_(max)
{
	histogram_.resize(bins,0);
	measurements_=0;
}

ThetaHistogram::~ThetaHistogram(void)
{
}

void ThetaHistogram::Measure()
{
	for(int i=0,end=thetamodel_->getSize();i<end;i++)
	{
		for(int j=0;j<3;j++)
		{
			int val = (int)(bins_*(thetamodel_->getRealTheta(i,j)-min_)/(max_-min_));
			if( val >= 0 && val < bins_ )
			{
				histogram_[val]++;
			}
		}
	}
	measurements_++;
}

std::string ThetaHistogram::OutputData() const
{
	std::ostringstream stream;
	stream << "thetahistogram -> {measurements -> " << measurements_ << ", histogram -> ";
	PrintToStream(stream,histogram_.begin(),histogram_.end());
	stream << "}";
	return stream.str();
}