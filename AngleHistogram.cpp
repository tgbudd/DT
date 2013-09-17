#include "AngleHistogram.h"



AngleHistogram::AngleHistogram(Embedding * const embedding, int bins, double min, double max)
	: embedding_(embedding), bins_(bins), min_(min), max_(max)
{
	histogram_.resize(bins,0);
	measurements_=0;
}


AngleHistogram::~AngleHistogram(void)
{
}

void AngleHistogram::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		for(int i=0,end=embedding_->getSize();i<end;i++)
		{
			for(int j=0;j<3;j++)
			{
				int val = (int)(bins_*(CalculateAngle(i,j)-min_)/(max_-min_));
				if( val >= 0 && val < bins_ )
				{
					histogram_[val]++;
				}
			}
		}
		measurements_++;
	}
}

std::string AngleHistogram::OutputData() const
{
	std::ostringstream stream;
	stream << "anglehistogram -> {measurements -> " << measurements_ << ", histogram -> ";
	PrintToStream(stream,histogram_.begin(),histogram_.end());
	stream << "}";
	return stream.str();
}

double AngleHistogram::CalculateAngle(int triangle, int edge) const
{
	Vector2D form1 = embedding_->getForm(triangle, (edge+2)%3);
	Vector2D form2 = NegateVector2D(embedding_->getForm(triangle, (edge+1)%3));
	return VectorAngle(form1, form2);
}
