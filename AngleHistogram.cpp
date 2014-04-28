#include "AngleHistogram.h"
#include "TriangulationProperties.h"


AngleHistogram::AngleHistogram(Embedding * const embedding, Triangulation * triangulation, int bins, double min, double max)
	: embedding_(embedding), triangulation_(triangulation), bins_(bins), min_(min), max_(max),
	min_shear_(0.0), max_shear_(12.5), shear_bins_(250),
	min_log_length_(-15.0), max_log_length_(2.0), log_length_bins_(170),
	min_log_area_(-25.0), max_log_area_(0.0), log_area_bins_(250)
{
	histogram_.resize(bins,0);
	measurements_=0;
	measure_other_ = (triangulation_ != NULL);
	if( measure_other_ )
	{
		shear_histogram_.resize(shear_bins_,0);
		theta_histogram_.resize(bins_,0);
		log_area_histogram_.resize(log_area_bins_,0);
		log_length_histogram_.resize(log_length_bins_,0);
	}
}


AngleHistogram::~AngleHistogram(void)
{
}

void AngleHistogram::Measure()
{
	if( embedding_->IsUpToDate() || embedding_->MakeUpToDate() )
	{
		modulus_ = embedding_->CalculateModuli();
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

		if( measure_other_ )
		{
			for(int i=0,endi=triangulation_->NumberOfTriangles();i<endi;i++)
			{
				Triangle * triangle = triangulation_->getTriangle(i);
				for(int j=0;j<3;j++)
				{
					Edge * edge = triangle->getEdge(j);
					if( edge->getAdjacent()->getParent()->getId() > i )
					{
						double angle1 = CalculateAngle(i,j);
						double angle2 = CalculateAngle(edge->getAdjacent()->getParent()->getId(),edge->getAdjacent()->getId());
						double theta = PI - angle1 - angle2;
						int val = (int)(bins_*(theta-min_)/(max_-min_));
						if( val >= 0 && val <bins_ )
						{
							theta_histogram_[val]++;
						}
						double shear = std::fabs(CalculateShear(edge));
						int val2 = (int)(shear_bins_*(shear-min_shear_)/(max_shear_-min_shear_));
						if( val2 >= 0 && val2 < shear_bins_ )
						{
							shear_histogram_[val2]++;
						}
						double loglength = CalculateLogLength(i,j);
						int val3 = (int)(log_length_bins_*(loglength-min_log_length_)/(max_log_length_-min_log_length_));
						if( val3 >= 0 && val3 < log_length_bins_ )
						{
							log_length_histogram_[val3]++;
						}
					}
				}
				double logarea = CalculateLogArea(i);
				int val = (int)(log_area_bins_*(logarea-min_log_area_)/(max_log_area_-min_log_area_));
				if( val >= 0 && val < log_area_bins_ )
				{
					log_area_histogram_[val]++;
				}
			}
			std::vector<int> degree;
			properties::DegreeList(triangulation_,degree);
			for(int i=0,endi=degree.size();i<endi;i++)
			{
				if( degree[i] > static_cast<int>(degree_histogram_.size()) )
				{
					degree_histogram_.resize(degree[i],0);
				}
				degree_histogram_[degree[i]-1]++;
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
	stream << ", thetahistogram -> ";
	PrintToStream(stream,theta_histogram_.begin(),theta_histogram_.end());
	stream << ", shearbins -> " << shear_bins_ << ", shearmin -> " << min_shear_ << ", shearmax -> " << max_shear_ << ", shearhistogram -> ";
	PrintToStream(stream,shear_histogram_.begin(),shear_histogram_.end());
	stream << ", loglengthbins -> " << log_length_bins_ << ", loglengthmin -> " << min_log_length_ << ", loglengthmax -> " << max_log_length_ << ", loglengthhistogram -> ";
	PrintToStream(stream,log_length_histogram_.begin(),log_length_histogram_.end());
	stream << ", logareabins -> " << log_area_bins_ << ", logareamin -> " << min_log_area_ << ", logareamax -> " << max_log_area_ << ", logareahistogram -> ";
	PrintToStream(stream,log_area_histogram_.begin(),log_area_histogram_.end());
	stream << ", degreehist -> ";
	PrintToStream(stream,degree_histogram_.begin(),degree_histogram_.end());
	stream << "}";
	return stream.str();
}

double AngleHistogram::CalculateAngle(int triangle, int edge) const
{
	Vector2D form1 = TransformByModulus(embedding_->getForm(triangle, (edge+2)%3),modulus_);
	Vector2D form2 = TransformByModulus(NegateVector2D(embedding_->getForm(triangle, (edge+1)%3)),modulus_);
	return std::fabs(VectorAngle(form1, form2));
}

double AngleHistogram::CalculateShear(const Edge * edge) const
{
	double x = NormSquaredTransformedByModulus(embedding_->getForm(edge->getNext()),modulus_);
	x /= NormSquaredTransformedByModulus(embedding_->getForm(edge->getPrevious()),modulus_);
	x *= NormSquaredTransformedByModulus(embedding_->getForm(edge->getAdjacent()->getNext()),modulus_);
	x /= NormSquaredTransformedByModulus(embedding_->getForm(edge->getAdjacent()->getPrevious()),modulus_);
	return 0.5*std::log(x);
}

double AngleHistogram::CalculateLogLength(int triangle, int edge) const
{
	return 0.5*std::log(NormSquaredTransformedByModulus(embedding_->getForm(triangle,edge),modulus_));
}

double AngleHistogram::CalculateLogArea(int triangle) const
{
	Vector2D form1 = TransformByModulus(embedding_->getForm(triangle, 0),modulus_);
	Vector2D form2 = TransformByModulus(NegateVector2D(embedding_->getForm(triangle, 1)),modulus_);
	return std::log(0.5*std::fabs(form1[0]*form2[1] - form1[1]*form2[0])/modulus_.second);
}