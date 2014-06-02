#include "Histogram.h"
#include "utilities.h"

// Explicitly instantiate frequently used templates.
template class Histogram<double>;
template class Histogram<float>;
template class Histogram<int>;
template class ExtendableHistogram<int>;
template class ExtendableHistogram<double>;

template<typename T>
Histogram<T>::Histogram(T min, T max, int bins) :
	bins_(bins),
	min_(min),
	max_(max),
	measurements_(0),
	total_(0.0),
	total_squared_(0.0)
{
	histogram_.resize(bins_,0);
};

template<typename T>
bool Histogram<T>::Insert(T x)
{
	total_ += static_cast<double>(x);
	total_squared_ += static_cast<double>(x)*static_cast<double>(x);
	measurements_++;

	int bin = static_cast<int>(((x - min_)*bins_)/(max_-min_));
	if( bin >= 0 && bin < bins_ )
	{
		histogram_[bin]++;
		return true;
	}
	return false;
};

template<typename T>
void Histogram<T>::Resize(int bins)
{
	bins_ = bins;
	histogram_.resize(bins_,0);
};

template<typename T>
void Histogram<T>::PrintTo(std::ostream & stream) const
{
	stream << "{min -> " << min_;
	stream << ", max -> " << max_;
	stream << ", bins -> " << bins_;
	stream << ", measurements -> " << measurements_;
	stream << ", average -> " << total_/measurements_;
	stream << ", variance -> " << (total_squared_ - total_*total_/measurements_)/measurements_;
	stream << ", histogram -> ";
	PrintToStream( stream, histogram_.begin(), histogram_.end() );
	stream << "}";
};

template<typename T>
ExtendableHistogram<T>::ExtendableHistogram(T min, T delta) :
	Histogram<T>(min,min,0),
	delta_(delta)
{

};

template<typename T>
bool ExtendableHistogram<T>::Insert(T x)
{
	if( x >= max_ )
	{
		int bin = static_cast<int>((x - this->min_)/this->delta_);
		Resize(bin+1);
		this->max_ = (bin+1) * delta_;
	}
	return Histogram<T>::Insert(x);
};