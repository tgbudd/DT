#include "Histogram.h"
#include "utilities.h"

// Explicitly instantiate frequently used templates.
template class Histogram<double>;
template class Histogram<float>;
template class Histogram<int>;
template class ExtendableHistogram<int>;
template class ExtendableHistogram<double>;

template<typename T>
Histogram<T>::Histogram(T min, T max, int bins, bool calculateaverage) :
	bins_(bins),
	min_(min),
	max_(max),
	measurements_(0),
	total_(0.0),
	total_squared_(0.0),
	below_min_(0),
	calculateaverage_(calculateaverage)
{
	histogram_.resize(bins_,0);
}

template<typename T>
T Histogram<T>::GetMin() const
{
	return min_;
}

template<typename T>
T Histogram<T>::GetMax() const
{
	return max_;
}

template<typename T>
int Histogram<T>::GetBins() const
{
	return bins_;
}

template<typename T>
int Histogram<T>::Bin(T x) const
{
	return static_cast<int>(((x - min_)*bins_)/(max_-min_));
}

template<typename T>
bool Histogram<T>::Insert(T x)
{
	if( calculateaverage_ )
	{
		total_ += static_cast<double>(x);
		total_squared_ += static_cast<double>(x)*static_cast<double>(x);
	}
	measurements_++;

	int bin = Bin(x);
	if( bin < 0 )
	{
		below_min_++;
	} else if( bin < bins_ )
	{
		histogram_[bin]++;
		return true;
	}
	return false;
}

template<typename T>
int Histogram<T>::Cumulative(T value) const
{
	// return the number of measurements in the bins up to but excluding the bin containing value
	return CumulativeAtBin(Bin(value));
}

template<typename T>
int Histogram<T>::CumulativeAtBin(int bin) const
{
	// return the number of measurements in the bins up to but excluding bin
	if( bin > bins_ )
	{
		bin = bins_;
	}
	int cumul = below_min_;
	for(int i=0;i<bin;++i)
	{
		cumul += histogram_[i];
	}
	return cumul;
}

template<typename T>
std::vector<int> Histogram<T>::CumulToBin(const std::vector<int> & cumul) const
{
	// Returns a vector bins, where bins[i] is the bin at which the CumulativeAtBin(bins[i]) first exceeds cumul[i].
	// Assuming cumul is sorted!	
	std::vector<int> bins;
	bins.reserve(cumul.size());
	int num = below_min_;
	int current_bin = 0;
	for(int i=0,endi=cumul.size();i<endi;++i)
	{
		while( current_bin < bins_ && num < cumul[i] )
		{
			num += histogram_[current_bin];
			current_bin++;
		}
		bins.push_back(current_bin);
	}
	return bins;
}

template<typename T>
void Histogram<T>::Resize(int bins)
{
	bins_ = bins;
	histogram_.resize(bins_,0);
}

template<typename T>
void Histogram<T>::Reset()
{
	std::fill(histogram_.begin(),histogram_.end(),0);
	measurements_ = 0;
	total_ = 0.0;
	total_squared_ = 0.0;
	below_min_ = 0;
}

template<typename T>
void Histogram<T>::PrintTo(std::ostream & stream) const
{
	stream << "{min -> " << min_;
	stream << ", max -> " << max_;
	stream << ", bins -> " << bins_;
	stream << ", measurements -> " << measurements_;
	stream << ", belowmin -> " << below_min_;
	if( calculateaverage_ )
	{
		stream << ", average -> " << total_/measurements_;
		stream << ", variance -> " << (total_squared_ - total_*total_/measurements_)/measurements_;
	}
	stream << ", histogram -> ";
	PrintToStream( stream, histogram_.begin(), histogram_.end() );
	stream << "}";
}

template<typename T>
ExtendableHistogram<T>::ExtendableHistogram(T min, T delta) :
	Histogram<T>(min,min,0),
	delta_(delta)
{

}

template<typename T>
bool ExtendableHistogram<T>::Insert(T x)
{
	if( x >= this->max_ )
	{
		int bin = static_cast<int>((x - this->min_)/this->delta_);
		this->Resize(bin+1);
		this->max_ = (bin+1) * delta_;
	}
	return Histogram<T>::Insert(x);
}