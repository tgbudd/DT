#ifndef DT_HISTOGRAM_H
#define DT_HISTOGRAM_H

#include <vector>
#include <ostream>

template <typename T>
class Histogram {
public:
	Histogram(T min, T max, int bins, bool calculateaverage = true);
	void Resize(int bins);
	virtual bool Insert(T x);
	void PrintTo(std::ostream & stream) const;
	void Reset();
	int Cumulative(T value) const;
	int CumulativeAtBin(int bin) const;
	T GetMin() const;
	T GetMax() const;
	int GetBins() const;
	std::vector<int> CumulToBin(const std::vector<int> & cumul) const;
protected:
	int bins_;
	T min_, max_;

private:
	int Bin(T x) const;

	int measurements_;
	bool calculateaverage_;
	double total_, total_squared_;
	std::vector<int> histogram_;
	int below_min_;
};

template <typename T>
class ExtendableHistogram : public Histogram<T> {
public:
	ExtendableHistogram(T min, T delta);
	bool Insert(T x);
private:
	T delta_;
};

#endif