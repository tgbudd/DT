#ifndef DT_HISTOGRAM_H
#define DT_HISTOGRAM_H

#include <vector>
#include <ostream>

template <typename T>
class Histogram {
public:
	Histogram(T min, T max, int bins);
	void Resize(int bins);
	virtual bool Insert(T x);
	void PrintTo(std::ostream & stream) const;

protected:
	int bins_;
	T min_, max_;

private:
	int measurements_;
	double total_, total_squared_;
	std::vector<int> histogram_;
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