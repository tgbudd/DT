#define PI 3.14159265358979

#include <iostream>

inline double properfmod(double x,double a)
{
	return fmod(fmod(x,a)+a,a);
}

template <typename Iter>
inline void PrintToCout(Iter it, Iter end) {
    std::cout << "{";
	bool first = true;
	for (; it!=end; ++it) 
	{
		std::cout << (first?"",",") << *it;
	}
	std::cout << "}\n";
}