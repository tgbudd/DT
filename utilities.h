#ifndef DT_UTILITIES_H
#define DT_UTILITIES_H


#define PI 3.14159265358979

#include <iostream>
#include <sstream>
#include <vector>

#include "boost/array.hpp"

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
		std::cout << (first?"":",") << *it;
		first = false;
	}
	std::cout << "}\n";
}

template <typename Iter>
inline void PrintToStream(std::ostringstream & stream,Iter it, Iter end) {
    stream << "{";
	bool first = true;
	for (; it!=end; ++it) 
	{
		stream << (first?"":",") << *it;
		first = false;
	}
	stream << "}";
}

typedef boost::array<double,2> Vector2D;

inline Vector2D AddVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D sum(v1);
	sum[0] += v2[0];
	sum[1] += v2[1];
	return sum;
}

inline Vector2D SubtractVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D diff(v1);
	diff[0] -= v2[0];
	diff[1] -= v2[1];
	return diff;
}
inline Vector2D NegateVector2D( const Vector2D & v )
{
	Vector2D neg;
	neg[0] = -v[0];
	neg[1] = -v[1];
	return neg;
}
inline double SignedPolygonArea( const std::vector<Vector2D> & v )
{
	double area = 0.0;
	for(int i=0;i<v.size();i++)
	{
		area += 0.5 * ( v[i][0] * v[(i+1)%v.size()][1] - v[i][1] * v[(i+1)%v.size()][0] );
	}
	return area;
}

inline void ResizeAndAdd(std::vector<int> & v, int index, int value)
{
	if(index >= v.size())
	{
		v.resize(index+1,0);
	}
	v[index] += value;
}

#endif