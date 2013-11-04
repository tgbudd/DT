#ifndef DT_UTILITIES_H
#define DT_UTILITIES_H


#define PI 3.14159265358979

#include <iostream>
#include <sstream>
#include <vector>
#include <string> 
#include <cmath>

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

template <typename Iter>
inline void PrintToStream2D(std::ostringstream & stream,Iter it, Iter end) {
    stream << "{";
	bool first = true;
	for (; it!=end; ++it) 
	{
		stream << (first?"":",");
		PrintToStream( stream, it->begin(), it->end() );
		first = false;
	}
	stream << "}";
}

inline std::string RemoveNewline(std::string str)
{
	str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
	return str;
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
	for(int i=0,end=v.size();i<end;i++)
	{
		area += 0.5 * ( v[i][0] * v[(i+1)%static_cast<int>(v.size())][1] - v[i][1] * v[(i+1)%static_cast<int>(v.size())][0] );
	}
	return area;
}
inline double NormSquared2D( const Vector2D & v )
{
	return v[0]*v[0] + v[1]*v[1];
}

inline void ResizeAndAdd(std::vector<int> & v, int index, int value)
{
	if(index >= static_cast<int>(v.size()))
	{
		v.resize(index+1,0);
	}
	v[index] += value;
}

inline double VectorAngle( const Vector2D & from, const Vector2D & to )
{
	// Returns the angle (between -PI and PI) that one has to rotate "from" 
	// in clockwise direction to align it with "to".
	return atan2( from[0] * to[1] - from[1] * to[0], from[0] * to[0] + from[1] * to[1] );
}


class ParameterStream {
public:
	ParameterStream(int argc, char** argv) : argc_(argc), argv_(argv), current_(1) {}
	template<class T> T Read( std::string name )
	{
		T t;
		if( current_ < argc_ )
		{
			std::istringstream is(argv_[current_]);
			is >> t;
			std::cout << name << " = " << t << "\n";
		} else
		{
			std::cout << name << " = ";
			std::cin >> t;
		}
		current_++;
		return t;
	}
	bool UserInput() {
		return current_ > argc_;
	}
private: 
	int current_;
	int argc_;
	char** argv_;
};

namespace math {
	template<class T> 
	T square(const T & t)
	{
		return t*t;
	}
}


class ReusableFlag
{
public:
	ReusableFlag(int n) : current_flag_(0) {
		flag_.resize(n,0);
	}
	void Reset() {
		current_flag_++;
		if( current_flag_ == 0 )
		{
			std::fill(flag_.begin(),flag_.end(),0);
			current_flag_=1;
		}
	}
	bool isSet(int n) {
		return flag_[n] == current_flag_;
	}
	void Set(int n, bool f=true) {
		flag_[n] = (f ? current_flag_ : current_flag_-1 );
	}
private:
	unsigned int current_flag_;
	std::vector<unsigned int> flag_;
};

#endif