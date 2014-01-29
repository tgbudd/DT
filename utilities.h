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

inline Vector2D MakeVector2D( const double & x, const double & y )
{
	Vector2D v;
	v[0] = x;
	v[1] = y;
	return v;
}

inline Vector2D AddVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D sum(v1);
	sum[0] += v2[0];
	sum[1] += v2[1];
	return sum;
}

inline Vector2D AddScaledVectors2D( double a1, const Vector2D & v1, double a2, const Vector2D & v2 )
{
	Vector2D sum;
	sum[0] = a1*v1[0]+a2*v2[0];
	sum[1] = a1*v1[1]+a2*v2[1];
	return sum;
}

inline Vector2D SubtractVectors2D( const Vector2D & v1, const Vector2D & v2 )
{
	Vector2D diff(v1);
	diff[0] -= v2[0];
	diff[1] -= v2[1];
	return diff;
}
inline Vector2D ScaleVector2D( const Vector2D & v, double a )
{
	Vector2D v2;
	v2[0] = a*v[0];
	v2[1] = a*v[1];
	return v2;
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
inline double Norm2D( const Vector2D & v )
{
	return std::sqrt(NormSquared2D(v));
}

inline void ResizeAndAdd(std::vector<int> & v, int index, int value)
{
	if(index >= static_cast<int>(v.size()))
	{
		v.resize(index+1,0);
	}
	v[index] += value;
}

inline double DotProduct( const Vector2D & v1, const Vector2D & v2 )
{
	return v1[0] * v2[0] + v1[1] * v2[1];
}

inline double VectorAngle( const Vector2D & from, const Vector2D & to )
{
	// Returns the angle (between -PI and PI) that one has to rotate "from" 
	// in clockwise direction to align it with "to".
	return std::atan2( from[0] * to[1] - from[1] * to[0], from[0] * to[0] + from[1] * to[1] );
}

inline double VectorAngle( const Vector2D & v)
{
	// Returns the angle (between -PI and PI) between x-axis and v in counter-clockwise direction
	return std::atan2( v[1], v[0] );
}

inline Vector2D MobiusTransform(const Vector2D & z, const Vector2D & z0)
{
	// Transformation of the unit disk that maps z=z0 to 0.
	double normsqz0 = NormSquared2D(z0), normsqz = NormSquared2D(z);
	double denominator = 1.0 - 2.0 * DotProduct(z,z0) + normsqz * normsqz0;
	Vector2D y;
	y[0] = (z[0] - z0[0] - normsqz * z0[0] + z[0] * z0[0] * z0[0] - z[0] * z0[1] * z0[1] + 2.0 * z[1] * z0[0] * z0[1])/denominator;
	y[1] = (z[1] - z0[1] - normsqz * z0[1] - z[1] * z0[0] * z0[0] + z[1] * z0[1] * z0[1] + 2.0 * z[0] * z0[0] * z0[1])/denominator;
	return y;
}

inline Vector2D RotateVector2D(const Vector2D & v, double angle)
{
	Vector2D x;
	double sin = std::sin(angle), cos = std::cos(angle);
	x[0] = cos * v[0] - sin * v[1];
	x[1] = sin * v[0] + cos * v[1];
	return x;
}

inline Vector2D MobiusTransform(const Vector2D & z, const Vector2D & z0, double angle)
{
	return RotateVector2D(MobiusTransform(z,z0),angle);
}

inline Vector2D MobiusTransform(const Vector2D & z, const Vector2D & z0, const Vector2D & z1)
{
	// Mobius transform that maps z0 to 0 and z1 to the positive real axis
	return RotateVector2D(MobiusTransform(z,z0),-VectorAngle(MobiusTransform(z1,z0)));
}

inline Vector2D InverseMobiusTransform(const Vector2D & z, const Vector2D & z0, const Vector2D & z1)
{
	// Mobius transform that maps 0 to z0 and some point on the positive real axis to z1
	return MobiusTransform(RotateVector2D(z,VectorAngle(MobiusTransform(z1,z0))),NegateVector2D(z0));
}
inline double PoincareDistance(const Vector2D & z1, const Vector2D & z0)
{
	double r = Norm2D(MobiusTransform(z1,z0));
	return std::log(1+r) - std::log(1-r);
}
inline std::pair<Vector2D,double> EuclideanCircle(const Vector2D & z, double exprad)
{
	// Give the Euclidean center and radius of the hyperbolic circle around z with hyperbolic radius -log(exprad)/2
	double normsqz = NormSquared2D(z);
	double sqrtexprad = std::sqrt(exprad);
	double factor = 1.0/((1.0+sqrtexprad)*(1.0+sqrtexprad)-(1.0-sqrtexprad)*(1.0-sqrtexprad) * normsqz);
	double radius = factor * (1.0 - exprad) * (1.0-normsqz);
	return std::pair<Vector2D,double>(ScaleVector2D(z,1.0 - factor*(1.0-sqrtexprad)*(1.0-sqrtexprad)*(1.0-normsqz)),radius);
}

inline Vector2D CenterOfCircle(const boost::array<Vector2D,3> & v)
{
	Vector2D c;
	double norm0 = NormSquared2D(v[0]), norm1 = NormSquared2D(v[1]), norm2 = NormSquared2D(v[2]);
	double a = -v[0][1]*v[1][0] + v[0][0]*v[1][1] + v[0][1]*v[2][0] - v[1][1] * v[2][0] - v[0][0] * v[2][1] + v[1][0] * v[2][1]; 
	c[0] = (-v[0][1]*norm1 + norm0*v[1][1] + v[0][1]*norm2 - v[1][1] * norm2 - norm0 * v[2][1] + norm1 * v[2][1])/(2.0*a);
	c[1] = (-norm0*v[1][0] + v[0][0]*norm1 + norm0*v[2][0] - norm1 * v[2][0] - v[0][0] * norm2 + v[1][0] * norm2)/(2.0*a);
	return c;
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