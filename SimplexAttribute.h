#ifndef SIMPLEX_ATTRIBUTE_H
#define SIMPLEX_ATTRIBUTE_H

#include <vector>

#include "Triangulation.h"
#include "Vertex.h"
#include "Triangle.h"
#include "Edge.h"

template <typename simplex_t, typename T>
class SimplexAttribute {
public:
	typedef typename std::vector<T>::iterator iterator;
	iterator begin() 
	{
		return data_.begin();
	}
	iterator end()
	{
		return data_.end();
	}
	virtual T& operator[](const simplex_t * simplex)
	{
		return data_[GetIndex(simplex)];
	}
	virtual const T& operator[](const simplex_t * simplex) const
	{
		return data_[GetIndex(simplex)];
	}
	std::vector<T> & data()
	{
		return data_;
	}
protected:
	SimplexAttribute(size_t size) : data_(size)	{}
	SimplexAttribute(size_t size, const T & value) : data_(size,value) {}
private:
	virtual unsigned int GetIndex(const simplex_t * simplex) const = 0;
	virtual const simplex_t * GetSimplex(unsigned int index) const = 0;

	std::vector<T> data_;
};

template <typename simplex_t>
class SimplexAttribute<simplex_t,bool> {
public:
	typedef bool T;
	typedef std::vector<char>::iterator iterator;
	iterator begin() 
	{
		return data_.begin();
	}
	iterator end()
	{
		return data_.end();
	}
	virtual char& operator[](const simplex_t * simplex)
	{
		return data_[GetIndex(simplex)];
	}
	virtual bool operator[](const simplex_t * simplex) const
	{
		return !!data_[GetIndex(simplex)];
	}
protected:
	SimplexAttribute(size_t size) : data_(size)	{}
	SimplexAttribute(size_t size, const T & value) : data_(size,static_cast<char>(value)) {}
private:
	virtual unsigned int GetIndex(const simplex_t * simplex) const = 0;
	virtual const simplex_t * GetSimplex(unsigned int index) const = 0;

	std::vector<char> data_;
};

template <typename T>
class VertexAttribute : public SimplexAttribute<Vertex,T> {
public:
	VertexAttribute(const Triangulation * triangulation) :
		triangulation_(triangulation), 
		SimplexAttribute<Vertex,T>(triangulation->NumberOfVertices())
	{}
	VertexAttribute(const Triangulation * triangulation, const T & value) : 
		triangulation_(triangulation), 
		SimplexAttribute<Vertex,T>(triangulation->NumberOfVertices(),value)
	{}
private:
	unsigned int GetIndex(const Vertex * simplex) const
	{
		return static_cast<unsigned int>(simplex->getId());
	}
	const Vertex * GetSimplex(unsigned int index) const
	{
		return triangulation_->getVertex(index);
	}
	const Triangulation * const triangulation_;
};

template <typename T>
class TriangleAttribute : public SimplexAttribute<Triangle,T> {
public:
	TriangleAttribute(const Triangulation * triangulation) :
		triangulation_(triangulation),
		SimplexAttribute<Triangle,T>(triangulation->NumberOfTriangles())
	{}
	TriangleAttribute(const Triangulation * triangulation, const T & value) :
		triangulation_(triangulation),
		SimplexAttribute<Triangle,T>(triangulation->NumberOfTriangles(),value)
	{}
private:
	unsigned int GetIndex(const Triangle * simplex) const
	{
		return static_cast<unsigned int>(simplex->getId());
	}
	const Triangle * GetSimplex(unsigned int index) const
	{
		return triangulation_->getTriangle(index);
	}

	const Triangulation * const triangulation_;
};

template <typename T>
class EdgeAttribute : public SimplexAttribute<Edge,T> {
public:
	EdgeAttribute(const Triangulation * triangulation) :
		triangulation_(triangulation),
		SimplexAttribute<Edge,T>(3*triangulation->NumberOfTriangles())
	{}
	EdgeAttribute(const Triangulation * triangulation, const T & value) :
		triangulation_(triangulation),
		SimplexAttribute<Edge,T>(3*triangulation->NumberOfTriangles(),value)
	{}
private:
	unsigned int GetIndex(const Edge * simplex) const
	{
		return static_cast<unsigned int>(3*simplex->getParent()->getId()+simplex->getId());
	}
	const Edge * GetSimplex(unsigned int index) const
	{
		return triangulation_->getTriangle(index/3)->getEdge(index%3);
	}
	const Triangulation * const triangulation_;
};



#endif