#pragma once

#include <vector>

#include "triangulation.h"
#include "Vertex.h"
#include "Edge.h"

typedef boost::array<double,2> Vector2D;

class Embedding
{
public:
	Embedding() {}
	~Embedding() {}

	virtual bool FindEmbedding() = 0;

	const Vector2D & getCoordinate(Vertex * const & vertex) const
	{
		return coordinate_[vertex->getId()];
	}
	const Vector2D & getCoordinate(int i) const
	{
		return coordinate_[i];
	}
	void setCoordinate(Vertex * const & vertex, Vector2D coor) 
	{
		coordinate_[vertex->getId()] = coor;
	}
	void setCoordinate(Vertex * const & vertex, int i, double coor) 
	{
		coordinate_[vertex->getId()][i] = coor;
	}
	void setCoordinate(int vertex, int i, double coor) 
	{
		coordinate_[vertex][i] = coor;
	}
	const Vector2D & getForm(Edge * const & edge) const
	{
		return form_[edge->getParent()->getId()][edge->getId()];
	}
	const Vector2D & getForm(int triangle, int edge) const
	{
		return form_[triangle][edge];
	}
	void setForm(Edge * const & edge, const Vector2D & coor )
	{
		 form_[edge->getParent()->getId()][edge->getId()] = coor;
	}
	void setForm(Edge * const & edge, int i, double coor )
	{
		 form_[edge->getParent()->getId()][edge->getId()][i] = coor;
	}
	void setForm(int triangle, int edge, int i, double coor )
	{
		 form_[triangle][edge][i] = coor;
	}
	void setSize(int NumberOfTriangles,int NumberOfVertices)
	{
		if( coordinate_.size() != NumberOfVertices )
		{
			coordinate_.resize(NumberOfVertices);
		}
		if( form_.size() != NumberOfTriangles )
		{
			form_.resize(NumberOfTriangles );
		}
	}
private:
	std::vector<Vector2D> coordinate_;
	std::vector<boost::array<Vector2D,3> > form_;
};
