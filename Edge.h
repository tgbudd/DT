#pragma once

#include "triangulation.h"
#include "Triangle.h"
#include "Vertex.h"

class Edge
{
public:
	Edge() : adjacent_(NULL) {}
	Edge(Triangle * parent, int id) : parent_(parent), adjacent_(NULL), id_(id) {}

	~Edge(void);

	Edge * const & getNext() const {
		return next_;
	}
	void setNext(Edge * const next) {
		next_ = next;
	}
	Edge * const & getPrevious()  const {
		return previous_;
	}
	void setPrevious(Edge * const previous) {
		previous_ = previous;
	}
	Edge * const & getAdjacent() const {
		return adjacent_;
	}
	void setAdjacent(Edge * const adjacent) {
		adjacent_ = adjacent;
	}
	void bindAdjacent(Edge * const adjacent) {
		setAdjacent(adjacent);
		adjacent->setAdjacent(this);
	}

	Triangle * const & getParent() const {
		return parent_;
	}

	Vertex * const & getOpposite() const {
		return opposite_;
	}
	void setOpposite(Vertex * const & opposite) {
		opposite_ = opposite;
	}

	const int & getId() const {
		return id_;
	}

	bool IsFlipMovePossible() const;
	void DoFlipMove();



private:
	Edge *next_,		// points towards next edge in the triangle (edges are directed anti-clockwise)
		 *previous_,	// points towards previous edge
		 *adjacent_;	// points towards the edge in the neighbouring triangle

	int id_;			// position within the triangle

	Triangle *parent_;

	Vertex *opposite_;  // the vertex of parent_ opposite this edge
};

