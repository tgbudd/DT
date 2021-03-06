#ifndef DT_VERTEX_H
#define DT_VERTEX_H

#include "Edge.h"

class Vertex
{
public:
	Vertex(int id) : id_(id) {}
	Vertex(int id, Edge* parent) : id_(id), parent_(parent) {}
	~Vertex(void);

	const int & getId() const {
		return id_;
	}
	Edge * const & getParent() const {
		return parent_;
	}
	void setParent(Edge * const & parent)
	{
		parent_ = parent;
	}

	int getDegree() const;
private:
	int id_;
	Edge* parent_;
};

#endif
