#pragma once

#include <boost/array.hpp>

class Edge;

class Triangle
{
public:
	Triangle();
	~Triangle();

	const int& getId() const {
		return id_;
	}
	void setId(const int& id) {
		id_ = id;
	}
	Edge * const & getEdge(int i) const {
		return edges_[i];
	}
private:
	boost::array<Edge *,3> edges_;
	int id_;
};

