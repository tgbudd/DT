#ifndef DT_TRIANGLE_H
#define DT_TRIANGLE_H

#include <boost/array.hpp>

class Edge;

class Triangle
{
public:
	Triangle();
	Triangle(const Triangle & triangle);
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

#endif
