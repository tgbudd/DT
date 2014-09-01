#ifndef TRIANGULATION_PROPERTIES_H
#define TRIANGULATION_PROPERTIES_H

#include <vector>

#include "Triangulation.h"
#include "SimplexAttribute.h"

namespace properties {

void DegreeList(const Triangulation * const triangulation, std::vector<int> & degree);
void VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, std::vector<int> & distance);
void VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, VertexAttribute<int> & distance);
int VertexDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex);
void TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, std::vector<int> & distance);
void TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, TriangleAttribute<int> & distance);
int TriangleDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle);
void VertexWeightedDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, const EdgeAttribute<double> & weight, VertexAttribute<std::pair<double,int> > & distance);
void TriangleWeightedDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, const EdgeAttribute<double> & weight, TriangleAttribute<std::pair<double,int> > & distance);
std::pair<double,int> VertexWeightedDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex, const EdgeAttribute<double> & weight);
std::pair<double,int> TriangleWeightedDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle, const EdgeAttribute<double> & weight);

}

#endif