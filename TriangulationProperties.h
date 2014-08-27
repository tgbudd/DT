#ifndef TRIANGULATION_PROPERTIES_H
#define TRIANGULATION_PROPERTIES_H

#include "Triangulation.h"

namespace properties {

void DegreeList(const Triangulation * const triangulation, std::vector<int> & degree);
void VertexDistanceList(const Triangulation * const triangulation, const Vertex * const startVertex, std::vector<int> & distance);
int VertexDistance(const Triangulation * const triangulation, const Vertex * const startVertex, const Vertex * const endVertex);
void TriangleDistanceList(const Triangulation * const triangulation, const Triangle * const startTriangle, std::vector<int> & distance);
int TriangleDistance(const Triangulation * const triangulation, const Triangle * const startTriangle, const Triangle * const endTriangle);

}

#endif