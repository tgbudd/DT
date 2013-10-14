#ifndef BABY_UNIVERSE_DETECTOR_H
#define BABY_UNIVERSE_DETECTOR_H

#include "utilities.h"
#include "Triangulation.h"  

class BabyUniverseDetector
{
public:
	BabyUniverseDetector(const Triangulation * triangulation);
	~BabyUniverseDetector();

	std::pair<int,bool> VolumeEnclosed(const std::list<const Edge *> & boundary) const;
	void EnclosedTriangles(const std::list<const Edge *> & boundary, std::vector<Triangle *> & triangles, bool thisSide) const;
private:
	std::pair<int,bool> VolumeEnclosedSpherical(const std::list<const Edge *> & boundary) const;
	const Triangulation * const triangulation_;
	mutable ReusableFlag visited_;
	int genus_;
};

#endif