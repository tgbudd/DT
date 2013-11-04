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
	int VolumeEnclosed(const std::list<const Edge *> & boundary, bool thisSide) const;
	int VolumeEnclosed(std::list<const Edge *>::const_iterator begin, std::list<const Edge *>::const_iterator end, bool thisSide) const;
	void EnclosedTriangles(const std::list<const Edge *> & boundary, std::vector<Triangle *> & triangles, bool thisSide) const;
private:
	std::pair<int,bool> VolumeEnclosedSpherical(const std::list<const Edge *> & boundary) const;
	const Triangulation * const triangulation_;
	mutable ReusableFlag visited_;
	int genus_;
};

#endif