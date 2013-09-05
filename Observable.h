#pragma once

#include <string>

#include "triangulation.h"

class Observable
{
public:
	Observable() {}
	~Observable() {}

	virtual void Measure() = 0;
	std::string OutputData() {}
private:
	const Triangulation * triangulation_;
};

