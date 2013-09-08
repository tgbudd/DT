#pragma once

#include <string>

#include "triangulation.h"

class Observable
{
public:
	Observable() {}
	~Observable() {}

	virtual void Measure() = 0;
	virtual std::string OutputData() { return ""; }
private:
	const Triangulation * triangulation_;
};

