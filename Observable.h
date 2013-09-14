#pragma once

#include <string>

#include "triangulation.h"

class Observable
{
public:
	Observable() {}
	~Observable() {}

	virtual void Measure() = 0;
	virtual std::string OutputData() const { return ""; }
};

