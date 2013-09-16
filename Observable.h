#ifndef DT_OBSERVABLE_H
#define DT_OBSERVABLE_H

#include <string>

#include "Triangulation.h"

class Observable
{
public:
	Observable() {}
	~Observable() {}

	virtual void Measure() = 0;
	virtual std::string OutputData() const { return ""; }
};

#endif
