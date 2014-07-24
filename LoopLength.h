#pragma once
#include "Observable.h"
class LoopLength :
	public Observable
{
public:
	LoopLength(void);
	~LoopLength(void);

	void Measure() {}
};

