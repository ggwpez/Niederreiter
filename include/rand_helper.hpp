#pragma once

#include <NTL/ZZ.h>
#include <random>

class MersenneTwisterStream : private NTL::RandomStream
{
public:
	MersenneTwisterStream();
};
