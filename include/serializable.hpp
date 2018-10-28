#pragma once

#include <istream>
#include <ostream>

struct ISerializable
{
	inline ISerializable()
	{

	}

	virtual inline ~ISerializable()
	{ }

	virtual void serialize(std::ostream&) const = 0;
	virtual void deserialize(std::istream&) = 0;
};
