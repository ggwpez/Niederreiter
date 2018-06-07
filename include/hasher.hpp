#pragma once

#include <functional>
#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>

namespace std
{
	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
		seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
	}

	template<>
	struct hash<NTL::WordVector>
	{
		size_t operator()(const NTL::WordVector & x) const
		{
			return std::hash<bool>()(NTL::IsZero(x) ? false : true);
		}
	};

	template<>
	struct hash<NTL::GF2>
	{
		size_t operator()(const NTL::GF2 & x) const
		{
			return std::hash<bool>()(NTL::IsZero(x) ? false : true);
		}
	};

	template<>
	struct hash<NTL::GF2X>
	{
		size_t operator()(const NTL::GF2X & x) const
		{
			return std::hash<bool>()(x.);
		}
	};

	template<>
	struct hash<NTL::GF2E>
	{
		size_t operator()(const NTL::GF2E & x) const
		{
			return std::hash<NTL::GF2X>()(const_cast<NTL::GF2E&>(x).LoopHole());	// Well, thanks NTL man
		}
	};
}
