#pragma once

#include <functional>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <NTL/GF2.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>

#define HASH_OF(x) std::hash<decltype(x)>()(x)

namespace std
{
	// From Boost
	template <class T>
	inline void hash_combine(std::size_t& seed, T const& v)
	{
		seed ^= hash<T>()(v) + 0x9e3779b9 +(seed <<6) +(seed >>2);
	}

	template<>
	struct hash<NTL::WordVector>
	{
		size_t operator()(NTL::WordVector const& x) const
		{
			size_t h = 0;

			for (long i = 0; i < x.length(); ++i)
				hash_combine(h, x[i]);

			return h;
		}
	};

	template<>
	struct hash<NTL::GF2>
	{
		size_t operator()(NTL::GF2 const& x) const
		{
			return std::hash<bool>()(NTL::IsZero(x) ? false : true);
		}
	};

	template<>
	struct hash<NTL::GF2X>
	{
		size_t operator()(NTL::GF2X const& x) const
		{
			return std::hash<NTL::WordVector>()(x.xrep);
		}
	};

	template<>
	struct hash<NTL::GF2E>
	{
		size_t operator()(NTL::GF2E const& x) const
		{
			return std::hash<NTL::GF2X>()(const_cast<NTL::GF2E&>(x).LoopHole());  // LoopHole() has no const overload, thanks NTL man
		}
	};

	template<>
	struct hash<NTL::GF2EX>
	{
		size_t operator()(NTL::GF2EX const& x) const
		{
			size_t h = 0;

			for (long i = 0; i <= deg(x); ++i)
				hash_combine(h, NTL::coeff(x, i));

			return h;
		}
	};

	template<typename T, template<typename> class V>
	struct hash<V<T>>
	{
		size_t operator()(V<T> const& x) const
		{
			size_t h = 0;

			for (T const& e : x)
				hash_combine(h, e);

			return h;
		}
	};
}
