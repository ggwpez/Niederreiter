#pragma once

#include <NTL/vec_GF2.h>

namespace Binom
{
	void encode(long const n, long const t, char const* data, NTL::vec_GF2& out);
	void decode(const long n, const long t, NTL::vec_GF2 const& in, char*& out, long& size_bytes);
}
