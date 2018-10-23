#pragma once

#include <NTL/vec_GF2.h>

namespace Binom
{
	NTL::ZZ coeff(NTL::ZZ const& n, NTL::ZZ k);

	void encode(long const n, long const t, NTL::ZZ x, NTL::vec_GF2& out);
	void encode(long const n, long const t, const char* data, NTL::vec_GF2& out);
	void decode(long const n, long const t, NTL::vec_GF2 const& in, NTL::ZZ& out);
	void decode(long const n, long const t, NTL::vec_GF2 const& in, char*& out);
}
