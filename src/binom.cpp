#include "binom.hpp"

#include <NTL/ZZ.h>

using namespace NTL;

static ZZ fac(long x)
{
	if (x < 2)
		return ZZ(1);
	else
		return x *fac(x -1);
}

void Binom::encode(long const n, long const t, char const* data, vec_GF2& out)
{

}

///
/// https://www.cayrel.net/IMG/pdf/hymes_cw_buescher_meub.pdf
/// page 10 (paragraph 3.1)
///
void Binom::decode(long const n, long const t, vec_GF2 const& out, char*& data, long& size)
{
	ZZ s;
}
