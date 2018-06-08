#include "binom.hpp"
#include "helper.hpp"

#include <cassert>
#include <NTL/ZZ.h>

using namespace NTL;

static ZZ fac(long x)
{
	if (x < 2)
		return ZZ(1);
	else
		return x *fac(x -1);
}

static ZZ binom(ZZ const& n, ZZ k)
{
	if (k == ZZ::zero())
		return ZZ(1);
	if (2 *k > n)
		sub(k, n, k); // k = n -k;

	ZZ res;

	for (ZZ i = ZZ(1); i <= k; ++i)
	{
		res *= (n -k +i) /i;
	}

	return res;
}

void Binom::encode(long const n, long const t, char const* data, vec_GF2& out)
{
	ZZ c, i = binom(ZZ(n), ZZ(t));
	ZZFromBytes(c, reinterpret_cast<unsigned char const*>(data), n);

	if (i >= c)
		throw std::invalid_argument("Value too large");

	out.SetLength(n);

	long nn = n,
		 tt = t;

	for (long j = 0; j < n; ++j)
	{
		mul(c, c, nn -tt);
		div(c, c, nn);
		--nn;

		if (c <= i)
		{
			out.put(j, 1);
			sub(i, i, c);
			--tt;

			if (nn == tt)
				c = ZZ(1);
			else
			{
				mul(c, c, tt +1);
				div(c, c, nn -tt);
			}
		}
	}
}

///
/// https://www.cayrel.net/IMG/pdf/hymes_cw_buescher_meub.pdf
/// page 10 (paragraph 3.1)
///
void Binom::decode(long const n, long const t, vec_GF2 const& data, char*& out, long& size_bytes)
{
	assert(n == data.length());
	assert(t == width(data));
	ZZ res;

	for (long i = 0; i < t; ++i)
		if (! IsZero(data[i]))
			res += binom(ZZ(i), ZZ(i +1));

	long size_bits = NumBits(res);
	size_bytes = size_bits /8;
	assert(size_bits % sizeof(char) == 0 && "Im not going to return sub byte data");
	out = new char[size_bytes];
	BytesFromZZ(reinterpret_cast<unsigned char*>(out), res, size_bytes);
}
