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

ZZ Binom::coeff(ZZ const& n, ZZ k)
{
	if (k > n)
		throw std::invalid_argument("Binomial coefficient (n k), k was greater than n");

	if (IsZero(k))
		return ZZ(1);
	if (2 *k > n)
		sub(k, n, k);

	ZZ res(1);

	for (ZZ i = ZZ(1); i <= k; ++i)
	{
		mul(res, res, n -k +i);
		div(res, res, i);
	}

	return res;
}

void Binom::encode(long const n, long const t, ZZ x, vec_GF2& out)
{
	if (n < t)
		std::invalid_argument("n < t");

	ZZ c = coeff(ZZ(n), ZZ(t));
	//ZZFromBytes(x, reinterpret_cast<unsigned char const*>(data), n);

	out = vec_GF2(INIT_SIZE, n);
	if (x >= c)
		throw std::invalid_argument("Value too large");

	long nn = n,
		 tt = t;

	for (long j = 0; j < n; ++j)
	{
		mul(c, c, nn -tt);
		div(c, c, nn);
		--nn;

		if (c <= x)
		{
			out.put(j, 1);
			sub(x, x, c);
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

	assert(IsZero(x));
	assert(! tt);
}

void Binom::decode(const long n, const long t, const vec_GF2& in, ZZ& out)
{
	char* enc;
	decode(n, t, in, enc);

	ZZFromBytes(out, reinterpret_cast<unsigned char*>(enc), n);
	free(enc);
}

///
/// https://www.cayrel.net/IMG/pdf/hymes_cw_buescher_meub.pdf
///	paragraph 3.1
///
void Binom::decode(long const n, long const t, vec_GF2 const& data, char*& out)
{
	assert(n == data.length());
	assert(width(data) == t);

	ZZ bc = coeff(ZZ(n), ZZ(t));
	ZZ d = ZZ(0);

	long nn = n,
		 tt = t;

	for (long i = 0; i < n; ++i)
	{
		bc = (bc *(nn -tt)) /nn;
		--nn;

		if (IsOne(data[i]))
		{
			d = d +bc;
			--tt;

			if (nn == tt)
				bc = ZZ(1);
			else
				bc = (bc *(tt +1)) /(nn -tt);
		}
	}

	out = new char[n]();
	BytesFromZZ(reinterpret_cast<unsigned char*>(out), d, n);
	/*ZZ res = ZZ::zero();

	long tt = 0;
	for (long i = 0; i < n; ++i)
		if (! IsZero(data[i]))
			res += coeff(ZZ(i), ZZ(++tt));

	assert(tt == t);
	out = new char[n];
	BytesFromZZ(reinterpret_cast<unsigned char*>(out), res, n);*/
}

void Binom::encode(const long n, const long t, char const* data, vec_GF2& out)
{
	ZZ x;
	ZZFromBytes(x, reinterpret_cast<unsigned char const*>(data), n);

	encode(n, t, x, out);
}

