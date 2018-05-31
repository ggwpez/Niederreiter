#include "helper.hpp"
#include "printer.hpp"

#include <cassert>

using namespace NTL;

// TODO poly by reference but normalize is member only…
void monice(GF2EX& ret, GF2EX poly)
{
	if (IsZero(poly))
	{
		ret = poly;
	}
	else
	{
		poly.normalize();
		if (IsOne(LeadCoeff(poly)))
		{
			ret = poly;
		}
		else
		{
			//GF2E ic = inv(LeadCoeff(poly));
			//mul(ret, poly, ic);
			div(ret, poly, LeadCoeff(poly));
		}
	}
}

GF2EX monice(const GF2EX& poly)
{
	GF2EX ret;
	monice(ret, poly);
	return ret;
}

std::vector<GF2E> find_roots(const GF2EX& poly)	// TODO only support L needed to be tested
{
	std::vector<GF2E> roots;
	GF2E x = GF2E::zero();
	assert(GF2E::modulus().n < 64);

	for (long i = 0; i < (1l << GF2E::modulus().n); ++i)
	{
		x = generate_GF2E(i);
		GF2E y = call(poly, x);
		if (IsZero(y))
			roots.push_back(x);

		//std::cout << "poly(" << print(x) << ") = " << print(y) << std::endl;
	}

	return roots;
}

GF2E call(GF2EX const& p, const GF2E& x)
{
	GF2E res = GF2E::zero();

	for (int i = 0; i <= deg(p); ++i)
	{
		GF2E tmp;
		GF2E c = coeff(p, i);

		if (! IsZero(c))
		{
			power(tmp, x, i);
			if (! IsOne(c))
				mul(tmp, tmp, c);
			add(res, res, tmp);
		}
	}

	return res;
}
