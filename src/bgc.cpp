#include "bgc.hpp"
#include "helper.hpp"
#include <NTL/mat_GF2E.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <numeric>
#include <iomanip>

using namespace NTL;

bgc_t::bgc_t(long dimenion, long non_root_points, long errors)
	: /*d, m*/ dimension(dimenion), /*n*/ non_root_points(non_root_points), /*t*/ errors(errors)
{		// n k t
	if (! ((2 <= errors) && (errors <= ((1l << dimenion) -1) /dimenion)))
		throw std::invalid_argument("errors");
	if (! ((dimenion *errors +1 <= non_root_points) && (non_root_points <= (1l << dimenion))))
		throw std::invalid_argument("non_root_points");

	// TODO check 2 <= errors <= (non_root_points -1) /dimension
	// TODO needed?
	f = std::make_unique<GF2X>(NTL::BuildIrred_GF2X(dimenion));	// irreducible polynome of degree 10 over GF(2)
	GF2E::init(*f);
	SetCoeff(id, 1);

	g = NTL::BuildIrred_GF2EX(errors);	// TODO not only binary Coefficients, otherwise attackable with Liodreau-Sendrier

	// Generate L
	{
		L.SetLength(non_root_points);
		GF2X f_fix = *f;
		SetCoeff(f_fix, deg(f_fix), GF2::zero());
		f_fix.normalize();
		L[0] = gen = conv<GF2E>(f_fix);

		for (size_t i = 1; i < L.length() -1; ++i)
			L[i] = L[i -1] *gen;
		L[L.length() -1] = 0;
	}

	{
		mat_GF2E YZ, h;
		YZ.SetDims(errors, non_root_points);
		h.SetDims(errors, non_root_points);

		for (long col = 0; col < non_root_points; ++col)
			YZ[0][col] = inv(call(g, L[col]));

		for (long row = 1; row < errors; ++row)
			for (long col = 0; col < non_root_points; ++col)
				YZ[row][col] = YZ[row -1][col] *L[col];

		for (int i = 0; i < errors; ++i)
			for (int j = 0; j < non_root_points; ++j)
				for (int k = 0; k <= i; ++k)
					h[i][j] += YZ[k][j] *coeff(g, errors +k -i);

		H = trace_construct(h);
	}
}

GF2EX bgc_t::calculate_sc(const vec_GF2& c) const
{
	if (c.length() % deg(*f))
		throw std::invalid_argument("(c.length() % L.length()) was not NULL");
	NTL::GF2EX sc = NTL::GF2EX::zero();

	long b = c.length() /deg(*f);
	for (long i = 0; i < b; ++i)
	{
		GF2E e = GF2E::zero();

		for (long j = 0; j < deg(*f); ++j)
			if (! IsZero(c[i *deg(*f) +j]))
				SetCoeff(e.LoopHole(), j);

		SetCoeff(sc, b -i -1, e);
	}

	return sc;
}

GF2EX bgc_t::calculate_vc(const GF2EX& sc) const
{
	// Return the identity function of there are no errors
	if (IsZero(sc))
		throw std::runtime_error("sc not invertible");
	else
	{
		GF2EX v;

		// v = 1 / sc
		InvMod(v, sc, g);
		if (v == sc)
			return sc;
		// v = 1 / sc -x
		add(v, v, id);		// TODO try an add
		//std::cout << "Trying to square-root: " << print(v) << '\n';
		return sqr_root(v);
	}
}

GF2EX bgc_t::calculate_sigma(vec_GF2 e, long spare_i) const
{
	assert(size_t(e.length()) == L.length());
	GF2EX ret = GF2EX::zero();

	for (long i = 0; i < e.length(); ++i)
	{
		if (! IsZero(e[i]) && (spare_i != i))
		{
			GF2EX tmp = id;
			sub(tmp, tmp, L[i]);
			if (IsZero(ret))
				ret	= tmp;
			else
				MulMod(ret, ret, tmp, g);
		}
	}

	return ret;
}

GF2EX bgc_t::calculate_small_omega(vec_GF2 e) const
{
	assert(size_t(e.length()) == L.length());
	GF2EX ret = GF2EX::zero();

	for (long i = 0; i < e.length(); ++i)
	{
		if (! IsZero(e[i]))
		{
			GF2EX tmp = calculate_sigma(e, i);
			add(ret, ret, tmp);
		}
	}

	return ret;
}

std::vector<size_t> bgc_t::get_L_indices(const std::vector<GF2E>& l) const
{
	std::vector<size_t> ret(l.size(), 0);

	for (size_t	i = 0; i < l.size(); ++i)
	{
		auto found = std::find(L.begin(), L.end(), l[i]);

		if (found == L.end())
			throw std::runtime_error("Could not find element " +print(l[i]) +" in set L");
		else
			ret[i] = (found -L.begin());
	}

	return ret;
}

long bgc_t::k() const
{
	return non_root_points -errors *dimension;
}

vec_GF2 bgc_t::encode(const vec_GF2& msg) const
{
	/*if (msg.length() != k())
		throw std::invalid_argument("msg.length() != k()");

	// As easy as it gets
	return msg *G_XYZ;*/
	throw std::runtime_error("Not implemented");
}

vec_GF2 bgc_t::syndrom_decode_2(const vec_GF2& c) const
{
	GF2EX syndrome = calculate_sc(c);
	auto t = InvMod(syndrome, g);

	auto tau = t +conv<GF2EX>("[[][1]]");
	tau = sqr_root(tau);

	auto a2plusX2b = monice(calc_sigma(tau, g, g));

	vec_GF2 e;
	e.SetLength(non_root_points);
	for (long i = 0; i < e.length(); ++i)
		if (IsZero(call(a2plusX2b, L[i])))
			e[i] = GF2(1);

	return e;
}

vec_GF2 bgc_t::syndrom_decode(const vec_GF2& c) const
{
	/*if (uint64_t(c.length()) != (uint64_t(1) << dimension) || c.length() != L.size())
		throw std::invalid_argument("c.length() != (1 << dimension) || c.length() != L.size() was " +std::to_string(c.length()));
	if (is_codeword(c))
		return c;*/

	auto sc = calculate_sc(c);	// Calculate syndrom function
	auto s = calculate_vc(sc);
	auto split_me = monice(calc_sigma(s, g, g));
	auto roots = find_roots(split_me);
	//std::cout << "roots " << print(roots) << '\n';
	std::vector<size_t> B = get_L_indices(roots);
	//std::cout << "errors: " << print(B) << '\n';
	vec_GF2 e = calculate_error_vector(B);
	//std::cout << "error_vector: " << print(e) << '\n';

	//vec_GF2 corrected = c +e;

	//if (! is_codeword(corrected))
		//throw std::runtime_error("Internal decoding error; Corrected message was not a codeword");

	return e;
}

GF2EX bgc_t::berlekamp_massey(const vec_GF2& c) const
{
	auto s = calculate_sc(c);

	if (IsZero(s))
		throw std::runtime_error("sc not invertible");

	GF2EX gcd, sigma, omega;
	XGCD(gcd, sigma, omega,s,g);
	return sigma;
}

// (x)^(1/p) = x^(p^m-1)
GF2EX bgc_t::sqr_root(const GF2EX& p) const
{
	ZZ e = (ZZ(1) << (dimension *errors -1));

	auto i = PowerMod(p, e, g);
	if (SqrMod(i, g) != p)
		throw std::runtime_error("Could not find root");

	return i;
}

vec_GF2 bgc_t::calculate_error_vector(const std::vector<size_t>& L_indices) const
{
	vec_GF2 ret;
	ret.SetLength(L.length());

	for (size_t i = 0; i < L_indices.size(); ++i)
		ret[L_indices[i]] = GF2(1);

	return ret;
}

// L[a] * L[b] = gen^(a+1) * gen^(b+1) = gen^(a+b+2) = L[a+b+1]
size_t bgc_t::mul_L_elements(size_t a, size_t b)
{
	return (a+b+1) % L.length();
}

// L[a] ^ exp = gen^(a+1) ^ exp = gen^(exp*a+exp) = L[exp*a+exp-1]
size_t bgc_t::power_L_elements(size_t a, long exp)
{
	return (exp *a +exp -1) % L.length();
}

bool bgc_t::is_codeword(const vec_GF2& c) const
{
	//return (c.length() == (uint64_t(1) << dimension)) && IsZero(c *transpose(XYZ_bin));
	throw std::runtime_error("Not implemented");
}

ZZ bgc_t::GF_order() const
{
	return (ZZ(1) << (dimension *errors));
}

std::string bgc_t::to_str() const
{
	std::ostringstream ss;

	ss << "[n,k,d]-Code = [" << (uint64_t(1) << dimension) << ',' << k() << ',' << (2*errors +1) << "]" << std::endl
	   << "F(2^" << dimension << ") = " << "F(2)[x]/" << print(*f) << std::endl
	   << "Goppa Polynomial: " << print(g) << std::endl
	   << "Generator: " << print(gen) << std::endl
	   << "Support={ " << print(L) << '}'
	   << "\nH (" << H.NumRows() << 'x' << H.NumCols() << ")\n";

	return ss.str();
}
