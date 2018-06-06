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

bgc_t::bgc_t(long dimenion, long non_root_points, long errors, NTL::GF2X* supposed_f, NTL::GF2EX* supposed_g)
	: /*d, m*/ dimension(dimenion), /*n*/ non_root_points(non_root_points), /*t*/ errors(errors)
{		// n k t
	if (! ((2 <= errors) && (errors <= ((1l << dimenion) -1) /dimenion)))
		throw std::invalid_argument("errors");
	if (! ((dimenion *errors +1 <= non_root_points) && (non_root_points <= (1l << dimenion))))
		throw std::invalid_argument("non_root_points");

	// TODO check 2 <= errors <= (non_root_points -1) /dimension
	// TODO needed?
	if (supposed_f)
		f = std::make_unique<GF2X>(*supposed_f);
	else
		f = std::make_unique<GF2X>(NTL::BuildIrred_GF2X(dimenion));	// irreducible polynome of degree 10 over GF(2)
	GF2E::init(*f);
	SetCoeff(id, 1);

	if (supposed_g)
		g = *supposed_g;
	else
		g = NTL::BuildIrred_GF2EX(errors);	// TODO not only binary Coefficients, otherwise attackable with Liodreau-Sendrier

	// Generate L
	{
		L.resize(non_root_points);
		GF2X f_fix = *f;
		SetCoeff(f_fix, deg(f_fix), GF2::zero());
		f_fix.normalize();
		gen = conv<GF2E>(f_fix);

		for (size_t i = 0; i < L.size() -1; ++i)
		{
			L[i] = power(gen, i +1);
		}
		L[L.size() -1] = 0;
	}

	// Calculate Parity check matrix H
	{
		NTL::mat_GF2E X, Y, Y2, Z, Z2;
		X.SetDims(errors, errors);
		Y.SetDims(errors, non_root_points);
		Y2.SetDims(2 *errors, non_root_points);
		Z.SetDims(non_root_points, non_root_points);
		Z2.SetDims(non_root_points, non_root_points);

		// Teoplitz matrix
		for (int row = 0; row < X.NumRows(); ++row)
			for (int col = 0; col < X.NumCols(); ++col)
				if (row >= col)
					X.put(row, col, coeff(g, errors -(row -col)));

		// Vandermonde matrix
		for (int row = 0; row < Y.NumRows(); ++row)
			for (int col = 0; col < Y.NumCols(); ++col)
				Y.put(row, col, NTL::power(L[col], row));

		for (int row = 0; row < Y2.NumRows(); ++row)
			for (int col = 0; col < Y2.NumCols(); ++col)
				Y2.put(row, col, NTL::power(L[col], row));

		mul(XY, X, Y);

		// Diagonal matrix
		for (int i = 0; i < Z.NumRows(); ++i)
			Z.put(i, i, inv(call(g, L[i])));
		for (int i = 0; i < Z2.NumRows(); ++i)
		{
			GF2E e = call(g, L[i]);
			Z2.put(i, i, inv(e *e));
		}

		mul(YZ, Y, Z);
		mul(XYZ, XY, Z);
		mul(YZ2, Y2, Z2);
	}

	// Calculate generator
	{
		XY_bin = trace_construct(XY);
		XYZ_bin = trace_construct(XYZ);
		YZ_bin = trace_construct(YZ);
		Y2Z2_bin = trace_construct(YZ2);

		kernel(G_XY, transpose(XY_bin));
		kernel(G_XYZ, transpose(XYZ_bin));
	}
}

///
/// \brief bgc_t::calculate_fi
/// \param i
/// \return
///
/// Needed to compute fi = 1 / (x - Li) with:
/// 1 / (x - Li) = -(1 / g(Li)) * ((g(x) - g(Li)) / (x - Li))
///
GF2EX bgc_t::calculate_fi(long i) const
{
	// fi = -(1 / g(ai))
	NTL::GF2EX fi = NTL::GF2EX::zero() -NTL::inv(call(g, L[i]));

	NTL::GF2EX zaehler = g -call(g, L[i]), nenner, pro;
	SetCoeff(nenner, 1);
	nenner -= L[i];

	div(pro, zaehler, nenner);
	// fi *= ((g(x) - g(ai)) / (x - ai))
	MulMod(fi, fi, pro, g);

	return fi;
}

GF2EX bgc_t::calculate_sc(const vec_GF2& c) const
{
	//if (size_t(c.length()) != L.size())
		//throw std::invalid_argument("c.length() != L.length(), c.length() was " +std::to_string(c.length()) +" and L.length() was " +std::to_string(L.size()));
	if (c.length() % deg(*f))
		throw std::invalid_argument("(c.length() % L.length()) was not NULL");
	NTL::GF2EX sc_1 = NTL::GF2EX::zero(),
			   sc_2 = NTL::GF2EX::zero(),
			   sc_3 = NTL::GF2EX::zero();

	for (long i = 0; i < c.length(); ++i)
	{
		if (! IsZero(c[i]))
			add(sc_1, sc_1, calculate_fi(i));
	}

	for (long i = 0; i < c.length(); ++i)
	{
		if (! IsZero(c[i]))
		{
			NTL::GF2EX pol;
			SetCoeff(pol, 1);
			sub(pol, pol, L[i]);
			InvMod(pol, pol, g);

			add(sc_2, sc_2, pol);
		}
	}

	long b = c.length() /deg(*f);
	for (long i = 0; i < b; ++i)
	{
		GF2E e = GF2E::zero();

		for (long j = 0; j < deg(*f); ++j)
			if (! IsZero(c[i *deg(*f) +j]))
				SetCoeff(e.LoopHole(), j);

		SetCoeff(sc_3, b -i -1, e);
	}

	if (((sc_1 != sc_2) || (sc_2 != sc_3)), 0)
	{
		std::cerr << "sc_1\n" << sc_1
				  << "\nsc_2\n" << sc_2
				  << "\nsc_3\n" << sc_3
				  << "\nc\n" << print(c) << std::endl;
		throw std::runtime_error("Assert that the two alternate forms of calculating sc produce the same result. >inb4 they dont");
	}

	return sc_3;
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
	assert(size_t(e.length()) == L.size());
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
	assert(size_t(e.length()) == L.size());
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
		auto found = std::find(L.cbegin(), L.cend(), l[i]);

		if (found == L.cend())
			throw std::runtime_error("Could not find element " +print(l[i]) +" in set L");
		else
			ret[i] = (found -L.cbegin());
	}

	return ret;
}

long bgc_t::k() const
{
	return non_root_points -errors *dimension;
}

vec_GF2 bgc_t::encode(const vec_GF2& msg) const
{
	if (msg.length() != k())
		throw std::invalid_argument("msg.length() != k()");

	// As easy as it gets
	return msg *G_XYZ;
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
	ret.SetLength(L.size());

	for (size_t i = 0; i < L_indices.size(); ++i)
		ret[L_indices[i]] = GF2(1);

	return ret;
}

bool bgc_t::is_codeword(const vec_GF2& c) const
{
	return (c.length() == (uint64_t(1) << dimension)) && IsZero(c *transpose(XYZ_bin));
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
	   << "\nH (" << XYZ_bin.NumRows() << 'x' << XYZ_bin.NumCols() << ")\n";

	return ss.str();
}
