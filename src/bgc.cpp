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

BGC::BGC(long m, long n, long t, GF2X const& f, GF2EX const& g, GF2E const& gen, support_t const& L, mat_GF2 const& H)
	: m(m), n(n), t(t), f(f), g(g), gen(gen), L(L), H(H)
{

}

BGC BGC::create(long m, long n, long t)
{
	NTL::GF2X f;
	NTL::GF2EX g;
	calculate_f(m, f);
	GF2E::init(f);
	calculate_g(t, g);

	return BGC::create_with_field_and_gp_poly(m, n, t, f, g);
}

BGC BGC::create_with_field(long m, long n, long t, GF2X const& f)
{
	GF2E::init(f);
	NTL::GF2EX g;
	calculate_g(t, g);

	return BGC::create_with_field_and_gp_poly(m, n, t, f, g);
}

BGC BGC::create_with_field_and_gp_poly(long m, long n, long t, GF2X const& f, GF2EX const& g)
{
	check_args(m, n, t, f, g);
	BGC bgc = BGC(m, n, t, f, g, GF2E::zero(), support_t(), mat_GF2());

	calculate_gen(f, bgc.gen);
	calculate_L(n, bgc.gen, bgc.L);
	calculate_H(t, bgc.L, bgc.g, bgc.H);

	return bgc;
}

void BGC::calculate_f(long m, GF2X& f)
{
	GF2X tmp = BuildIrred_GF2X(m);
	BuildRandomIrred(f, tmp);
}

void BGC::calculate_g(long t, GF2EX& g)
{
	GF2EX tmp = BuildIrred_GF2EX(t);
#if DEBUG
	long tries = 0;
#endif
	do
	{
		BuildRandomIrred(g, tmp);
		++tries;
	} while (zero_coefficients(g) > 1);
#if DEBUG
	if (tries > 1)
		std::clog << "Tries for g: " << tries << std::endl;
#endif
}

void BGC::calculate_gen(GF2X f, GF2E& gen)
{
	SetCoeff(f, deg(f), GF2::zero());
	f.normalize();
	gen = conv<GF2E>(f);
}

void BGC::calculate_H(long t, support_t L, GF2EX const& g, mat_GF2& H)
{
	long n = L.length();
	mat_GF2E YZ, h;
	YZ.SetDims(t, n);
	h.SetDims(t, n);

	for (long col = 0; col < n; ++col)
		YZ[0][col] = inv(call(g, L[col]));

	for (long row = 1; row < t; ++row)
		for (long col = 0; col < n; ++col)
			// TODO, if n==2^m the multiplication can be converted to an int addition and an array access
			YZ[row][col] = YZ[row -1][col] *L[col];

	for (int i = 0; i < t; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k <= i; ++k)
				h[i][j] += YZ[k][j] *coeff(g, t +k -i);

	trace_construct(h, H);
}

void BGC::calculate_L(long n, GF2E const& gen, support_t& L)
{
	L.SetLength(n);
	L[0] = gen;

	for (long i = 1; i < L.length() -1; ++i)
		L[i] = L[i -1] *gen;

	L[L.length() -1] = 0;
}

void BGC::check_args(long m, long n, long t, GF2X const& f, GF2EX const& g)
{
	if (n <= t)
		throw std::invalid_argument("t < n");
	if (! ((2 <= t) && (t <= ((1l << m) -1) /m)))
		throw std::invalid_argument("t");
	if (! ((m *t +1 <= n) && (n <= (1l << m))))
		throw std::invalid_argument("n");
	if (deg(f) != m)
		throw std::invalid_argument("f");
	if (deg(g) != t)
		throw std::invalid_argument("g");
}

void BGC::calculate_error(GF2EX const& poly, vec_GF2& e) const
{
	e.SetLength(n);

	for (long i = 0; i < e.length(); ++i)
	{
		GF2E y = call(poly, L[i]);
		if (IsZero(y))
			e[i] = GF2(1);
	}
}

void BGC::calculate_sc(vec_GF2 const& c, GF2EX& sc) const
{
	if (c.length() % deg(f))
		throw std::invalid_argument("(c.length() % L.length()) was not NULL");
	sc = NTL::GF2EX::zero();

	long b = c.length() /deg(f);
	for (long i = 0; i < b; ++i)
	{
		GF2E e = GF2E::zero();

		for (long j = 0; j < deg(f); ++j)
			if (! IsZero(c[i *deg(f) +j]))
				SetCoeff(e.LoopHole(), j);

		SetCoeff(sc, b -i -1, e);
	}
}

void BGC::calculate_vc(GF2EX const& sc, GF2EX& vc) const
{
	// Return the identity function of there are no errors
	if (IsZero(sc))
		throw std::runtime_error("sc not invertible");	 // FIXME
	else
	{
		GF2EX v, id = conv<GF2EX>("[[][1]]");

		// v = 1 / sc
		InvMod(v, sc, g);
		if (v == sc)
			vc = sc;
		else
		{
			// v = 1 / sc -x
			add(v, v, id);		// TODO try an add
			//std::cout << "Trying to square-root: " << print(v) << '\n';
			sqr_root(v, vc);
		}
	}
}

long BGC::k() const
{
	return n -t *m;
}

void BGC::syndrom_decode(vec_GF2 const& c, vec_GF2& e) const
{
	GF2EX sc, vc, elp;

	calculate_sc(c, sc);
	calculate_vc(sc, vc);
	calculate_sigma(vc, g, g, elp);
	monice(elp, elp);				// Optional

	calculate_error(elp, e);
}

// (x)^(1/p) = x^(p^m-1)
void BGC::sqr_root(GF2EX const& p, GF2EX& res) const
{
	ZZ e = (ZZ(1) << (m *t -1));

	res = PowerMod(p, e, g);
#ifdef DEBUG
	if (SqrMod(res, g) != p)
		throw std::runtime_error("Could not find root");
#endif
}

// L[a] * L[b] = gen^(a+1) * gen^(b+1) = gen^(a+b+2) = L[a+b+1]
size_t BGC::mul_L_elements(size_t a, size_t b) const
{
	return (a+b+1) % L.length();
}

// L[a] ^ exp = gen^(a+1) ^ exp = gen^(exp*a+exp) = L[exp*a+exp-1]
size_t BGC::power_L_elements(size_t a, long exp) const
{
	return (exp *a +exp -1) % L.length();
}

ZZ BGC::order() const
{
	return (ZZ(1) << (m *t));
}

std::string BGC::to_str() const
{
	std::ostringstream ss;

	ss << "[n,k,d]-Code = [" << (uint64_t(1) << m) << ',' << k() << ',' << (2*t +1) << "]" << std::endl
	   << "F(2^" << m << ") = " << "F(2)[x]/" << print(f) << std::endl
	   << "Goppa Polynomial: " << print(g) << std::endl
	   << "Generator: " << print(gen) << std::endl
	   << "Support={ " << print(L) << '}'
	   << "\nH (" << H.NumRows() << 'x' << H.NumCols() << ")\n";

	return ss.str();
}
