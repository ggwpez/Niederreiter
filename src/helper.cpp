#include "binom.hpp"
#include "helper.hpp"
#include "printer.hpp"

#include <NTL/RR.h>
#include <NTL/vector.h>
#include <NTL/vec_GF2.h>

#include <cassert>
#include <algorithm>

using namespace NTL;

void monic(GF2EX poly, GF2EX& ret)
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
			GF2E ic = inv(LeadCoeff(poly));
			mul(ret, poly, ic);
			//div(ret, poly, LeadCoeff(poly));
		}
	}
}

GF2EX monic(GF2EX const& poly)
{
	GF2EX ret;
	monic(poly, ret);
	return ret;
}

mat_GF2 create_rand_permutation(size_t s)
{
	mat_GF2 ret;

	ident(ret, s);
	std::random_shuffle(ret._mat__rep.begin(), ret._mat__rep.end());

	return ret;
}

// Horner schema
GF2E call(GF2EX const& p, GF2E const& x)
{
	GF2E result = LeadCoeff(p);

	for (long i = deg(p); i --> 0;)
	{
		//result = (result *x) +coeff(p, i);
		mul(result, result, x);
		add(result, result, coeff(p, i));
	}

	return result;
}

mat_GF2 getLeftSubMatrix(mat_GF2 const& mat)
{
	if (mat.NumCols() <= mat.NumRows())
		throw std::invalid_argument("Empty Submatrix");

	long n = mat.NumRows();
	mat_GF2 res(INIT_SIZE, n, n);

	for (long i = 0; i < n; ++i)
		for (long j = 0; j < n; ++j)
			res.put(i,j, mat[i][j]);

	return res;
}

mat_GF2 getRightSubMatrix(mat_GF2 const& mat)
{
	if (mat.NumCols() <= mat.NumRows())
		throw std::invalid_argument("Empty Submatrix");

	long n = mat.NumRows(),
		 m = mat.NumCols() -n;
	mat_GF2 res(INIT_SIZE, mat.NumRows(), m);

	for (long i = 0; i < n; ++i)
		for (long j = 0; j < m; ++j)
			res.put(i,j, mat[i][j +n]);

	return res;
}

template<typename T>
inline bool is_power_of_two(T const& x)
{
	return ! (x & (x -1));
}

GF2 fast_dot_product(vec_GF2 const& a, vec_GF2 const& b)
{
	GF2 ret = GF2::zero();

	for (long i = 0; i < a.rep.length(); ++i)
	{
		ret += __builtin_parityl(a.rep[i] & b.rep[i]);
	}

	return ret;
}

void mat_mul_right_compact(mat_GF2 const& mat, vec_GF2 const& vec, vec_GF2& out)
{
	assert(sizeof(_ntl_ulong) == 8);
	if (vec.length() & (sizeof(_ntl_ulong) -1))
		throw std::invalid_argument("vec.length() must be a multiple of 8");
	if (vec.length() != mat.NumRows() +mat.NumCols())
		throw std::invalid_argument("Bad Vector length");

	out.SetLength(mat.NumRows());

	vec_GF2 vec2(INIT_SIZE, mat.NumCols());
	for (int i = 0; i < vec2.length(); ++i)
		vec2.put(i, vec[i +mat.NumRows()]);

	for (long i	= 0; i < mat.NumRows(); ++i)
	{
		out[i] = fast_dot_product(mat[i], vec2);

		if (IsOne(vec[i]))
			out[i] = GF2(1) -out[i];
	}
}

mat_GF2 mat_merge_colls(mat_GF2 const& a, mat_GF2 const& b)
{
	if (a.NumRows() != b.NumRows())
		throw std::invalid_argument("a and b need the same number of Rows");

	mat_GF2 res(INIT_SIZE, a.NumRows(), a.NumCols() +b.NumCols());

	for (long i = 0; i < res.NumRows(); ++i)
	{
		vec_GF2& row = res[i];

		for (long j = 0; j < a.NumCols(); ++j)
			row[j] = a[i][j];
		for (long j = 0; j < b.NumCols(); ++j)
			row[j +a.NumCols()] = b[i][j];
	}

	return res;
}

mat_GF2 mat_merge_ID_left(mat_GF2 const& b)
{
	mat_GF2 ID = ident_mat_GF2(b.NumRows());
	return mat_merge_colls(ID, b);
}

mat_GF2 mat_merge_ID_right(mat_GF2 const& a)
{
	mat_GF2 ID = ident_mat_GF2(a.NumRows());
	return mat_merge_colls(a, ID);
}

void trace_construct(NTL::Mat<GF2E> const& mat, NTL::Mat<GF2>& H)
{
	long m = NTL::GF2E::modulus().n;
	H.SetDims(m *mat.NumRows(), mat.NumCols());

	for (long row = 0; row < mat.NumRows(); ++row)
		for (long col = 0; col < mat.NumCols(); ++col)
		{
			NTL::GF2X const& e = NTL::conv<NTL::GF2X>(mat[row][col]);

			for (long i = 0; i < m; ++i)	// i < deg(e) TODO
				H.put(m *row +i, col, NTL::coeff(e, i));
		}
}

void calculate_sigma(const GF2EX& a, const GF2EX& b, const GF2EX& g, GF2EX& sigma)
{
	NTL::GF2EX F = a, G = b, B = NTL::conv<NTL::GF2EX>("[[1]]"), C = NTL::GF2EX::zero();
	long t = NTL::deg(b);

	while (NTL::deg(G) > (t /2))
	{
		NTL::swap(F, G); NTL::swap(B, C);

		while (NTL::deg(F) >= NTL::deg(G))
		{
			long j = NTL::deg(F) -NTL::deg(G);
			auto h = NTL::LeadCoeff(F) /NTL::LeadCoeff(G);

			F = F -h *NTL::PowerXMod(j, g) *G;
			B = B -h *NTL::PowerXMod(j, g) *C;
		}
	}

	GF2EX G_sq;
	mul(G_sq, G, G);

	mul(sigma, C, C);
	mul(sigma, sigma, NTL::PowerXMod(1, g));
	add(sigma, sigma, G_sq);
}

long count_coefficients(GF2EX& p, const GF2E& e)
{
	long sum = 0;

	for (long i = 0; i < deg(p); ++i)
		if (coeff(p, i) == e)
		++sum;

	return sum;
}

long log2_coeff(const long n, const long t)
{
	RRPush push_guard;
	RR::SetPrecision(4);
	RR x = to_RR(Binom::coeff(ZZ(n), ZZ(t)));

	x = NTL::log(x) /NTL::log(RR(2));

	ZZ z = FloorToZZ(x);
	if (x >= std::numeric_limits<long>::max())
		throw std::runtime_error("Overlow");

	return conv<long>(z);
}
