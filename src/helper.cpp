#include "helper.hpp"
#include "printer.hpp"

#include <cassert>
#include <algorithm>

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
			GF2E ic = inv(LeadCoeff(poly));
			mul(ret, poly, ic);
			//div(ret, poly, LeadCoeff(poly));
		}
	}
}

GF2EX monice(const GF2EX& poly)
{
	GF2EX ret;
	monice(ret, poly);
	return ret;
}

std::vector<GF2E> find_roots(const GF2EX& poly)	// TODO only support L needs to be tested
{
	std::vector<GF2E> roots;
	GF2E x = GF2E::zero();
	assert(GF2E::modulus().n < 64);

	for (long i = 0; i < (1ull << GF2E::modulus().n); ++i)
	{
		x = generate_GF2E(i);
		GF2E y = call(poly, x);
		if (IsZero(y))
			roots.push_back(x);

		//std::cout << "poly(" << print(x) << ") = " << print(y) << std::endl;
	}

	return roots;
}

GF2E call_slow(GF2EX const& p, const GF2E& x)
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

mat_GF2 create_rand_permutation(size_t s)
{
	mat_GF2 ret;

	ident(ret, s);
	std::random_shuffle(ret._mat__rep.begin(), ret._mat__rep.end());

	return ret;
}

GF2E call(const GF2EX& p, const GF2E& x)
{
	GF2E result = LeadCoeff(p);

	for (long i = deg(p); i --> 0;)
		result = (result *x) +coeff(p, i);

	return result;
}

vec_GF2E to_ext_field_poly(const vec_GF2& vec, GF2X const& field)
{
	int m = deg(field);
	if (vec.length() % m != 0)
		throw std::runtime_error("Conversion impossible");

	long t = vec.length() / m;
	vec_GF2E result;
	result.SetLength(t);
	long count = 0;
	for (long i = t - 1; i >= 0; --i)
	{
		for (long j = deg(field) - 1; j >= 0; --j)
		{
			long q = count >> 5;
			long r = count & 0x1f;

			long e = (conv<long>(vec[q]) >> r) & 1;
			if (e == 1)
			{
				result[i] += generate_GF2E(long(i) << j);
			}
			count++;
		}
	}

	return result;
}

void compute_systematic_form(const mat_GF2& H, mat_GF2& sInv, mat_GF2& m, mat_GF2& p)
{
	long n = H.NumCols();
	mat_GF2 hp;

	do
	{
		p = create_rand_permutation(n);
		mul(hp, H, p);
		sInv = getLeftSubMatrix(hp);
	} while (IsZero(determinant(sInv)));

	mat_GF2 shp = inv(sInv) *hp;
	m = getRightSubMatrix(shp);
}

mat_GF2 getLeftSubMatrix(const mat_GF2& mat)
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

mat_GF2 getRightSubMatrix(const mat_GF2& mat)
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

mat_GF2 mat_merge_colls(const mat_GF2& a, const mat_GF2& b)
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

mat_GF2 mat_merge_ID_left(const mat_GF2& b)
{
	mat_GF2 ID = ident_mat_GF2(b.NumRows());
	return mat_merge_colls(ID, b);
}

mat_GF2 mat_merge_ID_right(const mat_GF2& a)
{
	mat_GF2 ID = ident_mat_GF2(a.NumRows());
	return mat_merge_colls(a, ID);
}

NTL::Mat<GF2> trace_construct(const NTL::Mat<GF2E>& mat)
{
	NTL::Mat<NTL::GF2> ret;
	long m = NTL::GF2E::modulus().n;
	ret.SetDims(m *mat.NumRows(), mat.NumCols());

	for (long row = 0; row < mat.NumRows(); ++row)
		for (long col = 0; col < mat.NumCols(); ++col)
		{
			NTL::GF2X const& e = NTL::conv<NTL::GF2X>(mat[row][col]);

			for (long i = 0; i < m; ++i)	// i < deg(e) TODO
				ret.put(m *row +i, col, NTL::coeff(e, i));
		}

	return ret;
}
