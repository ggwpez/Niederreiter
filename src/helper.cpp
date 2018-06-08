#include "helper.hpp"
#include "printer.hpp"

#include <cassert>
#include <algorithm>

using namespace NTL;

void monice(GF2EX poly, GF2EX& ret)
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

GF2EX monice(GF2EX const& poly)
{
	GF2EX ret;
	monice(poly, ret);
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
		result = (result *x) +coeff(p, i);

	return result;
}

vec_GF2E to_ext_field_poly(vec_GF2 const& vec, GF2X const& field)
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
