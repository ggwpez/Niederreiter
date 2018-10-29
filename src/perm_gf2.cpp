#include "perm_gf2.hpp"

#include <cassert>
#include <algorithm>
#include <numeric>
#include <limits>

perm_GF2 create_rand_permutation(size_t s)
{
	if (s > std::numeric_limits<perm_GF2::value_type>::max())
		throw std::invalid_argument("perm_GF2's value type is to small.\nGo to perm_GF2.hpp and change uint16_t to something bigger");

	perm_GF2 ret;
	ret.FixLength(long(s));
	// std is ez lf
	std::iota(ret.begin(), ret.end(), 0);
	std::random_shuffle(ret.begin(), ret.end());

	return ret;
}

void mul(NTL::mat_GF2& out, const perm_GF2& perm, const NTL::mat_GF2& mat)
{
	assert(mat.NumRows() == perm.length());
	for (int i = 0; i < mat.NumRows(); ++i)
		out[i] = mat[perm[i]];
}

void mul(NTL::mat_GF2& out, const NTL::mat_GF2& mat, const perm_GF2& perm)
{
	out.SetDims(mat.NumRows(), mat.NumCols());
	std::cerr << mat.NumRows() << "x" << mat.NumCols() << '\n'
			  << perm.length() << '\n';
	assert(mat.NumCols() == perm.length());
	for (int i = 0; i < mat.NumRows(); ++i)
	{
		for (int j = 0; j < mat.NumCols(); ++j)
		{
			//std::cerr << perm[j] << '\n';
			out[i][j] = mat[i][perm[j]];
		}
	}
}

void mul(NTL::vec_GF2& out, const NTL::vec_GF2& vec, const perm_GF2& perm)
{
	out.FixLength(vec.length());
	for (int i = 0; i < vec.length(); ++i)
		out[i] = vec[perm[i]];
}
