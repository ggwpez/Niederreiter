#include "tests.hpp"
#include "helper.hpp"
#include "printer.hpp"
#include "bgc.hpp"
#include "binom.hpp"

#include <algorithm>
#include <numeric>

#include <NTL/mat_GF2.h>

using namespace NTL;

int test_binom_val(long const n, long const t, ZZ x)
{
	vec_GF2 enc;
	char* dec,
		* msg = new char[n];

	std::memset(msg, 0, n);
	BytesFromZZ(reinterpret_cast<unsigned char*>(msg), x, n);

	Binom::encode(n, t, x, enc);
	Binom::decode(n, t, enc, dec);

	if(std::strncmp(msg, dec, n))
		throw std::runtime_error("cw decoding failed");

	return 0;
}

int test_binom()
{
	long const n = 8192,
			   t = 117;
	long k = log2_coeff(n, t);
	ZZ biggest = Binom::coeff(ZZ(n), ZZ(t)) -1;

	std::cout << "n " << n << " t " << t << " k " << k << std::endl;
	for (long i = 0; i < 1000; ++i)
	{
		ZZ x = RandomBits_ZZ(k);
		test_binom_val(n, t, x);
	}
	test_binom_val(n, t, ZZ(0));
	test_binom_val(n, t, ZZ(1));

	test_binom_val(n, t, biggest);
	try
	{
		test_binom_val(n, t, biggest +1);	// Should fail
	}
	catch (std::invalid_argument const& e)
	{
		if (std::strcmp(e.what(), "Value too large"))	// TODO Bad style
			throw e;
	}

	return 0;
}

int test_mat()
{
	long m = 8, n = 15;
	mat_GF2 mat = random_mat_GF2(m,n);
	mat_GF2 mat_l = getLeftSubMatrix(mat),
			mat_r = getRightSubMatrix(mat),
			merged = mat_merge_colls(mat_l, mat_r),
			id_l = mat_merge_ID_left(mat_r),
			id_r = mat_merge_ID_right(mat_l);

	assert(ident_mat_GF2(id_l.NumRows()) == getLeftSubMatrix(id_l));
	assert(ident_mat_GF2(id_r.NumRows()) == getRightSubMatrix(id_r));
	assert(merged == mat);

	mat_GF2 h = random_mat_GF2(1 << 4, 1 << 5);
	mat_GF2 H = mat_merge_ID_left(h);	// TODO
	vec_GF2 v = random_vec_GF2(h.NumRows() +h.NumCols());

	vec_GF2 r1 = H *v,
			r2;
	mat_mul_right_compact(h, v, r2);
	assert(r1 == r2);

	return 0;
}

int test()
{
	return test_mat() || test_binom();
}
