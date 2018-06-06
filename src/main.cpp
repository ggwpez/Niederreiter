#include "bgc.hpp"
#include "tests.hpp"
#include "helper.hpp"

#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/mat_GF2.h>

#include <csignal>
#include <iostream>
#include <ctime>
#include <iomanip>
using namespace NTL;

mat_GF2 canonical_check_matrix(bgc_t const& bgc)
{
	long n = bgc.non_root_points, t = bgc.errors;
	mat_GF2E YZ, H;
	YZ.SetDims(t, n);
	H.SetDims(t, n);

	for (long j = 0; j < n; ++j)
		YZ[0][j] = inv(call(bgc.g, bgc.L[j]));

	for (long i = 1; i < t; ++i)
		for (long j = 0; j < n; ++j)
			YZ[i][j] = YZ[i -1][j] *bgc.L[j];

	for (int i = 0; i < t; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k <= i; ++k)
				H[i][j] += YZ[k][j] *coeff(bgc.g, t +k -i);

	return trace_construct(H);
}

struct nr_pubkey_t
{
	mat_GF2 h;
};

struct nr_seckey_t
{
	mat_GF2 Si, p;
};

bool test_petterson(bgc_t const& bgc);

void niederreiter_keygen(bgc_t const& bgc, nr_pubkey_t& pub_key, nr_seckey_t& sec_key)
{
	mat_GF2 sInv, m, p;

	mat_GF2 H = canonical_check_matrix(bgc);
	compute_systematic_form(H, sInv, m, p);
	std::cout << "H\n" << print(H) << "\nCompact_H\n" << print(m) << std::endl;

	pub_key = nr_pubkey_t{ m };
	sec_key = nr_seckey_t{ sInv, p };
}

void niederreiter_keygen_2(bgc_t const& bgc, NTL::mat_GF2& pub_key, NTL::mat_GF2& H, NTL::mat_GF2& Si, NTL::mat_GF2& Pi)
{
	H = bgc.XYZ_bin;
	long nk = bgc.non_root_points -bgc.k(),
		 n = bgc.non_root_points;
	NTL::mat_GF2 S, P, Key;

	do
	{
		S = NTL::random_mat_GF2(nk, nk);
	} while (NTL::IsZero(NTL::determinant(S)));

	P = create_rand_permutation(n);
	std::printf("S %ldx%ld\nH %ldx%ld\nP %ldx%ld\n", S.NumRows(), S.NumCols(), H.NumRows(), H.NumCols(), P.NumRows(), P.NumCols());
	pub_key = S *(H *P);
	assert(pub_key.NumRows() == nk && pub_key.NumCols() == bgc.non_root_points);

	NTL::inv(Si, S);
	NTL::inv(Pi, P);
}

NTL::vec_GF2 niederreiter_encode(NTL::vec_GF2 const& msg, nr_pubkey_t const& pub)
{
	if (msg.length() != pub.h.NumRows() +pub.h.NumCols())
		throw "lel";

	mat_GF2 H = mat_merge_ID_left(pub.h);

	return H *msg;
}

NTL::vec_GF2 niederreiter_decode(bgc_t const& bgc, NTL::vec_GF2 const& c, nr_seckey_t const& sec)
{
	vec_GF2 syndrome = sec.Si *c;
	vec_GF2 pe = bgc.syndrom_decode_2(syndrome);
	vec_GF2 e = pe *sec.p;	// inv(p)
	return e;
}

int main(int, char**)
{
	assert(! test());
	// Length
	int const m = 4;		// 10
	// Diemnsion
	int const n = 1 << m;
	int const t = 2;		// 9

	bgc_t bgc(m, n, t);
	std::cout << bgc.to_str() << std::endl;

	/*if (! test_petterson(bgc))
	{
		std::cerr << "Petterson test failed, abort" << '\n';
		abort();
	}*/

	nr_pubkey_t pub;
	nr_seckey_t sec;
	niederreiter_keygen(bgc, pub, sec);

	NTL::vec_GF2 msg;
	msg.SetLength(n);
	msg[5] = GF2(1) -msg[5];
	msg[4] = GF2(1) -msg[4];

	auto c = niederreiter_encode(msg, pub);
	auto recovered = niederreiter_decode(bgc, c, sec);

	std::cout << "msg\n" << print(msg)
			  << "\nmsg'\n" << print(recovered) << '\n';
}

bool test_petterson(bgc_t const& bgc)
{
	NTL::vec_GF2 message;
	message.SetLength(bgc.k());

	for (int i = 0; i < message.length(); ++i)
		message[i] = NTL::random_GF2();

	NTL::vec_GF2 resistant = bgc.encode(message);
	NTL::vec_GF2 modded = resistant;

	for (int i = 0; i < bgc.errors; ++i)
		modded[i*2] = NTL::GF2(1) -resistant[i*2];	// TODO random

	auto e = bgc.syndrom_decode(modded);
	auto recovered = modded +e;

	auto worked = (recovered == resistant);
	if (! worked)
		std::cerr << "Message:   " << print(message) << std::endl
				  << "Resistant:  " << print(resistant) << std::endl
				  << "Modded:    " << print(modded) << std::endl
				  << "Recovered: " << print(recovered) << std::endl
				  << "(Recovered != Resistant)" << std::endl;

	return worked;
}
