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

struct nr_pubkey_t
{
	mat_GF2 h;
};

struct nr_seckey_t
{
	mat_GF2 Si, p;
};

void niederreiter_keygen(bgc_t const& bgc, nr_pubkey_t& pub_key, nr_seckey_t& sec_key)
{
	mat_GF2 sInv, m, p;

	mat_GF2 H = bgc.H;//canonical_check_matrix(bgc);
	compute_systematic_form(H, sInv, m, p);
	//std::cout << "H\n" << print(H) << "\nCompact_H\n" << print(m) << std::endl;

	pub_key = nr_pubkey_t{ m };
	sec_key = nr_seckey_t{ sInv, p };
}

NTL::vec_GF2 niederreiter_encode(NTL::vec_GF2 const& msg, nr_pubkey_t const& pub)
{
	/*if (msg.length() != pub.h.NumRows() +pub.h.NumCols())
		throw "lel";*/

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
	SetSeed(ZZ(123456789));
	assert(! test());
	// Length
	int const m = 13;		// 10
	// Diemnsion
	int const n = 1 << m;
	int const t = 117;		// 9

	bgc_t bgc(m, n, t);
	std::cout << bgc.to_str() << std::endl;

	nr_pubkey_t pub;
	nr_seckey_t sec;
	niederreiter_keygen(bgc, pub, sec);

	//for (int j = 0; j < 10; ++j)
	{
		NTL::vec_GF2 msg;
		msg.SetLength(n);
		for (int i = 0; i < t; ++i)
			msg[i +(1 << (m -1))] = GF2(1);

		auto cypher = niederreiter_encode(msg, pub);
		auto recovered = niederreiter_decode(bgc, cypher, sec);

		std::cout << "msg\n" << print(msg)
				  << "\ncypher\n" << print(cypher)
				  << "\nmsg'\n" << print(recovered) << '\n';

		if (msg == recovered)
			std::cout << "OK\n";
		else
			throw std::runtime_error("Decoding error");
	}
}
