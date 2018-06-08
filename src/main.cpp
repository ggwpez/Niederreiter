#include "bgc.hpp"
#include "ncs.hpp"
#include "hasher.hpp"

#include <csignal>
#include <iostream>
#include <ctime>
#include <iomanip>
using namespace NTL;

int main(int, char**)
{
	SetSeed(ZZ(123456789));

	int const m = 13,
			  n = 1 << m,
			  t = 117;

	BGC bgc = BGC::create(m,n,t);
	std::cout << bgc.to_str() << std::endl;
	NCS::KeyPair keys = NCS::keygen(bgc);

	NTL::vec_GF2 msg, cipher, recovered;
	msg.SetLength(n);
	for (int i = 0; i < t; ++i)
		msg[i +(1 << (m -1))] = GF2(1);

	NCS::encode(msg, keys.pub, cipher);
	NCS::decode(cipher, keys.sec, recovered);

	std::cout << "HASH_OF(msg)\n" << HASH_OF(msg)
			  << "\nHASH_OF(cypher)\n" << HASH_OF(cipher)
			  << "\nHASH_OF(msg')\n" << HASH_OF(recovered) << '\n';

	if (msg != recovered)
		throw std::runtime_error("Decoding error");
}
