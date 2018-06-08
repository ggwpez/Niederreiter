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
	std::cout << "BGC generated" << std::endl;
	NCS::KeyPair keys = NCS::keygen(bgc);
	std::cout << "Key generated" << std::endl;

	NTL::vec_GF2 msg, cipher, recovered;
	msg.SetLength(n);
	for (int i = 0; i < t; ++i)
		msg[i +(1 << (m -1))] = GF2(1);

	NCS::encode(msg, keys.pub, cipher);
	NCS::decode(cipher, keys.sec, recovered);

	std::cout << "HASH_OF(msg)       " << HASH_OF(msg) << std::endl
			  // Should be 7577682492700538189 if no parameters were changed
			  << "HASH_OF(cypher)    " << HASH_OF(cipher) << std::endl
			  << "HASH_OF(recovered) " << HASH_OF(recovered) << std::endl;

	if (msg != recovered)
		throw std::runtime_error("Decoding error");
}
