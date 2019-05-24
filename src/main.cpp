#include "arg_parse.hpp"
#include "bgc.hpp"
#include "binom.hpp"
#include "hasher.hpp"
#include "helper.hpp"
#include "ncs.hpp"
#include "rand_helper.hpp"
#include "serializer.hpp"
#include "tests.hpp"

#include <fstream>
#include <iostream>
#include <istream>

#include <NTL/vec_GF2.h>

#define IMPL(x) x##_classic

static state_t state;
using namespace NTL;

std::ifstream open_key();
int			  keygen();
void		  encode(NCS::PubKey const& key);
void		  decode(NCS::SecKey const& key);

int main(int argc, char** argv)
{
	SetSeed(ZZ(std::time(nullptr))); // FIXME use secure random stream
	state = parse_args(argc, argv);

	if (state.mode != mENC && state.mode != mDEC && state.mode != mCPK && state.mode != mPINF && state.mode != mSINF)
	{
		if (state.mode == mGEN) // KEYGEN
			return keygen();

		std::cerr << "No operation mode specified, use --help for help" << std::endl;
		return 1;
	}

	std::ifstream fs = open_key();
	NCS::KeyPair  keys;

	if (state.mode == mENC) // ENC
	{
		keys.m_pk.deserialize(fs);

		encode(keys.m_pk);
	}
	else if (state.mode == mDEC) // DEC
	{
		keys.m_sk.deserialize(fs);

		decode(keys.m_sk);
	}
	else if (state.mode == mSIG)
	{
		keys.m_sk.deserialize(fs);

		// sign(keys.m_sk);
	}
	else if (state.mode == mSFY)
	{
		keys.m_pk.deserialize(fs);

		//	sig_verify(keys.m_sk);
	}
	else if (state.mode == mSINF)
	{
		keys.m_sk.deserialize(fs);

		std::cout << keys.m_sk.bgc.to_str() << std::endl;
	}
	else if (state.mode == mPINF)
	{
		keys.m_pk.deserialize(fs);

		std::cout << "PK-info not yet implemented" << std::endl;
		std::cout << (keys.m_pk.bits / 8) << " byte of userdata per message" << std::endl;
	}
	else if (state.mode == mCPK) // CPK
	{
		keys.m_sk.deserialize(fs);
		keys.IMPL(reconstruct_pk)(); // spooky

		keys.m_pk.serialize(std::cout);
	}
	else
		throw std::runtime_error("Unreachable");

	std::cout.flush();
	return 0;
}

int keygen()
{
	if (state.m < 0 || state.n < 0 || state.t < 0)
	{
		std::cerr << "Error, -m, -n and -t must be specified in --gen use --help for help" << std::endl;
		return 1;
	}

	BGC			 bgc  = BGC::create(state.m, state.n, state.t);
	NCS::KeyPair keys = NCS::IMPL(keygen)(bgc);

	keys.m_sk.serialize(std::cout);
	return 0;
}

void encode(NCS::PubKey const& key)
{
	NTL::vec_GF2 err, enc_err;

	uint32_t ml  = (key.bits / 8);
	char*	msg = new char[ml + 1]();
	state.is->read(msg, ml);

	Binom::encode(key.n, key.t, msg, ml, err);
	NCS::IMPL(encode)(err, key, enc_err);

	::serialize(std::cout, enc_err);
	delete[] msg;
}

void decode(const NCS::SecKey& key)
{
	NTL::vec_GF2 err, enc_err;
	char*		 msg;

	::deserialize(*state.is, enc_err);

	NCS::decode(enc_err, key, err);
	Binom::decode(key.bgc.n, key.bgc.t, err, msg);

	std::cout << msg;
	delete[] msg;
}

void sign(const NCS::SecKey& key)
{
	/*NTL::vec_GF2 err, enc_err;

	uint32_t ml = (key.bits /8);
	char* msg = new char[ml +1]();
	state.is->read(msg, ml);*/
	abort();
}

void sig_verify(const NCS::PubKey& key) { abort(); }

std::ifstream open_key()
{
	std::ifstream ret(state.path_key, std::ios::binary);

	if (!ret.is_open())
	{
		std::cerr << "Could not find file '" << state.path_key << "'" << std::endl;
		exit(EXIT_FAILURE);
	}

	return ret;
}
