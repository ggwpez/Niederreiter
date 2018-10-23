#include "bgc.hpp"
#include "ncs.hpp"
#include "hasher.hpp"
#include "helper.hpp"
#include "tests.hpp"
#include "rand_helper.hpp"
#include "binom.hpp"

#include <iostream>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <istream>

#include <getopt.h>

#include <NTL/vec_GF2.h>

using namespace NTL;

static std::chrono::steady_clock::time_point start;
static void timer_start()
{
	start = std::chrono::steady_clock::now();
}

static int64_t timer_get_elapsed()
{
	return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() -start).count();
}

#define mERR 0
#define mGEN 1
#define mENC 2
#define mDEC 3
#define mCPK 4

static struct
{
	long m = -1,t = -1, n = -1;
	int mode = mERR;

	std::string path_key;

	NCS::KeyPair* key = nullptr;
} state;

static struct option long_options[] =
{
	/* Flags */
	{ "help",	    no_argument, 0, 'h'},
	//{ "create-key", no_argument, &state.mode, mCREATE_KEY },

	{ "gen",    no_argument, 0, 'g' },
	{ "enc",    no_argument, 0, 'e' },
	{ "dec",    no_argument, 0, 'd' },
	{ "cpk",   	no_argument, 0, 'c' },
	//{ "verbosity",  optional_argument, 0, 'v' },

	{ "key",	required_argument, 0, 'k'},

	/* NULL descriptor */
	{0, 0, 0, 0}
};

// Input from stdin
// Output to stdout
// --gen -m 13 -t 117		// Writes sk to stdout
// --enc --key=key.pk		// Encrypts stdin
// --dec --key=key.sk		// Decryptis stdin
// --cpk --key=key.sk		// Calculates pk for given sk

void parse_options(int argc, char** argv)
{
	int c;
	while (1)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long(argc, argv, "hgedck:m:n:t:", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
			{
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag)
					break;
			} break;
			case 'h':
			{
				std::cerr << "HELP" << std::endl;
				exit(0);
			} /// break;
			case 'g':
			{
				state.mode = mGEN;
			} break;
			case 'e':
			{
				state.mode = mENC;
			} break;
			case 'd':
			{
				state.mode = mDEC;
			} break;
			case 'c':
			{
				state.mode = mCPK;
			} break;
			case 'k':
			{
				state.path_key = std::string(optarg);
			} break;
			case 'm':
			{
				state.m = std::strtol(optarg, nullptr, 10);
			} break;
			case 'n':
			{
				state.n = std::strtol(optarg, nullptr, 10);
			} break;
			case 't':
			{
				state.t = std::strtol(optarg, nullptr, 10);
			} break;
			case '?': default:
			{
				/* getopt_long already printed an error message. */
				exit(EXIT_FAILURE);
			}
		}
	}
}

std::ifstream open_key()
{
	std::ifstream ret(state.path_key, std::ios::binary);

	if (! ret.is_open())
	{
		std::cerr << "Could not find file '" << state.path_key << "'" << std::endl;
		exit(EXIT_FAILURE);
	}

	return ret;
}

int main(int argc, char** argv)
{
	parse_options(argc, argv);

	if (state.mode == mGEN)	// KEYGEN
	{
		if (state.m < 0 || state.n < 0 || state.t < 0)
		{
			std::cerr << "Error, -m,-n and -t must be specified in --gen use --help for help" << std::endl;
			return 1;
		}

		BGC bgc = BGC::create(state.m,state.n,state.t);
		NCS::KeyPair keys = NCS::keygen(bgc);

		keys.m_sk.serialize(std::cout);
		std::flush(std::cout);
	}
	else if (state.mode == mENC ||
			 state.mode == mDEC ||
			 state.mode == mCPK)
	{
		state.key = new NCS::KeyPair();
		std::ifstream fs = open_key();

		/**/ if (state.mode == mENC)	// ENC
		{
			state.key->m_pk.deserialize(fs);


		}
		else if(state.mode == mDEC)		// DEC
		{
			state.key->m_sk.deserialize(fs);

			// Do decryption

		}
		else							// CPK
		{
			state.key->m_sk.deserialize(fs);
			state.key->reconstruct_pk();		// spooky
		}
	}
	else
	{
		std::cerr << "No operation mode specified, use --help for help" << std::endl;
		return 1;
	}
	return 0;

	SetSeed(ZZ(std::time(nullptr)));	// FIXME use secure random stream

	int const m = 13,
			  n = (1 << m),	// 8192
			  t = 117;

	timer_start();
	BGC bgc = BGC::create(m,n,t);
	std::cout << "BGC generated in " << timer_get_elapsed() << " µs" << std::endl
			  << bgc.to_str() << std::endl;
	timer_start();
	NCS::KeyPair keys = NCS::keygen(bgc);
	std::cout << "NCS generated in " << timer_get_elapsed() << " µs" << std::endl;

	char* msg,* dec_msg;
	// Create a random message that we want to send
	{
		msg = new char[n]();
		//ZZ x = ; /*conv<ZZ>(Binom::coeff(ZZ(n), ZZ(t)) -1);*///RandomBits_ZZ(879);
		//std::cout << "Testing with\n" << x << '\n';
		//BytesFromZZ(reinterpret_cast<unsigned char*>(msg), x, n);		// FIXME
		std::strcpy(msg, "TEST DATA GOES HERE");
	}

	NTL::vec_GF2 err, enc_err, dec_err;
	timer_start();
	Binom::encode(n, t, msg, err);
	std::cout << "Binom Encoded in " << timer_get_elapsed() << " µs" << std::endl;
	timer_start();
	NCS::encode(err, keys.m_pk, enc_err);
	std::cout << "NCS Encoded in   " << timer_get_elapsed() << " µs" << std::endl;
	timer_start();
	NCS::decode(enc_err, keys.m_sk, dec_err);
	std::cout << "NCS Decoded in   " << timer_get_elapsed() << " µs" << std::endl;
	timer_start();
	Binom::decode(n, t, dec_err, dec_msg);
	std::cout << "Binom Decoded in " << timer_get_elapsed() << " µs" << std::endl;

	std::cout << "msg            " << msg << std::endl
			  << "dec_msg        " << dec_msg << std::endl;
			  //<< "err             " << print(err) << std::endl
			  //<< "enc_err         " << print(enc_err) << std::endl
			  //<< "dec_err         " << print(dec_err) << std::endl
	//		  << "dec_msg         " << dec_msg << std::endl;
	/*std::cout << "HASH_OF(msg)       " << HASH_OF(msg) << std::endl
			  // HASH_OF(cypher) should be 9019764880683253672 if no parameters were changed
			  << "HASH_OF(cypher)    " << HASH_OF(enc_err) << std::endl
			  << "HASH_OF(recovered) " << HASH_OF(dec_err) << std::endl;*/

	if (std::strncmp(msg, dec_msg, n))
		throw std::runtime_error("Decoding error");

	delete msg;
	return 0;
}
