#include "bgc.hpp"
#include "ncs.hpp"
#include "hasher.hpp"
#include "helper.hpp"
#include "tests.hpp"
#include "rand_helper.hpp"
#include "serializer.hpp"
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

#define mERR  0
#define mGEN  1
#define mENC  2
#define mDEC  3
#define mCPK  4
#define mSINF 5
#define mPINF 6

static struct
{
	long m = -1,t = -1, n = -1;
	int mode = mERR;

	std::istream* is = &std::cin;
	std::string path_key;

	NCS::KeyPair* key = nullptr;
} state;

#define kSEED 1000
static struct option long_options[] =
{
	/* Flags */
	{ "help",	    no_argument, nullptr, 'h'},
	{ "version",	no_argument, nullptr, 'v'},

	{ "gen",		no_argument, nullptr, 'g' },
	{ "enc",		no_argument, nullptr, 'e' },
	{ "dec",		no_argument, nullptr, 'd' },
	{ "cpk",		no_argument, nullptr, 'c' },
	{ "sk-info",	no_argument, nullptr, 's' },
	{ "pk-info",	no_argument, nullptr, 'p' },

	{ "seed",		required_argument, nullptr, kSEED },

	{ "key",		required_argument, nullptr, 'k' },
	{ "input",		required_argument, nullptr, 'i' },

	/* NULL descriptor */
	{ nullptr, 0, nullptr, 0 }
};

void print_help();
void print_version();

void parse_options(int argc, char** argv)
{
	int c;
	while (1)
	{
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long(argc, argv, "hgedck:m:n:t:i:spv", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
			case 0:
			{
				/* If this option set a flag, do nothing else now. */
				if (long_options[option_index].flag)
					;
			} break;
			case 'h':
			{
				print_help();
				exit(0);
			} /// break;
			case 'v':
			{
				print_version();
				exit(0);
			}
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
			case 's':
			{
				state.mode = mSINF;
			} break;
			case 'p':
			{
				state.mode = mPINF;
			} break;
			case kSEED:
			{
				SetSeed(ZZ(std::strtol(optarg, nullptr, 10)));
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
			case 'i':
			{
				std::string path(optarg);

				if (path != "-")	// not stdin
				{
					std::ifstream* fs = new std::ifstream(path, std::ios::binary);

					if (! fs->is_open())
					{
						std::cerr << "Could not find file '" << path << "'" << std::endl;
						free(fs);
						exit(EXIT_FAILURE);
					}
					else
						state.is = fs;
				}
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
	SetSeed(ZZ(std::time(nullptr)));	// FIXME use secure random stream
	parse_options(argc, argv);

	if (state.mode == mGEN)	// KEYGEN
	{
		if (state.m < 0 || state.n < 0 || state.t < 0)
		{
			std::cerr << "Error, -m, -n and -t must be specified in --gen use --help for help" << std::endl;
			return 1;
		}

		BGC bgc = BGC::create(state.m, state.n, state.t);
		NCS::KeyPair keys = NCS::keygen(bgc);

		keys.m_sk.serialize(std::cout);
	}
	else if (state.mode >= mENC && state.mode <= mPINF)
	{
		state.key = new NCS::KeyPair();
		std::ifstream fs = open_key();

		if (state.mode == mENC)	// ENC
		{
			state.key->m_pk.deserialize(fs);

			NTL::vec_GF2 err, enc_err;
			char* msg = new char[state.key->m_pk.n +1]();
			state.is->readsome(msg, state.key->m_pk.n);

			Binom::encode(state.key->m_pk.n, state.key->m_pk.t, msg, err);
			NCS::encode(err, state.key->m_pk, enc_err);

			::serialize(std::cout, enc_err);
			delete [] msg;
		}
		else if(state.mode == mDEC)		// DEC
		{
			state.key->m_sk.deserialize(fs);

			char* msg;
			NTL::vec_GF2 err, enc_err;
			::deserialize(*state.is, enc_err);

			NCS::decode(enc_err, state.key->m_sk, err);
			Binom::decode(state.key->m_sk.bgc.n, state.key->m_sk.bgc.t, err, msg);

			std::cout << msg;
			delete [] msg;
		}
		else if (state.mode == mSINF)
		{
			state.key->m_sk.deserialize(fs);

			std::cout << state.key->m_sk.bgc.to_str() << std::endl
					  << "sha256sum "; std::cout.flush();
			system((std::string("sha256sum ") +state.path_key).c_str());
		}
		else if (state.mode == mPINF)
		{
			state.key->m_pk.deserialize(fs);

			std::cout << "PK-info not yet implemented" << std::endl
					  << "sha256sum "; std::cout.flush();
			system((std::string("sha256sum ") +state.path_key).c_str());
		}
		else if (state.mode == mCPK)							// CPK
		{
			state.key->m_sk.deserialize(fs);
			state.key->reconstruct_pk();		// spooky

			state.key->m_pk.serialize(std::cout);
		}
		else
			throw std::runtime_error("Unreachable");	//SetSeed(ZZ(123));	// FIXME use secure random stream

	}
	else
	{
		std::cerr << "No operation mode specified, use --help for help" << std::endl;
		return 1;
	}

	std::cout.flush();
	return 0;

	int const m = 10,
			  n = (1 << m),	// 8192
			  t = 30;

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
		std::strcpy(msg, "TEST DATA\n");
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
	std::cout << "enc_err\n" << enc_err << '\n';

	if (std::strncmp(msg, dec_msg, n))
		throw std::runtime_error("Decoding error");

	delete msg;
	return 0;
}

void print_help()
{
	std::cout
		<< "Example for file 'msg':" << std::endl
		<< " ./niederreiter --gen -m10 -n1024 -t30 > key.sk" << std::endl
		<< " ./niederreiter --cpk --key=key.sk > key.pk" << std::endl
		<< " ./niederreiter --enc --input=msg     --key=key.pk > msg.enc" << std::endl
		<< " ./niederreiter --dec --input=msg.enc --key=key.sk" << std::endl
		<< "Explanation:" << std::endl
		<< " 1: Generate secret key and write to key.sk" << std::endl
		<< "    The arguments m,n,t are the security parameter λ" << std::endl
		<< " 2: Calculate public key from key.sk and write to key.pk" << std::endl
		<< " 3: Encrypt the file 'msg' and write cipher to msg.enc" << std::endl
		<< " 4: Decrypt msg from msg.enc and write to cout" << std::endl
		<< "    The output should match msg, if it wasnt too long" << std::endl
	;
}

void print_version()
{
	std::cout << "Git hash " << std::hex << GIT_HASH << std::endl;
}
