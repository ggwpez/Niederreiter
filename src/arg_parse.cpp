#include "arg_parse.hpp"

#include <fstream>
#include <getopt.h>
#include <NTL/ZZ.h>

#define kSEED 1000
#define kMODE 1001

static struct option long_options[] =
{
	/* Flags */
	{ "help",			no_argument, nullptr, 'h' },
	{ "version",	no_argument, nullptr, 'v' },

	{ "gen",		no_argument, nullptr, mGEN },
	{ "enc",		no_argument, nullptr, mENC },
	{ "dec",		no_argument, nullptr, mDEC },
	{ "sig",		no_argument, nullptr, mSIG },
	{ "sfy",		no_argument, nullptr, mSFY },
	{ "cpk",		no_argument, nullptr, mCPK },
	{ "sk-info",	no_argument, nullptr, mSINF },
	{ "pk-info",	no_argument, nullptr, mPINF },

	{ "seed",		required_argument, nullptr, kSEED },
	{ "mode",		required_argument, nullptr, kMODE },

	{ "key",		required_argument, nullptr, 'k' },
	{ "input",		required_argument, nullptr, 'i' },

	/* NULL descriptor */
	{ nullptr, 0, nullptr, 0 }
};

state_t parse_args(int argc, char** argv)
{
	state_t state;
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
				{ }
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
			case mGEN: case mENC: case mDEC: case mCPK: case mSINF: case mPINF: case mSIG: case mSFY:
			{
				state.mode = c;
			} break;
			case kSEED:
			{
				NTL::SetSeed(NTL::ZZ(std::strtol(optarg, nullptr, 10)));
			} break;
			case kMODE:
			{
				state.crypto_mode = int(std::strtol(optarg, nullptr, 10));
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

	return state;
}

void print_help()
{
	std::cout
		<< "Example for file 'msg':" << std::endl
		<< " ./niederreiter --gen -m10 -n1024 -t30 > key.sk" << std::endl
		<< " ./niederreiter --cpk --key=key.sk > key.pk" << std::endl
		<< " ./niederreiter --enc --input=msg		 --key=key.pk > msg.enc" << std::endl
		<< " ./niederreiter --dec --input=msg.enc --key=key.sk" << std::endl
		<< "Explanation:" << std::endl
		<< " 1: Generate secret key and write to key.sk" << std::endl
		<< "		The arguments m,n,t are the security parameter λ" << std::endl
		<< " 2: Calculate public key from key.sk and write to key.pk" << std::endl
		<< " 3: Encrypt the file 'msg' and write cipher to msg.enc" << std::endl
		<< " 4: Decrypt msg from msg.enc and write to cout" << std::endl
		<< "		The output should match msg, if it wasnt too long" << std::endl
	;
}

void print_version()
{
	std::cout << "Git hash " << std::hex << GIT_HASH << std::endl;
}
