#pragma once

#include <iostream>

#define mERR  0
#define mGEN  'g'
#define mENC  'e'
#define mDEC  'd'
#define mCPK  'c'
#define mSINF 's'
#define mPINF 'p'

#define cBLK 0
#define cCRT 1
#define cCBC 2

namespace NCS {
	class KeyPair;
}

typedef struct
{
	long m = -1,t = -1, n = -1;
	int mode = mERR;
	int crypto_mode = cBLK;

	std::istream* is = &std::cin;
	std::string path_key;
} state_t;

state_t parse_args(int argc, char** argv);
void print_help();
void print_version();
