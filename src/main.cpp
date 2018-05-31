#include "bgc.hpp"
#include "helper.hpp"

#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/mat_GF2.h>

#include <csignal>
#include <iostream>
#include <ctime>
#include <iomanip>

int main(int, char**)
{
	NTL::SetSeed(NTL::ZZ(time(0)));
	// Length
	int const m = 10;		// 10
	// Diemnsion
	int const n = 1 << m;
	int const t = 50;		// 9

	// Create Goppa Binary Code
	bgc_t bgc(m, n, t);
	std::cout << bgc.to_str() << std::endl;

	NTL::vec_GF2 message;
	message.SetLength(bgc.k());

	for (int i = 0; i < message.length(); ++i)
		message[i] = NTL::random_GF2();

	NTL::vec_GF2 resistant = bgc.encode(message);
	NTL::vec_GF2 modded = resistant;

	for (int i = 0; i < t; ++i)
		modded[i*2] = NTL::GF2(1) -resistant[i*2];

	auto recovered = bgc.patterson_decode(modded);

	bool worked = (recovered == resistant);
	std::cout << "Message:   " << print(message) << std::endl
			  << "Resistant:  " << print(resistant) << std::endl
			  << "Modded:    " << print(modded) << std::endl
			  << "Recovered: " << print(recovered) << std::endl
			  << "(Recovered == Resistant) == " << worked << std::endl;
}
