#pragma once

#include "perm_gf2.hpp"
#include "bgc.hpp"
#include <NTL/mat_GF2.h>

/// Niederreiter Cryptosystem
namespace NCS
{
	struct SecKey
	{
		NTL::mat_GF2 Si;
		perm_GF2 p;
		BGC bgc;
	};

	struct PubKey
	{
		NTL::mat_GF2 h;
	};

	struct KeyPair
	{
		SecKey sec;
		PubKey pub;
	};

	KeyPair keygen(BGC const&);
	void encode(NTL::vec_GF2 const&, PubKey const&, NTL::vec_GF2&);
	void decode(NTL::vec_GF2 const&, SecKey const&, NTL::vec_GF2&);

	// You dont need to call this
	void compute_systematic_form(NTL::mat_GF2 const& H, NTL::mat_GF2& sInv, NTL::mat_GF2& m, NTL::mat_GF2& p);
}
