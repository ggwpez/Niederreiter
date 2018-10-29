#pragma once

#include "perm_gf2.hpp"
#include "bgc.hpp"
#include "serializable.hpp"
#include <NTL/mat_GF2.h>

/// Niederreiter Cryptosystem
namespace NCS
{
	struct SecKey : ISerializable
	{
		SecKey() = default;
		SecKey(NTL::mat_GF2 const& Si, perm_GF2 const& p, BGC const& bgc);

		NTL::mat_GF2 Si;
		perm_GF2 p;
		BGC bgc;

		virtual void serialize(std::ostream&) const override;
		virtual	void deserialize(std::istream&) override;
	};

	struct PubKey : ISerializable
	{
		PubKey() = default;
		PubKey(NTL::mat_GF2 const& h, uint32_t n, uint32_t t);

		NTL::mat_GF2 h;
		uint32_t n, t, bits;

		virtual void serialize(std::ostream&) const override;
		virtual void deserialize(std::istream&) override;
	};

	struct KeyPair : ISerializable
	{
		KeyPair() = default;
		KeyPair(SecKey const&, PubKey const&);

		SecKey m_sk;
		PubKey m_pk;

		void reconstruct_pk();

		virtual void serialize(std::ostream&) const override;
		virtual void deserialize(std::istream&) override;
	};

	KeyPair keygen(BGC const&);
	void encode(NTL::vec_GF2 const&, PubKey const&, NTL::vec_GF2&);
	void decode(NTL::vec_GF2 const&, SecKey const&, NTL::vec_GF2&);

	// You dont need to call this
	void compute_systematic_form(NTL::mat_GF2 const& H, NTL::mat_GF2& sInv, NTL::mat_GF2& m, perm_GF2& p);
}
