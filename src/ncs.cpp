#include "ncs.hpp"
#include "hasher.hpp"
#include "helper.hpp"
#include <iterator>

using namespace NTL;

void NCS::compute_systematic_form(mat_GF2 const& H, mat_GF2& sInv, mat_GF2& m, mat_GF2& p)
{
	size_t n = size_t(H.NumCols());
	mat_GF2 hp;

	do
	{
		p = create_rand_permutation(n);
		mul(hp, H, p);
		sInv = getLeftSubMatrix(hp);
	} while (IsZero(determinant(sInv)));

	mat_GF2 shp = inv(sInv) *hp;
	m = getRightSubMatrix(shp);
}

NCS::KeyPair NCS::keygen(BGC const& bgc)
{
	mat_GF2 sInv, m, p;

	//mat_GF2 H = bgc.H; // CHANGED
	NCS::compute_systematic_form(bgc.H, sInv, m, p);

	return KeyPair{ SecKey{ sInv, p, bgc }, PubKey{ m } };
}

void NCS::encode(NTL::vec_GF2 const& msg, NCS::PubKey const& key, NTL::vec_GF2& cipher)
{

	//mat_GF2 H = mat_merge_ID_left(key.h);	// TODO

	//mul(cipher, H, msg);
	mat_mul_right_compact(key.h, msg, cipher);
}

void NCS::decode(NTL::vec_GF2 const& cipher, NCS::SecKey const& key, NTL::vec_GF2& msg)
{
	vec_GF2 syndrome, pe;

	mul(syndrome, key.Si, cipher);
	key.bgc.syndrom_decode(syndrome, pe);
	mul(msg, pe, key.p);
}

NCS::SecKey::SecKey(const mat_GF2& Si, const perm_GF2& p, const BGC& bgc)
	: Si(Si), p(p), bgc(bgc)
{

}

void NCS::SecKey::serialize(std::ostream& out) const
{
	bgc.serialize(out);

	::serialize(out, Si);
	::serialize(out, p);

	std::cerr << "HASH_OF(sInv) " << HASH_OF(Si._mat__rep) << std::endl;
}

void NCS::SecKey::deserialize(std::istream& in)
{
	bgc.deserialize(in);

	::deserialize(in, Si);
	::deserialize(in, p);

	std::cerr << "HASH_OF(sInv) " << HASH_OF(Si._mat__rep) << std::endl;
}

NCS::PubKey::PubKey(const mat_GF2& h)
		  : h(h)
{

}

void NCS::PubKey::serialize(std::ostream& ) const
{

}

void NCS::PubKey::deserialize(std::istream& )
{

}

NCS::KeyPair::KeyPair(const NCS::SecKey& sk, const NCS::PubKey& pk)
	: m_sk(sk), m_pk(pk)
{

}

void NCS::KeyPair::reconstruct_pk()
{
	mat_GF2 hp;

	mul(hp, m_sk.bgc.H, m_sk.p);
	mat_GF2 shp = inv(m_sk.Si) *hp;

	this->m_pk = getRightSubMatrix(shp);
}

void NCS::KeyPair::serialize(std::ostream&) const
{

}

void NCS::KeyPair::deserialize(std::istream&)
{

}
