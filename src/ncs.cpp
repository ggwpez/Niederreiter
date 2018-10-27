#include "ncs.hpp"
#include "serializer.hpp"
#include "helper.hpp"
#include "hasher.hpp"
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

	return KeyPair{ SecKey{ sInv, p, bgc }, PubKey( m, bgc.n, bgc.t ) };
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

	::serialize(out, p);
	::serialize(out, Si);
}

void NCS::SecKey::deserialize(std::istream& in)
{
	bgc.deserialize(in);

	::deserialize(in, p);
	::deserialize(in, Si);
}

NCS::PubKey::PubKey(const mat_GF2& h, uint32_t n, uint32_t t)
		  : h(h), n(n), t(t)
{

}

void NCS::PubKey::serialize(std::ostream& out) const
{
	::serialize(out, n);
	::serialize(out, t);

	::serialize(out, h);
}

void NCS::PubKey::deserialize(std::istream& in)
{
	::deserialize(in, n);
	::deserialize(in, t);

	::deserialize(in, h);
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

	this->m_pk = PubKey(getRightSubMatrix(shp), m_sk.bgc.n, m_sk.bgc.t);
}

void NCS::KeyPair::serialize(std::ostream&) const
{
	throw std::runtime_error("Unimplemented");
}

void NCS::KeyPair::deserialize(std::istream&)
{
	throw std::runtime_error("Unimplemented");
}
