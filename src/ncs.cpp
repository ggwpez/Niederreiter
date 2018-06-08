#include "ncs.hpp"
#include "helper.hpp"
using namespace NTL;

void NCS::compute_systematic_form(mat_GF2 const& H, mat_GF2& sInv, mat_GF2& m, mat_GF2& p)
{
	long n = H.NumCols();
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

	mat_GF2 H = bgc.H;
	NCS::compute_systematic_form(H, sInv, m, p);

	return KeyPair{ SecKey{ sInv, p, bgc }, PubKey{ m } };
}

void NCS::encode(NTL::vec_GF2 const& msg, NCS::PubKey const& key, NTL::vec_GF2& cipher)
{
	mat_GF2 H = mat_merge_ID_left(key.h);	// TODO

	mul(cipher, H, msg);
}

void NCS::decode(NTL::vec_GF2 const& cipher, NCS::SecKey const& key, NTL::vec_GF2& msg)
{
	vec_GF2 syndrome, pe;

	mul(syndrome, key.Si, cipher);
	key.bgc.syndrom_decode(syndrome, pe);
	mul(msg, pe, key.p);
}
