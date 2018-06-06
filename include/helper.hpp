#pragma once
#include "printer.hpp"
#include <NTL/mat_GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <cassert>
#include <vector>

NTL::GF2E call(NTL::GF2EX const& p, NTL::GF2E const& x);
NTL::GF2E call_horner(NTL::GF2EX const& p, NTL::GF2E const& x);
// Make Polynomial monic
void monice(NTL::GF2EX& ret, NTL::GF2EX poly);
NTL::GF2EX monice(NTL::GF2EX const&poly);

NTL::mat_GF2 create_rand_permutation(size_t s);

void compute_systematic_form(NTL::mat_GF2 const& H, NTL::mat_GF2& sInv, NTL::mat_GF2& m, NTL::mat_GF2& p);
NTL::mat_GF2 getLeftSubMatrix(NTL::mat_GF2 const& mat);
NTL::mat_GF2 getRightSubMatrix(NTL::mat_GF2 const& mat);

NTL::mat_GF2 mat_merge_colls(NTL::mat_GF2 const& a, NTL::mat_GF2 const& b);
NTL::mat_GF2 mat_merge_ID_left(NTL::mat_GF2 const& b);
NTL::mat_GF2 mat_merge_ID_right(NTL::mat_GF2 const& a);

template<typename T>
inline NTL::ZZ width(NTL::Vec<T> const& vec)
{
	NTL::ZZ ret = NTL::ZZ::zero();

	for (long i = 0; i < vec.length(); ++i)
		if (! NTL::IsZero(vec[i]))
			++ret;

	return ret;
}

std::vector<NTL::GF2E> find_roots(NTL::GF2EX const&);

//template<typename NTL::GF2E>
NTL::Mat<NTL::GF2> trace_construct(NTL::Mat<NTL::GF2E> const& m);
NTL::Mat<NTL::GF2> trace_construct_2(NTL::Mat<NTL::GF2E> const& m);

// FIXME: Please someone tell me how to do this better, its horrible!
inline NTL::GF2E generate_GF2E(uint64_t i)
{
	NTL::GF2E x = NTL::GF2E::zero();
	// Please someone tell me how to do this better, its horrible!
	for (uint64_t j = 0; j < 64; ++j)
		if ((uint64_t(1) << j) & i)
			NTL::SetCoeff(x.LoopHole(), j);

	return x;
}

///
/// https://www.cdc.informatik.tu-darmstadt.de/fileadmin/user_upload/Group_CDC/Documents/Lehre/SS15/pqc/code_based_Paulo_Barreto.pdf
/// page 56
///
inline NTL::GF2EX calc_sigma(const NTL::GF2EX& a, const NTL::GF2EX& b, const NTL::GF2EX& g)
{
	NTL::GF2EX F = a, G = b, B = NTL::conv<NTL::GF2EX>("[[1]]"), C = NTL::GF2EX::zero();
	long t = NTL::deg(b);

	while (NTL::deg(G) > (t /2))
	{
		NTL::swap(F, G); NTL::swap(B, C);

		while (NTL::deg(F) >= NTL::deg(G))
		{
			long j = NTL::deg(F) -NTL::deg(G);
			auto h = NTL::LeadCoeff(F) /NTL::LeadCoeff(G);

			F = F -h *NTL::PowerXMod(j, g) *G;
			B = B -h *NTL::PowerXMod(j, g) *C;
		}
	}

	return (G*G) +NTL::PowerXMod(1, g) *(C*C);
}

NTL::vec_GF2E to_ext_field_poly(const NTL::vec_GF2& vec, NTL::GF2X const& field);
