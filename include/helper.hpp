#pragma once
#include "printer.hpp"
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <cassert>
#include <vector>

NTL::GF2E call(NTL::GF2EX const& p, NTL::GF2E const& x);
// Make Polynomial monic
void monice(NTL::GF2EX& ret, NTL::GF2EX poly);
NTL::GF2EX monice(NTL::GF2EX const&poly);

std::vector<NTL::GF2E> find_roots(NTL::GF2EX const&);

//template<typename NTL::GF2E>
inline NTL::Mat<NTL::GF2> trace_construct(NTL::Mat<NTL::GF2E> const& m)
{
	NTL::Mat<NTL::GF2> ret;
	long t = NTL::GF2E::modulus().n;
	ret.SetDims(t *m.NumRows(), m.NumCols());

	for (long row = 0; row < m.NumRows(); ++row)
		for (long col = 0; col < m.NumCols(); ++col)
		{
			NTL::GF2X const& e = NTL::conv<NTL::GF2X>(m[row][col]);

			for (long i = 0; i < t; ++i)	// < deg(e) TODO
				ret.put(t *row +i, col, NTL::coeff(e, i));
		}

	return ret;
}

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
