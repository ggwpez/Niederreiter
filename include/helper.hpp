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
// Make Polynomial monic
void monice(NTL::GF2EX poly, NTL::GF2EX& ret);
NTL::GF2EX monice(NTL::GF2EX const&poly);

NTL::mat_GF2 create_rand_permutation(size_t s);

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

//template<typename NTL::GF2E>
void trace_construct(NTL::Mat<NTL::GF2E> const& m, NTL::Mat<NTL::GF2>& H);


///
/// https://www.cdc.informatik.tu-darmstadt.de/fileadmin/user_upload/Group_CDC/Documents/Lehre/SS15/pqc/code_based_Paulo_Barreto.pdf
/// page 56
///
void calculate_sigma(NTL::GF2EX const& a, NTL::GF2EX const& b, NTL::GF2EX const& g, NTL::GF2EX& sigma);
