#pragma once

#include <vector>
#include <NTL/mat_GF2.h>
#include <NTL/vector.h>

//typedef NTL::mat_GF2 perm_GF2;
// About to implement a more efficent representation of a mat_GF2 permutation matrix
typedef NTL::Vec<uint16_t> perm_GF2;

perm_GF2 create_rand_permutation(size_t s);

void mul(NTL::mat_GF2&out, const perm_GF2&perm, const NTL::mat_GF2&mat);
void mul(NTL::mat_GF2&out, NTL::mat_GF2 const&mat, perm_GF2 const&perm);
void mul(NTL::vec_GF2&out, NTL::vec_GF2 const&vec, perm_GF2 const&perm);
