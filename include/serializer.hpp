#pragma once

#include "printer.hpp"
#include <NTL/mat_GF2.h>
#include <NTL/GF2.h>
#include <NTL/GF2X.h>
#include <NTL/GF2EX.h>
#include <NTL/vec_GF2.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>

#include <cassert>
#include <vector>

void serialize(std::ostream&, uint32_t const&);
void deserialize(std::istream&, uint32_t&);

void serialize(std::ostream&, NTL::ZZ const&);
void deserialize(std::istream&, NTL::ZZ&);

template<typename T>
void serialize(std::ostream& out, NTL::Vec<T> const&);
template<typename T>
void deserialize(std::istream& in, NTL::Vec<T>& vec);

template<typename T>
void serialize(std::ostream&, NTL::Mat<T> const&);
template<typename T>
void deserialize(std::istream&, NTL::Mat<T>&);

void serialize(std::ostream&, NTL::GF2 const&);
void deserialize(std::istream&, NTL::GF2&);

void serialize(std::ostream&, NTL::GF2E const&);
void deserialize(std::istream&, NTL::GF2E&);

void serialize(std::ostream&, NTL::GF2X const&);
void deserialize(std::istream&, NTL::GF2X&);

void serialize(std::ostream&, NTL::GF2EX const&);
void deserialize(std::istream&, NTL::GF2EX&);
