#pragma once

#include "printer.hpp"
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2E.h>

#include <vector>
#include <memory>

typedef NTL::vec_GF2E support_t;

// Binary Goppa Code
class BGC
{
public:
	BGC() = default;
	BGC(long m, long n, long t, NTL::GF2X const& f, NTL::GF2EX const& g, const NTL::GF2E& gen, support_t const& L, const NTL::mat_GF2& H);

	static BGC create(long m, long n, long t);
	static BGC create_with_field(long m, long n, long t, NTL::GF2X const& f);
	static BGC create_with_gp_poly(long m, long n, long t, NTL::GF2EX const& g);
	static BGC create_with_field_and_gp_poly(long m, long n, long t, NTL::GF2X const& f, NTL::GF2EX const& g);

	static void calculate_f(long m, NTL::GF2X& f);
	static void calculate_g(long t, NTL::GF2EX& g);
	static void calculate_gen(NTL::GF2X f, NTL::GF2E& gen);

	static void calculate_H(long t, support_t L, NTL::GF2EX const& g, NTL::mat_GF2& H);
	static void calculate_L(long n, NTL::GF2E const& gen, support_t& L);

	static void check_args(long m, long n, long t, NTL::GF2X const& f, NTL::GF2EX const& g);
	long k() const;
	// Returns the number of elements of the Galois field
	NTL::ZZ order() const;

	void syndrom_decode(NTL::vec_GF2 const& c, NTL::vec_GF2& e) const;
	//NTL::GF2EX berlekamp_massey(NTL::vec_GF2 const& c) const;

	std::string to_str() const;

	void sqr_root(NTL::GF2EX const& p, NTL::GF2EX& res) const;

	void calculate_sc(NTL::vec_GF2 const& c, NTL::GF2EX& sc) const;
	void calculate_vc(NTL::GF2EX const& sc, NTL::GF2EX& vc) const;

	void calculate_error(NTL::GF2EX const& poly, NTL::vec_GF2& e) const;

	// L[a] * L[b]
	size_t mul_L_elements(size_t a, size_t b) const;
	// L[a] ^ b
	size_t power_L_elements(size_t a, long exp) const;

	long m, n, t;
	NTL::GF2X f;
	NTL::GF2EX g;
	NTL::GF2E gen;
	NTL::vec_GF2E L;
	NTL::mat_GF2 H;
};
