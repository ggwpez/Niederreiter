#pragma once

#include "printer.hpp"
#include <NTL/GF2E.h>
#include <NTL/GF2EX.h>
#include <NTL/mat_GF2E.h>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2E.h>

#include <vector>
#include <memory>

// Binary Goppa Code
class bgc_t
{
public:
	bgc_t(long dimenion, long non_root_points, long errors, NTL::GF2X* supposed_f = nullptr, NTL::GF2EX* supposed_g = nullptr);

	// Vectors passed to encode() need to be k() in length
	long k() const;
	bool is_codeword(NTL::vec_GF2 const&);
	// Returns the number of elements of the Galois field
	NTL::ZZ GF_order();

	NTL::vec_GF2 encode(NTL::vec_GF2 const& msg);
	NTL::vec_GF2 patterson_decode(const NTL::vec_GF2& c);

	std::string to_str() const;

	NTL::GF2EX sqr_root(NTL::GF2EX const& p);

	// Helper
	NTL::GF2EX calculate_fi(long i);
	// Syndrom function for code word c
	NTL::GF2EX calculate_sc(const NTL::vec_GF2& c);
	// Feed with syndrom function
	NTL::GF2EX calculate_vc(const NTL::GF2EX& sc);
	NTL::GF2E calculate_sc_val(const NTL::vec_GF2& c);
	// Error locator polynom
	NTL::GF2EX calculate_sigma(NTL::vec_GF2 e, long spare_i = -1);
	NTL::GF2EX calculate_small_omega(NTL::vec_GF2 e);

	std::vector<size_t> get_L_indices(std::vector<NTL::GF2E> const& l) const;
	NTL::vec_GF2 calculate_error_vector(const std::vector<size_t>& L_indices);


	//NTL::GF2EX sth_root(NTL::GF2EX const& p, NTL::ZZ const& s);

	long dimension, non_root_points, errors;
	NTL::GF2EX g;
	NTL::GF2E gen;
	std::vector<NTL::GF2E> L;
	// H* is the H from the amsterdam paper
	NTL::mat_GF2E H, H_star;
	NTL::mat_GF2 H_bin, H_star_bin,
				 G,     G_star;
	std::unique_ptr<NTL::GF2X> f;
	NTL::GF2EX id;
};
