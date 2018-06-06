#include "printer.hpp"
#include <NTL/GF2EX.h>
#include <iomanip>
#include <sstream>

std::string print_x_power(long i)
{
	if (! i)
		return std::string();
	else if (i == 1)
		return std::string("x");
	else
		return std::string("x^" +std::to_string(i));
}

// Poly over GF(2)
std::string print(NTL::GF2X const& p)
{
	if (NTL::IsZero(p))
		return "0";
	if (NTL::IsOne(p))
		return "1";

	std::ostringstream ss;
	bool place_plus = false;

	for (int i = NTL::deg(p); i --> 0;)
	{
		NTL::GF2 c = NTL::coeff(p, i);
		if (place_plus && ! NTL::IsZero(c))
		{
			ss << " +";
			place_plus = false;
		}

		if (! IsZero(c))
		{
			if (! IsOne(c))
				ss << print(c);
			ss << print_x_power(i);
			place_plus = true;
		}
	}

	return ss.str();
}

// Extension Ring of GF(2)
std::string print(NTL::GF2E const& _p)
{
	NTL::GF2X p = NTL::conv<NTL::GF2X>(_p);
	NTL::ZZ res = NTL::ZZ(0);

	for (int i = 0; i <= deg(p); ++i)
	{
		NTL::GF2 c = NTL::coeff(p, i);

		if (! NTL::IsZero(c))
			res += NTL::ZZ(1) << i;
	}

	std::ostringstream ss;
	ss << res;
	return ss.str();
}

std::string print(NTL::vec_GF2E const& l)
{
	std::ostringstream ss;

	for (int i = 0; i < l.length(); ++i)
		ss << print(l[i]) << (i == l.length() -1 ? "" : ", ");

	return ss.str();
}

std::string print(NTL::vec_pair_GF2EX_long const& vec)
{
	std::stringstream ss;

	for (long i = 0; i < vec.length(); ++i)
	{
		ss << (vec[i].b) << " times " << print(vec[i].a) << (i == vec.length() -1 ? "" : ", ");
	}

	return ss.str();
}

// Poly over GF2E
std::string print(NTL::GF2EX const& p)
{
	if (NTL::IsZero(p))
		return "0";
	if (NTL::IsOne(p))
		return "1";

	std::ostringstream ss;
	bool place_plus = false;

	for (int i = NTL::deg(p); i --> 0;)
	{
		NTL::GF2E c = NTL::coeff(p, i);
		if (place_plus && ! NTL::IsZero(c))
		{
			ss << " +";
			place_plus = false;
		}

		if (! IsZero(c))
		{
			if (! IsOne(c))
				ss << print(c);
			ss << print_x_power(i);
			place_plus = true;
		}
	}

	return ss.str();
}

std::string print(const NTL::GF2& p)
{
	return (NTL::IsZero(p) ? "0" : "1");
}
