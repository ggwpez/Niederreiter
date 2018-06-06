#pragma once

#include <NTL/mat_GF2E.h>
#include <NTL/pair_GF2EX_long.h>
#include <string>
#include <sstream>
#include <iomanip>

std::string print_x_power(long i);
std::string print(NTL::GF2 const&);
std::string print(NTL::GF2X const&);
std::string print(NTL::GF2E const&);
inline std::string print(size_t s)
{
	return std::to_string(s);
}
template<typename T>
inline std::string print(NTL::Vec<T> const& vec)
{
	std::stringstream ss;

	for (long i = 0; i < vec.length(); ++i)
		ss << print(vec[i]) << (i == vec.length() -1 ? "" : ", ");

	return ss.str();
}
template<typename T, template<typename> class Vec>
inline std::string print(Vec<T> const& vec)
{
	std::stringstream ss;

	for (auto it = vec.begin(); it != vec.end(); ++it)
		ss << print(*it) << ' ';// << (it == std::prev(vec.size()) ? "" : ", ");

	return ss.str();
}

std::string print(NTL::GF2EX const& p);
template<typename T>
inline std::string print(NTL::Mat<T> const& m)
{
	std::stringstream ss;

	for (int i = 0; i < m.NumRows(); ++i)
	{
		for (int j = 0; j < m.NumCols(); ++j)
			ss << std::setw(1) << print(m.get(i, j)) << ' ';

		ss << std::endl;
	}

	return ss.str();
}
std::string print(NTL::vec_pair_GF2EX_long const& vec);
