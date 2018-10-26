#include "serializer.hpp"

using namespace NTL;

void serialize(std::ostream& out, const uint32_t& val)
{
	uint32_t cpy = val;
	out.write(reinterpret_cast<char const*>(&cpy), 4);
}

void deserialize(std::istream& in, uint32_t& val)
{
	uint32_t cpy;
	in.read(reinterpret_cast<char*>(&cpy), 4);
	val = cpy;
}

void serialize(std::ostream& out, const ZZ& val)
{
	uint32_t len  = uint32_t(NumBytes(val)),
			 size = uint32_t(val.size());
	serialize(out, size);
	serialize(out, len);
	uint8_t* data = new uint8_t[len +1]();

	BytesFromZZ(data, val, len);
	out.write(reinterpret_cast<char*>(data), len);

	delete[] data;
}

void deserialize(std::istream& in, ZZ& val)
{
	uint32_t len, size;
	deserialize(in, size);
	deserialize(in, len);
	uint8_t* data = new uint8_t[len +1]();
	val.SetSize(size);

	in.read(reinterpret_cast<char*>(data), len);
	ZZFromBytes(val, data, len);

	delete[] data;
}

template<typename T>
void serialize(std::ostream& out, const Vec<T>& vec)
{
	serialize(out, uint32_t(vec.length()));

	for (uint32_t i = 0; i < vec.length(); ++i)
		serialize(out, vec[i]);
}

template<typename T>
void deserialize(std::istream& in, Vec<T>& vec)
{
	uint32_t l;
	deserialize(in, l);
	vec.FixLength(l);

	for (uint32_t i = 0; i < l; ++i)
	{
		T tmp;
		deserialize(in, tmp);
		vec.put(i, tmp);
		// deserialize(in, vec.at(i)); // why not
	}
}

template<typename T>
void serialize(std::ostream& out, const NTL::Mat<T>& mat)
{
	uint32_t r = mat.NumRows(),
			 c = mat.NumCols();
	serialize(out, r);
	serialize(out, c);

	for (long i = 0; i < r; ++i)
		for (long j = 0; j < c; ++j)
			serialize(out, mat[i][j]);
}

template<typename T>
void deserialize(std::istream& in, NTL::Mat<T>& mat)
{
	uint32_t r,c;
	deserialize(in, r);
	deserialize(in, c);

	mat.SetDims(r,c);
	for (long i = 0; i < r; ++i)
		for (long j = 0; j < c; ++j)
		{
			T tmp;
			deserialize(in, tmp);
			mat.put(i, j, tmp);
			// deserialize(in, mat[i][j]); // why not
		}
}

void serialize(std::ostream& out, const GF2& val)
{
	serialize(out, conv<ZZ>(val));
}

void deserialize(std::istream& in, GF2& val)
{
	ZZ tmp;
	deserialize(in, tmp);
	val = conv<GF2>(tmp);
}

void serialize(std::ostream& out, const GF2E& val)
{
	serialize(out, conv<GF2X>(val));
}

void deserialize(std::istream& in, GF2E& val)
{
	GF2X tmp;
	deserialize(in, tmp);
	val = conv<GF2E>(tmp);
}

//
void serialize(std::ostream& out, const GF2X& val)
{
	uint32_t len  = uint32_t(NumBytes(val)),
			 size = uint32_t(deg(val));
	serialize(out, size);
	serialize(out, len);
	uint8_t* data = new uint8_t[len +1]();

	BytesFromGF2X(data, val, len);
	out.write(reinterpret_cast<char*>(data), len);

	delete[] data;
}

void deserialize(std::istream& in, GF2X& val)
{
	uint32_t len, size;
	deserialize(in, size);
	deserialize(in, len);
	uint8_t* data = new uint8_t[len +1]();
	val.SetLength(size);

	in.read(reinterpret_cast<char*>(data), len);
	GF2XFromBytes(val, data, len);

	delete[] data;
}

//
void serialize(std::ostream& out, const GF2EX& val)
{
	serialize(out, deg(val));

	for (int i = 0; i <= deg(val); ++i)
		serialize(out, val[i]);
}

void deserialize(std::istream& in, GF2EX& val)
{
	uint32_t l;
	deserialize(in, l);
	val.SetLength(l +1);

	for (int i = 0; i <= l; ++i)
		deserialize(in, val[i]);
}

template void serialize<NTL::GF2>(std::ostream& out, NTL::Vec<NTL::GF2> const&);
template void serialize<NTL::GF2E>(std::ostream& out, NTL::Vec<NTL::GF2E> const&);
template void serialize<NTL::GF2EX>(std::ostream& out, NTL::Vec<NTL::GF2EX> const&);

template void deserialize<NTL::GF2>(std::istream&, NTL::Vec<NTL::GF2>&);
template void deserialize<NTL::GF2E>(std::istream&, NTL::Vec<NTL::GF2E>&);
template void deserialize<NTL::GF2EX>(std::istream&, NTL::Vec<NTL::GF2EX>&);

template void serialize<NTL::GF2>(std::ostream& out, NTL::Mat<NTL::GF2> const&);
template void serialize<NTL::GF2E>(std::ostream& out, NTL::Mat<NTL::GF2E> const&);
template void serialize<NTL::GF2EX>(std::ostream& out, NTL::Mat<NTL::GF2EX> const&);

template void deserialize<NTL::GF2>(std::istream&, NTL::Mat<NTL::GF2>&);
template void deserialize<NTL::GF2E>(std::istream&, NTL::Mat<NTL::GF2E>&);
template void deserialize<NTL::GF2EX>(std::istream&, NTL::Mat<NTL::GF2EX>&);
