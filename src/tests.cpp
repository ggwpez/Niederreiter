#include "tests.hpp"
#include "helper.hpp"
#include "printer.hpp"
#include "bgc.hpp"

#include <NTL/mat_GF2.h>

using namespace NTL;

int test()
{
	long m = 8, n = 15;
	mat_GF2 mat = random_mat_GF2(m,n);
	mat_GF2 mat_l = getLeftSubMatrix(mat),
			mat_r = getRightSubMatrix(mat),
			merged = mat_merge_colls(mat_l, mat_r),
			id_l = mat_merge_ID_left(mat_r),
			id_r = mat_merge_ID_right(mat_l);

	/*std::clog << "mat\n" << print(mat)
			  << "\nmat_l\n" << print(mat_l)
			  << "\nmat_r\n" << print(mat_r)
			  << "\nmerged\n" << print(merged)
			  << "\nid_l\n" << print(id_l)
			  << "\nid_r\n" << print(id_r)
			  << "\n";*/

	assert(ident_mat_GF2(id_l.NumRows()) == getLeftSubMatrix(id_l));
	assert(ident_mat_GF2(id_r.NumRows()) == getRightSubMatrix(id_r));
	assert(merged == mat);
	return 0;
}

bool test_petterson(bgc_t const& bgc)
{
	NTL::vec_GF2 message;
	message.SetLength(bgc.k());

	for (int i = 0; i < message.length(); ++i)
		message[i] = NTL::random_GF2();

	NTL::vec_GF2 resistant = bgc.encode(message);
	NTL::vec_GF2 modded = resistant;

	for (int i = 0; i < bgc.errors; ++i)
		modded[i*2] = NTL::GF2(1) -resistant[i*2];	// TODO random

	auto e = bgc.syndrom_decode_2(modded);
	auto recovered = modded +e;

	auto worked = (recovered == resistant);
	if (! worked)
		std::cerr << "Message:   " << print(message) << std::endl
				  << "Resistant:  " << print(resistant) << std::endl
				  << "Modded:    " << print(modded) << std::endl
				  << "Recovered: " << print(recovered) << std::endl
				  << "(Recovered != Resistant)" << std::endl;

	return worked;
}

mat_GF2 canonical_check_matrix(bgc_t const& bgc)
{

}
