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
