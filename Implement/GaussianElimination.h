#ifndef GaussianElimination_H
#define GaussianElimination_H

#include "NTL/ZZ_p.h"
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>

using namespace NTL;

void GE_ColDivideBy(Mat<ZZ_p>& A, long colNum, ZZ_p div);

void GE_RowMultConstAndSub(Mat<ZZ_p>& A, long rowNumTo, long rowNumFrom, ZZ_p mult);

void GE_RowMultConstAndSub(Mat<ZZ>& A, long rowNumTo, long rowNumFrom, ZZ mult);

void GE_ColMultConstAndSub(Mat<ZZ_p>& A, long colNumTo, long colNumFrom, ZZ_p mult);

void GE_ColSwitch(Mat<ZZ_p>& A, long colNumLeft, long colNumRight);

void findKernel(Mat<ZZ_p>& K, Mat<ZZ_p>& A);

#endif