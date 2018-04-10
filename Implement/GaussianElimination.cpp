#include "GaussianElimination.h"

void GE_ColDivideBy(Mat<ZZ_p>& A, long colNum, ZZ_p div)
{
    for (long i = 0; i < A.NumRows(); i++) {
        A[i][colNum - 1] /= div;
    }
}

void GE_RowMultConstAndSub(Mat<ZZ_p>& A, long rowNumTo, long rowNumFrom, ZZ_p mult)
{
    for (long i = 0; i < A.NumCols(); i++) {
        A[rowNumTo - 1][i] -= A[rowNumFrom - 1][i] * mult;
    }
}

void GE_RowMultConstAndSub(Mat<ZZ>& A, long rowNumTo, long rowNumFrom, ZZ mult)
{
    for (long i = 0; i < A.NumCols(); i++) {
        A[rowNumTo - 1][i] -= A[rowNumFrom - 1][i] * mult;
    }
}

void GE_ColMultConstAndSub(Mat<ZZ_p>& A, long colNumTo, long colNumFrom, ZZ_p mult)
{
    for (long i = 0; i < A.NumRows(); i++) {
        A[i][colNumTo - 1] += A[i][colNumFrom - 1] * mult;
    }
}

void GE_ColSwitch(Mat<ZZ_p>& A, long colNumLeft, long colNumRight)
{
    for (long i = 0; i < A.NumRows(); i++) {
        ZZ_p dummy = A[i][colNumLeft - 1];
        A[i][colNumLeft - 1] = A[i][colNumRight - 1];
        A[i][colNumRight - 1] = dummy;
    }
}

void findKernel(Mat<ZZ_p>& K, Mat<ZZ_p>& A)
{
    long m = A.NumRows();
    long n = A.NumCols();
    long rank = 0;

    Mat<ZZ_p> AI;
    AI.SetDims(m + n, n);
    for (long j = 0; j < n; j++) {
        for (long i = 0; i < m; i++) {
            AI[i][j] = A[i][j];
        }
    }
    for (long i = 0; i < n; i++) {
        AI[m + i][i] = 1;
    }

    long currentRow = 1;
    while (1) {
        long currentCol = currentRow;
        while (currentCol <= n && AI(currentRow, currentCol) == 0) {
            currentCol += 1;
        }

        if (currentCol > n) {
            break;
        }
        if (currentCol != rank + 1) {
            GE_ColSwitch(AI, rank + 1, currentCol);
            currentCol = rank + 1;
        }

        GE_ColDivideBy(AI, currentCol, AI(currentRow, currentCol));

        for (long i = currentRow + 1; i <= m + n; i++) {
            GE_RowMultConstAndSub(AI, i, currentRow, AI(i, currentCol));
        }

        for (long j = currentCol + 1; j <= n; j++) {
            AI(currentRow, j) = 0;
            // if (AI(currentRow, j) != 0) {
            // GE_ColMultConstAndSub(AI, j, currentCol, AI(currentRow, j));
            // }
        }
        currentRow++;
        rank++;
        if (currentRow > m) {
            break;
        }
    }
    K.SetDims(n, n - rank);
    for (long i = 0; i < K.NumRows(); i++) {
        for (long j = 0; j < K.NumCols(); j++) {
            K[i][j] = AI[m + i][rank + j];
        }
    }
}