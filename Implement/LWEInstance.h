#ifndef LWEInstance_H
#define LWEInstance_H

#include "NTL/ZZ.h"
#include <NTL/matrix.h>
#include <NTL/LLL.h>
#include <NTL/vec_vec_ZZ.h>

#include "stdlib.h"

using namespace std;
using namespace NTL;

class LWEInstance
{
  public:
    long m;
    long n;
    Mat<ZZ> A;
    Vec<ZZ> b;
    Vec<ZZ> s;
    Vec<ZZ> e;

    void print();

    LWEInstance(long q, long m_, long n_, long h);

    void generateA(long q, long m, long n);
    void generateSecret(long dim, long HW);
    void generateError(long dim);
    void generateB(long q, long dim);
};
#endif