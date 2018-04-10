
#include "LWEInstance.h"

#include "iostream"

void LWEInstance::print()
{
    cout << "A = " << endl
         << A << endl
         << "s = " << endl
         << s << endl
         << "e = " << endl
         << e << endl
         << "b = " << endl
         << b << endl;
}

LWEInstance::LWEInstance(long q, long m_, long n_, long h)
{
    srand((unsigned int)time(NULL));
    ZZ_p::init(to_ZZ(q));
    m = m_;
    n = n_;

    generateA(q, m, n); // gen A
    generateSecret(n, h); // gen s, HW = h
    generateError(m); // gen e
    generateB(q, m); // gen b = As + e (mod q)
}

void LWEInstance::generateA(long q, long m, long n)
{
    srand((unsigned int)time(NULL));
    A.SetDims(m, n);
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < n; j++) {
            A[i][j] = rand() % q - q / 2;
        }
    }
}

void LWEInstance::generateSecret(long dim, long HW)
{
    s.SetLength(dim);
    long numOfNonzero = 0;
    while (numOfNonzero < HW) {
        long loc = rand() % dim;
        if (s[loc] == 0) {
            if (rand() % 2 == 0) {
                s[loc] = 1;
            } else {
                s[loc] = -1;
            }
            numOfNonzero++;
        }
    }
}

void LWEInstance::generateError(long dim)
{
    e.SetLength(dim);
    for (long i = 0; i < dim; i += 1) {
        e[i] = (rand() % 3) - 1;
    }
}
void LWEInstance::generateB(long q, long dim)
{
    b.SetLength(m);
    b = A * s + e;
    for (long i = 0; i < m; i++) {
        b[i] = b[i] % q - q / 2;
    }
}