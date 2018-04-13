#include "iostream"

#include "stdlib.h"

#include "NTL/ZZ_p.h"
#include <NTL/HNF.h>
#include <NTL/LLL.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>
#include <NTL/mat_ZZ.h>

#include "LWEInstance.h"
#include "GaussianElimination.h"

using namespace std;
using namespace NTL;



int main()
{
    srand((unsigned int)time(NULL));
    LWEInstance Inst(11, 30, 50, 2);
    Inst.print();
    // ZZ p = to_ZZ(11);
    // ZZ_p::init(p);
    // cout << "A = " << endl
    //      << Inst.A << endl;

    ZZ detA;

    // Mat<ZZ> AA = Inst.A * transpose(Inst.A);
    // cout << AA << endl;

    // determinant(det, Inst.A, 0);
    // cout << det << endl;

    // HNF(W, Inst.A, det);

    // cout << "HNF " << endl;
    // cout << W << endl;

    
    Mat<ZZ> AT;
    transpose(AT, Inst.A);
    Mat<ZZ> U;
    long rankA = image(detA, AT, U, 0);
    // cout << "rank = " << rankA << endl;
    // cout << "A = " << endl << Inst.A << endl;
    // cout << "U = " << endl
    //      << U << endl;
    
    Mat<ZZ> K;
    K.SetDims(U.NumRows() - rankA, U.NumCols());
    for(long i = 0 ; i < K.NumRows() ; i++){
        for(long j = 0 ; j < K.NumCols() ; j++){
            K[i][j] = U[i][j];
        }
    }

    // cout << "K = " << endl << K << endl;
    
    Mat<ZZ> Y;
    long tau = 10;
    Y.SetDims(K.NumCols(), tau);
    
    for(long i = 0 ; i < tau ; i++){
        ZZ detK;
        Mat<ZZ> UK;
        cout << "before" << endl;
        long ran_num = rand() % (K.NumRows() - 1) + 1;
        cout << "ran_num = " << ran_num << endl;
        GE_RowMultConstAndSub(K, 1, ran_num, to_ZZ(ran_num + rand() % 5));
        cout << "after" << endl;
        long rankK = LLL(detK, K, UK, 0);
        for(long j = 0; j < Y.NumRows(); j++){
            Y[j][i] = K[0][j];
            // cout << Y[j][i] << ' ';
        }
        // cout << endl;
    }

    cout << "Y = " << endl << Y << endl;
    Mat<ZZ> YT;
    transpose(YT, Y);
    cerr << "YT size = " << YT.NumRows() << " " << YT.NumCols() << endl;
    
    
    // Mat<ZZ> C = AT * YT;
    mat_ZZ C;
    mul(C, YT, AT);
    cout << "C = " << endl
        << C << endl;

    cerr << "C size = " << C.NumRows() << " " << C.NumCols() << endl;
    cerr << "s length = " << Inst.s.length() << endl;
    vec_ZZ ans;
    ans = C * Inst.s;
    cout << ans << endl;

    // LLL(det, W, Inst.A);

    // cout << "LLL" << endl;
    // cout << W << endl;

    // Vec<ZZ> w;
    // Vec<ZZ> zeroVector;
    // zeroVector.SetLength(10);
    // NearVector(w, Inst.A, Inst.b);

    // cout << w << endl;
}