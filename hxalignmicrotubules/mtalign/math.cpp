#include <hxalignmicrotubules/mtalign/math.h>

#include <mclib/McAssert.h>
#include <mcla/f77_blas.h>
#include <mcla/f77_lapack.h>

namespace ma = mtalign;

void ma::multMatrices(const McDMatrix<double>& A, const McDMatrix<double>& B,
                      const McDMatrix<double>& C, const double alpha,
                      const double beta, McDMatrix<double>& result) {
    result = C;
    mcassert((A.nRows() == C.nRows()) || (beta == 0));
    mcassert((B.nCols() == C.nCols()) || (beta == 0));
    mcassert(A.nCols() == B.nRows());
    if (beta == 0.0) {
        result.resize(B.nCols(), A.nRows());
        result.fill(0.0);
    } else {
        result.transpose();
    }
    const int M = A.nRows();
    const int N = B.nCols();
    const int K = A.nCols();
    const int LDA = A.nCols();
    const int LDB = B.nCols();
    const int LDC = result.nCols();
    f77_dgemm('T', 'T', M, N, K, alpha, A.dataPtr(), LDA, B.dataPtr(), LDB,
              beta, result.dataPtr(), LDC);
    result.transpose();
}

void ma::multThreeMatrices(const McDMatrix<double>& A,
                           const McDMatrix<double>& B,
                           const McDMatrix<double>& C,
                           McDMatrix<double>& result) {
    McDMatrix<double> tempResult(A.nRows(), B.nCols());
    McDMatrix<double> dummy;
    multMatrices(A, B, dummy, 1.0, 0.0, tempResult);
    result.resize(tempResult.nRows(), C.nCols());
    multMatrices(tempResult, C, dummy, 1.0, 0.0, result);
}

double ma::trace(const McDMatrix<double>& mat) {
    double trace = 0.0;
    for (int i = 0; i < mat.nCols(); i++)
        trace += mat[i][i];
    return trace;
}

// I(spr) think the two LAPACK calls could be replaced by a single call to the
// driver routine SPOSV, since A is symmetric.
void ma::solve(const McDMatrix<double>& A, const McDMatrix<double>& B,
               McDMatrix<double>& X) {
    // Require square matrices.
    mcassert(A.nCols() == A.nRows());
    mcassert(A.nCols() != 0);

    // LU factorization of A.
    McDMatrix<double> AT = A;
    AT.transpose();
    int M = A.nRows();
    int N = A.nCols();
    int* piv = new int[M];
    int info;
    f77_dgetrf(M, N, AT.dataPtr(), M, piv, info);
    mcassert(info == 0);

    // Solve.
    X = B;
    X.transpose();
    N = A.nCols();
    int LDA = A.nRows();
    int LDB = B.nRows();
    int NRHS = B.nCols();
    f77_dgetrs('N', N, NRHS, AT.dataPtr(), LDA, piv, X.dataPtr(), LDB, info);
    mcassert(info == 0);
    X.transpose();
    delete[] piv;
}

double ma::gauss(const McDVector<double>& mean, const double sigmaSquare,
                 const McDVector<double> x) {
    McDVector<double> distvec = x;
    distvec -= mean;
    const double dist = distvec.length2();
    return 1.0 /
           (pow((2.0 * sigmaSquare * M_PI), (double)(distvec.size()) / 2.0)) *
           exp(-1 * dist / (2.0 * sigmaSquare));
}

// cf Eq. 6 'Fisher-Mises distribution'.
double ma::fisherMises(const McDVector<double>& mean, const double kappa,
                       const McDVector<double> x) {
    const double dotProd = mean[0] * x[0] + mean[1] * x[1] + mean[2] * x[2];
    return kappa / (2.0 * M_PI * (exp(kappa) - exp(-1 * kappa))) *
           exp(kappa * dotProd);
}
