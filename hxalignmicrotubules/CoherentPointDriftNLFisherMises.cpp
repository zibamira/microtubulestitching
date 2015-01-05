#include <hxalignmicrotubules/CoherentPointDriftNLFisherMises.h>

#include <QString>

#include <mclib/McVec2d.h>
#include <mcla/f77_blas.h>
#include <mcla/f77_lapack.h>
#include <mclib/McWatch.h>

#ifdef _OPENMP
#include <omp.h>
#include <hxcore/HxSettingsMgr.h>
#endif

// I(spr) think the two LAPACK calls could be replaced by a single call to the
// driver routine SPOSV, since A is symmetric.
static void solve(const McDMatrix<double>& A, const McDMatrix<double>& B,
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

static void multMatrices(const McDMatrix<double>& A, const McDMatrix<double>& B,
                         const McDMatrix<double>& C, const double alpha,
                         const double beta, McDMatrix<double>& result) {
    result = C;
    mcassert((A.nRows() == C.nRows()) || (beta == 0));
    mcassert((B.nCols() == C.nCols()) || (beta == 0));
    mcassert(A.nCols() == B.nRows());
    if (beta == 0.0) {
        result.resize(B.nCols(), A.nRows());
        result.fill(0.0);
    } else
        result.transpose();
    int M = A.nRows();
    int N = B.nCols();
    int K = A.nCols();
    int LDA = A.nCols();
    int LDB = B.nCols();
    int LDC = result.nCols();

    f77_dgemm('T', 'T', M, N, K, alpha, A.dataPtr(), LDA, B.dataPtr(), LDB,
              beta, result.dataPtr(), LDC);
    result.transpose();
}

static void multThreeMatrices(const McDMatrix<double>& A,
                              const McDMatrix<double>& B,
                              const McDMatrix<double>& C,
                              McDMatrix<double>& result) {
    McDMatrix<double> tempResult(A.nRows(), B.nCols());
    McDMatrix<double> dummy;
    multMatrices(A, B, dummy, 1.0, 0.0, tempResult);
    result.resize(tempResult.nRows(), C.nCols());
    multMatrices(tempResult, C, dummy, 1.0, 0.0, result);
}

static double gauss(const McDVector<double>& mean, const double sigmaSquare,
                    const McDVector<double> x) {
    McDVector<double> distvec = x;
    distvec -= mean;
    const double dist = distvec.length2();
    return 1.0 /
           (pow((2.0 * sigmaSquare * M_PI), (double)(distvec.size()) / 2.0)) *
           exp(-1 * dist / (2.0 * sigmaSquare));
}

// cf Eq. 6 'Fisher-Mises distribution'.
static double fisherMises(const McDVector<double>& mean, const double kappa,
                          const McDVector<double> x) {
    const double dotProd = mean[0] * x[0] + mean[1] * x[1] + mean[2] * x[2];
    return kappa / (2.0 * M_PI * (exp(kappa) - exp(-1 * kappa))) *
           exp(kappa * dotProd);
}

static double trace(const McDMatrix<double>& mat) {
    double trace = 0.0;
    for (int i = 0; i < mat.nCols(); i++)
        trace += mat[i][i];
    return trace;
}

static void printStdout(const char* msg) {
    printf("%s", msg);
}

CoherentPointDriftNLFisherMises::CoherentPointDriftNLFisherMises(void) {
    mMeansAndStds.std = 1.0;
    params.maxIterations = 100;
    params.eDiffRelStop = 1.e-5;
    params.sigmaSquareStop = 1.e-7;
    params.useDirections = false;
    mPrint = &printStdout;
}

CoherentPointDriftNLFisherMises::~CoherentPointDriftNLFisherMises(void) {}

void CoherentPointDriftNLFisherMises::setPrint(print_t print) {
    mPrint = print;
}

void CoherentPointDriftNLFisherMises::print(QString msg) {
    mPrint(qPrintable(msg + "\n"));
}

// This implements Figure S3 'Algorithm for elastic transformation'.  The
// correspondence between variables here and in the paper should be obvious
// unless noted otherwise.
mtalign::AlignInfo
CoherentPointDriftNLFisherMises::align(McDMatrix<double>& G,
                                       McDMatrix<double>& W,
                                       McDArray<McVec2i>& correspondences) {
    const double beta = params.beta;
    const double lambda = params.lambda;
    const double w = params.w;
    const int maxIter = params.maxIterations;
    const double eDiffRelStop = params.eDiffRelStop;
    const double sigmaSquareStop = params.sigmaSquareStop;
    const bool useDirections = params.useDirections;

    McWatch watch;
    watch.start();
    print("Align CPD NL.");
    print(QString("beta: %1").arg(beta));
    print(QString("lambda: %1").arg(lambda));
    print(QString("w: %1").arg(w));
    print(QString("maxIter: %1").arg(maxIter));
    print(QString("sigmaSquareStop: %1").arg(sigmaSquareStop));
    print(QString("eDiffRelStop: %1").arg(eDiffRelStop));

    normalize();

    const int M = ys.nRows();
    const int N = xs.nRows();
    const int D = xs.nCols();
    print(QString("N: %1").arg(N));
    print(QString("M: %1").arg(M));
    print(QString("D: %1").arg(D));

    const McDMatrix<double> X = xs;
    const McDMatrix<double> Y = ys;

    print("Initializing W...");
    W.resize(M, D);
    W.fill(0.0);

    print("Computing initial sigma...");

    // D = 2 in the paper.
    double sigmaSquare = (1.0 / (D * M * N)) * sumSquaredDistances(xs, ys);

    // The paper states "Initialization: \kappa \ll 1.".  We observed kappa in
    // the range > 10-100 for convergence.  So initializing with 1.0 should be
    // ok, since it at least one order smaller than the converged kappa.
    double kappa = 1.0;

    print("..init G..");
    initializeG(G, ys);

    // Create P, will be filled below.
    McDMatrix<double> P;

    // E-step: Compute the P matrix, conditional P(x|m) and compute log
    // likelihood E for convergence check.
    double E;
    // ... shiftedYs is (y_m + G(m, .) * W) in the paper.
    McDMatrix<double> shiftedYs = ys;
    print("..shifting ys G..");
    shiftYs(ys, G, W, shiftedYs);
    print("..computing P..");
    // ... computeP expects shifted ys.
    computeP(xs, shiftedYs, sigmaSquare, kappa, w, P, E);

    double eDiffRel = FLT_MAX;
    int i;
    for (i = 0; (i < maxIter) && (sigmaSquare > sigmaSquareStop) &&
                    (eDiffRel > eDiffRelStop);
         i++) {
        ////////////////////
        // M-Step.

        // item Solve (G...) W = ...
        // Get diagonal of P with identity vector multiplied.
        McDMatrix<double> dP;
        print("..getdP..");
        getDP(P, dP);

        // Create the second term in the brackets.
        McDMatrix<double> lambdaSigmaEye(M, M);
        lambdaSigmaEye.makeIdentity();
        lambdaSigmaEye.multScal(lambda * sigmaSquare);

        // get W.
        //
        // Different from the manuscript, dP, which is diag(P1) in the
        // manuscript, is left multiplied.
        //
        // Note that the left-hand side is symmetric, positive definite, so we
        // could use Cholesky instead of general LU decomposition.
        print("..preparing matrices for solve..");
        McDMatrix<double> denominatorMatrix;
        multMatrices(dP, G, lambdaSigmaEye, 1.0, 1.0, denominatorMatrix);
        McDMatrix<double> numeratorPart, numerator;
        McDMatrix<double> emptyMatrix(0, 0);
        multMatrices(dP, Y, emptyMatrix, 1.0, 0.0, numeratorPart);
        multMatrices(P, X, numeratorPart, 1.0, -1.0, numerator);

        print("..solving..");
        solve(denominatorMatrix, numerator, W);

        // item N_P = ....
        // Sum Entries in P.
        double Np = 0.0;
        for (int k = 0; k < P.nRows(); k++)
            for (int j = 0; j < P.nCols(); j++)
                Np += P[k][j];

        // Compute T.
        McDMatrix<double> shiftedYs2(ys.nRows(), ys.nCols());
        shiftYs(ys, G, W, shiftedYs2);
        McDMatrix<double> T = shiftedYs;

        // item \sigma^2 = ...
        // Create first term in brackets.
        print("..compute new sigma..");
        McDMatrix<double> dPT;
        McDMatrix<double> temp = P;
        getDP(temp.transpose(), dPT);
        temp = X;
        McDMatrix<double> firstTermMatPart, firstTermMat;
        temp.transpose();
        multMatrices(temp, dPT, emptyMatrix, 1.0, 0.0, firstTermMatPart);
        multMatrices(firstTermMatPart, X, emptyMatrix, 1.0, 0.0, firstTermMat);

        // Second term.
        getDP(P, dP);
        McDMatrix<double> secondTermMatPart, secondTermMat;
        multMatrices(P, X, emptyMatrix, 1.0, 0.0, secondTermMatPart);
        secondTermMatPart.transpose();
        multMatrices(secondTermMatPart, T, emptyMatrix, 1.0, 0.0,
                     secondTermMat);

        // And third term.
        temp = T;
        temp.transpose();
        McDMatrix<double> thirdTermMat, thirdTermMatPart;
        multMatrices(temp, dP, emptyMatrix, 1.0, 0.0, thirdTermMatPart);
        multMatrices(thirdTermMatPart, T, emptyMatrix, 1.0, 0.0, thirdTermMat);

        // D = 2 in paper.
        sigmaSquare =
            1.0 / (Np * D) * (trace(firstTermMat) - 2.0 * trace(secondTermMat) +
                              trace(thirdTermMat));

        print(QString("Sigma: %1").arg(sigmaSquare));

        // item iterate \kappa = ...
        if (useDirections)
            computeKappa(P, Np, kappa);

        //////////
        // E-step: Compute the P matrix, conditional P(x|m) and compute log
        // likelihood for convergence check.

        // ... shiftedYs is (y_m + G(m, .) * W) in the paper.
        print("..shift ys G..");
        shiftYs(ys, G, W, shiftedYs);
        print("..compute P..");
        const double oldE = E;

        // ... computeP expects shifted ys.
        computeP(xs, shiftedYs, sigmaSquare, kappa, w, P, E);

        // For checking convergence on log likelihood.
        eDiffRel = (oldE - E) / fabs(oldE);

        print(QString("E: %1").arg(E));
        print(QString("EDiff: %1").arg(eDiffRel));
        print(QString("This was EM iteration: %1").arg(i));
    }

    mtalign::AlignInfo info;
    info.timeInSec = watch.stop();
    print(QString("This took %1 seconds.").arg(info.timeInSec));
    info.sigmaSquare = sigmaSquare;
    info.kappa = kappa;
    info.numIterations = i;
    info.eDiffRel = eDiffRel;
    info.e = E;
    return info;
}

void CoherentPointDriftNLFisherMises::convertDirectionsToMatrix(
    const McDArray<McVec3f>& directions, McDMatrix<double>& matrix) {
    matrix.resize(directions.size(), 3);
    for (int i = 0; i < directions.size(); i++) {
        McVec3f dir = directions[i];
        dir.normalize();
        matrix[i][0] = dir.x;
        matrix[i][1] = dir.y;
        matrix[i][2] = dir.z;
    }
}

void CoherentPointDriftNLFisherMises::convertCoordsToMatrix(
    const McDArray<McVec3f>& points, McDMatrix<double>& matrix) {
    matrix.resize(points.size(), 2);
    for (int i = 0; i < points.size(); i++) {
        matrix[i][0] = points[i].x;
        matrix[i][1] = points[i].y;
    }
}
void CoherentPointDriftNLFisherMises::convertMatrixToCoords(
    McDArray<McVec3f>& points, const McDMatrix<double>& matrix) {
    points.resize(matrix.nRows());
    for (int i = 0; i < points.size(); i++) {
        points[i].x = matrix[i][0];
        points[i].y = matrix[i][1];
    }
}

void CoherentPointDriftNLFisherMises::normalize() {
    // Compute mean.
    mMeansAndStds.mean = McDVector<double>(xs.nCols());
    mMeansAndStds.mean.fill(0.0);
    for (int i = 0; i < xs.nRows(); i++) {
        mMeansAndStds.mean += xs.getRowVector(i);
    }
    for (int i = 0; i < ys.nRows(); i++) {
        mMeansAndStds.mean += ys.getRowVector(i);
    }
    mMeansAndStds.mean /= (double)(xs.nRows() + ys.nRows());

    // Compute std dev, assuming a scalar.
    mMeansAndStds.std = 0.0;
    for (int i = 0; i < xs.nRows(); i++) {
        McDVector<double> theRow = xs.getRowVector(i);
        theRow -= mMeansAndStds.mean;
        xs.setRowVector(theRow, i);
        mMeansAndStds.std += xs.getRowVector(i).length2();
    }
    for (int i = 0; i < ys.nRows(); i++) {
        McDVector<double> theRow = ys.getRowVector(i);
        theRow -= mMeansAndStds.mean;
        ys.setRowVector(theRow, i);
        mMeansAndStds.std += ys.getRowVector(i).length2();
    }
    mMeansAndStds.std /= (double)(xs.nRows() + ys.nRows());

    // Normalize.
    mMeansAndStds.std = sqrt(mMeansAndStds.std);
    for (int i = 0; i < xs.nRows(); i++) {
        McDVector<double> cur = xs.getRowVector(i);
        cur *= 1.0 / mMeansAndStds.std;
        xs.setRowVector(cur, i);
    }
    for (int i = 0; i < ys.nRows(); i++) {
        McDVector<double> cur = ys.getRowVector(i);
        cur *= 1.0 / mMeansAndStds.std;
        ys.setRowVector(cur, i);
    }
}

bool CoherentPointDriftNLFisherMises::computeKappa(const McDMatrix<double>& P,
                                                   const double Np,
                                                   double& kappa) {
    print("computeKappa.");

    McDMatrix<double> XdT = xDirs;
    XdT.transpose();
    McDMatrix<double> PT = P;
    PT.transpose();
    McDMatrix<double> XdT_x_PT_x_Yd;
    multThreeMatrices(XdT, PT, yDirs, XdT_x_PT_x_Yd);

    const double c = trace(XdT_x_PT_x_Yd);
    const bool ret = newtonsMethodForKappa(kappa, c, Np);
    if (kappa < 1.e-6)
        kappa = 1.e-6;

    return ret;
}

bool CoherentPointDriftNLFisherMises::newtonsMethodForKappa(double& kappa,
                                                            const double c,
                                                            const double Np) {
    double diff = FLT_MAX;
    int counter = 0;
    if (c > Np - 1.e-4)
        return false;

    while (fabs(diff) > 1.e-5 && counter < 200) {
        diff = kappaFirstDerivative(kappa, c, Np) /
               kappaSecondDerivative(kappa, Np);
        print(QString("kappa: %1").arg(kappa));
        kappa = kappa - diff;
        print(QString("kappa-diff: %1").arg(kappa));
        counter++;
    }
    return kappa > 0.0;
}

double CoherentPointDriftNLFisherMises::kappaFirstDerivative(const double kappa,
                                                             const double c,
                                                             const double Np) {
    const double first = c + Np * (1.0 / kappa - 1.0 / tanh(kappa));
    print(QString("first: %1").arg(first));
    print(QString("c: %1").arg(c));
    return first;
}

double
CoherentPointDriftNLFisherMises::kappaSecondDerivative(const double kappa,
                                                       const double Np) {
    const double second =
        +Np * (-1.0 / (kappa * kappa) + pow(1.0 / sinh(kappa), 2.0));
    print(QString("second: %1").arg(second));
    return second;
}

void CoherentPointDriftNLFisherMises::computeP(
    const McDMatrix<double>& X, const McDMatrix<double>& Y,
    const double sigmaSquare, const double kappa, const double w,
    McDMatrix<double>& P, double& LL) {
    const bool useDirections = params.useDirections;
    const int M = Y.nRows();
    const int N = X.nRows();
    P.resize(M, N);
    McDVector<double> denominators;
    denominators.resize(N);
    denominators.fill(0.0);

    // Compute the P matrix, conditional P(x|m).

    // Here, the numerator and the first term in the denominator are both
    // multiplied by the Gauss factor (1 / (2 pi sigma^2)) and the Fisher Mises
    // factor (2 pi (exp(kappa) - exp(-kappa)) / kappa), while in the paper the
    // Fisher Mises factor is in the last term in the denominator.
#pragma omp parallel for
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            double numerator =
                gauss(Y.getRowVector(i), sigmaSquare, X.getRowVector(j));
            if (useDirections)
                numerator *= fisherMises(yDirs.getRowVector(i), kappa,
                                         xDirs.getRowVector(j));
            denominators[j] += numerator;
            P[i][j] = numerator;
        }
    }

    // Compute log likelihood for convergence check.
    LL = 0.0;
    for (int j = 0; j < N; j++) {
        LL -= log(1.0 / (double)M * (1.0 - w) * denominators[j] +
                  w * 1.0 / (double)N);
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            P[i][j] /=
                (denominators[j] + w / (1.0 - w) * (double)M / (double)N);
        }
    }
}

void CoherentPointDriftNLFisherMises::getDP(const McDMatrix<double>& P,
                                            McDMatrix<double>& dP) {
    dP.resize(P.nRows(), P.nRows());
    dP.fill(0.0);
    for (int i = 0; i < P.nRows(); i++) {
        double diagEntry = 0.0;
        for (int j = 0; j < P.nCols(); j++) {
            diagEntry += P[i][j];
        }
        dP[i][i] = diagEntry;
    }
}

void CoherentPointDriftNLFisherMises::shiftYs(const McDMatrix<double>& oldYs,
                                              const McDMatrix<double>& G,
                                              const McDMatrix<double>& W,
                                              McDMatrix<double>& shiftedYs) {
    shiftedYs = oldYs;
    for (int i = 0; i < oldYs.nRows(); i++) {
        McDVector<double> shift(oldYs.nCols());
        // multVec() from left.  Therefore transpose (different from paper).
        McDMatrix<double> WT = W;
        WT = WT.transpose();

        // G is symmetric MxM matrix, so we can take row or column.
        WT.multVec(G[i], shift);

        McDVector<double> newY = oldYs.getRowVector(i);
        newY += shift;

        shiftedYs.setRowVector(newY, i);
    }
}

void
CoherentPointDriftNLFisherMises::rescaleYs(McDMatrix<double>& oldYs) const {
    for (int i = 0; i < oldYs.nRows(); i++) {
        McDVector<double> newY = oldYs.getRowVector(i);
        newY *= mMeansAndStds.std;
        newY += mMeansAndStds.mean;
        oldYs.setRowVector(newY, i);
    }
}

void CoherentPointDriftNLFisherMises::initializeG(McDMatrix<double>& G,
                                                  const McDMatrix<double>& ps) {
    const double beta = params.beta;
    G.resize(ps.nRows(), ps.nRows());
    for (int i = 0; i < ps.nRows(); i++)
        for (int j = i; j < ps.nRows(); j++) {
            McDVector<double> row = ps.getRowVector(i);
            row -= ps.getRowVector(j);
            double dist = exp(-1.0 / (2.0 * beta * beta) * row.length2());
            G[i][j] = dist;
            G[j][i] = dist;
        }
}

double CoherentPointDriftNLFisherMises::sumSquaredDistances(
    const McDMatrix<double>& p1, const McDMatrix<double>& p2) {
    double sumOfSquaredDistances = 0.0;
    for (int i = 0; i < p1.nRows(); i++)
        for (int j = 0; j < p2.nRows(); j++) {
            McDVector<double> row = p1.getRowVector(i);
            row -= p2.getRowVector(j);
            sumOfSquaredDistances += (row).length2();
        }
    return sumOfSquaredDistances;
}
