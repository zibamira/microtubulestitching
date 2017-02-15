#include <hxalignmicrotubules/mtalign/CPDElasticAligner.h>

#include <QString>

#include <mcla/internal/f77_blas.h>
#include <mcla/internal/f77_lapack.h>
#include <mclib/internal/McWatch.h>

#ifdef _OPENMP
#include <omp.h>
#include <hxcore/HxSettingsMgr.h>
#endif

#include <hxalignmicrotubules/mtalign/math.h>

namespace ma = mtalign;

ma::CPDElasticAligner::CPDElasticAligner() {
    mMeansAndStds.std = 1.0;
    params.maxIterations = 100;
    params.eDiffRelStop = 1.e-5;
    params.sigmaSquareStop = 1.e-7;
    params.useDirections = false;
    mContext = &defaultContext();
}

void ma::CPDElasticAligner::print(QString msg) {
    mContext->print(msg);
}

void ma::CPDElasticAligner::setContext(ma::Context* ctx) {
    mContext = ctx;
}

static void convertDirectionsToMatrix(const McDArray<McVec3f>& directions,
                                      McDMatrix<double>& matrix) {
    matrix.resize(directions.size(), 3);
    for (int i = 0; i < directions.size(); i++) {
        McVec3f dir = directions[i];
        dir.normalize();
        matrix[i][0] = dir.x;
        matrix[i][1] = dir.y;
        matrix[i][2] = dir.z;
    }
}

void ma::CPDElasticAligner::setPoints(const FacingPointSets& points) {
    mcrequire(!params.useDirections ||
              points.ref.directions.size() == points.ref.positions.size());
    mcrequire(!params.useDirections ||
              points.trans.directions.size() == points.trans.positions.size());
    convertCoordsToMatrix(points.ref.positions, xs);
    convertDirectionsToMatrix(points.ref.directions, xDirs);
    convertCoordsToMatrix(points.trans.positions, ys);
    convertDirectionsToMatrix(points.trans.directions, yDirs);
    normalize();
}

static void convertMatrixToCoords(McDArray<McVec3f>& points,
                                  const McDMatrix<double>& matrix) {
    points.resize(matrix.nRows());
    for (int i = 0; i < points.size(); i++) {
        points[i].x = matrix[i][0];
        points[i].y = matrix[i][1];
        points[i].z = 0;
    }
}

McDArray<McVec3f> ma::CPDElasticAligner::align(AlignInfo& info) {
    McDMatrix<double> G;
    McDMatrix<double> W;
    info = align(G, W);
    McDMatrix<double> transCoordsShiftedM;
    shiftYs(ys, G, W, transCoordsShiftedM);
    rescaleYs(transCoordsShiftedM);
    McDArray<McVec3f> transCoords;
    convertMatrixToCoords(transCoords, transCoordsShiftedM);
    return transCoords;
}

static McDMatrix<double> makeG(const McDMatrix<double>& ps, const double beta) {
    McDMatrix<double> G;
    G.resize(ps.nRows(), ps.nRows());
    for (int i = 0; i < ps.nRows(); i++) {
        for (int j = i; j < ps.nRows(); j++) {
            McDVector<double> row = ps.getRowVector(i);
            row -= ps.getRowVector(j);
            const double dist = exp(-1.0 / (2.0 * beta * beta) * row.length2());
            G[i][j] = dist;
            G[j][i] = dist;
        }
    }
    return G;
}

static McDMatrix<double> getDP(const McDMatrix<double>& P) {
    McDMatrix<double> dP;
    dP.resize(P.nRows(), P.nRows());
    dP.fill(0.0);
    for (int i = 0; i < P.nRows(); i++) {
        double diagEntry = 0.0;
        for (int j = 0; j < P.nCols(); j++) {
            diagEntry += P[i][j];
        }
        dP[i][i] = diagEntry;
    }
    return dP;
}

// This implements Figure S3 'Algorithm for elastic transformation'.  The
// correspondence between variables here and in the paper should be obvious
// unless noted otherwise.
ma::AlignInfo ma::CPDElasticAligner::align(McDMatrix<double>& G,
                                           McDMatrix<double>& W) {
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
    G = makeG(ys, params.beta);

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
        dP = getDP(P);

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
        dPT = getDP(temp.transpose());
        temp = X;
        McDMatrix<double> firstTermMatPart, firstTermMat;
        temp.transpose();
        multMatrices(temp, dPT, emptyMatrix, 1.0, 0.0, firstTermMatPart);
        multMatrices(firstTermMatPart, X, emptyMatrix, 1.0, 0.0, firstTermMat);

        // Second term.
        dP = getDP(P);
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

    AlignInfo info;
    info.timeInSec = watch.stop();
    print(QString("This took %1 seconds.").arg(info.timeInSec));
    info.sigmaSquare = sigmaSquare;
    info.kappa = kappa;
    info.numIterations = i;
    info.eDiffRel = eDiffRel;
    info.e = E;
    return info;
}

void
ma::CPDElasticAligner::convertCoordsToMatrix(const McDArray<McVec3f>& points,
                                             McDMatrix<double>& matrix) {
    matrix.resize(points.size(), 2);
    for (int i = 0; i < points.size(); i++) {
        matrix[i][0] = points[i].x;
        matrix[i][1] = points[i].y;
    }
}

void ma::CPDElasticAligner::normalize() {
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

bool ma::CPDElasticAligner::computeKappa(const McDMatrix<double>& P,
                                         const double Np, double& kappa) {
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

bool ma::CPDElasticAligner::newtonsMethodForKappa(double& kappa, const double c,
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

double ma::CPDElasticAligner::kappaFirstDerivative(const double kappa,
                                                   const double c,
                                                   const double Np) {
    const double first = c + Np * (1.0 / kappa - 1.0 / tanh(kappa));
    print(QString("first: %1").arg(first));
    print(QString("c: %1").arg(c));
    return first;
}

double ma::CPDElasticAligner::kappaSecondDerivative(const double kappa,
                                                    const double Np) {
    const double second =
        +Np * (-1.0 / (kappa * kappa) + pow(1.0 / sinh(kappa), 2.0));
    print(QString("second: %1").arg(second));
    return second;
}

void ma::CPDElasticAligner::computeP(const McDMatrix<double>& X,
                                     const McDMatrix<double>& Y,
                                     const double sigmaSquare,
                                     const double kappa, const double w,
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

void ma::CPDElasticAligner::shiftYs(const McDMatrix<double>& oldYs,
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

void ma::CPDElasticAligner::rescaleYs(McDMatrix<double>& oldYs) const {
    for (int i = 0; i < oldYs.nRows(); i++) {
        McDVector<double> newY = oldYs.getRowVector(i);
        newY *= mMeansAndStds.std;
        newY += mMeansAndStds.mean;
        oldYs.setRowVector(newY, i);
    }
}

double ma::CPDElasticAligner::sumSquaredDistances(const McDMatrix<double>& p1,
                                                  const McDMatrix<double>& p2) {
    double sumOfSquaredDistances = 0.0;
    for (int i = 0; i < p1.nRows(); i++)
        for (int j = 0; j < p2.nRows(); j++) {
            McDVector<double> row = p1.getRowVector(i);
            row -= p2.getRowVector(j);
            sumOfSquaredDistances += (row).length2();
        }
    return sumOfSquaredDistances;
}
