#include <hxalignmicrotubules/mtalign/CPDLinearAligner.h>

#include <QString>

#include <mcla/internal/f77_blas.h>
#include <mcla/internal/f77_lapack.h>
#include <mclib/McVec2.h>
#include <mclib/internal/McWatch.h>

#include <coin/IpIpoptApplication.hpp>
#include <hxalignmicrotubules/mtalign/IPOPTForCPD.h>
#include <hxalignmicrotubules/mtalign/IPOPTForCPDLinearAligner.h>
#include <hxalignmicrotubules/mtalign/IPOPTForCPDLinearAlignerWithoutScale.h>
#include <hxalignmicrotubules/mtalign/data.h>
#include <hxalignmicrotubules/mtalign/rotation.h>
#include <hxalignmicrotubules/mtalign/math.h>

/****************************************************************************/
// This class implements Agorithm 1 and 2 (Box 1,2) of the Paper "Automated
// stitching of microtubule networks in serial sections of electron tomograms"
// For Algorithm 1, set
// useDirections=true;
// usePositions=false;

// For Algorithm 2, set
// useDirections=true;
// usePositions=true;

/****************************************************************************/

namespace ma = mtalign;
namespace ip = Ipopt;

ma::CPDLinearAligner::CPDLinearAligner() {
    mMeansAndStds.stdC = 1.0;
    // switch this off, if you do not like the data to be scaled
    params.withScaling = true;
    params.maxIterations = 100;
    params.eDiffRelStop = 1.e-5;
    params.sigmaSquareStop = 1.e-5;
    params.useDirections = true;
    params.usePositions = true;
    mContext = &defaultContext();
}

void ma::CPDLinearAligner::print(QString msg) {
    mContext->print(msg);
}

void ma::CPDLinearAligner::setContext(ma::Context* ctx) {
    mContext = ctx;
}

static void convertMcDArraysToMatrix(const McDArray<McVec3f>& coords,
                                     const McDArray<McVec3f>& directions,
                                     McDMatrix<double>& CMat,
                                     McDMatrix<double>& DMat) {
    mcassert(coords.size() == directions.size());
    CMat.resize(coords.size(), 2);
    DMat.resize(coords.size(), 3);

    for (int i = 0; i < coords.size(); i++) {
        CMat[i][0] = coords[i].x;
        CMat[i][1] = coords[i].y;
        McVec3f normalized = directions[i];
        normalized.normalize();
        DMat[i][0] = normalized.x;
        DMat[i][1] = normalized.y;
        DMat[i][2] = normalized.z;
    }
}

void ma::CPDLinearAligner::setPoints(const mtalign::FacingPointSets& points) {
    convertMcDArraysToMatrix(points.ref.positions, points.ref.directions, Xc,
                             Xd);
    convertMcDArraysToMatrix(points.trans.positions, points.trans.directions,
                             Yc, Yd);
}

static double sumSquaredDistances(const McDMatrix<double>& p1,
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

static void shiftYDirections(const McDMatrix<double>& oldYd,
                             const McDMatrix<double>& Rd,
                             McDMatrix<double>& shiftedYd) {
    shiftedYd = oldYd;
    for (int i = 0; i < oldYd.nRows(); i++) {
        McDVector<double> shiftedDNew(3);
        McDVector<double> curOldYd = oldYd.getRowVector(i);
        mcassert(Rd.nCols() == 2);
        McDVector<double> yOldxy(2), rotatedYdxy(2);
        yOldxy[0] = curOldYd[0];
        yOldxy[1] = curOldYd[1];
        Rd.multVec(yOldxy, rotatedYdxy);
        shiftedDNew[0] = rotatedYdxy[0];
        shiftedDNew[1] = rotatedYdxy[1];
        shiftedDNew[2] = curOldYd[2];
        shiftedYd.setRowVector(shiftedDNew, i);
    }
}

mtalign::AlignInfo ma::CPDLinearAligner::align(McDMatrix<double>& Rc, double& s,
                                               McDVector<double>& t,
                                               McDMatrix<double>& Rd) {

    const int maxIter = params.maxIterations;
    const double w = params.w;
    const bool withScaling = params.withScaling;
    const double sigmaTol = params.sigmaSquareStop;
    const double eDiffTol = params.eDiffRelStop;
    const bool usePositions = params.usePositions;
    const bool useDirections = params.useDirections;

    double kappa, sigmaSquare;

    McWatch watch;
    watch.start();

    normalize();

    const int M = Yc.nRows();
    const double Md = M;
    const int N = Xc.nRows();
    print(QString("N: %1").arg(N));
    print(QString("M: %1").arg(M));
    print(QString("w: %1").arg(w));
    const double Nd = N;
    const int Dc = Xc.nCols();

    // Initialize T
    t.resize(Dc);
    t.fill(0.0);
    // initialize sigma
    print("Compute initial sigma...");
    sigmaSquare = (1.0 / (2.0 * Md * Nd)) * sumSquaredDistances(Xc, Yc);

    print(QString("Initial sigma^2 is %1").arg(sigmaSquare));
    s = 1.0;
    Rc.resize(Dc, Dc);
    Rc.makeIdentity();
    Rd.resize(3, 3);
    Rd.makeIdentity();
    // initialize P, just create it
    McDMatrix<double> P;
    P.resize(M, N);

    McDMatrix<double> shiftedYc = Yc;
    McDMatrix<double> shiftedYd = Yd;
    print(QString("maxIter is: %1").arg(maxIter));
    double E;  //=computeLogLikelihood(Rc,s,t,sigmaSquare,Rd,kappa);
    kappa = 1.e-4;
    print(QString("Initial kappa is %1").arg(kappa));

    computeP(Xc, Xd, shiftedYc, shiftedYd, sigmaSquare, kappa, w, P, E);
    print(QString("E: %1").arg(E));
    double EDiffPerc = FLT_MAX;

    int i;

    for (i = 0; (i < maxIter) && (sigmaSquare > sigmaTol) &&
                    (EDiffPerc > eDiffTol) && (kappa < 600.0);
         i++) {

        /*******************************************************************************/
        // M-Step
        // compute Np
        double NP = getNP(P);
        print(QString("...NP is %1").arg(NP));
        // compue muX and muY
        McDMatrix<double> PT = P;
        PT.transpose();
        McDMatrix<double> muX, muY;
        getMu(Xc, PT, NP, muX);
        getMu(Yc, P, NP, muY);
        McDMatrix<double> XcHat, YcHat;
        getHat(Xc, muX, XcHat);
        print(QString("muX is %1 %2").arg(muX[0][0]).arg(muX[1][0]));
        getHat(Yc, muY, YcHat);
        print(QString("muY is %1 %2").arg(muY[0][0]).arg(muY[1][0]));
        // Initialize all the long expressions
        print("...initializing terms...");
        CPDLinearAlignerTerms terms(XcHat, Xd, YcHat, Yd, P);
        if (usePositions && !useDirections) {
            print("...computing R...");
            terms.computeRc(Rc);
            print(QString("Rc is : %1 %2 %3 %4")
                      .arg(Rc[0][0])
                      .arg(Rc[0][1])
                      .arg(Rc[1][0])
                      .arg(Rc[1][1]));
            if (withScaling)
                terms.computeS(Rc, s);

            computeT(muX, s, Rc, muY, t);
            sigmaSquare = terms.computeSigmaSquare(Rc, s, NP);
        } else if (!usePositions && useDirections) {
            print("...computing R...");
            terms.computeRd2d(Rd);
            print(QString("Rd is : %1 %2 %3 %4")
                      .arg(Rd[0][0])
                      .arg(Rd[0][1])
                      .arg(Rd[1][0])
                      .arg(Rd[1][1]));
            terms.computeKappa(Rd, NP, kappa, mContext->print);
            Rc[0][0] = Rd[0][0];
            Rc[1][0] = Rd[1][0];
            Rc[0][1] = Rd[0][1];
            Rc[1][1] = Rd[1][1];
            // we compute t so that the result looks OK, but is will not
            // influence the outcome
            computeT(muX, s, Rc, muY, t);
        } else {
            print("..optimize parameters with IPOpt...");
            terms.optimizeParameters(Rc, s, sigmaSquare, kappa, NP,
                                     params.withScaling, mContext->print);
            computeT(muX, s, Rc, muY, t);
        }
        print(QString("...new s is %1 ...").arg(s));
        print(QString("...new t is %1 %2 ...").arg(t[0]).arg(t[1]));
        print(QString("...new sigma^2 is %1 ...").arg(sigmaSquare));
        print(QString("...new kappa is %1 ...").arg(kappa));
        print(QString("...new rho is %1 ...").arg(ma::rotationAngle2d(Rc)));
        /*******************************************************************************/
        // perform E-step
        // 1. Compute new ys
        print(QString("..shift ys .."));

        // compute the P matrix, conditional P(x|m)
        print(QString("..compute P.."));
        double oldE = E;
        shiftYCoords(Yc, Rc, s, t, shiftedYc);
        shiftYDirections(Yd, Rc, shiftedYd);

        computeP(Xc, Xd, shiftedYc, shiftedYd, sigmaSquare, kappa, w, P, E);
        print(QString("E: %1").arg(E));
        EDiffPerc = (oldE - E) / fabs(oldE);
        print(QString("EDiff: %1").arg(EDiffPerc));
        if (sigmaSquare > sigmaTol && kappa < 300)
            mcassert(EDiffPerc > -1.e-6);
        print(QString("This was EM iteration %1").arg(i));
        print("");
        if (i == maxIter - 1)
            print("EM max number of iterations reached.");
    }

    mtalign::AlignInfo info;
    info.timeInSec = watch.stop();
    print(QString("This took %1 seconds.").arg(info.timeInSec));
    info.sigmaSquare = sigmaSquare;
    info.kappa = kappa;
    info.numIterations = i;
    info.eDiffRel = EDiffPerc;
    info.e = E;
    return info;
}

void ma::CPDLinearAligner::computeT(const McDMatrix<double>& muX,
                                    const double s, const McDMatrix<double>& R,
                                    const McDMatrix<double>& muY,
                                    McDVector<double>& t) {
    McDVector<double> rotatedMuY(2);
    R.multVec(muY.getColVector(0), rotatedMuY);
    print(QString("...rotatedMuY is %1 %2")
              .arg(rotatedMuY[0])
              .arg(rotatedMuY[1]));
    rotatedMuY *= s;
    print(QString("...rotatedMuY*s is %1 %2")
              .arg(rotatedMuY[0])
              .arg(rotatedMuY[1]));
    print(QString("...muX.getColVector(0) is %1 %2")
              .arg(muX.getColVector(0)[0])
              .arg(muX.getColVector(0)[1]));
    t = muX.getColVector(0);
    t -= rotatedMuY;
}

void ma::CPDLinearAligner::getHat(const McDMatrix<double>& X,
                                  const McDMatrix<double>& muX,
                                  McDMatrix<double>& XHat) {

    McDMatrix<double> muXT = muX;
    muXT.transpose();
    XHat.resize(X.nRows(), X.nCols());

    McDMatrix<double> ones(X.nRows(), 1);
    ones.fill(1.0);
    McDMatrix<double> dummyC;
    ma::multMatrices(ones, muXT, dummyC, 1.0, 0.0, XHat);
    XHat = X - XHat;
}

void ma::CPDLinearAligner::getMu(const McDMatrix<double>& X,
                                 const McDMatrix<double>& P, const double NP,
                                 McDMatrix<double>& muX) {
    McDMatrix<double> XT = X;
    XT.transpose();
    McDMatrix<double> ones(P.nCols(), 1);
    ones.fill(1.0);
    ma::multThreeMatrices(XT, P, ones, muX);
    muX.multScal(1.0 / NP);
}

double ma::CPDLinearAligner::getNP(const McDMatrix<double>& P) {

    // compute Np -> Kahan sum Entries in P
    double Np = 0.0;
    double c = 0.0;
    for (int i = 0; i < P.nRows(); i++) {
        for (int j = 0; j < P.nCols(); j++) {
            double y = P[i][j] - c;
            double t = Np + y;
            c = (t - Np) - y;
            Np = t;
        }
    }
    return Np;
}

void ma::CPDLinearAligner::normalize() {
    // get Mean
    mMeansAndStds.meanC = McDVector<double>(Xc.nCols());
    mMeansAndStds.meanC.fill(0.0);

    for (int i = 0; i < Xc.nRows(); i++) {
        mMeansAndStds.meanC += Xc.getRowVector(i);
    }
    for (int i = 0; i < Yc.nRows(); i++) {
        mMeansAndStds.meanC += Yc.getRowVector(i);
    }
    mMeansAndStds.meanC /= (double)(Xc.nRows() + Yc.nRows());
    // get std, we assume a scalar
    mMeansAndStds.stdC = 0.0;
    // we do not want to shift direction with the mean

    for (int i = 0; i < Xc.nRows(); i++) {
        McDVector<double> tmpVec = Xc.getRowVector(i);
        tmpVec -= mMeansAndStds.meanC;
        Xc.setRowVector(tmpVec, i);
        mMeansAndStds.stdC += tmpVec.length2();
    }
    for (int i = 0; i < Yc.nRows(); i++) {
        McDVector<double> tmpVec = Yc.getRowVector(i);
        tmpVec -= mMeansAndStds.meanC;
        Yc.setRowVector(tmpVec, i);
        mMeansAndStds.stdC += tmpVec.length2();
    }
    mMeansAndStds.stdC /= (double)(Xc.nRows() + Yc.nRows());

    mMeansAndStds.stdC = sqrt(mMeansAndStds.stdC);
    for (int i = 0; i < Xc.nRows(); i++) {
        McDVector<double> cur = Xc.getRowVector(i);
        cur *= 1.0 / mMeansAndStds.stdC;
        Xc.setRowVector(cur, i);
    }
    for (int i = 0; i < Yc.nRows(); i++) {
        McDVector<double> cur = Yc.getRowVector(i);
        cur *= 1.0 / mMeansAndStds.stdC;
        Yc.setRowVector(cur, i);
    }
    print(QString("MeanC is: %1 %2")
              .arg(mMeansAndStds.meanC[0])
              .arg(mMeansAndStds.meanC[1]));
    print(QString("stdC is: %1").arg(mMeansAndStds.stdC));
}

void ma::CPDLinearAligner::computeP(const McDMatrix<double>& Xc,
                                    const McDMatrix<double>& Xd,
                                    const McDMatrix<double>& Yc,
                                    const McDMatrix<double>& Yd,
                                    const double sigmaSquare,
                                    const double kappa, const double w,
                                    McDMatrix<double>& P, double& LL) {
    const bool useDirections = params.useDirections;
    const bool usePositions = params.usePositions;

    int M = Yc.nRows();
    int N = Xc.nRows();
    P.resize(M, N);
    McDVector<double> denominators;
    denominators.resize(N);
    denominators.fill(0.0);
    // compute the P matrix, conditional P(x|m)
    LL = 0.0;

    for (int j = 0; j < N; j++) {
        // Kahan sum the denominator
        double sum = 0.0;
        double c = 0.0;
        for (int i = 0; i < M; i++) {
            // compute numerator
            double numerator = 1.0;
            if (usePositions)
                numerator *= ma::gauss(Yc.getRowVector(i), sigmaSquare,
                                       Xc.getRowVector(j));
            if (useDirections)
                numerator *= ma::fisherMises(Yd.getRowVector(i), kappa,
                                             Xd.getRowVector(j));
            double y = numerator - c;
            double t = sum + y;
            c = t - sum - y;
            sum = t;

            P[i][j] = numerator;
        }
        denominators[j] = sum;
    }
    for (int j = 0; j < N; j++) {
        LL -= log(1.0 / ((double)M) * (1.0 - w) * (denominators[j]) +
                  w * 1.0 / (double)N);
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            P[i][j] /=
                (denominators[j] + w / (1.0 - w) * (double)M / (double)N);
        }
    }
}

void ma::CPDLinearAligner::shiftYCoords(const McDMatrix<double>& oldYc,
                                        const McDMatrix<double>& Rc,
                                        const double& s,
                                        const McDVector<double>& t,
                                        McDMatrix<double>& shiftedYc) {
    shiftedYc = oldYc;
    for (int i = 0; i < oldYc.nRows(); i++) {
        McDVector<double> shiftedCNew(2);
        McDVector<double> curOldYc = oldYc.getRowVector(i);
        Rc.multVec(curOldYc, shiftedCNew);
        shiftedCNew *= s;
        shiftedCNew += t;
        shiftedYc.setRowVector(shiftedCNew, i);
    }
}

McVec3f ma::CPDLinearAligner::warpPoint(const McVec3f& point,
                                        const McDMatrix<double>& Rc,
                                        const double s,
                                        const McDVector<double>& t) const {

    // scale
    McDVector<double> pointcStdMean(2);
    pointcStdMean.dataPtr()[0] = point.x;
    pointcStdMean.dataPtr()[1] = point.y;
    pointcStdMean.dataPtr()[0] -= mMeansAndStds.meanC[0];
    pointcStdMean.dataPtr()[1] -= mMeansAndStds.meanC[1];
    pointcStdMean /= mMeansAndStds.stdC;

    // transform
    McDVector<double> pointcStdMeanShifted(pointcStdMean.size());
    Rc.multVec(pointcStdMean, pointcStdMeanShifted);
    pointcStdMeanShifted *= s;
    pointcStdMeanShifted += t;

    // rescale
    pointcStdMeanShifted *= mMeansAndStds.stdC;
    pointcStdMeanShifted.dataPtr()[0] += mMeansAndStds.meanC[0];
    pointcStdMeanShifted.dataPtr()[1] += mMeansAndStds.meanC[1];

    return McVec3f(pointcStdMeanShifted[0], pointcStdMeanShifted[1], point.z);
}

McMat4f ma::CPDLinearAligner::getTransformMat4f(McDMatrix<double> R, double s,
                                                McDVector<double> t) {
    McMat4f mat = McMat4f::IDENTITY;
    t *= mMeansAndStds.stdC;
    t[0] += mMeansAndStds.meanC[0];
    t[1] += mMeansAndStds.meanC[1];
    McDVector<double> muRotScaled(2), mu(2);
    mu[0] = mMeansAndStds.meanC[0];
    mu[1] = mMeansAndStds.meanC[1];
    R.multVec(mu, muRotScaled);
    muRotScaled *= s;
    t -= muRotScaled;
    R.transpose();
    mat[0][0] = R[0][0] * s;
    mat[0][1] = R[0][1] * s;
    mat[1][0] = R[1][0] * s;
    mat[1][1] = R[1][1] * s;
    mat[3][0] = t[0];
    mat[3][1] = t[1];
    return mat;
}

static void getDP(const McDMatrix<double>& P, McDMatrix<double>& dP) {
    dP.resize(P.nRows(), P.nRows());
    dP.fill(0.0);
    for (int i = 0; i < P.nRows(); i++) {
        // Kahan sum
        double sum = 0.0;
        double c = 0.0;
        for (int j = 0; j < P.nCols(); j++) {
            double y = P[i][j] - c;
            double t = sum + y;
            c = t - sum - y;
            sum = t;
        }
        dP[i][i] = sum;
    }
}

ma::CPDLinearAlignerTerms::CPDLinearAlignerTerms(const McDMatrix<double>& XcHat,
                                                 const McDMatrix<double>& Xd,
                                                 const McDMatrix<double>& YcHat,
                                                 const McDMatrix<double>& Yd,
                                                 const McDMatrix<double>& P) {
    McDMatrix<double> XcHatT = XcHat;
    XcHatT.transpose();
    McDMatrix<double> YcHatT = YcHat;
    YcHatT.transpose();
    McDMatrix<double> PT = P;
    PT.transpose();
    // get the strange diagonal matrix dP1
    McDMatrix<double> dP1, dPT1;
    getDP(P, dP1);
    getDP(PT, dPT1);
    ma::multThreeMatrices(XcHatT, dPT1, XcHat, XcHatT_x_dPT1_x_XcHat);
    ma::multThreeMatrices(YcHatT, dP1, YcHat, YcHatT_x_dP1_x_YcHat);
    ma::multThreeMatrices(XcHatT, PT, YcHat, XcHatT_x_PT_x_YcHat);

    McDMatrix<double> XdT = Xd;
    XdT.transpose();
    ma::multThreeMatrices(XdT, PT, Yd, XdT_x_PT_x_Yd);
}

void ma::CPDLinearAlignerTerms::computeS(const McDMatrix<double>& R,
                                         double& s) const {
    McDMatrix<double> numeratorMatrix = XcHatT_x_PT_x_YcHat;
    numeratorMatrix.transpose();
    McDMatrix<double> numeratorMatrixRotated, dummyMat;
    ma::multMatrices(numeratorMatrix, R, dummyMat, 1.0, 0.0,
                     numeratorMatrixRotated);
    double numerator = ma::trace(numeratorMatrixRotated);
    McDMatrix<double> denominatorMatrix = YcHatT_x_dP1_x_YcHat;
    double denominator = ma::trace(denominatorMatrix);
    s = numerator / denominator;
}

void ma::CPDLinearAlignerTerms::computeRc(McDMatrix<double>& Rc) const {
    McDMatrix<double> U, V, SVD, dummyMat(2, 2);

    SVD = XcHatT_x_PT_x_YcHat;
    double dummyd[2];
    SVD.SVD(U, &dummyd[0], V);
    McDMatrix<double> VT = V;
    VT.transpose();

    McDMatrix<double> lastDiagEntryMat;
    ma::multMatrices(U, VT, dummyMat, 1.0, 0.0, lastDiagEntryMat);
    double lastDiagEntry = lastDiagEntryMat.det();
    McDMatrix<double> C(2, 2);
    C.makeIdentity();
    C[1][1] = lastDiagEntry;
    ma::multThreeMatrices(U, C, VT, Rc);
}

void ma::CPDLinearAlignerTerms::computeRd2d(McDMatrix<double>& Rd) const {
    McDMatrix<double> U, V, SVD(2, 2), dummyMat(2, 2), R2d;

    SVD[0][0] = XdT_x_PT_x_Yd[0][0];
    SVD[1][0] = XdT_x_PT_x_Yd[1][0];
    SVD[0][1] = XdT_x_PT_x_Yd[0][1];
    SVD[1][1] = XdT_x_PT_x_Yd[1][1];
    double dummyd[2];
    SVD.SVD(U, &dummyd[0], V);
    McDMatrix<double> VT = V;
    VT.transpose();

    McDMatrix<double> lastDiagEntryMat;
    ma::multMatrices(U, VT, dummyMat, 1.0, 0.0, lastDiagEntryMat);
    double lastDiagEntry = lastDiagEntryMat.det();
    McDMatrix<double> C(2, 2);
    C.makeIdentity();
    C[1][1] = lastDiagEntry;
    ma::multThreeMatrices(U, C, VT, R2d);

    Rd.makeIdentity();
    Rd[0][0] = R2d[0][0];
    Rd[1][0] = R2d[1][0];
    Rd[0][1] = R2d[0][1];
    Rd[1][1] = R2d[1][1];
}

double ma::CPDLinearAlignerTerms::computeSigmaSquare(
    const McDMatrix<double>& Rc, const double s, const double NP) const {
    float D = 2.0;

    // compute first term
    McDMatrix<double> firstTermMat = XcHatT_x_dPT1_x_XcHat;

    McDMatrix<double> secondTermMat = XcHatT_x_PT_x_YcHat;
    secondTermMat.transpose();
    McDMatrix<double> secondTermMatRotated, dummyMat;
    ma::multMatrices(secondTermMat, Rc, dummyMat, 1.0, 0.0,
                     secondTermMatRotated);
    double firstTerm = ma::trace(firstTermMat);
    double secondTerm = ma::trace(secondTermMatRotated);

    return 1.0 / (NP * D) * (firstTerm - s * secondTerm);
}

static double kappaFirstDerivative(const double kappa, const double c,
                                   const double Np) {
    return c + Np * (1.0 / kappa - 1.0 / tanh(kappa));
}

static double kappaSecondDerivative(const double kappa, const double Np) {
    return +Np * (-1.0 / (kappa * kappa) + pow(1.0 / sinh(kappa), 2.0));
}

static bool newtonsMethodForKappa(double& kappa, const double c,
                                  const double Np,
                                  ma::Context::print_t& print) {
    double diff = FLT_MAX;
    int counter = 0;
    if (c > Np - 1.e-4)
        return false;
    while (fabs(diff) > 1.e-10 && counter < 200) {
        const double first = kappaFirstDerivative(kappa, c, Np);
        const double second = kappaSecondDerivative(kappa, Np);
        diff = first / second;
        print(QString("first: %1").arg(first));
        print(QString("c: %1").arg(c));
        print(QString("second: %1").arg(second));
        print(QString("kappa: %1").arg(kappa));
        kappa = kappa - diff;
        print(QString("kappa-diff: %1").arg(kappa));
        counter++;
    }
    return true;
}

bool
ma::CPDLinearAlignerTerms::computeKappa(const McDMatrix<double>& Rd,
                                        const double Np, double& kappa,
                                        ma::Context::print_t& print) const {
    print("computeKappa:");
    // compute first term
    McDMatrix<double> cMat = XdT_x_PT_x_Yd;
    cMat.transpose();
    McDMatrix<double> cMatRotated, dummyMat;
    ma::multMatrices(cMat, Rd, dummyMat, 1.0, 0.0, cMatRotated);
    double c = ma::trace(cMatRotated);
    return newtonsMethodForKappa(kappa, c, Np, print);
}

bool ma::CPDLinearAlignerTerms::optimizeParameters(
    McDMatrix<double>& R, double& s, double& sigmaSquare, double& kappa,
    const double Np, const bool withScaling,
    ma::Context::print_t& print) const {

    CPDLinearAlignerDerivatives func;
    func.Np = Np;
    func.A = XcHatT_x_dPT1_x_XcHat;
    func.BT = XcHatT_x_PT_x_YcHat;
    func.BT.transpose();
    func.C = YcHatT_x_dP1_x_YcHat;
    func.DT.resize(2, 2);
    func.DT[0][0] = XdT_x_PT_x_Yd[0][0];
    func.DT[1][0] = XdT_x_PT_x_Yd[1][0];
    func.DT[0][1] = XdT_x_PT_x_Yd[0][1];
    func.DT[1][1] = XdT_x_PT_x_Yd[1][1];
    func.DT.transpose();
    func.E.resize(1, 1);
    func.E[0][0] = XdT_x_PT_x_Yd[2][2];

    double rho = mtalign::rotationAngle2d(R);

    IPOPTForCPD* optimizer;
    if (withScaling) {
        optimizer = createIPOPTForCPDLinearAligner();
    } else {
        optimizer = createIPOPTForCPDLinearAlignerWithoutScale();
    }
    optimizer->gradAndHessAndFunc = func;
    optimizer->s = s;
    McDVector<double> startVector(4);
    startVector[0] = s;

    optimizer->rho = rho;
    startVector[1] = rho;
    optimizer->kappa = kappa;
    startVector[2] = kappa;
    startVector[3] = sigmaSquare;
    func.initCurParams(startVector);
    print(QString("sigma^2 before init for opt: %1").arg(sigmaSquare));
    func.initStartValues(startVector);
    print(QString("sigma^2 after init for opt: %1").arg(startVector[3]));
    optimizer->resultValues.resize(4);
    optimizer->sigmaSquare = startVector[3];
    double* result = optimizer->resultValues.dataPtr();

    // Create a new instance of your nlp
    //  (use a SmartPtr, not raw)
    ip::SmartPtr<ip::TNLP> mynlp(optimizer);

    // Create a new instance of IpoptApplication
    //  (use a SmartPtr, not raw)
    // We are using the factory, since this allows us to compile this
    // example with an Ipopt Windows DLL
    ip::SmartPtr<ip::IpoptApplication> app = IpoptApplicationFactory();

    // Change some options
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("derivative_test", "second-order");
    app->Options()->SetIntegerValue("max_iter", 1000);
    app->Options()->SetNumericValue("max_cpu_time", 30);

    // It can be useful to increase the `print_level` to 4 or higher for
    // debugging.
    app->Options()->SetIntegerValue("print_level", 0);

    // Intialize the IpoptApplication and process the options
    ip::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != ip::Solve_Succeeded) {
        print("*** Error during initialization!");
        print(QString("*** Errorcode is %1").arg(status));
        return (int)status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == ip::Solve_Succeeded) {
        print("*** The problem solved.");
    } else {
        print(
            QString("*** The problem did not solve.  Code is %1").arg(status));
    }

    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.

    s = result[0];
    sigmaSquare = result[3];
    kappa = result[2];
    R[0][0] = cos(result[1]);
    R[0][1] = -sin(result[1]);
    R[1][0] = sin(result[1]);
    R[1][1] = cos(result[1]);
    return 1;
}
