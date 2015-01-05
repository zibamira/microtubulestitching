#include <hxalignmicrotubules/CoherentPointDriftRigidFisherMises.h>

#include <mcla/f77_blas.h>
#include <mcla/f77_lapack.h>
#include <mclib/McVec2d.h>
#include <mclib/McWatch.h>

#include <coin/IpIpoptApplication.hpp>
#include <hxalignmicrotubules/IPOPTForCPD.h>
#include <hxalignmicrotubules/IPOPTForRigidAlign.h>
#include <hxalignmicrotubules/IPOPTForRigidAlignWithoutScale.h>
#include <hxalignmicrotubules/mtalign/rotation.h>

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

namespace ip = Ipopt;

CoherentPointDriftRigidFisherMises::CoherentPointDriftRigidFisherMises() {
    mMeansAndStds.stdC = 1.0;
    // switch this off, if you do not like the data to be scaled
    params.withScaling = true;
    params.maxIterations = 100;
    params.eDiffRelStop = 1.e-5;
    params.sigmaSquareStop = 1.e-5;
    params.useDirections = true;
    params.usePositions = true;
}

CoherentPointDriftRigidFisherMises::~CoherentPointDriftRigidFisherMises() {}

mtalign::AlignInfo CoherentPointDriftRigidFisherMises::align(
    McDMatrix<double>& Rc, double& s, McDVector<double>& t,
    McDMatrix<double>& Rd, McDArray<McVec2i>& correspondences) {
    double kappa, sigmaSquare;
    return align(Rc, s, t, Rd, correspondences, sigmaSquare, kappa);
}

mtalign::AlignInfo CoherentPointDriftRigidFisherMises::align(
    McDMatrix<double>& Rc, double& s, McDVector<double>& t,
    McDMatrix<double>& Rd, McDArray<McVec2i>& correspondences,
    double& sigmaSquare, double& kappa) {
    const int maxIter = params.maxIterations;
    const double w = params.w;
    const bool withScaling = params.withScaling;
    const double sigmaTol = params.sigmaSquareStop;
    const double eDiffTol = params.eDiffRelStop;
    const bool usePositions = params.usePositions;
    const bool useDirections = params.useDirections;

    McWatch watch;
    watch.start();

    normalize();

    const int M = Yc.nRows();
    const double Md = M;
    const int N = Xc.nRows();
    std::cout << "\nN= " << N;
    std::cout << "\nM= " << M;
    std::cout << "\nw is:" << w;
    const double Nd = N;
    const int Dc = Xc.nCols();

    // Initialize T
    t.resize(Dc);
    t.fill(0.0);
    // initialize sigma
    std::cout << "\nCompute initial sigma...";
    sigmaSquare = (1.0 / (2.0 * Md * Nd)) * sumSquaredDistances(Xc, Yc);

    std::cout << "\nInitial sigma^2 is" << sigmaSquare;
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
    std::cout << "\n maxIter is: " << maxIter;
    double E;  //=computeLogLikelihood(Rc,s,t,sigmaSquare,Rd,kappa);
    allTerms terms;
    kappa = 1.e-4;
    std::cout << "\nInitial kappa is" << kappa;

    computeP(Xc, Xd, shiftedYc, shiftedYd, sigmaSquare, kappa, w, P, E);
    std::cout << "\nE: " << E << "\n";
    double EDiffPerc = FLT_MAX;

    int i;

    for (i = 0; (i < maxIter) && (sigmaSquare > sigmaTol) &&
                    (EDiffPerc > eDiffTol) && (kappa < 600.0);
         i++) {

        /*******************************************************************************/
        // M-Step
        // compute Np
        double NP = getNP(P);
        std::cout << "...NP is" << NP;
        // compue muX and muY
        McDMatrix<double> PT = P;
        PT.transpose();
        McDMatrix<double> muX, muY;
        getMu(Xc, PT, NP, muX);
        getMu(Yc, P, NP, muY);
        McDMatrix<double> XcHat, YcHat;
        getHat(Xc, muX, XcHat);
        std::cout << "\nmuX is " << muX[0][0] << " " << muX[1][0];
        getHat(Yc, muY, YcHat);
        std::cout << "\nmuY is " << muY[0][0] << " " << muY[1][0];
        // Initialize all the long expressions
        initTerms(terms, XcHat, Xd, YcHat, Yd, P);
        if (usePositions && !useDirections) {
            computeRc(terms, Rc);
            if (withScaling)
                computeS(terms, Rc, s);

            computeT(muX, s, Rc, muY, t);
            sigmaSquare = computeSigmaSquare(terms, Rc, s, NP);
        } else if (!usePositions && useDirections) {
            computeRd2d(terms, Rd);
            computeKappa(terms, Rd, NP, kappa);
            Rc[0][0] = Rd[0][0];
            Rc[1][0] = Rd[1][0];
            Rc[0][1] = Rd[0][1];
            Rc[1][1] = Rd[1][1];
            // we compute t so that the result looks OK, but is will not
            // influence the outcome
            computeT(muX, s, Rc, muY, t);
        } else {
            std::cout << "..optimize parameters with IPOpt...";
            optimizeParameters(terms, Rc, s, sigmaSquare, kappa, NP);
            computeT(muX, s, Rc, muY, t);
        }
        std::cout << "...new s is " << s << "...";
        std::cout << "...new t is " << t[0] << " " << t[1] << "...";
        std::cout << "...new sigma^2 is " << sigmaSquare << "...";
        std::cout << "...new kappa is " << kappa << "...";
        std::cout << "...new rho is " << mtalign::rotationAngle2d(Rc) << "...";
        /*******************************************************************************/
        // perform E-step
        // 1. Compute new ys
        std::cout << "\n..shift ys ..";

        // compute the P matrix, conditional P(x|m)
        std::cout << "\n..compute P..";
        double oldE = E;
        shiftYCoords(Yc, Rc, s, t, shiftedYc);
        shiftYDirections(Yd, Rc, shiftedYd);

        computeP(Xc, Xd, shiftedYc, shiftedYd, sigmaSquare, kappa, w, P, E);
        std::cout << "\nE: " << E << "\n";
        EDiffPerc = (oldE - E) / fabs(oldE);
        std::cout << "EDiff: " << EDiffPerc << "\n";
        if (sigmaSquare > sigmaTol && kappa < 300)
            mcassert(EDiffPerc > -1.e-6);
        std::cout << "\n This was EM iteration " << i;
        if (i == maxIter - 1)
            std::cout << "\n EM max number of iterations reached.";
    }

    mtalign::AlignInfo info;
    info.timeInSec = watch.stop();
    std::cout << "\nThis took " << info.timeInSec << " seconds.\n";
    info.sigmaSquare = sigmaSquare;
    info.kappa = kappa;
    info.numIterations = i;
    info.eDiffRel = EDiffPerc;
    info.e = E;
    return info;
}

bool CoherentPointDriftRigidFisherMises::optimizeParameters(
    const allTerms& terms, McDMatrix<double>& R, double& s, double& sigmaSquare,
    double& kappa, const double Np) {
    const bool withScaling = params.withScaling;

    DerivativesForRigidRegistration func;
    func.Np = Np;
    func.A = terms.XcHatT_x_dPT1_x_XcHat;
    func.BT = terms.XcHatT_x_PT_x_YcHat;
    func.BT.transpose();
    func.C = terms.YcHatT_x_dP1_x_YcHat;
    func.DT.resize(2, 2);
    func.DT[0][0] = terms.XdT_x_PT_x_Yd[0][0];
    func.DT[1][0] = terms.XdT_x_PT_x_Yd[1][0];
    func.DT[0][1] = terms.XdT_x_PT_x_Yd[0][1];
    func.DT[1][1] = terms.XdT_x_PT_x_Yd[1][1];
    func.DT.transpose();
    func.E.resize(1, 1);
    func.E[0][0] = terms.XdT_x_PT_x_Yd[2][2];

    double rho = mtalign::rotationAngle2d(R);

    IPOPTForCPD* optimizer;
    if (withScaling) {
        optimizer = createIPOPTForRigidAlign();
    } else {
        optimizer = createIPOPTForRigidAlignWithoutScale();
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
    std::cout << "\n sigma^2 before init for opt: " << sigmaSquare;
    func.initStartValues(startVector);
    std::cout << "\n sigma^2 after init for opt: " << startVector[3];
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

    app->Options()->SetIntegerValue("print_level", 4);

    // Intialize the IpoptApplication and process the options
    ip::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != ip::Solve_Succeeded) {
        std::cout << std::endl << std::endl
                  << "*** Error during initialization!" << std::endl
                  << "Errorcode is " << status << std::endl;
        return (int)status;
    }

    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);

    if (status == ip::Solve_Succeeded) {
        std::cout << std::endl << std::endl << "*** The problem solved."
                  << std::endl;
    } else {
        std::cout << std::endl << std::endl << "*** The problem did not solve."
                  << std::endl << "Code is " << status << std::endl;
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

bool CoherentPointDriftRigidFisherMises::computeKappa(
    const allTerms& terms, const McDMatrix<double>& Rd, const double Np,
    double& kappa) {
    std::cout << "\n computeKappa:";

    // compute first term
    McDMatrix<double> cMat = terms.XdT_x_PT_x_Yd;
    cMat.transpose();
    McDMatrix<double> cMatRotated, dummyMat;
    multMatrices(cMat, Rd, dummyMat, 1.0, 0.0, cMatRotated);
    double c = trace(cMatRotated);
    return newtonsMethodForKappa(kappa, c, Np);
}

bool CoherentPointDriftRigidFisherMises::newtonsMethodForKappa(
    double& kappa, const double c, const double Np) {
    double diff = FLT_MAX;
    int counter = 0;
    if (c > Np - 1.e-4)
        return false;
    while (fabs(diff) > 1.e-10 && counter < 200) {
        diff = kappaFirstDerivative(kappa, c, Np) /
               kappaSecondDerivative(kappa, Np);
        std::cout << "kappa:" << kappa;
        kappa = kappa - diff;
        std::cout << "kappa-diff:" << kappa;
        counter++;
    }
    return true;
}

double CoherentPointDriftRigidFisherMises::kappaFirstDerivative(
    const double kappa, const double c, const double Np) {
    std::cout << "\n first:" << c + Np*(1.0 / kappa - 1.0 / tan(kappa));
    std::cout << "\n c is:" << c;
    return c + Np * (1.0 / kappa - 1.0 / tanh(kappa));
}

double
CoherentPointDriftRigidFisherMises::kappaSecondDerivative(const double kappa,
                                                          const double Np) {
    std::cout << "\n second:"
              << +Np*(-1.0 / (kappa * kappa) + pow(1.0 / sinh(kappa), 2.0));
    std::cout << "\n cur kappa: " << kappa;
    return +Np * (-1.0 / (kappa * kappa) + pow(1.0 / sinh(kappa), 2.0));
}

double CoherentPointDriftRigidFisherMises::computeSigmaSquare(
    const allTerms& terms, const McDMatrix<double>& Rc, const double s,
    const double NP) {
    float D = 2.0;

    // compute first term
    McDMatrix<double> firstTermMat = terms.XcHatT_x_dPT1_x_XcHat;

    McDMatrix<double> secondTermMat = terms.XcHatT_x_PT_x_YcHat;
    secondTermMat.transpose();
    McDMatrix<double> secondTermMatRotated, dummyMat;
    multMatrices(secondTermMat, Rc, dummyMat, 1.0, 0.0, secondTermMatRotated);
    double firstTerm = trace(firstTermMat);
    double secondTerm = trace(secondTermMatRotated);

    return 1.0 / (NP * D) * (firstTerm - s * secondTerm);
}

void CoherentPointDriftRigidFisherMises::computeT(const McDMatrix<double>& muX,
                                                  const double s,
                                                  const McDMatrix<double>& R,
                                                  const McDMatrix<double>& muY,
                                                  McDVector<double>& t) {
    McDVector<double> rotatedMuY(2);
    R.multVec(muY.getColVector(0), rotatedMuY);
    std::cout << "...rotatedMuY is " << rotatedMuY[0] << " " << rotatedMuY[1];
    rotatedMuY *= s;
    std::cout << "...rotatedMuY*s is " << rotatedMuY[0] << " " << rotatedMuY[1];
    std::cout << "...muX.getColVector(0) is " << muX.getColVector(0)[0] << " "
              << muX.getColVector(0)[1];
    t = muX.getColVector(0);
    t -= rotatedMuY;
}

void CoherentPointDriftRigidFisherMises::computeS(const allTerms& terms,
                                                  const McDMatrix<double>& R,
                                                  double& s) {
    McDMatrix<double> numeratorMatrix = terms.XcHatT_x_PT_x_YcHat;
    numeratorMatrix.transpose();
    McDMatrix<double> numeratorMatrixRotated, dummyMat;
    multMatrices(numeratorMatrix, R, dummyMat, 1.0, 0.0,
                 numeratorMatrixRotated);
    double numerator = trace(numeratorMatrixRotated);
    std::cout << "\n numerator for s: " << numerator;

    McDMatrix<double> denominatorMatrix = terms.YcHatT_x_dP1_x_YcHat;
    double denominator = trace(denominatorMatrix);
    std::cout << "\n denominator for s: " << denominator;
    s = numerator / denominator;
}

void CoherentPointDriftRigidFisherMises::initTerms(
    allTerms& terms, const McDMatrix<double>& XcHat,
    const McDMatrix<double>& Xd, const McDMatrix<double>& YcHat,
    const McDMatrix<double>& Yd, const McDMatrix<double>& P) {
    std::cout << "...initializing terms...";
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
    multThreeMatrices(XcHatT, dPT1, XcHat, terms.XcHatT_x_dPT1_x_XcHat);
    multThreeMatrices(YcHatT, dP1, YcHat, terms.YcHatT_x_dP1_x_YcHat);
    multThreeMatrices(XcHatT, PT, YcHat, terms.XcHatT_x_PT_x_YcHat);

    McDMatrix<double> XdT = Xd;
    XdT.transpose();
    multThreeMatrices(XdT, PT, Yd, terms.XdT_x_PT_x_Yd);
}

void CoherentPointDriftRigidFisherMises::getHat(const McDMatrix<double>& X,
                                                const McDMatrix<double>& muX,
                                                McDMatrix<double>& XHat) {

    McDMatrix<double> muXT = muX;
    muXT.transpose();
    XHat.resize(X.nRows(), X.nCols());

    McDMatrix<double> ones(X.nRows(), 1);
    ones.fill(1.0);
    McDMatrix<double> dummyC;
    multMatrices(ones, muXT, dummyC, 1.0, 0.0, XHat);
    XHat = X - XHat;
}

void CoherentPointDriftRigidFisherMises::getMu(const McDMatrix<double>& X,
                                               const McDMatrix<double>& P,
                                               const double NP,
                                               McDMatrix<double>& muX) {
    McDMatrix<double> XT = X;
    XT.transpose();
    McDMatrix<double> ones(P.nCols(), 1);
    ones.fill(1.0);
    multThreeMatrices(XT, P, ones, muX);
    muX.multScal(1.0 / NP);
}

void CoherentPointDriftRigidFisherMises::computeRc(const allTerms& terms,
                                                   McDMatrix<double>& Rc) {
    std::cout << "\n...computing R...";
    McDMatrix<double> U, V, SVD, dummyMat(2, 2);

    SVD = terms.XcHatT_x_PT_x_YcHat;
    double dummyd[2];
    SVD.SVD(U, &dummyd[0], V);
    McDMatrix<double> VT = V;
    VT.transpose();

    McDMatrix<double> lastDiagEntryMat;
    multMatrices(U, VT, dummyMat, 1.0, 0.0, lastDiagEntryMat);
    double lastDiagEntry = lastDiagEntryMat.det();
    McDMatrix<double> C(2, 2);
    C.makeIdentity();
    C[1][1] = lastDiagEntry;
    multThreeMatrices(U, C, VT, Rc);
    std::cout << "\nRc is :\n" << Rc[0][0] << " " << Rc[0][1] << "\n"
              << Rc[1][0] << " " << Rc[1][1] << "\n";
}

void CoherentPointDriftRigidFisherMises::computeRd2d(const allTerms& terms,
                                                     McDMatrix<double>& Rd) {
    std::cout << "\n...computing R...";
    McDMatrix<double> U, V, SVD(2, 2), dummyMat(2, 2), R2d;

    SVD[0][0] = terms.XdT_x_PT_x_Yd[0][0];
    SVD[1][0] = terms.XdT_x_PT_x_Yd[1][0];
    SVD[0][1] = terms.XdT_x_PT_x_Yd[0][1];
    SVD[1][1] = terms.XdT_x_PT_x_Yd[1][1];
    double dummyd[2];
    SVD.SVD(U, &dummyd[0], V);
    McDMatrix<double> VT = V;
    VT.transpose();

    McDMatrix<double> lastDiagEntryMat;
    std::cout << "one\n";
    multMatrices(U, VT, dummyMat, 1.0, 0.0, lastDiagEntryMat);
    double lastDiagEntry = lastDiagEntryMat.det();
    McDMatrix<double> C(2, 2);
    C.makeIdentity();
    C[1][1] = lastDiagEntry;
    std::cout << "two\n";
    multThreeMatrices(U, C, VT, R2d);

    Rd.makeIdentity();
    Rd[0][0] = R2d[0][0];
    Rd[1][0] = R2d[1][0];
    Rd[0][1] = R2d[0][1];
    Rd[1][1] = R2d[1][1];

    std::cout << "\nRd is :\n" << Rd[0][0] << " " << Rd[0][1] << "\n"
              << Rd[1][0] << " " << Rd[1][1] << "\n";
}
void CoherentPointDriftRigidFisherMises::computeRd(const allTerms& terms,
                                                   McDMatrix<double>& Rd) {
    std::cout << "\n...computing R...";
    McDMatrix<double> U, V, SVD, dummyMat(3, 3);

    SVD = terms.XdT_x_PT_x_Yd;
    double dummyd[3];
    SVD.SVD(U, &dummyd[0], V);
    McDMatrix<double> VT = V;
    VT.transpose();

    McDMatrix<double> lastDiagEntryMat;
    std::cout << "one\n";
    multMatrices(U, VT, dummyMat, 1.0, 0.0, lastDiagEntryMat);
    double lastDiagEntry = lastDiagEntryMat.det();
    McDMatrix<double> C(3, 3);
    C.makeIdentity();
    C[2][2] = lastDiagEntry;
    std::cout << "two\n";
    multThreeMatrices(U, C, VT, Rd);
    std::cout << "\nRd is :\n" << Rd[0][0] << " " << Rd[0][1] << "\n"
              << Rd[1][0] << " " << Rd[1][1] << "\n";
}
void CoherentPointDriftRigidFisherMises::multThreeMatrices(
    const McDMatrix<double>& A, const McDMatrix<double>& B,
    const McDMatrix<double>& C, McDMatrix<double>& result) {
    McDMatrix<double> tempResult(A.nRows(), B.nCols());
    McDMatrix<double> dummy;
    multMatrices(A, B, dummy, 1.0, 0.0, tempResult);
    result.resize(tempResult.nRows(), C.nCols());
    multMatrices(tempResult, C, dummy, 1.0, 0.0, result);
}

double CoherentPointDriftRigidFisherMises::getNP(const McDMatrix<double>& P) {

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

void CoherentPointDriftRigidFisherMises::convertMcDArraysToMatrix(
    const McDArray<McVec3f>& coords, const McDArray<McVec3f>& directions,
    McDMatrix<double>& CMat, McDMatrix<double>& DMat) {
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

void CoherentPointDriftRigidFisherMises::normalize() {
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
    std::cout << "\nMeanC is: " << mMeansAndStds.meanC[0] << " "
              << mMeansAndStds.meanC[1];
    std::cout << "\nstdC is: " << mMeansAndStds.stdC << " ";
}

void CoherentPointDriftRigidFisherMises::solve(const McDMatrix<double>& A,
                                               const McDMatrix<double>& B,
                                               McDMatrix<double>& X)

{
    // can only do square matrices
    mcassert(A.nCols() == A.nRows());
    mcassert(A.nCols() != 0);

    // LU factorization of A;
    McDMatrix<double> AT = A;
    AT.transpose();
    int M = A.nRows();
    int N = A.nCols();
    int* piv = new int[M];
    int info;
    std::cout << "..LU decomp..";
    f77_dgetrf(M, N, AT.dataPtr(), M, piv, info);
    mcassert(info == 0);

    // solve
    X = B;
    X.transpose();
    N = A.nCols();
    int LDA = A.nRows();
    int LDB = B.nRows();
    int NRHS = B.nCols();
    std::cout << "..compute result..";
    f77_dgetrs('N', N, NRHS, AT.dataPtr(), LDA, piv, X.dataPtr(), LDB, info);
    mcassert(info == 0);
    X.transpose();
    delete[] piv;
}

void CoherentPointDriftRigidFisherMises::multMatrices(
    const McDMatrix<double>& A, const McDMatrix<double>& B,
    const McDMatrix<double>& C, const double alpha, const double beta,
    McDMatrix<double>& result) {
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

double CoherentPointDriftRigidFisherMises::gauss(const McDVector<double>& mean,
                                                 const double sigmaSquare,
                                                 const McDVector<double> x) {

    McDVector<double> distvec = x;
    distvec -= mean;
    double dist = distvec.length2();
    return 1.0 /
           (pow((2.0 * sigmaSquare * M_PI), (double)(distvec.size()) / 2.0)) *
           exp(-1 * dist / (2.0 * sigmaSquare));
}

double
CoherentPointDriftRigidFisherMises::fisherMises(const McDVector<double>& mean,
                                                const double kappa,
                                                const McDVector<double> x) {

    double dotProd = mean[0] * x[0] + mean[1] * x[1] + mean[2] * x[2];

    return kappa / (4.0 * M_PI * (exp(kappa) - exp(-1 * kappa))) *
           exp(kappa * dotProd);
}

void CoherentPointDriftRigidFisherMises::computeP(
    const McDMatrix<double>& Xc, const McDMatrix<double>& Xd,
    const McDMatrix<double>& Yc, const McDMatrix<double>& Yd,
    const double sigmaSquare, const double kappa, const double w,
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
                numerator *=
                    gauss(Yc.getRowVector(i), sigmaSquare, Xc.getRowVector(j));
            if (useDirections)
                numerator *=
                    fisherMises(Yd.getRowVector(i), kappa, Xd.getRowVector(j));
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

double CoherentPointDriftRigidFisherMises::trace(const McDMatrix<double>& mat) {
    double trace = 0.0;
    for (int i = 0; i < mat.nCols(); i++)
        trace += mat[i][i];
    return trace;
}

void CoherentPointDriftRigidFisherMises::getDP(const McDMatrix<double>& P,
                                               McDMatrix<double>& dP) {
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

void CoherentPointDriftRigidFisherMises::shiftYDirections(
    const McDMatrix<double>& oldYd, const McDMatrix<double>& Rd,
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

void CoherentPointDriftRigidFisherMises::shiftYCoords(
    const McDMatrix<double>& oldYc, const McDMatrix<double>& Rc,
    const double& s, const McDVector<double>& t, McDMatrix<double>& shiftedYc) {
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

void CoherentPointDriftRigidFisherMises::warpPoint(const McVec3f& point,
                                                   const McDMatrix<double>& Rc,
                                                   const double s,
                                                   const McDVector<double>& t,
                                                   McVec3f& warped) const {

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
    warped = McVec3f(pointcStdMeanShifted[0], pointcStdMeanShifted[1], point.z);
}

void
CoherentPointDriftRigidFisherMises::rescaleYs(McDMatrix<double>& oldYs) const {
    for (int i = 0; i < oldYs.nRows(); i++) {
        McDVector<double> newY = oldYs.getRowVector(i);
        newY *= mMeansAndStds.stdC;
        newY += mMeansAndStds.meanC;
        oldYs.setRowVector(newY, i);
    }
}

double CoherentPointDriftRigidFisherMises::sumSquaredDistances(
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
