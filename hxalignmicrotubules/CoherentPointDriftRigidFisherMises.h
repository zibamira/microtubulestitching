#pragma once

#include <mclib/McVec2i.h>
#include <mclib/McVec3f.h>
#include <mclib/McDArray.h>
#include <mclib/McDMatrix.h>

#include <hxalignmicrotubules/mtalign/cpd.h>
#include <hxalignmicrotubules/api.h>

/// This class implements CDP as described in the paper "Point Set
/// registration: Coherent Point Drift", Myronenko et al., 2010.
class HXALIGNMICROTUBULES_API CoherentPointDriftRigidFisherMises {
  public:
    CoherentPointDriftRigidFisherMises();
    ~CoherentPointDriftRigidFisherMises();

    struct meanAndStd {
        McDVector<double> meanC;
        double stdC;
    };

    struct allTerms {
        // these are the terms from chapter 5
        McDMatrix<double> XcHatT_x_dPT1_x_XcHat;
        McDMatrix<double> YcHatT_x_dP1_x_YcHat;
        McDMatrix<double> XcHatT_x_PT_x_YcHat;
        McDMatrix<double> XdT_x_PT_x_Yd;
    };

    McDMatrix<double> Xc;
    McDMatrix<double> Yc;
    McDMatrix<double> Xd;
    McDMatrix<double> Yd;

    meanAndStd mMeansAndStds;

    mtalign::AlignParamsRigidFisherMises params;

    mtalign::AlignInfo align(McDMatrix<double>& Rc, double& s,
                             McDVector<double>& t, McDMatrix<double>& Rd,
                             McDArray<McVec2i>& correspondences);

    mtalign::AlignInfo align(McDMatrix<double>& Rc, double& s,
                             McDVector<double>& t, McDMatrix<double>& Rd,
                             McDArray<McVec2i>& correspondences,
                             double& sigmaSquare, double& kappa);

    double sumSquaredDistances(const McDMatrix<double>& p1,
                               const McDMatrix<double>& p2);

    void shiftYDirections(const McDMatrix<double>& oldYd,
                          const McDMatrix<double>& Rd,
                          McDMatrix<double>& shiftedYd);
    static void shiftYCoords(const McDMatrix<double>& oldYc,
                             const McDMatrix<double>& Rc, const double& s,
                             const McDVector<double>& t,
                             McDMatrix<double>& shiftedYc);

    void getDP(const McDMatrix<double>& P, McDMatrix<double>& dP);
    static double trace(const McDMatrix<double>& mat);
    void computeP(const McDMatrix<double>& Xc, const McDMatrix<double>& Xd,
                  const McDMatrix<double>& Yc, const McDMatrix<double>& Yd,
                  const double sigmaSquare, const double kappa, const double w,
                  McDMatrix<double>& P, double& LL);

    static void multMatrices(const McDMatrix<double>& A,
                             const McDMatrix<double>& B,
                             const McDMatrix<double>& C, const double alpha,
                             const double beta, McDMatrix<double>& result);
    static void solve(const McDMatrix<double>& A, const McDMatrix<double>& B,
                      McDMatrix<double>& X);
    void normalize();
    void rescaleYs(McDMatrix<double>& oldYs) const;

    static double gauss(const McDVector<double>& mean, const double sigma,
                        const McDVector<double> x);
    static double fisherMises(const McDVector<double>& mean, const double kappa,
                              const McDVector<double> x);

    static void convertMcDArraysToMatrix(const McDArray<McVec3f>& coords,
                                         const McDArray<McVec3f>& directions,
                                         McDMatrix<double>& CMat,
                                         McDMatrix<double>& DMat);
    static void convertMatrixTo2McDArrays(McDArray<McVec3f>& points,
                                          McDArray<McVec3f>& directionPoints,
                                          const McDMatrix<double>& matrix);

    double getNP(const McDMatrix<double>& P);
    void computeRc(const allTerms& terms, McDMatrix<double>& Rc);
    void computeRd(const allTerms& terms, McDMatrix<double>& Rd);
    void computeRd2d(const allTerms& terms, McDMatrix<double>& Rd);
    void getMu(const McDMatrix<double>& X, const McDMatrix<double>& P,
               const double NP, McDMatrix<double>& muX);
    void getHat(const McDMatrix<double>& X, const McDMatrix<double>& muX,
                McDMatrix<double>& XHat);
    void splitVector(const McDVector<double>& tmpVec, McDVector<double>& xc,
                     McDVector<double>& xd);
    static void multThreeMatrices(const McDMatrix<double>& A,
                                  const McDMatrix<double>& B,
                                  const McDMatrix<double>& C,
                                  McDMatrix<double>& result);
    void initTerms(allTerms& terms, const McDMatrix<double>& XcHat,
                   const McDMatrix<double>& Xd, const McDMatrix<double>& YcHat,
                   const McDMatrix<double>& Yd, const McDMatrix<double>& P);
    void computeS(const allTerms& terms, const McDMatrix<double>& R, double& s);
    void computeT(const McDMatrix<double>& muX, const double s,
                  const McDMatrix<double>& R, const McDMatrix<double>& muY,
                  McDVector<double>&);
    double computeSigmaSquare(const allTerms& terms,
                              const McDMatrix<double>& Rc, const double s,
                              const double NP);
    bool computeKappa(const allTerms& terms, const McDMatrix<double>& Rd,
                      const double NP, double& kappa);
    void warpPoint(const McVec3f& point, const McDMatrix<double>& R,
                   const double s, const McDVector<double>& t,
                   McVec3f& warped) const;

    double kappaFirstDerivative(const double kappa, const double c,
                                const double Np);
    double kappaSecondDerivative(const double kappa, const double Np);
    bool newtonsMethodForKappa(double& kappa, const double c, const double Np);
    bool optimizeParameters(const allTerms& terms, McDMatrix<double>& R,
                            double& s, double& sigmaSquare, double& kappa,
                            const double Np);
};
