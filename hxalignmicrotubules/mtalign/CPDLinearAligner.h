#pragma once

#include <mclib/McVec2i.h>
#include <mclib/McVec3f.h>
#include <mclib/McDArray.h>
#include <mclib/McDMatrix.h>

#include <hxalignmicrotubules/mtalign/cpd.h>

namespace mtalign {

/// `CPDLinearAligner` implements the linear alignment algorithms from [Weber
/// 2014].
class CPDLinearAligner {
  public:
    /// `CPDLinearAligner` constructs an aligner with default parameters.  It
    /// is unclear whether they are useful.
    CPDLinearAligner();

    /// `setPoints()` sets the points with directions to be aligned.
    void setPoints(const mtalign::FacingPointSets& points);

    /// `params` control details of the algorithm.
    mtalign::AlignParamsLinear params;

    /// `align()` computes the alignment between the `points`.  It returns
    /// information about convergence and stores the alignment result in the
    /// out params `Rc`, `s`, `t`, and `Rd`.  The results are valid in a
    /// normalized coordinate system, as computed with `normalize()`.  Use
    /// `getTransformMat4f()` or `warpPoint()` to apply the transformation to
    /// the original points.
    mtalign::AlignInfo align(McDMatrix<double>& Rc, double& s,
                             McDVector<double>& t, McDMatrix<double>& Rd);

    /// `warpPoint()` can be used after calling `align()` to apply the
    /// transformation to `point`.  `warpPoint()` correctly considers the
    /// rescaling that has been computed by `align()`.
    McVec3f warpPoint(const McVec3f& point, const McDMatrix<double>& R,
                      const double s, const McDVector<double>& t) const;

    /// `getTransformMat4f()` can be used after calling `align()` to compute a
    /// 4d transformation matrix that can be applied to the original points.
    McMat4f getTransformMat4f(McDMatrix<double> R, double s,
                              McDVector<double> t);

  protected:  // Protected members are used in tests.
    McDMatrix<double> Xc;
    McDMatrix<double> Yc;
    McDMatrix<double> Xd;
    McDMatrix<double> Yd;

    /// `MeanAndStd` is used by `normalize()` to store statistics.
    struct MeanAndStd {
        McDVector<double> meanC;
        double stdC;
    };

    /// `mMeansAndStds` is filled by `normalize()`.
    MeanAndStd mMeansAndStds;

    static void shiftYCoords(const McDMatrix<double>& oldYc,
                             const McDMatrix<double>& Rc, const double& s,
                             const McDVector<double>& t,
                             McDMatrix<double>& shiftedYc);

    void computeP(const McDMatrix<double>& Xc, const McDMatrix<double>& Xd,
                  const McDMatrix<double>& Yc, const McDMatrix<double>& Yd,
                  const double sigmaSquare, const double kappa, const double w,
                  McDMatrix<double>& P, double& LL);

    static double getNP(const McDMatrix<double>& P);

    static void getMu(const McDMatrix<double>& X, const McDMatrix<double>& P,
                      const double NP, McDMatrix<double>& muX);

    static void getHat(const McDMatrix<double>& X, const McDMatrix<double>& muX,
                       McDMatrix<double>& XHat);

    static void computeT(const McDMatrix<double>& muX, const double s,
                         const McDMatrix<double>& R,
                         const McDMatrix<double>& muY, McDVector<double>&);

  private:
    /// `normalize()` transforms the `Xc` and `Yc` into a standard coordinate
    /// system.
    void normalize();
};

/// `CPDLinearAlignerTerms` is used internally.  It implements some more
/// complex terms used for linear alignment (see supporting information [Weber
/// 2014]).
class CPDLinearAlignerTerms {
  public:
    CPDLinearAlignerTerms(const McDMatrix<double>& XcHat,
                          const McDMatrix<double>& Xd,
                          const McDMatrix<double>& YcHat,
                          const McDMatrix<double>& Yd,
                          const McDMatrix<double>& P);

    void computeS(const McDMatrix<double>& R, double& s) const;

    void computeRc(McDMatrix<double>& Rc) const;

    void computeRd2d(McDMatrix<double>& Rd) const;

    double computeSigmaSquare(const McDMatrix<double>& Rc, const double s,
                              const double NP) const;

    bool computeKappa(const McDMatrix<double>& Rd, const double NP,
                      double& kappa) const;

    bool optimizeParameters(McDMatrix<double>& R, double& s,
                            double& sigmaSquare, double& kappa, const double Np,
                            const bool withScaling) const;

  private:
    McDMatrix<double> XcHatT_x_dPT1_x_XcHat;
    McDMatrix<double> YcHatT_x_dP1_x_YcHat;
    McDMatrix<double> XcHatT_x_PT_x_YcHat;
    McDMatrix<double> XdT_x_PT_x_Yd;
};

}  // namespace mtalign
