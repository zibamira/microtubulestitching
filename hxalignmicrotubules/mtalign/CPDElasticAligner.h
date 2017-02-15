#pragma once

#include <boost/function.hpp>

#include <mclib/McDArray.h>
#include <mclib/internal/McDMatrix.h>

#include <hxalignmicrotubules/mtalign/Context.h>
#include <hxalignmicrotubules/mtalign/cpd.h>
#include <hxalignmicrotubules/mtalign/data.h>

class QString;

namespace mtalign {

/// `CPDElasticAligner` implements 'Algorithm 3: Elastic transformation' of
/// [Weber 2014].  It is used internally to implement `cpd()`.  It is also used
/// directly by `HxCPDSpatialGraphWarp`, which is kept for backward
/// compatibility with the supplementary data to [Weber 2014].  Future clients
/// should not use `CPDElasticAligner` but instead call `cpd()`.
///
/// A `Context` can be used to configure the run-time environment used for
/// printing.  The default is to print to stdout.
class CPDElasticAligner {
  public:
    ///
    CPDElasticAligner();

    /// `setContext()` configures the run-time environment.
    void setContext(Context* ctx);

    /// `setPoints()` set the point sets that are used in `align()`.
    void setPoints(const mtalign::FacingPointSets& points);

    /// `align()` computes the alignment.  It returns the shifted positions of
    /// the original `points.trans.positions` and stores additional info into
    /// the out parameter `info`.
    McDArray<McVec3f> align(mtalign::AlignInfo& info);

    /// `params` controls details of the algorithm.
    mtalign::AlignParamsElastic params;

  protected:  // Members in this block are used by tests.
    static void shiftYs(const McDMatrix<double>& oldYs,
                        const McDMatrix<double>& G, const McDMatrix<double>& W,
                        McDMatrix<double>& shiftedYs);

    /// `convertCoordsToMatrix` converts a vector of positions to its matrix
    /// representation.
    static void convertCoordsToMatrix(const McDArray<McVec3f>& points,
                                      McDMatrix<double>& matrix);

    /// `align()` computes the alignment without applying it to compute shifted
    /// positions.
    mtalign::AlignInfo align(McDMatrix<double>& G, McDMatrix<double>& W);

    /// `sumSquaredDistances()` is a helper function.
    double sumSquaredDistances(const McDMatrix<double>& p1,
                               const McDMatrix<double>& p2);

    /// `computeP()` computes the `P` matrix (see algorithm in [Weber 2014]).
    void computeP(const McDMatrix<double>& X, const McDMatrix<double>& Y,
                  const double sigmaSquare, const double kappa, const double w,
                  McDMatrix<double>& P, double& LL);

    /// `xs` stores the normalized `points.ref.positions` in matrix
    /// representation.
    McDMatrix<double> xs;

    /// `ys` stores the normalized `points.trans.positions` in matrix
    /// representation.
    McDMatrix<double> ys;

    /// `xDirs` stores the `points.ref.directions` in matrix representation.
    McDMatrix<double> xDirs;

    /// `yDirs` stores the `points.trans.directions` in matrix representation.
    McDMatrix<double> yDirs;

  private:  // These are truly private.

    /// `MeanAndStd` is used by `normalize()` to store statistics.
    struct MeanAndStd {
        McDVector<double> mean;
        double std;
    };

    /// `mMeansAndStds` are filled by `normalize()`.
    MeanAndStd mMeansAndStds;

    /// `normalize()` transforms the `xs` and `ys` into a standard coordinate
    /// system.
    void normalize();

    /// `rescaleYs()` reverses the effect of `normalize()`, which is called by
    /// `align()` to transform the `xs` and `ys` into a standard coordinate
    /// system.  The `oldYs` are modified in place.
    void rescaleYs(McDMatrix<double>& oldYs) const;

    /// `print()` is used for printing via the `Context`.
    void print(QString msg);

    Context* mContext;

    /// `computeKappa()` calls Newton's method for computing kappa.  It is a
    /// members as all related functions, because they use `print()`.
    bool computeKappa(const McDMatrix<double>& P, const double NP,
                      double& kappa);

    /// `newtonsMethodForKappa()` searches the optimal kappa.
    bool newtonsMethodForKappa(double& kappa, const double c, const double Np);

    /// `kappaFirstDerivative()` implements the first derivative for the
    /// Newton's method (see algorithm in [Weber 2014]).
    double kappaFirstDerivative(const double kappa, const double c,
                                const double Np);

    /// `kappaSecondDerivative()` implements the second derivative for the
    /// Newton's method (see algorithm in [Weber 2014]).
    double kappaSecondDerivative(const double kappa, const double Np);
};

}  // namespace mtalign
