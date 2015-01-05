#pragma once

#include <mclib/McVec2d.h>
#include <mclib/McVec2i.h>
#include <mclib/McDArray.h>
#include <mclib/McMat2d.h>
#include <mclib/McDMatrix.h>
#include <vector>
#include <boost/function.hpp>

#include <hxalignmicrotubules/mtalign/cpd.h>
#include <hxalignmicrotubules/api.h>

class QString;

/// Implementation of 'Algorithm 3: Elastic transformation' of the paper Weber
/// et al. "Automated stitching of microtubule centerlines across serial
/// electron tomograms".
///
/// To disable orientation, set `mUseDirection=false`.
class HXALIGNMICROTUBULES_API CoherentPointDriftNLFisherMises {
  public:
    static void shiftYs(const McDMatrix<double>& oldYs,
                        const McDMatrix<double>& G, const McDMatrix<double>& W,
                        McDMatrix<double>& shiftedYs);

    static void convertCoordsToMatrix(const McDArray<McVec3f>& points,
                                      McDMatrix<double>& matrix);

    static void convertDirectionsToMatrix(const McDArray<McVec3f>& directions,
                                          McDMatrix<double>& matrix);

    static void convertMatrixToCoords(McDArray<McVec3f>& points,
                                      const McDMatrix<double>& matrix);

  public:
    CoherentPointDriftNLFisherMises();
    ~CoherentPointDriftNLFisherMises();

    typedef boost::function<void (const char*)> print_t;

    /// Set function used for printing.
    void setPrint(print_t print);

    mtalign::AlignInfo align(McDMatrix<double>& G, McDMatrix<double>& W,
                             McDArray<McVec2i>& correspondences);

    double sumSquaredDistances(const McDMatrix<double>& p1,
                               const McDMatrix<double>& p2);

    void computeP(const McDMatrix<double>& X, const McDMatrix<double>& Y,
                  const double sigmaSquare, const double kappa, const double w,
                  McDMatrix<double>& P, double& LL);

    void rescaleYs(McDMatrix<double>& oldYs) const;

    McDMatrix<double> xs;
    McDMatrix<double> ys;
    McDMatrix<double> xDirs;
    McDMatrix<double> yDirs;

    mtalign::AlignParamsElastic params;

  private:
    // Print message.  Function that is used for printing can be configured
    // with `setPrint()`.
    void print(QString msg);

    print_t mPrint;

    struct MeanAndStd {
        McDVector<double> mean;
        double std;
    };

    MeanAndStd mMeansAndStds;

  private:
    void initializeG(McDMatrix<double>& G, const McDMatrix<double>& ps);

    void getDP(const McDMatrix<double>& P, McDMatrix<double>& dP);

    void normalize();

    double kappaFirstDerivative(const double kappa, const double c,
                                const double Np);

    double kappaSecondDerivative(const double kappa, const double Np);

    bool newtonsMethodForKappa(double& kappa, const double c, const double Np);

    bool computeKappa(const McDMatrix<double>& P, const double NP,
                      double& kappa);
};
