#pragma once

#include <mclib/McMat4f.h>
#include <mclib/McVec3f.h>
#include <mclib/McDArray.h>

#include <hxalignmicrotubules/MovingLeastSquares.h>

namespace mtalign {

struct DirectionalPoints;
struct FacingPointSets;

/// `CPDType` enumerates coherent point drift variants.
enum CPDType {
    // The first index is 1, because `CPD_RIGID` has been removed.

    /// `CPD_RIGID_FISHER_MISES` selects algorithm 1 'rotation from line
    /// orientation' or algorithm 2 'rotation, scale, and translation from
    /// orientation and endpoint position' (see Weber 2014), depending on the
    /// value of `CPDParams.useDirections` and `CPDParams.usePositions`.
    CPD_RIGID_FISHER_MISES = 1,

    /// `CPD_ELASTIC` selects algorithm 3 'elastic transformation' (see
    /// Weber 2014).
    CPD_ELASTIC = 2
};

/// `AlignParams` contains the parameters that are common to the coherent point
/// drift algorithms.
struct AlignParams {
    double w;
    int maxIterations;
    float eDiffRelStop;
    float sigmaSquareStop;
    bool useDirections;
};

/// `AlignParamsElastic` contains additional parameters for `cpdElastic()`.
struct AlignParamsElastic : public AlignParams {
    double lambda;
    double beta;
    float sampleDistForWarpingLandmarks;
};

/// `AlignParamsRigidFisherMises` contains additional parameters for
/// `cpdRigidFisherMises()`.
struct AlignParamsRigidFisherMises : public AlignParams {
    bool usePositions;
    bool withScaling;
};

/// `CPDParams` controls the parameters of the coherent point drift
/// algorithms (functions `cpd*()`).
struct CPDParams {
    CPDType type;

    union {
        /// When `type` is `CPD_RIGID_FISHER_MISES`.
        AlignParamsRigidFisherMises rigid;

        /// When `type` is `CPD_ELASTIC`.
        AlignParamsElastic elastic;
    };

    /// `alphaForMLS` controls the moving least squares warp for the elastic
    /// alignment.  It would better be a member of `AlignParamsElastic`, but
    /// `MicrotubuleSpatialGraphAligner` uses it in an unconditional code path.
    float alphaForMLS;

    /// `maxAcceptableSigmaSquare` controls when the result should be
    /// considered unreliable.  If `WarpResult::sigmaSquare` is greater than
    /// the threshold, the transformation should be ignored.
    float maxAcceptableSigmaSquare;
};

/// `AlignInfo` is used to return information about the convergence of the
/// coherent point drift algorithms.
struct AlignInfo {
    /// `sigma^2` upon convergence.
    float sigmaSquare;

    /// `kappa` upon convergence.
    float kappa;

    /// Number of iterations until convergence.
    int numIterations;

    /// Relative expectation difference upon convergence.
    double eDiffRel;

    /// Expectation upon convergence.
    double e;

    /// Time until convergence.
    double timeInSec;
};

/// `WarpResult` is returned by the coherent point drift functions (`cpd()` and
/// variants).
struct WarpResult {
    int type;  // 0 if transform matrix, 1 if mls interpolation.
    MovingLeastSquares mls;
    McMat4f transformMatrix;
    int refSlice;
    int transSlice;
    AlignInfo alignInfo;
};

/// `cpdRigidFisherMises()` applies the coherent point drift algorithm with a
/// rigid deformation model.  It takes the position of the `points` and their
/// direction into account.  The direction distribution is modeled as a
/// Fisher-Mises distribution.
void cpdRigidFisherMises(const FacingPointSets& points, const CPDParams& params,
                         WarpResult& warpResult);

/// `cpdElastic()` applies the coherent point drift algorithm with an elastic
/// deformation model.  It takes the position of the `points` and their
/// direction into account.
void cpdElastic(const FacingPointSets& points, const CPDParams& params,
                WarpResult& warpResult);

/// `cpd()` dispatches to the appropriate coherent point drift algorithm based
/// on the `CPDType` stored in `params.type`.
void cpd(const FacingPointSets& points, const CPDParams& params,
         WarpResult& warpResult);

}  // namespace mtalign
