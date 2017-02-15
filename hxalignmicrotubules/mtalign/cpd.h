#pragma once


#include <hxalignmicrotubules/api.h>
#include <mclib/McDArray.h>
#include <mclib/McMat4.h>
#include <mclib/McVec2.h>
#include <mclib/McVec3.h>

namespace mtalign {

struct Context;

struct DirectionalPoints;
struct FacingPointSets;

/// `CPDType` enumerates coherent point drift variants.
enum CPDType {
    // The first index is 1, because `CPD_RIGID` has been removed.

    /// `CPD_LINEAR` selects algorithm 1 'rotation from line
    /// orientation' or algorithm 2 'rotation, scale, and translation from
    /// orientation and endpoint position' (see Weber 2014), depending on the
    /// value of `AlignParamsLinear.useDirections` and
    /// `AlignParamsLinear.usePositions`.
    CPD_LINEAR = 1,

    /// `CPD_ELASTIC` selects algorithm 3 'elastic transformation' (see
    /// Weber 2014).
    CPD_ELASTIC = 2
};

/// `AlignParams` contains the parameters that are common to the coherent point
/// drift algorithms.
///
/// See supplementary figures 1 to 3 of [Weber 2014].  `w` controls the weight of
/// the outlier term.   `useDirections` controls whether the directions are
/// used.  The remaining parameters `maxIterations`, `eDiffRelStop`, and
/// `sigmaSquareStop` control the termination of the EM algorithms.
struct AlignParams {
    // No ctor, because it would be incompatible with the union in `CPDParams`.
    double w;
    bool useDirections;
    int maxIterations;
    float eDiffRelStop;
    float sigmaSquareStop;
};

/// `AlignParamsLinear` contains additional parameters for `cpdLinear()`.
///
/// `usePositions` controls whether the endpoint positions are used.
/// `withScaling` controls whether the transformation includes uniform scaling.
struct AlignParamsLinear : public AlignParams {
    bool usePositions;
    bool withScaling;
};

/// `AlignParamsElastic` contains additional parameters for `cpdElastic()`.
///
/// `lambda` and `beta` control the CPD algorithm (see supplementary figure 3
/// of [Weber 2014]).
///
/// `sampleDistForWarpingLandmarks` specifies the minimal distance between
/// landmarks in the output.  If two landmarks are closer than the distance,
/// one of the two will be removed.  Use `sampleDistForWarpingLandmarks=0` to
/// disable subsampling.
struct AlignParamsElastic : public AlignParams {
    double lambda;
    double beta;
    float sampleDistForWarpingLandmarks;
};

/// `CPDParams` controls the parameters of the coherent point drift algorithms
/// (functions `cpd*()`).
///
/// Use `defaultsLinear()` or `defaultsElastic()` to get reasonable defaults as
/// used for [Weber 2014].
struct CPDParams {
    CPDType type;

    union {
        /// When `type` is `CPD_LINEAR`.
        AlignParamsLinear linear;

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

    static CPDParams defaultsLinear() {
        CPDParams p;
        p.type = CPD_LINEAR;
        p.linear.w = 0.1;
        p.linear.maxIterations = 200;
        p.linear.eDiffRelStop = 1.e-5f;
        p.linear.sigmaSquareStop = 1.e-7f;
        p.linear.useDirections = true;
        p.linear.usePositions = true;
        p.linear.withScaling = true;
        p.alphaForMLS = 2;
        p.maxAcceptableSigmaSquare = 0;
        return p;
    }

    static CPDParams defaultsElastic() {
        CPDParams p;
        p.type = CPD_ELASTIC;
        p.elastic.w = 0.1;
        p.elastic.maxIterations = 200;
        p.elastic.eDiffRelStop = 1.e-5f;
        p.elastic.sigmaSquareStop = 1.e-7f;
        p.elastic.useDirections = true;
        p.elastic.lambda = 1;
        p.elastic.beta = 10;
        p.elastic.sampleDistForWarpingLandmarks = 0;
        p.alphaForMLS = 2;
        p.maxAcceptableSigmaSquare = 0;
        return p;
    }
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

/// `MLSParams` are pairs of corresponding landmarks `(p, q)` that are used to
/// define a moving least squares warp together with the `alpha` parameter.
struct MLSParams {
    float alpha;
    McDArray<McVec2d> ps;
    McDArray<McVec2d> qs;
};

/// `WarpType` enumerates the possible transformation types in `WarpResult`.
/// The integer values should be kept for backward compatibility.
enum WarpType {
    WT_LINEAR = 0,
    WT_ELASTIC = 1
};

/// `WarpResult` is returned by the coherent point drift functions (`cpd()` and
/// variants).
struct WarpResult {
    WarpType type;

    /// Only used when `type == WT_LINEAR`.
    McMat4f transformMatrix;

    /// Only used when `type == WT_ELASTIC`.
    MLSParams mlsParams;

    int refSlice;
    int transSlice;
    AlignInfo alignInfo;
};

/// `cpdLinear()` applies the coherent point drift algorithm with a
/// rigid deformation model.  It takes the position of the `points` and their
/// direction into account.  The direction distribution is modeled as a
/// Fisher-Mises distribution.
void cpdLinear(const FacingPointSets& points, WarpResult& warpResult,
               const CPDParams& params, Context* ctx = 0);

/// `cpdElastic()` applies the coherent point drift algorithm with an elastic
/// deformation model.  It takes the position of the `points` and their
/// direction into account.
void cpdElastic(const FacingPointSets& points, WarpResult& warpResult,
                const CPDParams& params, Context* ctx = 0);

/// `cpd()` dispatches to the appropriate coherent point drift algorithm based
/// on the `CPDType` stored in `params.type`.
HXALIGNMICROTUBULES_API void cpd(const FacingPointSets& points, WarpResult& warpResult,
                                 const CPDParams& params, Context* ctx = 0);

}  // namespace mtalign
