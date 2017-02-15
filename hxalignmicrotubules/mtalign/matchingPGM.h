#pragma once

#include <hxalignmicrotubules/api.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <mclib/McDArray.h>

class McVec2i;

namespace mtalign {

struct Context;
struct FacingPointSets;
struct MatchingPGM;

/// `MatchingPGMParams` contains the parameters for `matchingPGM()`.
///
/// `defaultsWeber2014()` returns the parameters used in [Weber 2014].
struct MatchingPGMParams {
    /// `evidence` contains fixed user-defined assignments as pairs of 'ref'
    /// and 'trans' indices into the `FacingPointSets`.
    McDArray<McVec2i> evidence;

    /// `weightConfig` specifies the params of the `PGMPairWeights` used for
    /// the `PGMMatcher`.
    PGMPairWeightsParams weightConfig;

    /// `pairFactorParam` is `lambda_s^-1` from [Weber 2014].  It is used for
    /// the `PGMMatcher`.
    double pairFactorParam;

    static MatchingPGMParams defaultsWeber2014() {
        MatchingPGMParams p;
        p.weightConfig.weightType = PGMPairWeightsParams::EXPONENTIAL;
        p.weightConfig.dist3dParam = 243;  // `lambda_c^-1`.
        p.weightConfig.distProjectedParam = 170;  // `lambda_p^-1`.
        p.weightConfig.angleWeightParam = 5.8;  // `lambda_alpha^-1`.
        p.pairFactorParam = 293;  // `lambda_s^1`.
        p.weightConfig.dummySignificance = 0.01;  // `r`.

        // `-log(r)` times params from above; default values from microtubule
        // align filament editor toolbox.
        p.weightConfig.distanceThreshold3d = 1119;
        p.weightConfig.distanceThresholdProjected = 783;
        p.weightConfig.angleThreshold = 26.7;

        // Activate all factors.
        p.weightConfig.useDist3dWeight = true;
        p.weightConfig.useProjectedDistWeight = true;
        p.weightConfig.useAngleWeight = true;
        p.weightConfig.useDistanceThreshold3d = true;
        p.weightConfig.useDistanceThresholdProjected = true;
        p.weightConfig.useAngleThreshold = true;

        return p;
    }
};

/// `matchingPGM()` computes a matching between the `FacingPointSets` as
/// specified by `MatchingPGMParams`.  `params.evidence` contains fixed
/// user-defined assignments as pairs of 'ref' and 'trans' indices into the
/// `FacingPointSets`.
HXALIGNMICROTUBULES_API MatchingPGM matchingPGM(const FacingPointSets& pts,
                                                const MatchingPGMParams& params, Context* ctx = 0);

}  // namespace mtalign.
