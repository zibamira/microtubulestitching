#pragma once

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <mclib/McDArray.h>

class McVec2i;

namespace mtalign {

struct FacingPointSets;
struct MatchingPGM;

/// `MatchingPGMParams` contains the parameters for `matchingPGM()`.
struct MatchingPGMParams {
    /// `evidence` contains fixed user-defined assignments as pairs of 'ref'
    /// and 'trans' indices into the `FacingPointSets`.
    McDArray<McVec2i> evidence;

    /// `weightConfig` specifies the params of the `PGMPairWeights` used for
    /// the `PGMMatcher`.
    PGMPairWeightsParams weightConfig;

    /// `pairFactorParam` is used for the `PGMMatcher`.
    double pairFactorParam;
};

/// `matchingPGM()` computes a matching between the `FacingPointSets` as
/// specified by `MatchingPGMParams`.  `params.evidence` contains fixed
/// user-defined assignments as pairs of 'ref' and 'trans' indices into the
/// `FacingPointSets`.
MatchingPGM matchingPGM(const FacingPointSets& pts,
                        const MatchingPGMParams& params);

}  // namespace mtalign.
