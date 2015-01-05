#pragma once

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <mclib/McDArray.h>

class McVec2i;

namespace mtalign {

struct FacingPointSets;
struct MatchingPGM;

struct MatchingPGMParams {
    /// `evidence` contains fixed user-defined assignments as pairs of 'ref'
    /// and 'trans' indices into the two `FacingPointSets`.
    McDArray<McVec2i> evidence;

    /// `weightConfig` specifies the `PGMPairWeights` used for the
    /// `PGMMatcher`.
    PGMPairWeightsParams weightConfig;

    /// `pairFactorParam` is used in the `PGMMatcher`.
    double pairFactorParam;
};

/// `matchingPGM()` computes a matching between the two `FacingPointSets` as
/// specified by `MatchingPGMParams`.  `params.evidence` contains fixed
/// user-defined assignments as pairs of 'ref' and 'trans' indices into the two
/// `FacingPointSets`.
MatchingPGM matchingPGM(const FacingPointSets& pts,
                        const MatchingPGMParams& params);

}  // namespace mtalign.
