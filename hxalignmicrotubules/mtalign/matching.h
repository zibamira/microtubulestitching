#pragma once

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <hxalignmicrotubules/mtalign/data.h>

#include <mclib/McMat4.h>

namespace mtalign {

struct FacingPointSets;
struct Matching;
struct MatchingParams;

/// `matchingDirect()` calls the matching algorithm specified in
/// `MatchingAlgorithm`, i.e. `matchingGreedy()` or `matchingExact()`.
Matching matchingDirect(const FacingPointSets& pts, McMat4f startTransformation,
                        const PGMPairWeightsParams& weightConfig,
                        const MatchingAlgorithm algo);

/// `matchingTwoStepBest()` computes a matching in two stages:  It first uses
/// `matchingCliqueTransforms()` to compute a set of initial transformations.
/// It then uses the algorithm specified in `MatchingAlgorithm` to compute a
/// matching for each transformation and returns the largest matching.
Matching matchingTwoStepBest(const FacingPointSets& pts,
                             const MatchingParams& params,
                             const PGMPairWeightsParams& weightConfig,
                             const MatchingAlgorithm algo);

}  // namespace mtalign

