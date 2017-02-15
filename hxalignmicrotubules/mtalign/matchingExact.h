#pragma once

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

#include <mclib/McMat4.h>

class PointMatchingDataStruct;

namespace mtalign {

struct FacingPointSets;
struct Matching;

/// `matchingExact()` uses `pointmatching/ExactPointMatchingAlgorithm` to
/// compute a refined matching of two facing point sets assuming the initial
/// transformation `startTransformation` for the second 'trans' set of
/// `FacingPointSets`.
Matching matchingExact(const FacingPointSets& pts, McMat4f startTransformation,
                       const PGMPairWeightsParams& weightConfig);

}  // namespace mtalign.

