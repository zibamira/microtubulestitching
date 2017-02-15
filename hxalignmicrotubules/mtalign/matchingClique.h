#pragma once

#include <mclib/McDArray.h>
#include <hxalignmicrotubules/mtalign/data.h>

#include <mclib/McMat4.h>

namespace mtalign {

struct FacingPointSets;

/// `MatchingParams` configures details of `matchingCliqueTransforms()`.
struct MatchingParams {
    /// `maxDistanceForGraphConstruction` limits which points are considered by
    /// `matchingCliqueTransforms()`.  A pair of points is only considered if
    /// they are closer than the threshold.
    float maxDistanceForGraphConstruction;

    /// `maxAngleDiffForInitMatching` limits which points are considered by
    /// `matchingCliqueTransforms()`.  A pair of points is only considered if
    /// the difference in direction at the points is below the threshold.
    float maxAngleDiffForInitMatching;

    /// `minCliqueSizeFraction` configures `matchingCliqueTransforms()` to
    /// expect cliques that are at least as large as the specified fraction of
    /// points of the smaller endpoint set.
    float minCliqueSizeFraction;

    /// `maxNumCliques` limits the maximum number of cliques that
    /// `matchingCliqueTransforms()` uses.
    int maxNumCliques;

    /// `transformType` specifies the transform used by
    /// `matchingCliqueTransforms()`.  It must be `TF_RIGID` or
    /// `TF_RIGID_ISO_SCALE`.
    TransformType transformType;
};

/// `matchingCliqueTransforms()` computes initial transformations for two
/// facing point sets.  Algorithm details are configured via `MatchingParams`.
McDArray<McMat4f> matchingCliqueTransforms(const FacingPointSets& pts,
                                           const MatchingParams& params);

}  // namespace mtalign.

