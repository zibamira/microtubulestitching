#pragma once

#include <hxalignmicrotubules/mtalign/data.h>

#include <mclib/McMat4.h>

namespace mtalign {

struct FacingPointSets;
struct Matching;

/// `fitTransform()` fits an transform as specified by `TransformType` that
/// transforms the second `trans` set onto the first `ref` set of
/// `FacingPointSets` assuming the correspondences in `Matching`.
McMat4f fitTransform(const FacingPointSets& pts, const Matching& matching,
                     TransformType tfType);

/// `fitTransformRigid()` fits a transform that is either completely rigid or
/// rigid + iso-scale as specified by `TransformType`.
McMat4f fitTransformRigid(const FacingPointSets& pts, const Matching& matching,
                          TransformType tfType);

/// `fitTransformAffine()` fits an affine transform.
McMat4f fitTransformAffine(const FacingPointSets& pts,
                           const Matching& matching);

}  // namespace mtalign

