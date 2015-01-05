#pragma once

#include <mclib/McHandle.h>
#include <mclib/McString.h>

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

class HxSpatialGraph;
class SpatialGraphSelection;

namespace mtalign {

struct FacingPointSets;

/// `EndPointParams` specifies how endpoints on sections boundaries are
/// extracted from a spatial graph stack.
struct EndPointParams {
    int refSliceNum;
    int transSliceNum;
    float endPointRegion;
    McString projectionType;
    int transformType;
    float projectionPlane;
    float angleToPlaneFilter;
    float maxDistForAngle;
    bool useAbsouteValueForEndPointRegion;
    int numMaxPointsForInitTransform;
};

/// `projectEndPoints()` computes two sets of endpoints for the two facing
/// boundaries of two consecutive sections, according to the `EndPointParams`.
/// Sections are defined by different values of the sptial graph attribute
/// `TransformInfo`.
FacingPointSets projectEndPoints(const HxSpatialGraph* mGraph,
                                 const EndPointParams& params);

/// This variant of `projectEndPoints()` stores the subsets of the spatial
/// graph that were used for the two sections in `refSelection` and
/// `transSelection`.
FacingPointSets projectEndPoints(const HxSpatialGraph* mGraph,
                                 SpatialGraphSelection& refSelection,
                                 SpatialGraphSelection& transSelection,
                                 const EndPointParams& params);

/// `projectEndPointsSubset()` is identical to `projectEndPoints()` but limits
/// the number of projected points to
/// `EndPointParams::numMaxPointsForInitTransform`.
FacingPointSets projectEndPointsSubset(const HxSpatialGraph* mGraph,
                                       SpatialGraphSelection& refSelection,
                                       SpatialGraphSelection& transSelection,
                                       const EndPointParams& params);

}  // namespace mtalign
