#pragma once

#include <mclib/McHandle.h>
#include <mclib/McString.h>

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

class HxSpatialGraph;
class SpatialGraphSelection;

namespace mtalign {

struct FacingPointSets;

/// `ProjectionType` selects the algorithm to project points in
/// `projectEndPoints()`.  Use `P_ORTHOGONAL` to select the projection
/// described in [Weber 2014].
enum ProjectionType {
    // Uses same indices as MicrotubuleSpatialGraphAligner::ProjectionTypes.

    /// Orthogonal projection as described in [Weber 2014].
    P_ORTHOGONAL,

    /// See source.
    P_LINEAR,

    /// See source.
    P_TANGENT,

    /// See source.
    P_FIT_0,

    /// See source.
    P_FIT_1,

    /// See source.
    P_FIT_2,

    /// See source.
    P_FIT_3,

    /// See source.
    P_APPROX_TANGENT,

    /// See source.
    P_NONE
};

/// `EndPointParams` specifies how endpoints on sections boundaries are
/// extracted from a spatial graph stack.
struct EndPointParams {
    int refSliceNum;
    int transSliceNum;
    float endPointRegion;
    ProjectionType projectionType;
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
