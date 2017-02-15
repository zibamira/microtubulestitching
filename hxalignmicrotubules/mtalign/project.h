#pragma once

#include <mclib/McHandle.h>
#include <mclib/internal/McString.h>

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
///
/// The slice numbers `refSliceNum` and `transSliceNum` are zero-based.
/// `projectionPlane` specifies the z-value of the target xy-plane.  It is
/// usually computed with `SliceSelector::computeMidPlane()`.  `endPointRegion`
/// specifies the region to be projected either as a percentage of the slice
/// thickness (if `useAbsoluteValueForEndPointRegion` is `false`) or in absolute
/// physical units (if `useAbsoluteValueForEndPointRegion` is `true`).  Use
/// `projectionType = P_ORTHOGONAL` to get the projection as described in the
/// [Weber 2014].
///
/// `maxDistForAngle` limits the distance (in physical units) between the
/// boundary point and a second point on the line that is used in
/// `projectEndPoints()` for computing the `DirectionalPoints::directions`.
///
/// `angleToPlaneFilter` controls which lines are rejected during the
/// projection.  Lines are rejected if their direction points towards the
/// midplane instead of the section center; or if the direction points towards
/// the section center, but the angle with the boundary plane is less than
/// `angleToPlaneFilter * 180 degrees`.  A larger `angleToPlaneFilter` rejects
/// more lines.
///
/// `numMaxPointsForInitTransform` limits the number of points returned by
/// `projectEndPointsSubset()`.  Points whose associated direction is more
/// perpendicular to the boundary are preferred over points whose associated
/// direction is more parallel.
struct EndPointParams {
    EndPointParams() {
        refSliceNum = 0;
        transSliceNum = 1;
        projectionPlane = 0.0;
        useAbsoluteValueForEndPointRegion = false;
        endPointRegion = 40;
        projectionType = P_ORTHOGONAL;
        angleToPlaneFilter = 0.01;
        maxDistForAngle = 2000;
        numMaxPointsForInitTransform = 50;
    }

    int refSliceNum;
    int transSliceNum;
    float projectionPlane;
    bool useAbsoluteValueForEndPointRegion;
    float endPointRegion;
    ProjectionType projectionType;
    float angleToPlaneFilter;
    float maxDistForAngle;
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
