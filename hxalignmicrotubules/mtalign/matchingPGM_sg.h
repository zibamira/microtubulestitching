#pragma once

class SpatialGraphSelection;
class HxSpatialGraph;

namespace mtalign {

struct FacingPointSets;
struct MatchingPGM;
struct MatchingPGMParams;

/// This variant of `matchingPGM()` uses the fixed user-defined assignments
/// from the attribute "Evidence" on `sg`; it ignores `params.evidence`.
/// Vertices that are explicitly unassigned in "Evidence" are removed from
/// `refVertexSelection` and `transVertexSelection`.  The `MatchingPGM` results
/// are stored as attributes onto `sg` after clearing the attributes for
/// vertices in `aroundMidPlaneSelection`.
MatchingPGM matchingPGM(const FacingPointSets& pts, MatchingPGMParams params,
                        HxSpatialGraph* sg,
                        SpatialGraphSelection& refVertexSelection,
                        SpatialGraphSelection& transVertexSelection,
                        const SpatialGraphSelection& aroundMidPlaneSelection);

}  // namespace mtalign