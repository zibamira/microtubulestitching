#pragma once

#include <mclib/McDArray.h>
#include <mclib/McVec3.h>
#include <mclib/McVec2i.h>

namespace mtalign {

/// `TransformType` enumerates which transform to use, e.g. to compute from a
/// matching in `fitTransform()`.
enum TransformType { TF_RIGID, TF_RIGID_ISO_SCALE, TF_AFFINE, TF_NONE };

/// `MatchingAlgorithm` enumerates matching algorithms.
enum MatchingAlgorithm {

    /// `MA_EXACT` indicates to use `matchingExact()`.
    MA_EXACT,

    /// `MA_GREEDY` indicates to use `matchingGreedy()`.
    MA_GREEDY
};

/// `DirectionalPoints` represents points with `positions` and `directions` on
/// a section boundary.  The `directions` point from the boundary point towards
/// the section center.  The direction is reversed in `cpd()` to match the
/// figure in [Weber 2014] that illustrates the model at the section boundary.
struct DirectionalPoints {
    McDArray<McVec3f> positions;
    McDArray<McVec3f> directions;
};

/// `FacingPointSets` represents points on each side of a section boundary.
/// `ref` are the fixed points.  `trans` are the moving points.
struct FacingPointSets {
    DirectionalPoints ref;
    DirectionalPoints trans;
};

/// `Matching` represents a matching computed for `FacingPointSets`.  The
/// integers are indices into the arrays of the `DirectionalPoints` of
/// `FacingPointSets`.
struct Matching {
    McDArray<int> matchedRefPointIds;
    McDArray<int> matchedTransPointIds;
};

/// `MatchingPGM` represents a matching computed for `FacingPointSets` with a
/// probabilistic graphical model (see `matchingPGM()`).  The integers of
/// `matchedRefPointIds` and `matchedTransPointIds` are indices into the arrays
/// of the `DirectionalPoints` of `FacingPointSets`.  The other `int` arrays
/// contain special ref points. `queerEvidence` contains pairs of ref and trans
/// points.  `maxDiff` contains the maximum message difference for each ref
/// point; `entropie` contains the belief entropies.  See implementation of
/// `matchingPGM()` for details.
///
/// The typo in 'entropie' (instead of 'entropy') is kept for backward
/// compatibility with previously computed spatial graphs.
struct MatchingPGM {
    McDArray<int> matchedRefPointIds;
    McDArray<int> matchedTransPointIds;
    McDArray<int> ambiguities;
    McDArray<int> ambiguitiesFromMarginalDiff;
    McDArray<int> assignToMakeTree;
    McDArray<McVec2i> queerEvidence;
    McDArray<int> criticalNodes;
    McDArray<int> criticalNodesWithAdaptiveThreshold;
    McDArray<int> maxEntropieNodes;
    McDArray<float> maxDiff;
    McDArray<float> entropie;
};

}  // namespace mtalign
