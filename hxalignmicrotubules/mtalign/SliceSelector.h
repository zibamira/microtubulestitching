#pragma once

#include <map>
#include <string>

class HxSpatialGraph;
class SpatialGraphSelection;

namespace mtalign {

/// `MinMax` is used by `SliceSelector` to return the z range of a slice.
struct MinMax {
    float min;
    float max;
};

/// `SliceSelector` partitions a spatial graph into slices based on an
/// attribute and provides functions to compute properties of the slices.
///
/// The attribute is expected to fulfill `isOrderedAttribute()`, that is nodes
/// are grouped by attribute value if visited in index order.  Slices indexes
/// start at 0.  Slices consist of vertices and complete edges.  `getSlice()`
/// will never select individual points.
class SliceSelector {
  public:
    /// `SliceSelector()` initializes an invalid object.
    SliceSelector();

    /// `SliceSelector()` computes the two-way mapping from attribute values to
    /// slice indices.  The `graph` must remain valid as long as the
    /// `SliceSelector` instance is used.
    SliceSelector(const HxSpatialGraph* graph, std::string sliceAttributeName);

    /// `isOrderedAttribute()` tests whether values of the node int attribute
    /// are grouped when nodes are visited in index order.  The attribute
    /// `TransformInfo` on a spatial graph stack fulfills this property.
    static bool isOrderedAttribute(const HxSpatialGraph* sg,
                                   const char* attrName);

    /// `getNumSlices()` does what it suggests.
    int getNumSlices() const;

    /// `getSliceAttributeValueFromIndex()` returns the attribute value that
    /// has been associated with the zero-based `sliceIndex`.
    int getSliceAttributeValueFromIndex(const int sliceIndex) const;

    /// `getSliceIdxOfVertex()` returns the zero-based slice index for a node.
    int getSliceIdxOfVertex(const int vertexIdx) const;

    /// `getSlice()` selects the nodes and edges for `attrValue` in the output
    /// selection `slice`.  The output selection contains no individual points:
    /// `slice.getNumSelectedPoints() == 0`.
    void getSlice(const int attrValue, SpatialGraphSelection& slice) const;

    /// `getZRange()` computes the z-range of the nodes selected in `sel`.
    MinMax getZRange(const SpatialGraphSelection& sel) const;

    /// `getZRangeByIndex()` computes the z-range for the given slice index.
    MinMax getZRangeByIndex(const int sliceIndex) const;

    /// `selectAdjacentHalfSlices()` computes two selections that each contain
    /// one half-slice.  The slices are specified by attribute value.
    void selectAdjacentHalfSlices(const int slice1AttVal,
                                  const int slice2AttVal,
                                  SpatialGraphSelection& halfSlice1,
                                  SpatialGraphSelection& halfSlice2) const;

    /// `computeMidPlane()` returns the position of a z-mid-plane between two
    /// slices specified by zero-based slice index.
    float computeMidPlane(const int slice1Idx, const int slice2Idx) const;

    /// `selectCloseToMidplane()` returns a selection around the midplane of
    /// the two slices for indices `refNum` and `transNum`.  The result is the
    /// union of the two selections computed by `selectAdjacentHalfSlices()`.
    SpatialGraphSelection selectCloseToMidplane(const int refNum,
                                                const int transNum) const;

  private:
    int getSliceIndexFromAttributeValue(const int sliceAttribute) const;

    int getSliceIdxOfEdge(const int edgeIdx) const;

    void
    deselectVerticesFartherAwayFromMidPlaneThan(SpatialGraphSelection& slice,
                                                const float zMidPlane,
                                                const float maxDist) const;

    bool isSliceAboveOther(const int slice1Idx, const int slice2Idx) const;

  private:
    std::string mSliceAttributeName;
    const HxSpatialGraph* mGraph;

    // Maps from the transform attribute to the slice index.
    std::map<int, int> mSliceAttributeToIndexMap;
    std::map<int, int> mSliceIndexToAttributeMap;

    void fillSliceAttributeMap();
};

}  // namespace mtalign
