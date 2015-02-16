#include <hxalignmicrotubules/mtalign/SliceSelector.h>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <hxalignmicrotubules/hxtesting.h>
#include <hxspatialgraph/HxSpatialGraph.h>

namespace ht = hxtesting;
namespace ma = mtalign;

namespace {

class mtalign__SliceSelector_withSpatialGraphThreeStackedSections
    : public ::testing::Test {
  protected:
    virtual void SetUp() {
        mSpatialGraph = ht::makeSpatialGraphThreeStackedSections();
        mSliceSelector = ma::SliceSelector(mSpatialGraph, "section");
    }

    McHandle<HxSpatialGraph> mSpatialGraph;
    ma::SliceSelector mSliceSelector;
};

}  // namespace.

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldDistinguishOrderedAndUnorderedAttributes) {
    EXPECT_TRUE(
        ma::SliceSelector::isOrderedAttribute(mSpatialGraph, "section"));
    EXPECT_FALSE(
        ma::SliceSelector::isOrderedAttribute(mSpatialGraph, "unordered"));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldDetermineExpectedNumberOfSections) {
    EXPECT_EQ(3, mSliceSelector.getNumSlices());
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldMapFromSliceIndexToAttribute) {
    const int expectedAttribute = 1;  // of 1, 2, 3.
    EXPECT_EQ(expectedAttribute,
              mSliceSelector.getSliceAttributeValueFromIndex(0));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldMapFromVertexToSliceIndex) {
    const int expectedIndex = 2;  // of 0, 1, 2.
    EXPECT_EQ(expectedIndex, mSliceSelector.getSliceIdxOfVertex(5));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldComputeZRange) {
    SpatialGraphSelection sel(mSpatialGraph);
    sel.selectVertex(2);
    sel.selectVertex(3);
    const ma::MinMax mm = mSliceSelector.getZRange(sel);
    EXPECT_FLOAT_EQ(1.0, mm.min);
    EXPECT_FLOAT_EQ(1.9, mm.max);
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldComputeMidPlane) {
    const int sliceIdxFirst = 0;  // of 0, 1, 2.
    const int sliceIdxSecond = 1;
    EXPECT_FLOAT_EQ(
        0.95, mSliceSelector.computeMidPlane(sliceIdxFirst, sliceIdxSecond));
    EXPECT_FLOAT_EQ(
        0.95, mSliceSelector.computeMidPlane(sliceIdxSecond, sliceIdxFirst));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldSelectSlice) {
    const int attrVal = 2;  // of 1, 2, 3.
    SpatialGraphSelection sel;
    mSliceSelector.getSlice(attrVal, sel);
    // Should select the nodes {2, 3}.
    EXPECT_EQ(2, sel.getNumSelectedVertices());
    EXPECT_EQ(2, sel.getSelectedVertex(0));
    EXPECT_EQ(3, sel.getSelectedVertex(1));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       shouldSelectAdjacentHalfSlices) {
    const int attrValLower = 1;  // of 1, 2, 3.
    const int attrValHigher = 2;
    SpatialGraphSelection selLower;
    SpatialGraphSelection selHigher;
    mSliceSelector.selectAdjacentHalfSlices(attrValLower, attrValHigher,
                                            selLower, selHigher);
    // Should select two node groups: {1} and {2}
    EXPECT_EQ(1, selLower.getNumSelectedVertices());
    EXPECT_EQ(1, selLower.getSelectedVertex(0));
    EXPECT_EQ(1, selHigher.getNumSelectedVertices());
    EXPECT_EQ(2, selHigher.getSelectedVertex(0));
}

TEST_F(mtalign__SliceSelector_withSpatialGraphThreeStackedSections,
       selectsCloseToMidplane) {
    const int idxLower = 0;
    const int idxHigher = 1;
    const SpatialGraphSelection sel =
        mSliceSelector.selectCloseToMidplane(idxLower, idxHigher);
    // Expected to select {1, 2}.
    EXPECT_EQ(2, sel.getNumSelectedVertices());
    EXPECT_EQ(1, sel.getSelectedVertex(0));
    EXPECT_EQ(2, sel.getSelectedVertex(1));
}
