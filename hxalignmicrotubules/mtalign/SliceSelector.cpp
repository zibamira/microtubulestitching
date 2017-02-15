#include <hxalignmicrotubules/mtalign/SliceSelector.h>

#include <algorithm>
#include <limits>
#include <vector>

#include <mclib/McMath.h>

#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/SpatialGraphSelection.h>

namespace ma = mtalign;

ma::SliceSelector::SliceSelector() : mGraph(0) {}

ma::SliceSelector::SliceSelector(const HxSpatialGraph* graph,
                                 std::string sliceAttributeName)
    : mSliceAttributeName(sliceAttributeName), mGraph(graph) {
    fillSliceAttributeMap();
}

int ma::SliceSelector::getNumSlices(void) const {
    return mSliceIndexToAttributeMap.size();
}

void ma::SliceSelector::fillSliceAttributeMap() {
    EdgeVertexAttribute* sliceAtt =
        dynamic_cast<EdgeVertexAttribute*>(mGraph->findAttribute(
            HxSpatialGraph::VERTEX, McString(mSliceAttributeName.c_str())));
    if (!sliceAtt)
        return;
    mcrequire(isOrderedAttribute(mGraph, mSliceAttributeName.c_str()));
    int sliceCounter = -1;
    for (int i = 0; i < mGraph->getNumVertices(); i++) {
        int sliceVal = sliceAtt->getIntDataAtIdx(i);
        if (mSliceAttributeToIndexMap.find(sliceVal) ==
            mSliceAttributeToIndexMap.end())
            sliceCounter++;
        mSliceAttributeToIndexMap.insert(
            std::pair<int, int>(sliceVal, sliceCounter));
        mSliceIndexToAttributeMap.insert(
            std::pair<int, int>(sliceCounter, sliceVal));
    }
}

bool ma::SliceSelector::isOrderedAttribute(const HxSpatialGraph* sg,
                                           const char* attrName) {
    EdgeVertexAttribute* attr = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, McString(attrName)));
    mcrequire(attr);

    // If node attr is ordered in groups, unique will compress each group to a
    // single element, and the number of such groups will match the number of
    // unique elements, which is determined after sorting.
    const int nNodes = sg->getNumVertices();
    std::vector<int> vec;
    vec.resize(nNodes);
    for (int i = 0; i < nNodes; i++) {
        vec[i] = attr->getIntDataAtIdx(i);
    }
    const std::vector<int>::iterator b = vec.begin();
    const std::vector<int>::iterator e = vec.end();
    const std::vector<int>::iterator u = std::unique(b, e);
    const unsigned int nUniqUnsorted = u - b;
    std::sort(b, u);
    const unsigned int nUniqSorted = std::unique(b, u) - b;
    return nUniqUnsorted == nUniqSorted;
}

int
ma::SliceSelector::getSliceAttributeValueFromIndex(const int sliceIndex) const {
    return mSliceIndexToAttributeMap.find(sliceIndex)->second;
}

int ma::SliceSelector::getSliceIndexFromAttributeValue(
    const int sliceAttributeValue) const {
    std::map<int, int>::const_iterator indexIterator =
        mSliceAttributeToIndexMap.find(sliceAttributeValue);
    if (indexIterator != mSliceAttributeToIndexMap.end())
        return indexIterator->second;
    else
        return -1;
}

void ma::SliceSelector::getSlice(const int attValue,
                                 SpatialGraphSelection& slice) const {

    slice.resize(mGraph->getNumVertices(), mGraph->getNumEdges());
    int sliceIndex = getSliceIndexFromAttributeValue(attValue);

    for (int i = 0; i < mGraph->getNumVertices(); i++) {
        if (getSliceIdxOfVertex(i) == sliceIndex) {
            slice.selectVertex(i);
        }
    }
    for (int i = 0; i < mGraph->getNumEdges(); i++) {
        if (getSliceIdxOfEdge(i) == sliceIndex) {
            slice.selectEdge(i);
        }
    }
}

int ma::SliceSelector::getSliceIdxOfVertex(const int vertexIdx) const {
    EdgeVertexAttribute* vAtt =
        dynamic_cast<EdgeVertexAttribute*>(mGraph->findAttribute(
            HxSpatialGraph::VERTEX, McString(mSliceAttributeName.c_str())));
    if (!vAtt)
        return -1;
    int vertexAttVal = vAtt->getIntDataAtIdx(vertexIdx);
    return getSliceIndexFromAttributeValue(vertexAttVal);
}

int ma::SliceSelector::getSliceIdxOfEdge(const int edgeIdx) const {
    EdgeVertexAttribute* eAtt =
        dynamic_cast<EdgeVertexAttribute*>(mGraph->findAttribute(
            HxSpatialGraph::EDGE, McString(mSliceAttributeName.c_str())));
    if (!eAtt)
        return -1;
    if (!eAtt)
        return -1;
    int edgeAttVal = eAtt->getIntDataAtIdx(edgeIdx);
    return getSliceIndexFromAttributeValue(edgeAttVal);
}

ma::MinMax
ma::SliceSelector::getZRange(const SpatialGraphSelection& sel) const {
    float min = std::numeric_limits<float>::infinity();
    float max = -std::numeric_limits<float>::infinity();

    SpatialGraphSelection::Iterator iter(sel);
    iter.edges.reset();
    for (int edgeNum = iter.edges.nextSelected(); edgeNum != -1;
         edgeNum = iter.edges.nextSelected()) {
        int numPoints;
        const McVec3f* pts = mGraph->getEdgePoints(edgeNum, numPoints);
        for (int i = 0; i < numPoints; ++i) {
            if (pts[i].z < min) {
                min = pts[i].z;
            }
            if (pts[i].z > max) {
                max = pts[i].z;
            }
        }
    }

    for (int i = 0; i < sel.getNumSelectedVertices(); i++) {
        const McVec3f v = mGraph->getVertexCoords(sel.getSelectedVertex(i));
        if (v.z < min) {
            min = v.z;
        }
        if (v.z > max) {
            max = v.z;
        }
    }

    ma::MinMax mm = { min, max };
    return mm;
}

ma::MinMax ma::SliceSelector::getZRangeByIndex(const int sliceIndex) const {
    SpatialGraphSelection sel;
    getSlice(getSliceAttributeValueFromIndex(sliceIndex), sel);
    return getZRange(sel);
}

static float midPlane(const ma::MinMax mm1, const ma::MinMax mm2) {
    if (fabs(mm2.max - mm1.min) > fabs(mm1.max - mm2.min)) {
        // Slice 2 above 1.
        return (0.5 * (mm2.min + mm1.max));
    } else {
        // Slice 1 above 2.
        return (0.5 * (mm1.min + mm2.max));
    }
}

float ma::SliceSelector::computeMidPlane(const int slice1Idx,
                                         const int slice2Idx) const {
    const ma::MinMax mm1 = getZRangeByIndex(slice1Idx);
    const ma::MinMax mm2 = getZRangeByIndex(slice2Idx);
    return midPlane(mm1, mm2);
}

bool ma::SliceSelector::isSliceAboveOther(const int slice1Idx,
                                          const int slice2Idx) const {
    float mid = computeMidPlane(slice1Idx, slice2Idx);
    const ma::MinMax mm = getZRangeByIndex(slice1Idx);
    return (mm.min + mm.max) / 2.0 > mid;
}

// From mid plane half way to slice boundary that is farther away.
static float halfwayDistance(const ma::MinMax mm, const float midPlane) {
    return 0.5 * std::max(fabs(midPlane - mm.min), fabs(midPlane - mm.max));
}

void ma::SliceSelector::selectAdjacentHalfSlices(
    const int slice1Attribute, const int slice2Attribute,
    SpatialGraphSelection& halfSlice1,
    SpatialGraphSelection& halfSlice2) const {
    getSlice(slice1Attribute, halfSlice1);
    getSlice(slice2Attribute, halfSlice2);
    const ma::MinMax mm1 = getZRange(halfSlice1);
    const ma::MinMax mm2 = getZRange(halfSlice2);
    const float mid = midPlane(mm1, mm2);
    const float maxDist1 = halfwayDistance(mm1, mid);
    deselectVerticesFartherAwayFromMidPlaneThan(halfSlice1, mid, maxDist1);
    const float maxDist2 = halfwayDistance(mm2, mid);
    deselectVerticesFartherAwayFromMidPlaneThan(halfSlice2, mid, maxDist2);
}

void ma::SliceSelector::deselectVerticesFartherAwayFromMidPlaneThan(
    SpatialGraphSelection& slice, const float zMidPlane,
    const float maxDist) const {
    for (int i = slice.getNumSelectedVertices() - 1; i >= 0; i--) {
        int vertex = slice.getSelectedVertex(i);
        if (fabs(mGraph->getVertexCoords(vertex).z - zMidPlane) > maxDist)
            slice.deselectVertex(slice.getSelectedVertex(i));
    }
}

static SpatialGraphSelection unionSel(SpatialGraphSelection first,
                                      const SpatialGraphSelection& second) {
    first.addSelection(second);
    return first;
}

SpatialGraphSelection
ma::SliceSelector::selectCloseToMidplane(const int refNum,
                                         const int transNum) const {
    const int refAttr = getSliceAttributeValueFromIndex(refNum);
    const int transAttr = getSliceAttributeValueFromIndex(transNum);
    SpatialGraphSelection ref;
    SpatialGraphSelection trans;
    selectAdjacentHalfSlices(refAttr, transAttr, ref, trans);
    return unionSel(ref, trans);
}
