#include <hxalignmicrotubules/MovingLeastSquares.h>

#include <hxfield/HxUniformVectorField3.h>
#include <hxspatialgraph/HxSpatialGraph.h>

#include <hxalignmicrotubules/mtalign/SliceSelector.h>
#include <hxalignmicrotubules/HxMovingLeastSquaresSpatialGraphWarp.h>

HX_INIT_CLASS(HxMovingLeastSquaresSpatialGraphWarp, HxCompModule);

HxMovingLeastSquaresSpatialGraphWarp::HxMovingLeastSquaresSpatialGraphWarp()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portMethod(this, "method", 1),
      portAlpha(this, "parameter"),
      portAction(this, "action") {
    portAction.setAliasName("doIt");
    portMethod.setLabel(0, "Rigid");
    portAlpha.setLabel("Alpha:");
    portAlpha.setMinMax(0, 10);
    portAlpha.setValue(4);
}

HxMovingLeastSquaresSpatialGraphWarp::~HxMovingLeastSquaresSpatialGraphWarp() {}

void HxMovingLeastSquaresSpatialGraphWarp::update() {}

HxUniformVectorField3*
HxMovingLeastSquaresSpatialGraphWarp::createOutputVectorDataSet() {
    HxSpatialGraph* sg = dynamic_cast<HxSpatialGraph*>(portData.source());
    HxUniformVectorField3* output =
        dynamic_cast<HxUniformVectorField3*>(getResult(1));

    if (!sg)
        return (0);

    int dims[3];
    float bbox[6];

    if (sg) {
        sg->getBoundingBox(bbox);
        dims[0] = (bbox[1] - bbox[0]) / 400.0;
        dims[1] = (bbox[3] - bbox[2]) / 400.0;
        dims[2] = 1;
    }
    if (!output || output->lattice.dimsInt()[0] != dims[0] ||
        output->lattice.dimsInt()[1] != dims[1] ||
        output->lattice.dimsInt()[2] != dims[2])
        output = 0;

    if (!output) {
        output = new HxUniformVectorField3(dims, McPrimType::mc_float);
    }
    output->lattice.setBoundingBox(bbox);
    // memcpy(output->bbox(),bbox,6*sizeof(float));
    output->composeLabel(sg->getLabel().getString(), "Displacement");
    setResult(1, output);
    return (output);
}

HxSpatialGraph* HxMovingLeastSquaresSpatialGraphWarp::createOutputDataSet() {
    HxSpatialGraph* inputSG = dynamic_cast<HxSpatialGraph*>(portData.source());
    HxSpatialGraph* warpedSG = dynamic_cast<HxSpatialGraph*>(getResult(0));

    if (!inputSG)
        return (0);
    if (!warpedSG)
        warpedSG = new HxSpatialGraph();
    else
        warpedSG->clear();

    warpedSG->composeLabel(inputSG->getLabel().getString(), "Warped");
    setResult(0, warpedSG);
    return (warpedSG);
}

void HxMovingLeastSquaresSpatialGraphWarp::prepareLandmarks(
    McDArray<McVec2d>& p1, McDArray<McVec2d>& p2,
    const HxSpatialGraph* spatialGraph, const int slice1Index,
    const int slice2Index) {

    mtalign::SliceSelector selectionHelper(spatialGraph, "TransformInfo");
    // find all existing pairs
    std::map<int, std::vector<int> > sameValueGroups;
    getSameValueVertices(sameValueGroups, spatialGraph, "WarpPairs");

    p1.resize(0);
    p2.resize(0);
    std::map<int, std::vector<int> >::iterator groupIterator =
        sameValueGroups.begin();

    for (; groupIterator != sameValueGroups.end(); groupIterator++) {
        if (isValidPair(groupIterator->second, selectionHelper)) {
            McVec3f p1Coord =
                spatialGraph->getVertexCoords(groupIterator->second.front());
            McVec3f p2Coord =
                spatialGraph->getVertexCoords(groupIterator->second.back());
            int trafoP1 = selectionHelper.getSliceIdxOfVertex(
                groupIterator->second.front());
            int trafoP2 = selectionHelper.getSliceIdxOfVertex(
                groupIterator->second.back());
            if (trafoP1 == slice1Index && trafoP2 == slice2Index) {
                p1.append(McVec2d(p1Coord.x, p1Coord.y));
                p2.append(McVec2d(p2Coord.x, p2Coord.y));
            } else if (trafoP1 == slice2Index && trafoP2 == slice1Index) {
                p2.append(McVec2d(p1Coord.x, p1Coord.y));
                p1.append(McVec2d(p2Coord.x, p2Coord.y));
            }
        }
    }
}

void HxMovingLeastSquaresSpatialGraphWarp::getSameValueVertices(
    std::map<int, std::vector<int> >& sameValueVertices,
    const HxSpatialGraph* spatialGraph, const std::string& attributeName) {
    EdgeVertexAttribute* attribute = dynamic_cast<EdgeVertexAttribute*>(
        spatialGraph->findVertexAttribute(McString(attributeName.c_str())));
    if (attribute == 0)
        return;
    // find mapping from unsorted labels to new labels
    for (int i = 0; i < spatialGraph->getNumVertices(); ++i) {
        int pairNum = attribute->getIntDataAtIdx(i);
        if (sameValueVertices.find(pairNum) != sameValueVertices.end()) {
            sameValueVertices.find(pairNum)->second.push_back(i);
        } else {
            std::pair<int, std::vector<int> > newGroup;
            newGroup.first = pairNum;
            newGroup.second.push_back(i);
            sameValueVertices.insert(newGroup);
        }
    }
}

bool HxMovingLeastSquaresSpatialGraphWarp::isValidPair(
    const std::vector<int>& vertexIndices,
    const mtalign::SliceSelector& selectionHelper) {
    if (vertexIndices.size() != 2)
        return false;
    int sliceIdx1 = selectionHelper.getSliceIdxOfVertex(vertexIndices.front());
    int sliceIdx2 = selectionHelper.getSliceIdxOfVertex(vertexIndices.back());
    return abs(sliceIdx1 - sliceIdx2) == 1;
}

void HxMovingLeastSquaresSpatialGraphWarp::compute() {
    if (!portAction.wasHit())
        return;

    // Der Volumendatensatz
    HxSpatialGraph* inputSG = (HxSpatialGraph*)portData.source();

    if (!inputSG)
        return;

    mtalign::SliceSelector selectionHelper(inputSG, "TransformInfo");
    std::vector<TransformationOrderInfo> orderOfDeformationComputation;
    generateOrderOfDeformationSimple(selectionHelper,
                                     orderOfDeformationComputation);

    HxSpatialGraph* resultSG = createOutputDataSet();
    resultSG->copyFrom(inputSG);
    std::vector<TransformationOrderInfo>::iterator transformationIterator =
        orderOfDeformationComputation.begin();
    for (; transformationIterator != orderOfDeformationComputation.end();
         transformationIterator++) {
        applyDeformation(*transformationIterator, resultSG);
    }

    resultSG->touch();
    resultSG->fire();
}

void HxMovingLeastSquaresSpatialGraphWarp::applyDeformation(
    TransformationOrderInfo& transformation, HxSpatialGraph* spatialGraph) {
    mtalign::SliceSelector selectionHelper(spatialGraph, "TransformInfo");
    MovingLeastSquares lowerTransform, upperTransform;
    generateTransformsForSlicesSimple(
        spatialGraph, transformation.lowerSliceIdx,
        transformation.upperSliceIdx, lowerTransform, upperTransform,
        portAlpha.getValue());
    for (std::vector<int>::iterator sliceIterator =
             transformation.lowerTransformSlices.begin();
         sliceIterator != transformation.lowerTransformSlices.end();
         sliceIterator++) {
        int currentSliceIndex = *sliceIterator;
        SpatialGraphSelection currentSlice;
        selectionHelper.getSlice(
            selectionHelper.getSliceAttributeValueFromIndex(currentSliceIndex),
            currentSlice);
        applyDeformationToSlice(currentSlice, spatialGraph, lowerTransform);
    }
    for (std::vector<int>::iterator sliceIterator =
             transformation.upperTransformSlices.begin();
         sliceIterator != transformation.upperTransformSlices.end();
         sliceIterator++) {
        int currentSliceIndex = *sliceIterator;
        SpatialGraphSelection currentSlice;
        selectionHelper.getSlice(
            selectionHelper.getSliceAttributeValueFromIndex(currentSliceIndex),
            currentSlice);
        applyDeformationToSlice(currentSlice, spatialGraph, upperTransform);
    }
}

void HxMovingLeastSquaresSpatialGraphWarp::applyDeformationToSlice(
    SpatialGraphSelection& slice, HxSpatialGraph* spatialGraph,
    MovingLeastSquares& deformation) {
    for (int i = 0; i < slice.getNumSelectedEdges(); i++) {
        int currentEdge = slice.getSelectedEdge(i);
        McDArray<McVec3f> currentEdgePoints =
            spatialGraph->getEdgePoints(currentEdge);
        McDArray<McVec3f> warpedEdgePoints = currentEdgePoints;
        for (int j = 0; j < warpedEdgePoints.size(); j++) {
            warpPoint(currentEdgePoints[j], warpedEdgePoints[j], deformation);
        }
        spatialGraph->setEdgePoints(currentEdge, warpedEdgePoints);
        // move the vertices as well
        int sourceIdx = spatialGraph->getEdgeSource(currentEdge);
        McVec3f sourceCoord = spatialGraph->getVertexCoords(sourceIdx);
        int targetIdx = spatialGraph->getEdgeTarget(currentEdge);
        McVec3f targetCoord = spatialGraph->getVertexCoords(targetIdx);
        McVec3f warpedSource, warpedTarget;
        warpPoint(sourceCoord, warpedSource, deformation);
        warpPoint(targetCoord, warpedTarget, deformation);
        spatialGraph->setVertexCoords(sourceIdx, warpedSource);
        spatialGraph->setVertexCoords(targetIdx, warpedTarget);
    }
}

void HxMovingLeastSquaresSpatialGraphWarp::generateTransformsForSlicesSimple(
    const HxSpatialGraph* spatialGraph, const int slice1Index,
    const int slice2Index, MovingLeastSquares& lowerTransform,
    MovingLeastSquares& upperTransform, int alpha) {
    McDArray<McVec2d> lm1, lm2;
    prepareLandmarks(lm1, lm2, spatialGraph, slice1Index, slice2Index);
    // Clear lower transform.
    lowerTransform = MovingLeastSquares();
    lowerTransform.setAlpha(alpha);
    upperTransform.setAlpha(alpha);
    upperTransform.setLandmarks(lm2, lm1);
}

void HxMovingLeastSquaresSpatialGraphWarp::generateOrderOfDeformationSimple(
    const mtalign::SliceSelector& selectionHelper,
    std::vector<TransformationOrderInfo>& orderOfDeformation) {

    int numSlices = selectionHelper.getNumSlices();
    for (int i = 0; i < numSlices - 1; i++) {
        TransformationOrderInfo tranformationInfo;
        tranformationInfo.lowerSliceIdx = i;
        tranformationInfo.upperSliceIdx = i + 1;
        tranformationInfo.lowerTransformSlices.push_back(i);
        for (int j = i + 1; j < numSlices; j++) {
            tranformationInfo.upperTransformSlices.push_back(j);
        }

        orderOfDeformation.push_back(tranformationInfo);
    }
}

void HxMovingLeastSquaresSpatialGraphWarp::warpPoint(
    const McVec3f& source, McVec3f& target,
    MovingLeastSquares& mlsInterpolator) {

    McVec2d curPoint = McVec2d(source.x, source.y);
    McVec2d warpedPoint = mlsInterpolator.interpolate(curPoint);
    target = McVec3f(warpedPoint.x, warpedPoint.y, source.z);
}
