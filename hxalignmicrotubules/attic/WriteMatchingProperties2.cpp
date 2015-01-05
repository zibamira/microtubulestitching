#include <hxalignmicrotubules/attic/WriteMatchingProperties2.h>

#include <QString>

#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxspatialgraph/SpatialGraphSelection.h>
#include <mclib/McException.h>

#include <hxalignmicrotubules/mtalign/PGMMatcher.h>

#include <hxalignmicrotubules/mtalign.h>

namespace ma = mtalign;

HX_INIT_CLASS(WriteMatchingProperties2, HxCompModule);

WriteMatchingProperties2::WriteMatchingProperties2()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portDistanceThreshold(this, "dist"),
      portProjectedDistanceThreshold(this, "projDist"),
      portAngleThreshold(this, "angle"),
      portAction(this, "action") {
    mSelectionHelper = NULL;
    mOutputNamesMap.insert(std::pair<std::string, int>("3dDistance", 0));
    mOutputNamesMap.insert(std::pair<std::string, int>("ProjectedDistance", 1));
    mOutputNamesMap.insert(std::pair<std::string, int>("MutualAngles", 2));
    mOutputNamesMap.insert(std::pair<std::string, int>("PairwiseShiftDiff", 3));
    portDistanceThreshold.setValue(1000);
    portProjectedDistanceThreshold.setValue(600);
    portAngleThreshold.setValue(30);
}

WriteMatchingProperties2::~WriteMatchingProperties2() {
    delete mSelectionHelper;
}

void WriteMatchingProperties2::initPairMap() {
    const char* attr = "UserDefinedMatchings";
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    mPairMap.resize(graph->getNumVertices());
    for (unsigned int i = 0; i < mPairMap.size(); i++)
        mPairMap[i] = -1;

    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, attr));
    if (!att) {
        const QString msg =
            QString("Could not find node attribute '%1'.").arg(attr);
#ifdef HX_AMIRA5_COMPAT
        mcthrow(qPrintable(msg));
#else
        mcthrow(msg);
#endif
    }
    std::map<int, int> labelToVertexMap;

    for (int i = 0; i < graph->getNumVertices(); i++) {
        int pairValue = att->getIntDataAtIdx(i);
        if (pairValue > 0) {
            std::map<int, int>::iterator other =
                labelToVertexMap.find(pairValue);
            if (other != labelToVertexMap.end()) {
                mPairMap[i] = other->second;
                mPairMap[other->second] = i;
            } else {
                labelToVertexMap.insert(std::pair<int, int>(pairValue, i));
            }
        }
    }
}

void WriteMatchingProperties2::compute() {
    mtalign::PGMPairWeightsParams config;
    createWeightConfig(config);

    if (!portAction.wasHit())
        return;

    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    if (graph == NULL)
        return;

    initPairMap();
    McDArray<McVec3f> refSliceCoordsToMatch(1);
    McDArray<McVec3f> transSliceCoordsToMatch(1);
    McDArray<McVec3f> refDirections;
    McDArray<McVec3f> transDirections;
    SpatialGraphSelection selection(graph);
    for (unsigned int i = 0; i < mPairMap.size(); i++) {
        if (mPairMap[i] != -1) {
            std::cout << "\nPair " << i << " " << mPairMap[i];

            selection.clear();
            selection.selectVertex(i);
            getDirections(selection, refDirections);
            selection.clear();
            selection.selectVertex(mPairMap[i]);
            getDirections(selection, transDirections);
            refSliceCoordsToMatch[0] = graph->getVertexCoords(i);
            refSliceCoordsToMatch[0].z = 0;  // ortho projection
            transSliceCoordsToMatch[0] = graph->getVertexCoords(mPairMap[i]);
            transSliceCoordsToMatch[0].z = 0;  // ortho projection
            ma::PGMPairWeights weighter(refSliceCoordsToMatch,
                                        transSliceCoordsToMatch, refDirections,
                                        transDirections, config);
            writeSingleProbParameters(0, 0, weighter);
        }
    }
    mSelectionHelper = new ma::SliceSelector(
        hxconnection_cast<HxSpatialGraph>(portData), "TransformInfo");

    writeParameters();
}

void WriteMatchingProperties2::writeParametersForTwoConsecutiveSlices(
    const int slice1Idx, const int slice2Idx) {

    float midplane = mSelectionHelper->computeMidPlane(slice1Idx, slice2Idx);
    const HxSpatialGraph* sg = hxconnection_cast<HxSpatialGraph>(portData);
    ma::EndPointParams params;
    params.refSliceNum = slice1Idx;
    params.transSliceNum = slice2Idx;
    params.endPointRegion = 50;
    params.projectionType = "Orthogonal";
    params.projectionPlane = midplane;
    params.useAbsouteValueForEndPointRegion = false;
    params.maxDistForAngle = 2000;
    params.transformType = 0;
    SpatialGraphSelection slice1, slice2;
    ma::FacingPointSets pr = ma::projectEndPoints(sg, slice1, slice2, params);

    simulateAlignPairFunction(pr, slice1, slice2);
}

void WriteMatchingProperties2::writeSingleProbParameters(
    const int variableIdx, const int assignmentIdx,
    const ma::PGMPairWeights& weighter) {
    float threeDDistance = weighter.get3dDistance(variableIdx, assignmentIdx);
    addResultEntry("3dDistance", threeDDistance);
    float angle = weighter.getAngle(variableIdx, assignmentIdx);
    addResultEntry("MutualAngles", angle);
    float projectedDistance =
        weighter.getProjectedDistance(variableIdx, assignmentIdx);
    addResultEntry("ProjectedDistance", projectedDistance);
}

void WriteMatchingProperties2::addResultEntry(const std::string& whichTable,
                                              const float value) {
    SpreadSheetWrapper* spreadSheet = getResultSpreadSheet(whichTable);
    McString row;
    row.printf("%d", spreadSheet->getNumRows(whichTable.c_str()));
    spreadSheet->addEntry(whichTable.c_str(), whichTable.c_str(), row, value);
}

void WriteMatchingProperties2::simulateAlignPairFunction(
    const ma::FacingPointSets& pr, SpatialGraphSelection& refSelection,
    SpatialGraphSelection& transSelection) {

    McDArray<int> refSlicePointsToMatch;
    McDArray<int> transSlicePointsToMatch;
    McDArray<McVec3f> refSliceCoordsToMatch;
    McDArray<McVec3f> transSliceCoordsToMatch;
    McDArray<McVec3f> refDirections, transDirections;

    // get the indexes of the vertices
    for (int i = 0; i < refSelection.getNumSelectedVertices(); i++)
        refSlicePointsToMatch.append(refSelection.getSelectedVertex(i));
    for (int i = 0; i < transSelection.getNumSelectedVertices(); i++)
        transSlicePointsToMatch.append(transSelection.getSelectedVertex(i));
    refSliceCoordsToMatch = pr.ref.positions;
    transSliceCoordsToMatch = pr.trans.positions;

    // get the directions
    refDirections = pr.ref.directions;
    transDirections = pr.trans.directions;
    McDArray<McVec2i> evidenceDummy;

    mtalign::PGMPairWeightsParams config;
    createWeightConfig(config);
    ma::PGMPairWeights weighter(refSliceCoordsToMatch, transSliceCoordsToMatch,
                                refDirections, transDirections, config);
    mtalign::PGMMatcher bfom(refSliceCoordsToMatch, transSliceCoordsToMatch,
                             refDirections, transDirections, evidenceDummy,
                             weighter, 200);

    bfom.initConnectedComponentIteration();
    mtalign::ConnectedFactorGraph graph;
    while (bfom.getNextConnectedComponent(graph)) {
        theMsg->printf("conncomp");
        for (int i = 0; i < graph.variables.size(); i++) {
            int curVar = graph.variables[i];
            // first, get the assigned node
            int otherVertex = getMatchedVertex(refSlicePointsToMatch[curVar]);
            if (otherVertex < 0)
                continue;

            // find the index for the weight computer.
            int otherVertexAssignmentIndex = -1;
            for (int k = 0; k < transSlicePointsToMatch.size(); k++) {
                if (transSlicePointsToMatch[k] == otherVertex)
                    otherVertexAssignmentIndex = k;
            }
            // if we could not find the other vertex, the user defined matchings
            // might contain pairs that do not fit the
            // parameters.
            if (otherVertexAssignmentIndex == -1) {
                continue;
            }
            // writeSingleProbParameters(curVar, otherVertexAssignmentIndex,
            // weighter);

            for (int j = i + 1; j < graph.variables.size(); j++) {
                if (graph.adjacencyMatrix[graph.variables[i]]
                                         [graph.variables[j]]) {
                    // get the vertex indices of the two variables and their
                    // assignments
                    // we need them to get the coordinates
                    int vertexIdxVar1 = refSlicePointsToMatch[curVar];
                    int vertexIdxVar2 =
                        refSlicePointsToMatch[graph.variables[j]];
                    int matchedIdx1 = getMatchedVertex(vertexIdxVar1);
                    int matchedIdx2 = getMatchedVertex(vertexIdxVar2);
                    if (matchedIdx1 != -1 && matchedIdx2 != -1) {
                        // find index of the assignments in the array of
                        // assignments of the point matching function
                        int matchedIdx1InArray = -1;
                        int matchedIdx2InArray = -1;
                        for (int k = 0; k < transSliceCoordsToMatch.size();
                             k++) {
                            if (transSlicePointsToMatch[k] == matchedIdx2)
                                matchedIdx2InArray = k;
                            if (transSlicePointsToMatch[k] == matchedIdx1)
                                matchedIdx1InArray = k;
                        }
                        mcassert(matchedIdx1InArray > -1);
                        mcassert(matchedIdx2InArray > -1);
                        if (matchedIdx2InArray < 0 || matchedIdx1InArray < 0)
                            continue;
                        McVec3f var1Coord = weighter.mCoords1[curVar];
                        McVec3f var2Coord =
                            weighter.mCoords1[graph.variables[j]];
                        McVec3f matchedVar2Coord =
                            weighter.mCoords2[matchedIdx2InArray];
                        McVec3f matchedVar1Coord =
                            weighter.mCoords2[matchedIdx1InArray];
                        float shiftDistance =
                            mtalign::PGMMatcher::computeShiftDistance(
                                var1Coord, var2Coord, matchedVar1Coord,
                                matchedVar2Coord);
                        addResultEntry("PairwiseShiftDiff", shiftDistance);
                    }
                }
            }
        }
    }
}

void WriteMatchingProperties2::writeParameters() {
    for (int i = 0; i < mSelectionHelper->getNumSlices() - 1; i++) {
        theMsg->printf("Write for slices %d %d", i, i + 1);
        writeParametersForTwoConsecutiveSlices(i, i + 1);
    }
}
void WriteMatchingProperties2::createWeightConfig(
    mtalign::PGMPairWeightsParams& config) {
    config.distanceThreshold3d = portDistanceThreshold.getValue();
    config.useDistanceThreshold3d = true;
    config.distanceThresholdProjected =
        portProjectedDistanceThreshold.getValue();
    config.useDistanceThresholdProjected = true;

    config.angleThreshold = portAngleThreshold.getValue();
    config.useAngleThreshold = true;

    config.angleWeightParam = 10;
    config.useAngleWeight = true;
    config.dist3dParam = 600;
    config.useDist3dWeight = true;
    config.distProjectedParam = 200;
    config.useProjectedDistWeight = true;
    config.distProjectedParam = 200;
    config.dummySignificance = 0.01;
    config.weightType = mtalign::PGMPairWeightsParams::EXPONENTIAL;
}

int WriteMatchingProperties2::getMatchedVertex(const int vertex) {
    return mPairMap[vertex];
}

SpreadSheetWrapper*
WriteMatchingProperties2::getResultSpreadSheet(const std::string& outputName) {
    SpreadSheetWrapper* resultSpreadSheet = dynamic_cast<SpreadSheetWrapper*>(
        getResult(mOutputNamesMap[outputName]));
    if (resultSpreadSheet == NULL) {
        resultSpreadSheet = new SpreadSheetWrapper();
        setResult(mOutputNamesMap[outputName], resultSpreadSheet);
        resultSpreadSheet->setLabel(outputName.c_str());
    }
    return resultSpreadSheet;
}

McVec3f
WriteMatchingProperties2::getNextPointAtDist(const McDArray<McVec3f> edgePoints,
                                             const float maxDist) {

    float dist = 0;
    int actIndex = 1;

    while ((dist < maxDist) && (actIndex < edgePoints.size())) {
        dist += (edgePoints[actIndex - 1] - edgePoints[actIndex]).length();
        actIndex++;
    }
    // if the max Dist was not reached at all
    if (dist <= maxDist) {
        return edgePoints[actIndex - 1];
    }
    // we interpolate the last point
    else {
        // get difference
        float diffDist = dist - maxDist;
        McVec3f x2 = edgePoints[actIndex - 1];
        McVec3f x1 = edgePoints[actIndex - 2];
        float distLastPoints = (x2 - x1).length();
        float frac = diffDist / distLastPoints;
        mcassert(frac < 1);

        return (1 - frac) * x2 + frac * x1;
    }
}

void WriteMatchingProperties2::getDirections(
    const SpatialGraphSelection& selectedVertices,
    McDArray<McVec3f>& directions) {

    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    directions.resize(selectedVertices.getNumSelectedVertices());
    for (int i = 0; i < selectedVertices.getNumSelectedVertices(); i++) {
        int selectedVertex = selectedVertices.getSelectedVertex(i);
        IncidenceList connectedEdges = graph->getIncidentEdges(selectedVertex);

        McVec3f point1 = graph->getVertexCoords(selectedVertex);
        McVec3f point2;

        int connectedEdge = connectedEdges[0];
        int source = graph->getEdgeSource(connectedEdge);
        McDArray<McVec3f> edgePoints = graph->getEdgePoints(connectedEdge);
        if (source != selectedVertex) {
            edgePoints.reverse();
        }

        point2 = getNextPointAtDist(edgePoints, 2000);
        McVec3f direction = point2 - point1;
        direction.normalize();
        directions[i] = direction;
    }
}
