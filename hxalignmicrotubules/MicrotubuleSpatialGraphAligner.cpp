#include <hxalignmicrotubules/MicrotubuleSpatialGraphAligner.h>

#include <limits>

#include <QTime>

#include <hxcore/HxObjectPool.h>
#include <hxcore/HxWorkArea.h>
#include <hxneuroneditor/HxNeuronEditorSubApp.h>
#include <hxspatialgraph/EdgeVertexAttribute.h>
#include <hxspatialgraph/HierarchicalLabels.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxspatialgraph/SpatialGraphOperationSet.h>
#include <mclib/McDArray.h>
#include <mclib/McString.h>

#include <hxalignmicrotubules/MicrotubuleTransformOperation.h>
#include <hxalignmicrotubules/MovingLeastSquares.h>
#include <hxalignmicrotubules/SpreadSheetWrapper.h>
#include <hxalignmicrotubules/mtalign.h>

namespace ma = mtalign;

static mtalign::PGMPairWeightsParams
makeWeightConfigForInitMatching(const ma::MatchingParams& params) {
    mtalign::PGMPairWeightsParams cfg;
    cfg.useDistanceThreshold3d = true;
    cfg.distanceThreshold3d = params.maxDistanceForGraphConstruction;
    cfg.useDistanceThresholdProjected = false;
    cfg.angleThreshold = params.maxAngleDiffForInitMatching;
    cfg.useAngleThreshold = true;
    cfg.useAngleWeight = false;
    cfg.dist3dParam = params.maxDistanceForGraphConstruction;
    cfg.useDist3dWeight = true;
    cfg.useProjectedDistWeight = false;
    cfg.weightType = mtalign::PGMPairWeightsParams::EXPONENTIAL;
    return cfg;
}

const int MicrotubuleSpatialGraphAligner::NumProjectionTypes = 9;
const McString MicrotubuleSpatialGraphAligner::ProjectionTypes[] = {
    McString("Orthogonal"),   McString("Linear"),
    McString("Tangent"),      McString("Fit degree 0"),
    McString("Fit degree 1"), McString("Fit degree 2"),
    McString("Fit degree 3"), McString("ApproxTangent theshold"),
    McString("None")
};

static ma::ProjectionType asEnumProjectionType(const McString& s) {
    const int nTypes = MicrotubuleSpatialGraphAligner::NumProjectionTypes;
    for (int i = 0; i < nTypes; i++) {
        if (s == MicrotubuleSpatialGraphAligner::ProjectionTypes[i]) {
            return ma::ProjectionType(i);
        }
    }
    return ma::P_NONE;
}

const int MicrotubuleSpatialGraphAligner::NumTransformTypes = 4;
const McString MicrotubuleSpatialGraphAligner::TransformTypes[] = {
    McString("Rigid"),  McString("Rigid + iso-scale"),
    McString("Affine"), McString("None"),
};

static ma::TransformType asTransformType(McString t) {
    if (t == McString("None")) {
        return ma::TF_NONE;
    } else if (t == McString("Rigid")) {
        return ma::TF_RIGID;
    } else if (t == McString("Rigid + iso-scale")) {
        return ma::TF_RIGID_ISO_SCALE;
    } else if (t == McString("Affine")) {
        return ma::TF_AFFINE;
    } else {
        mcerror("Invalid transform type.");
    }
    return ma::TF_NONE;
}

const int MicrotubuleSpatialGraphAligner::NumPointMatchingAlgorithmTypes = 3;
const McString MicrotubuleSpatialGraphAligner::PointMatchingAlgorithmTypes[] = {
    McString("Greedy"), McString("Exact"), McString("PGM")
};

MicrotubuleSpatialGraphAligner::MicrotubuleSpatialGraphAligner(
    HxSpatialGraph* graph)
    : mSelectionHelper(graph, "TransformInfo"),
      mGraph(graph),
      mMaxDistanceDifference(0.0f),
      mMinCliqueSizeFraction(0.0f),
      mAlpha(0.0001f),
      mNumMaxPointsForInitialTransform(50),
      mEndPointRegion(0.0f),
      mCreateMatchingLabels(true),
      mAutoEstimateGapSize(0) {
    mMaxAngleDiffForInitMatching = 180;
    mNotUsedLabelValue = 1;
    mNotMatchedLabelValue = 2;
}

MicrotubuleSpatialGraphAligner::~MicrotubuleSpatialGraphAligner() {}

int MicrotubuleSpatialGraphAligner::getNumSlices() const {
    return mSelectionHelper.getNumSlices();
}

McDArray<McString>
MicrotubuleSpatialGraphAligner::getPossibleAlignmentAttributeNames() const {
    McDArray<McString> result(0);
    if (hasTransformInfo()) {
        result.append(getTransformInfoAttribute());
    } else {
        result.append(McString("Error: no alignment attribute found"));
    }
    return result;
}

McDArray<McString>
MicrotubuleSpatialGraphAligner::getPossiblePointMatchingAlgorithms() {
    McDArray<McString> result(0);
    for (int i = 0; i < NumPointMatchingAlgorithmTypes; ++i) {
        result.append(PointMatchingAlgorithmTypes[i]);
    }
    return result;
}

McDArray<McString>
MicrotubuleSpatialGraphAligner::getPossibleProjectionTypes() {
    McDArray<McString> result(0);
    for (int i = 0; i < NumProjectionTypes; ++i) {
        result.append(ProjectionTypes[i]);
    }
    return result;
}

McDArray<McString> MicrotubuleSpatialGraphAligner::getPossibleTransformTypes() {
    McDArray<McString> result(0);
    for (int i = 0; i < NumTransformTypes; ++i) {
        result.append(TransformTypes[i]);
    }
    return result;
}

void
MicrotubuleSpatialGraphAligner::setAutoEstimateGapSize(const int autoEstimate) {
    mAutoEstimateGapSize = autoEstimate;
}

int MicrotubuleSpatialGraphAligner::getAutoEstimateGapSize() const {
    return mAutoEstimateGapSize;
}

void
MicrotubuleSpatialGraphAligner::setEndPointRegion(const float slicePercentage) {
    mEndPointRegion = slicePercentage;
}

void MicrotubuleSpatialGraphAligner::setMaxDistanceDifference(
    const float maxDistDiff) {
    mMaxDistanceDifference = maxDistDiff;
}

void MicrotubuleSpatialGraphAligner::setMinCliqueSizeFraction(
    const float minCliqueSizeFraction) {
    mMinCliqueSizeFraction = minCliqueSizeFraction;
}

void MicrotubuleSpatialGraphAligner::setMaxNumCliques(const int maxCliques) {
    mMaxNumCliques = maxCliques;
}

void MicrotubuleSpatialGraphAligner::setAlpha(const float alpha) {
    mAlpha = alpha;
}

void MicrotubuleSpatialGraphAligner::setPointMatchingAlgorithm(
    const McString pmAlg) {
    mPMAlgorithm = pmAlg;
}

void MicrotubuleSpatialGraphAligner::setGapSize(const float gapSize) {
    mGapSize = gapSize;
}

void MicrotubuleSpatialGraphAligner::setDMin(const float dMin) { mDMin = dMin; }

void MicrotubuleSpatialGraphAligner::setDMax(const float dMax) { mDMax = dMax; }

void MicrotubuleSpatialGraphAligner::setDeltaH(const float deltaH) {
    mDeltaH = deltaH;
}

void MicrotubuleSpatialGraphAligner::setProjectionType(const McString pType) {
    mProjectionType = pType;
}

void MicrotubuleSpatialGraphAligner::setTransformType(const McString tType) {
    mTransformType = tType;
}

void
MicrotubuleSpatialGraphAligner::setCreateMatchingLabels(const bool enabled) {
    mCreateMatchingLabels = enabled;
}

void MicrotubuleSpatialGraphAligner::setScaleTestMin(const float scaleTestMin) {
    mScaleTestMin = scaleTestMin;
}

void MicrotubuleSpatialGraphAligner::setScaleTestMax(const float scaleTestMax) {
    mScaleTestMax = scaleTestMax;
}

void MicrotubuleSpatialGraphAligner::setScaleTestIncrement(
    const float scaleTestIncrement) {
    mScaleTestIncrement = scaleTestIncrement;
}

void MicrotubuleSpatialGraphAligner::setMaxAngleDiffForInitMatching(
    const float maxAngleDiffForInitMatching) {
    mMaxAngleDiffForInitMatching = maxAngleDiffForInitMatching;
}

void MicrotubuleSpatialGraphAligner::setMaxDistForAngle(
    const float maxDistForAngle) {
    mMaxDistForAngle = maxDistForAngle;
}

void MicrotubuleSpatialGraphAligner::setUseAbsouteValueForEndPointRegion(
    const bool useAbsoluteValueForEndPointRegion) {
    mUseAbsoluteValueForEndPointRegion = useAbsoluteValueForEndPointRegion;
}

void
MicrotubuleSpatialGraphAligner::warpAll(HxNeuronEditorSubApp* editor,
                                        const ma::CPDParams& cpdParameters) {
    theWorkArea->startWorking("");
    McDArray<ma::WarpResult> warpResults;
    McHandle<SpreadSheetWrapper> spreadSheet = getEvaluationSpreadsheet();
    warpResults.resize(mSelectionHelper.getNumSlices() - 1);
    for (int i = 1; i < mSelectionHelper.getNumSlices(); ++i) {
        const QString info = QString("Warping slices %1-%2")
            .arg(mSelectionHelper.getSliceAttributeValueFromIndex(i - 1))
            .arg(mSelectionHelper.getSliceAttributeValueFromIndex(i));

#ifdef HX_AMIRA5_COMPAT
        theWorkArea->setProgressInfo("%s", qPrintable(info));
#else
        theWorkArea->setProgressInfo(info);
#endif
        theWorkArea->setProgressValue(
            float(i) / float(mSelectionHelper.getNumSlices() - 1));

        ma::WarpResult& result = warpResults[i - 1];
        ma::AlignInfo& alignInfo = result.alignInfo;
        warpSlices(i - 1, i, cpdParameters, result);
        McString row;
        row.printf("slices-%d-%d", i, i + 1);
        spreadSheet->addEntry("CPDResults", "sigmaSquare", row.dataPtr(),
                              alignInfo.sigmaSquare);
        spreadSheet->addEntry("CPDResults", "kappa", row.dataPtr(),
                              alignInfo.kappa);
        spreadSheet->addEntry("CPDResults", "numIter", row.dataPtr(),
                              alignInfo.numIterations);
        spreadSheet->addEntry("CPDResults", "E", row.dataPtr(), alignInfo.e);
        spreadSheet->addEntry("CPDResults", "EDiff", row.dataPtr(),
                              alignInfo.eDiffRel);
        spreadSheet->addEntry("CPDResults", "time", row.dataPtr(),
                              alignInfo.timeInSec);

        // Ignore transform if sigma^2 is too large.  Such transformations have
        // never been useful in practice.
        if (alignInfo.sigmaSquare <= cpdParameters.maxAcceptableSigmaSquare) {
            spreadSheet->addEntry("CPDResults", "note", row.dataPtr(),
                                  "transform applied");
        } else {
            McString reason(0, "transform ignored, sigmaSquare > %f",
                            cpdParameters.maxAcceptableSigmaSquare);
            spreadSheet->addEntry("CPDResults", "note", row.dataPtr(),
                                  reason.getString());

            result.transformMatrix = McMat4f::identity();

            // Create a moving-least-squares identity transform by mapping the
            // ps to itself, so that the downstream pipeline should work
            // without changes, since it looks like a regular MLS.
            result.mlsParams.qs = result.mlsParams.ps;
        }
    }

    if (cpdParameters.type != ma::CPD_ELASTIC) {
        for (int i = mSelectionHelper.getNumSlices() - 1; i >= 1; --i) {
            mcassert(warpResults[i - 1].type == ma::WT_LINEAR);
            applyTransform(warpResults[i - 1].transformMatrix, i, -1, editor);
        }
    } else {
        HxParamBundle* bundle =
            mGraph->parameters.bundle("CPDTransformLandmarks", false);
        if (bundle != 0) {
            theMsg->printf("No nl transform applied! Was transformed already.");
            theWorkArea->stopWorking();
            return;
        }
        HxParamBundle* pb = new HxParamBundle("CPDTransformLandmarks");
        pb->set("alpha", cpdParameters.alphaForMLS);
        mGraph->parameters.insert(pb);
        for (int i = mSelectionHelper.getNumSlices() - 1; i >= 1; --i) {
            mcassert(warpResults[i - 1].type == ma::WT_ELASTIC);
            addWarpPointsToParams(warpResults[i - 1].mlsParams, i);
            MovingLeastSquares mls;
            mls.setAlpha(warpResults[i - 1].mlsParams.alpha);
            mls.setLandmarks(warpResults[i - 1].mlsParams.ps,
                             warpResults[i - 1].mlsParams.qs);
            applyMLS(mls, i, editor);
        }
    }

    mGraph->touch();
    mGraph->fire();
    theWorkArea->stopWorking();
}

void MicrotubuleSpatialGraphAligner::addWarpPointsToParams(
    const ma::MLSParams& mlsParams, const int sliceNum) {

    HxParamBundle* bundle =
        mGraph->parameters.bundle("CPDTransformLandmarks", false);
    mcassert(bundle);
    McString sliceNums;
    sliceNums.printf("Slices%d-%d", sliceNum, sliceNum + 1);
    HxParamBundle* sliceLandmarks = new HxParamBundle(sliceNums.dataPtr());
    McString label;

    for (int i = 0; i < mlsParams.ps.size(); i++) {
        label.printf("PQ%d", i);
        float coord[4];
        coord[0] = mlsParams.ps[i].x;
        coord[1] = mlsParams.ps[i].y;
        coord[2] = mlsParams.qs[i].x;
        coord[3] = mlsParams.qs[i].y;
        sliceLandmarks->set(label.dataPtr(), 4, coord);
    }
    bundle->insert(sliceLandmarks);
}

void MicrotubuleSpatialGraphAligner::warpSlices(const int refSliceNum,
                                                const int transSliceNum,
                                                const ma::CPDParams& cpdParams,
                                                ma::WarpResult& deformation) {
    const float midPlane =
        mSelectionHelper.computeMidPlane(refSliceNum, transSliceNum);
    ma::EndPointParams params =
        makePointRepresentationParams(refSliceNum, transSliceNum, midPlane);
    const ma::FacingPointSets points = ma::projectEndPoints(mGraph, params);
    ma::cpd(points, deformation, cpdParams);
}

void MicrotubuleSpatialGraphAligner::applyMLS(MovingLeastSquares& mls,
                                              const int start,
                                              HxNeuronEditorSubApp* editor) {
    SpatialGraphSelection sel(mGraph);
    sel.clear();
    for (int i = start; i <= mSelectionHelper.getNumSlices() - 1; ++i) {
        SpatialGraphSelection tmpsel;
        mSelectionHelper.getSlice(
            mSelectionHelper.getSliceAttributeValueFromIndex(i), tmpsel);
        sel.addSelection(tmpsel);
    }
    applyMLSToSelection(mls, mGraph, sel);
}

void MicrotubuleSpatialGraphAligner::applyMLSToSelection(
    MovingLeastSquares& mls, HxSpatialGraph* sg,
    const SpatialGraphSelection& sel) {
    for (int i = 0; i < sel.getNumSelectedVertices(); i++) {
        McVec3f curCoord = sg->getVertexCoords(sel.getSelectedVertex(i));
        McVec2d warpedCoord = mls.interpolate(McVec2d(curCoord.x, curCoord.y));
        sg->setVertexCoords(sel.getSelectedVertex(i),
                            McVec3f(warpedCoord.x, warpedCoord.y, curCoord.z));
    }
    for (int i = 0; i < sel.getNumSelectedEdges(); i++) {
        int edge = sel.getSelectedEdge(i);
        McDArray<McVec3f> newEdgePoints;
        for (int j = 0; j < sg->getNumEdgePoints(edge); j++) {
            McVec3f curCoord = sg->getEdgePoint(edge, j);
            McVec2d warpedCoord =
                mls.interpolate(McVec2d(curCoord.x, curCoord.y));
            newEdgePoints.append(
                McVec3f(warpedCoord.x, warpedCoord.y, curCoord.z));
        }
        sg->setEdgePoints(edge, newEdgePoints);
    }
}

void
MicrotubuleSpatialGraphAligner::alignAllOrPair(HxNeuronEditorSubApp* editor,
                                               const int refSliceNum,
                                               const int transSliceNum) {

    theWorkArea->startWorking("");

    McDArray<McMat4f> matrices;
    McHandle<SpreadSheetWrapper> spreadSheet = getEvaluationSpreadsheet();

    for (int i = 1; i < mSelectionHelper.getNumSlices(); ++i) {
        McMat4f transMat = McMat4f::identity();
        if (transSliceNum == i || transSliceNum == -1) {
            const QString info = QString("Aligning slices %1-%2")
                .arg(mSelectionHelper.getSliceAttributeValueFromIndex(i - 1))
                .arg(mSelectionHelper.getSliceAttributeValueFromIndex(i));
#ifdef HX_AMIRA5_COMPAT
            theWorkArea->setProgressInfo("%s", qPrintable(info));
#else
            theWorkArea->setProgressInfo(info);
#endif
            theWorkArea->setProgressValue(
                float(i) / float(mSelectionHelper.getNumSlices() - 1));
            alignPairAndTestScaling(i - 1, i, transMat, editor, spreadSheet);
        }
        matrices.append(transMat);
    }
    editor->updateLabelTreeView();

    for (int i = mSelectionHelper.getNumSlices() - 1; i >= 1; --i) {
        applyTransform(matrices[i - 1], i, -1, editor);
    }

    if (editor) {
        HxNeuronEditorSubApp::SpatialGraphChanges changes =
            HxNeuronEditorSubApp::SpatialGraphLabelTreeChange |
            HxNeuronEditorSubApp::SpatialGraphAttributeValueChange;
        editor->updateAfterSpatialGraphChange(changes);
    }
    theWorkArea->stopWorking();
}

void MicrotubuleSpatialGraphAligner::alignPairAndTestScaling(
    const int refSliceNum, const int transSliceNum, McMat4f& transMat,
    HxNeuronEditorSubApp* editor, McHandle<SpreadSheetWrapper> spreadSheet) {

    float scaleFactor = mScaleTestMin;
    McMat4f bestMatrix = McMat4f::identity();
    float bestScaleFactor = 1;
    AlignmentResults bestAlignmentResults;
    bestAlignmentResults.score = 0;
    float bestGapVariation = 0;
    float midPlane =
        mSelectionHelper.computeMidPlane(refSliceNum, transSliceNum);

    for (; scaleFactor <= mScaleTestMax; scaleFactor += mScaleTestIncrement) {
        // apply x,y scale
        McMat4f scaleMatrix(scaleFactor, 0, 0, 0, 0, scaleFactor, 0, 0, 0, 0,
                            /*zScaling=*/1.0, 0, 0, 0, 0, 1.0);
        applyTransform(scaleMatrix, transSliceNum, -1);

        for (float varyGap = mDMin; varyGap <= mDMax; varyGap += mDeltaH) {
            QTime startTime = QTime::currentTime();

            // apply z-transform
            McMat4f transMatrix = McMat4f::identity();
            transMatrix.setTranslate(McVec3f(0, 0, varyGap));
            applyTransform(transMatrix, transSliceNum, 1);

            AlignmentInputParams inputParams(midPlane + 0.5 * varyGap,
                                             refSliceNum, transSliceNum);
            AlignmentResults alignmentResults;

            McMat4f resultMat = align(inputParams, alignmentResults, editor);

            theMsg->printf(" projectionPlane %f,score: %f",
                           midPlane + 0.5 * varyGap, alignmentResults.score);
            if ((alignmentResults.score > bestAlignmentResults.score) &&
                (alignmentResults.matchedRefPointIds.size() > 2)) {
                bestAlignmentResults = alignmentResults;
                bestMatrix = resultMat;
                bestGapVariation = varyGap;
                bestScaleFactor = scaleFactor;
            }

            // undo the translation
            transMatrix = McMat4f::identity();
            transMatrix.setTranslate(McVec3f(0, 0, -1 * varyGap));
            applyTransform(transMatrix, transSliceNum, 1);
            QTime endTime = QTime::currentTime();

            alignmentResults.duration_s = endTime.msecsTo(startTime) / 1000;
            writeResultsToSpreadSheet(inputParams, alignmentResults, varyGap,
                                      scaleFactor, spreadSheet);
        }
        // undo the scale
        McMat4f scaleMatrixInv = scaleMatrix.inverse();
        applyTransform(scaleMatrixInv, transSliceNum, -1);
    }
    McMat4f finalMatrix(bestScaleFactor, 0, 0, 0, 0, bestScaleFactor, 0, 0, 0,
                        0, /*zScaling=*/1.0, 0, 0, 0, 0, 1.0);
    McMat4f transMatrix = McMat4f::identity();
    transMatrix.setTranslate(McVec3f(0, 0, bestGapVariation));
    finalMatrix.multRight(transMatrix);
    finalMatrix.multRight(bestMatrix);
    transMat = finalMatrix;

    clearPointMatchingAttributeForSlices(refSliceNum, transSliceNum, editor);

    setPointMatchingAttribute(bestAlignmentResults.refVertexSelection,
                              bestAlignmentResults.transVertexSelection,
                              bestAlignmentResults.matchedRefPointIds,
                              bestAlignmentResults.matchedTransPointIds,
                              refSliceNum, transSliceNum, editor);
    setUsedForAlignmentLabel(bestAlignmentResults.refVertexSelection,
                             bestAlignmentResults.transVertexSelection);
}

void MicrotubuleSpatialGraphAligner::setUsedForAlignmentLabel(
    const SpatialGraphSelection& refSelection,
    const SpatialGraphSelection& transSelection) {
    rewriteLabel(refSelection, "UsedForMatching");
    rewriteLabel(transSelection, "UsedForMatching");
}

void MicrotubuleSpatialGraphAligner::rewriteLabel(
    const SpatialGraphSelection& selection, const McString& labelName) {
    int numPoints = selection.getNumSelectedVertices();
    McDArray<int> points(numPoints);
    for (int i = 0; i < numPoints; i++) {
        points[i] = selection.getSelectedVertex(i);
    }
    rewriteBinaryLabel(points, labelName);
}

void MicrotubuleSpatialGraphAligner::rewriteBinaryLabel(
    const McDArray<int>& newPositives, McString labelName) {

    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(mGraph->addAttribute(
            labelName, HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    mGraph->addNewLabelGroup(labelName, false, true);

    HierarchicalLabels* labelGroup = mGraph->getLabelGroup(labelName);
    if (!labelGroup) {
        return;
    }

    McDArray<int> existingPositives;

    // find existing ambiguities
    for (int i = 0; i < mGraph->getNumVertices(); ++i) {
        int ambNum = att->getIntDataAtIdx(i);
        if (ambNum > 0) {
            existingPositives.append(i);
        }
    }

    existingPositives.appendArray(newPositives);

    // rewrite labels
    labelGroup->removeChildLabels();

    SbColor color;
    color[0] = 1;
    color[1] = 0;
    color[2] = 0;
    int val = mGraph->addLabel(labelName, 0, "P", color);
    for (int i = 0; i < existingPositives.size(); ++i) {
        att->setIntDataAtIdx(existingPositives[i], val);
    }
}

void MicrotubuleSpatialGraphAligner::writeResultsToSpreadSheet(
    const AlignmentInputParams& inputParams,
    const AlignmentResults& alignmentResults, const float gap,
    const float scale, McHandle<SpreadSheetWrapper> spreadSheet) {

    McString refName = getNameOfIthSlice(inputParams.refSliceNum);
    McString transName = getNameOfIthSlice(inputParams.transSliceNum);
    McString rowN, dummy;
    char* rowName = rowN.printf("%d-%d-gap%f-scale%f", inputParams.refSliceNum,
                                inputParams.transSliceNum, gap, scale);
    spreadSheet->addEntry("AlignmentResults", "Slice1", rowName,
                          refName.dataPtr());
    spreadSheet->addEntry("AlignmentResults", "Slice2", rowName,
                          transName.dataPtr());

    char* numPointsUsed1 = dummy.printf(
        "%d", alignmentResults.refVertexSelection.getNumSelectedVertices());
    spreadSheet->addEntry("AlignmentResults", "numUsed1", rowName,
                          numPointsUsed1);
    char* numPointsUsed2 = dummy.printf(
        "%d", alignmentResults.transVertexSelection.getNumSelectedVertices());
    spreadSheet->addEntry("AlignmentResults", "numUsed2", rowName,
                          numPointsUsed2);

    char* numMatched =
        dummy.printf("%d", alignmentResults.matchedRefPointIds.size());
    spreadSheet->addEntry("AlignmentResults", "numMatched", rowName,
                          numMatched);
    spreadSheet->addEntry("AlignmentResults", "ratio", rowName,
                          (float)(alignmentResults.matchedRefPointIds.size()) /
                              (float)(alignmentResults.refVertexSelection
                                          .getNumSelectedVertices()));

    char* time = dummy.printf("%f sec", alignmentResults.duration_s);
    spreadSheet->addEntry("AlignmentResults", "time", rowName, time);

    char* scalec = dummy.printf("%f", scale);
    spreadSheet->addEntry("AlignmentResults", "scale", rowName, scalec);

    char* gapc = dummy.printf("%f", gap);
    spreadSheet->addEntry("AlignmentResults", "gap", rowName, gapc);
}

ma::EndPointParams
MicrotubuleSpatialGraphAligner::makePointRepresentationParams(
    int refSliceNum, int transSliceNum, float projectionPlane) {
    ma::EndPointParams params;
    params.refSliceNum = refSliceNum;
    params.transSliceNum = transSliceNum;
    params.projectionPlane = projectionPlane;
    params.endPointRegion = mEndPointRegion;
    params.useAbsoluteValueForEndPointRegion =
        mUseAbsoluteValueForEndPointRegion;
    params.projectionType = asEnumProjectionType(mProjectionType);
    params.numMaxPointsForInitTransform = mNumMaxPointsForInitialTransform;
    params.maxDistForAngle = mMaxDistForAngle;
    params.angleToPlaneFilter = mAngleToPlaneFilter;
    return params;
}

// Use non-PGM matching to compute an alignment.
//
// For `mAlignType==false`, an initial transform is computed: only a subset of
// points is used and the matching is computed in two stages (clique-detection
// followed by greedy or exact).
//
// For `mAlignType==true`, all points are projected and a greedy or exact
// matching is computed directly without clique-based initialization.
McMat4f MicrotubuleSpatialGraphAligner::alignNonPGM(
    const AlignmentInputParams& inputParams,
    AlignmentResults& alignmentResults) {

    const bool useProjectSubset = mAlignType ? false : true;
    const bool useBestMatching = mAlignType ? false : true;
    const bool useInitWeights = mAlignType ? false : true;

    ma::MatchingAlgorithm algo;
    if (mPMAlgorithm == PointMatchingAlgorithmTypes[1]) {
        algo = ma::MA_EXACT;
    } else if (mPMAlgorithm == PointMatchingAlgorithmTypes[0]) {
        algo = ma::MA_GREEDY;
    } else {
        // Set to fallback value to avoid uninitialized value.
        algo = ma::MA_GREEDY;
        mcerror("Invalid matching algorithm.");
    }

    const ma::MatchingParams cliqueParams = makeMatchingParams();
    const mtalign::PGMPairWeightsParams weightConfig =
        useInitWeights ? makeWeightConfigForInitMatching(cliqueParams)
                       : mWeightConfig;

    const ma::EndPointParams params = makePointRepresentationParams(
        inputParams.refSliceNum, inputParams.transSliceNum,
        inputParams.projectionPlane);

    const ma::FacingPointSets pts =
        useProjectSubset
            ? ma::projectEndPointsSubset(
                  mGraph, alignmentResults.refVertexSelection,
                  alignmentResults.transVertexSelection, params)
            : ma::projectEndPoints(mGraph, alignmentResults.refVertexSelection,
                                   alignmentResults.transVertexSelection,
                                   params);

    theMsg->printf(
        "start\t ref- points:\t %d",
        alignmentResults.refVertexSelection.getNumSelectedVertices());

    theMsg->printf(
        "start\t trans- points:\t %d",
        alignmentResults.transVertexSelection.getNumSelectedVertices());

    const ma::Matching matching =
        useBestMatching
            ? matchingTwoStepBest(pts, cliqueParams, weightConfig, algo)
            : matchingDirect(pts, McMat4f::identity(), weightConfig, algo);

    alignmentResults.matchedRefPointIds = matching.matchedRefPointIds;
    alignmentResults.matchedTransPointIds = matching.matchedTransPointIds;
    alignmentResults.score =
        (float)(alignmentResults.matchedRefPointIds.size()) /
        (float)(alignmentResults.refVertexSelection.getNumSelectedVertices());

    if (matching.matchedRefPointIds.size() == 0) {
        return McMat4f::identity();
    }
    return fitTransform(pts, matching, asTransformType(mTransformType));
}

// Apply PGM matching to all points and compute an alignment for the matching.
McMat4f MicrotubuleSpatialGraphAligner::alignPGM(
    const AlignmentInputParams& inputParams,
    AlignmentResults& alignmentResults) {

    const ma::EndPointParams params = makePointRepresentationParams(
        inputParams.refSliceNum, inputParams.transSliceNum,
        inputParams.projectionPlane);
    const ma::FacingPointSets pts =
        ma::projectEndPoints(mGraph, alignmentResults.refVertexSelection,
                             alignmentResults.transVertexSelection, params);

    theMsg->printf(
        "opt\t ref- points:\t %d",
        alignmentResults.refVertexSelection.getNumSelectedVertices());

    theMsg->printf(
        "opt\t trans- points:\t %d",
        alignmentResults.transVertexSelection.getNumSelectedVertices());

    const SpatialGraphSelection aroundMidPlaneSelection =
        mSelectionHelper.selectCloseToMidplane(inputParams.refSliceNum,
                                               inputParams.transSliceNum);
    ma::MatchingPGMParams mparams;
    mparams.weightConfig = mWeightConfig;
    mparams.pairFactorParam = mPairFactorParam;
    const ma::MatchingPGM matching = ma::matchingPGM(
        pts, mparams, mGraph, alignmentResults.refVertexSelection,
        alignmentResults.transVertexSelection, aroundMidPlaneSelection);

    alignmentResults.matchedRefPointIds = matching.matchedRefPointIds;
    alignmentResults.matchedTransPointIds = matching.matchedTransPointIds;
    alignmentResults.score =
        (float)(alignmentResults.matchedRefPointIds.size()) /
        (float)(alignmentResults.refVertexSelection.getNumSelectedVertices());

    if (matching.matchedRefPointIds.size() == 0) {
        return McMat4f::identity();
    }
    const ma::Matching m = { matching.matchedRefPointIds,
                             matching.matchedTransPointIds };
    return fitTransform(pts, m, asTransformType(mTransformType));
}

ma::MatchingParams MicrotubuleSpatialGraphAligner::makeMatchingParams() const {
    ma::MatchingParams params;
    params.maxDistanceForGraphConstruction = mMaxDistanceDifference;
    params.maxAngleDiffForInitMatching = mMaxAngleDiffForInitMatching;
    params.minCliqueSizeFraction = mMinCliqueSizeFraction;
    params.maxNumCliques = mMaxNumCliques;
    if (mTransformType == TransformTypes[1]) {
        params.transformType = ma::TF_RIGID_ISO_SCALE;
    } else {
        params.transformType = ma::TF_RIGID;
    }
    return params;
}

void MicrotubuleSpatialGraphAligner::createPointMatchingAttribute(
    HxNeuronEditorSubApp* editor) {
    McString curPMAttName;
    getPointMatchingString(curPMAttName);
    const char* attName = curPMAttName.dataPtr();
    HierarchicalLabels* labelGroup = mGraph->getLabelGroup(attName);
    if (labelGroup) {
        mNotUsedLabelValue =
            labelGroup->getLabelIdFromName(McString("NotUsed"));
        mNotMatchedLabelValue =
            labelGroup->getLabelIdFromName(McString("NotMatched"));
        if (mNotUsedLabelValue < 1 || mNotMatchedLabelValue < 1)
            mNotMatchedLabelValue = 2;
        mNotUsedLabelValue = 1;
    } else {
        if (mCreateMatchingLabels) {
            mGraph->addNewLabelGroup(attName, false, true);
            mNotUsedLabelValue = mGraph->addLabel(attName, 0, "NotUsed",
                                                  SbColor(1.0f, 1.0f, 1.0f));
            mNotMatchedLabelValue = mGraph->addLabel(attName, 0, "NotMatched",
                                                     SbColor(0.0f, 0.0f, 0.0f));
            EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
                mGraph->findAttribute(HxSpatialGraph::VERTEX, attName));
            // Set all vertices to unused.
            for (int i = 0; i < att->size(); ++i) {
                att->setIntDataAtIdx(i, mNotUsedLabelValue);
            }
            mGraph->touch(HxData::NEW_PARAMETERS);
            HxNeuronEditorSubApp::SpatialGraphChanges changes =
                HxNeuronEditorSubApp::SpatialGraphAttributeChange |
                HxNeuronEditorSubApp::SpatialGraphLabelTreeChange;
            editor->updateAfterSpatialGraphChange(changes);
        } else {
            EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
                mGraph->findAttribute(HxSpatialGraph::VERTEX, attName));
            if (!att) {
                att = dynamic_cast<EdgeVertexAttribute*>(mGraph->addAttribute(
                    attName, HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
                mcassert(att);
                // Set all vertices to unused.
                for (int i = 0; i < att->size(); ++i) {
                    att->setIntDataAtIdx(i, mNotUsedLabelValue);
                }
                mGraph->touch(HxData::NEW_PARAMETERS);
                HxNeuronEditorSubApp::SpatialGraphChanges changes =
                    HxNeuronEditorSubApp::SpatialGraphAttributeChange;
                editor->updateAfterSpatialGraphChange(changes);
            }
        }
    }
}

void MicrotubuleSpatialGraphAligner::setPointMatchingAttribute(
    const SpatialGraphSelection& refVertexSel,
    const SpatialGraphSelection& transVertexSel,
    const McDArray<int>& matchedRefPoints,
    const McDArray<int>& matchedTransPoints, const int refSliceNum,
    const int transSliceNum, HxNeuronEditorSubApp* editor) {

    McDArray<int> selectedRefIds;
    McDArray<int> selectedTransIds;

    // Find vertex ids of matched points.
    int vertexNum;
    SpatialGraphSelection::Iterator refIter(refVertexSel);
    refIter.vertices.reset();
    vertexNum = refIter.vertices.nextSelected();
    while (vertexNum != -1) {
        selectedRefIds.append(vertexNum);
        vertexNum = refIter.vertices.nextSelected();
    }

    SpatialGraphSelection::Iterator transIter(transVertexSel);
    transIter.vertices.reset();
    vertexNum = transIter.vertices.nextSelected();
    while (vertexNum != -1) {
        selectedTransIds.append(vertexNum);
        vertexNum = transIter.vertices.nextSelected();
    }

    McString pmAttName;
    getPointMatchingString(pmAttName);
    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        mGraph->findAttribute(HxSpatialGraph::VERTEX, pmAttName.dataPtr()));

    // Reset half reference and half transformed slice.
    SpatialGraphSelection refVertices = refVertexSel;
    SpatialGraphSelection transVertices = transVertexSel;

    SpatialGraphSelection::Iterator it(refVertices);
    it.vertices.reset();
    vertexNum = it.vertices.nextSelected();
    while (vertexNum != -1) {
        att->setIntDataAtIdx(vertexNum, mNotUsedLabelValue);
        vertexNum = it.vertices.nextSelected();
    }

    SpatialGraphSelection::Iterator it2(transVertices);
    it2.vertices.reset();
    vertexNum = it2.vertices.nextSelected();
    while (vertexNum != -1) {
        att->setIntDataAtIdx(vertexNum, mNotUsedLabelValue);
        vertexNum = it2.vertices.nextSelected();
    }

    // Set all vertices used for matching to not matched.
    for (int i = 0; i < selectedRefIds.size(); ++i) {
        att->setIntDataAtIdx(selectedRefIds[i], mNotMatchedLabelValue);
    }
    for (int i = 0; i < selectedTransIds.size(); ++i) {
        att->setIntDataAtIdx(selectedTransIds[i], mNotMatchedLabelValue);
    }

    // Get next label value, find the actual highest labels val.
    int maxAttVal = 0;
    for (int i = 0; i < mGraph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (pairNum > maxAttVal)
            maxAttVal = pairNum;
    }
    int nextLabel = maxAttVal + 1;

    for (int i = 0; i < matchedRefPoints.size(); ++i) {
        att->setIntDataAtIdx(selectedRefIds[matchedRefPoints[i]], nextLabel);
        att->setIntDataAtIdx(selectedTransIds[matchedTransPoints[i]],
                             nextLabel);
        nextLabel++;
    }

    rewriteLabels(editor);
}

void
MicrotubuleSpatialGraphAligner::rewriteLabels(HxNeuronEditorSubApp* editor) {
    McString curPMAttName;
    getPointMatchingString(curPMAttName);
    const char* attName = curPMAttName.dataPtr();
    HierarchicalLabels* labelGroup = mGraph->getLabelGroup(attName);
    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        mGraph->findAttribute(HxSpatialGraph::VERTEX, attName));

    // Get maximum att val.
    int maxAttVal = 0;
    for (int i = 0; i < mGraph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (pairNum > maxAttVal)
            maxAttVal = pairNum;
    }

    McDArray<int> attLabelMap(maxAttVal + 2);
    attLabelMap.fill(-1);
    attLabelMap[mNotUsedLabelValue] = 1;
    attLabelMap[mNotMatchedLabelValue] = 2;
    int countNewLabels = 3;
    // Find mapping from unsorted labels to new labels.
    for (int i = 0; i < mGraph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (attLabelMap[pairNum] == -1) {
            attLabelMap[pairNum] = countNewLabels;
            countNewLabels++;
        }
    }

    // Set the new attribute values.
    for (int i = 0; i < mGraph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        att->setIntDataAtIdx(i, attLabelMap[pairNum]);
    }

    // Rewrite labels.
    if (labelGroup) {
        labelGroup->removeChildLabels();
        mGraph->addLabel(attName, 0, "NotUsed", SbColor(1.0f, 1.0f, 1.0f));
        mGraph->addLabel(attName, 0, "NotMatched", SbColor(0.0f, 0.0f, 0.0f));
        for (int i = 0; i < mGraph->getNumVertices(); ++i) {
            int pairNum = att->getIntDataAtIdx(i);
            if ((pairNum == mNotUsedLabelValue) ||
                (pairNum == mNotMatchedLabelValue))
                continue;
            McString s;
            s.printf("Pair%d", pairNum);
            SbColor color;
            color[0] = float(rand()) / float(RAND_MAX);
            color[1] = float(rand()) / float(RAND_MAX);
            color[2] = float(rand()) / float(RAND_MAX);
            mGraph->addLabel(attName, 0, s.getString(), color);
        }
    }
    mGraph->touch(HxData::NEW_PARAMETERS);
}

void
MicrotubuleSpatialGraphAligner::applyTransform(const McMat4f& transform,
                                               const int start, const int num,
                                               HxNeuronEditorSubApp* editor) {
    int last;
    if (num < 0) {
        last = mSelectionHelper.getNumSlices() - 1;
    } else {
        last = start + num - 1;
        if (last >= mSelectionHelper.getNumSlices()) {
            last = mSelectionHelper.getNumSlices() - 1;
        }
    }
    SpatialGraphSelection sel(mGraph->getNumVertices(), mGraph->getNumEdges());
    sel.clear();
    McDArray<HxParameter*> transParams(0);
    for (int i = start; i <= last; ++i) {
        SpatialGraphSelection tmpsel;
        mSelectionHelper.getSlice(
            mSelectionHelper.getSliceAttributeValueFromIndex(i), tmpsel);
        sel.addSelection(tmpsel);
        transParams.append(getTransformParameter(i));
    }

    McMat4f dummyTransformForCast = transform;
    SbMatrix transformMat =
        SbMatrix(*(reinterpret_cast<SbMatrix*>(&dummyTransformForCast)));
    if (editor) {

        MicrotubuleTransformOperation* op = new MicrotubuleTransformOperation(
            mGraph, sel, editor->getVisibleElements(), transformMat);
        op->setTransformParameterList(transParams);
        if (!editor->execNewOp(op)) {
            theMsg->printf("Error performing transformation");
            return;
        }
    } else {
        SpatialGraphSelection visSel(mGraph->getNumVertices(),
                                     mGraph->getNumEdges());
        MicrotubuleTransformOperation* op = new MicrotubuleTransformOperation(
            mGraph, sel, visSel, transformMat);
        op->setTransformParameterList(transParams);
        op->exec();
    }
}

HxParameter*
MicrotubuleSpatialGraphAligner::getTransformParameter(const int slice) const {

    HxParamBundle* transformInfoPB =
        mGraph->parameters.bundle("TransformInfo", 0);
    if (!transformInfoPB) {
        theMsg->printf(
            "Warning: could not find TransformInfo parameter bundle");
        return 0;
    }

    const char* attStr;
    if (!transformInfoPB->findString("AlignAttribute", attStr)) {
        theMsg->printf("Warning: could not find AlignAttribute");
        return 0;
    }

    HxParamBundle* slicePB = transformInfoPB->bundle(slice);
    if (slicePB) {
        int v;
        if (slicePB->findNum("Value", v)) {
            if (v != mSelectionHelper.getSliceAttributeValueFromIndex(slice)) {
                theMsg->printf("Warning: param value for slice %d differs from "
                               "graph value",
                               slice);
            }
        }
        return slicePB->find("Transform");
    } else {
        theMsg->printf("Warning: No transform parameter for slice %d", slice);
    }
    return 0;
}

McMat4f MicrotubuleSpatialGraphAligner::getTransform(const int sliceNum) const {
    HxParameter* param = getTransformParameter(sliceNum);
    McMat4f mat;
    if (!param) {
        mat.makeIdentity();
        theMsg->printf("Error: no transformation found for slice %d.",
                       sliceNum);
    } else if (param->dim() == 16) {
        double res[16];
        param->getReal(res);
        float* ptr = &mat[0][0];
        for (int i = 0; i < 16; ++i) {
            ptr[i] = float(res[i]);
        }
    } else {
        mat.makeIdentity();
        theMsg->printf("Error: no transformation found for slice %d.",
                       sliceNum);
    }
    return mat;
}

float
MicrotubuleSpatialGraphAligner::getSliceZPosition(const int sliceNum) const {
    McMat4f mat = getTransform(sliceNum);
    return mat[3][2];
}

void MicrotubuleSpatialGraphAligner::setSliceZPosition(
    const int sliceNum, const float newZ, HxNeuronEditorSubApp* editor) {
    float oldZ = getSliceZPosition(sliceNum);
    McMat4f transform = McMat4f::identity();
    McVec3f translation(0.0f, 0.0f, newZ - oldZ);
    transform.setTranslate(translation);
    applyTransform(transform, sliceNum, -1, editor);
}

bool MicrotubuleSpatialGraphAligner::hasTransformInfo() const {
    return mGraph->parameters.bundle("TransformInfo", 0);
}

McString MicrotubuleSpatialGraphAligner::getTransformInfoAttribute() const {
    HxParamBundle* transformInfoPB =
        mGraph->parameters.bundle("TransformInfo", 0);
    if (!transformInfoPB) {
        theMsg->printf(
            "Warning: could not find TransformInfo parameter bundle");
        return 0;
    }

    const char* attStr;
    if (!transformInfoPB->findString("AlignAttribute", attStr)) {
        theMsg->printf("Warning: could not find AlignAttribute");
        return 0;
    }
    return McString(attStr);
}

McString
MicrotubuleSpatialGraphAligner::getNameOfIthSlice(const int sliceNum) const {
    HxParamBundle* transformInfoPB =
        mGraph->parameters.bundle("TransformInfo", 0);
    if (!transformInfoPB) {
        theMsg->printf(
            "Warning: could not find TransformInfo parameter bundle");
        return 0;
    }
    if (transformInfoPB->nBundles() <= sliceNum) {
        theMsg->printf("Warning: Slice %d does not exist", sliceNum);
        return 0;
    }
    HxParamBundle* sliceParameters = transformInfoPB->bundle(sliceNum);

    const char* nameString = NULL;
    if (!sliceParameters->findString("FileName", nameString)) {
        theMsg->printf("Warning: Cannot find slice name");
        return 0;
    }

    McString filename(nameString);
    McDArray<McString> substrings;
    char slash = '/';
    filename.explode(slash, substrings);
    McString finalName = substrings[substrings.size() - 1];

    return finalName;
}

void MicrotubuleSpatialGraphAligner::getSliceSelection(
    const int sliceNum, SpatialGraphSelection& slice) const {
    mSelectionHelper.getSlice(
        mSelectionHelper.getSliceAttributeValueFromIndex(sliceNum), slice);
}

McHandle<SpreadSheetWrapper>
MicrotubuleSpatialGraphAligner::getEvaluationSpreadsheet() {

    McHandle<SpreadSheetWrapper> resultSpreadSheet = new SpreadSheetWrapper();
    theObjectPool->addObject(resultSpreadSheet);
    McString label(mGraph->getLabel());
    label += "-alignResults";
    resultSpreadSheet->setLabel(label);
    return resultSpreadSheet;
}

void MicrotubuleSpatialGraphAligner::clearPointMatchingAttributeForSlices(
    const int refSliceNum, const int transSliceNum,
    HxNeuronEditorSubApp* editor) {

    SpatialGraphSelection refSliceSelection, transSliceSelection;
    mSelectionHelper.selectAdjacentHalfSlices(
        mSelectionHelper.getSliceAttributeValueFromIndex(refSliceNum),
        mSelectionHelper.getSliceAttributeValueFromIndex(transSliceNum),
        refSliceSelection, transSliceSelection);
    createPointMatchingAttribute(editor);
    McString curPMAttName;
    getPointMatchingString(curPMAttName);
    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        mGraph->findAttribute(HxSpatialGraph::VERTEX, curPMAttName.dataPtr()));

    SpatialGraphSelection::Iterator it(refSliceSelection);
    it.vertices.reset();
    int vertexNum = it.vertices.nextSelected();
    while (vertexNum != -1) {
        att->setIntDataAtIdx(vertexNum, mNotUsedLabelValue);
        vertexNum = it.vertices.nextSelected();
    }

    SpatialGraphSelection::Iterator it2(transSliceSelection);
    it2.vertices.reset();
    vertexNum = it2.vertices.nextSelected();
    while (vertexNum != -1) {
        att->setIntDataAtIdx(vertexNum, mNotUsedLabelValue);
        vertexNum = it2.vertices.nextSelected();
    }
}

McMat4f
MicrotubuleSpatialGraphAligner::align(const AlignmentInputParams& inputParams,
                                      AlignmentResults& alignmentResults,
                                      HxNeuronEditorSubApp* editor) {
    if (mPMAlgorithm == PointMatchingAlgorithmTypes[2]) {
        return alignPGM(inputParams, alignmentResults);
    } else {
        return alignNonPGM(inputParams, alignmentResults);
    }
}

void MicrotubuleSpatialGraphAligner::getPointMatchingString(McString& curName) {
    if (mPMAlgorithm == PointMatchingAlgorithmTypes[2]) {
        curName = McString("PGMAssignedPairs");
    } else if (mPMAlgorithm == PointMatchingAlgorithmTypes[1]) {
        curName = McString("ExactAssignedPairs");
    } else if (mPMAlgorithm == PointMatchingAlgorithmTypes[0]) {
        curName = McString("GreedyAssignedPairs");
    } else {
        theMsg->printf("No valid algorithm selected");
    }
}

void
MicrotubuleSpatialGraphAligner::joinMatchedEnds(const char* matchingLabel) {
    EdgeVertexAttribute* pointAtt = mGraph->findVertexAttribute(matchingLabel);
    if (pointAtt == NULL) {
        theMsg->printf("Attribute %s not found", matchingLabel);
        return;
    }

    // Join the matched ends.
    for (int i = 0; i < mGraph->getNumVertices(); i++) {
        int* curlabel = pointAtt->intDataAtIdx(i);
        if (curlabel[0] > 2) {
            int otherIdx = -1;
            for (int j = i + 1; j < mGraph->getNumVertices(); j++) {
                int* tmplabel = pointAtt->intDataAtIdx(j);
                if (tmplabel[0] == curlabel[0]) {
                    otherIdx = j;
                    break;
                }
            }
            // Join and reset label.
            if (otherIdx > -1) {
                mGraph->addEdge(i, otherIdx);
            }
        }
    }
    addLabelForUnmatched(matchingLabel);
}

void MicrotubuleSpatialGraphAligner::addLabelForUnmatched(
    const char* matchingLabel) {
    // Give all ends and connected lines that were not matched a different
    // label.
    EdgeVertexAttribute* pointAtt = mGraph->findVertexAttribute(matchingLabel);
    if (pointAtt == NULL) {
        theMsg->printf("Attribute %s not found", matchingLabel);
        return;
    }

    EdgeVertexAttribute* unmatchedVertex = dynamic_cast<EdgeVertexAttribute*>(
        mGraph->findAttribute(HxSpatialGraph::VERTEX, "Unmatched"));
    if (!unmatchedVertex)
        unmatchedVertex =
            dynamic_cast<EdgeVertexAttribute*>(mGraph->addAttribute(
                "Unmatched", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    unmatchedVertex->clearMemory(0);

    EdgeVertexAttribute* unmatchedEdge = dynamic_cast<EdgeVertexAttribute*>(
        mGraph->findAttribute(HxSpatialGraph::EDGE, "Unmatched"));
    if (!unmatchedEdge)
        unmatchedEdge = dynamic_cast<EdgeVertexAttribute*>(mGraph->addAttribute(
            "Unmatched", HxSpatialGraph::EDGE, McPrimType::mc_int32, 1));
    unmatchedVertex->clearMemory(0);

    mGraph->addNewLabelGroup("Unmatched", false, true);
    SbColor color;
    color[0] = 1;
    color[1] = 0;
    color[2] = 0;

    int valUnmatched =
        mGraph->addLabel("Unmatched", 0, "unmatchedElements", color);
    color[1] = 1;
    color[0] = 0;
    int valMatched = mGraph->addLabel("Unmatched", 0, "matchedElements", color);

    for (int e = 0; e < mGraph->getNumEdges(); ++e) {

        int source = mGraph->getEdgeSource(e);
        unmatchedVertex->setIntDataAtIdx(source, valMatched);

        int target = mGraph->getEdgeTarget(e);
        unmatchedVertex->setIntDataAtIdx(target, valMatched);
        unmatchedEdge->setIntDataAtIdx(e, valMatched);

        if (pointAtt->getIntDataAtIdx(source) == 2) {
            unmatchedVertex->setIntDataAtIdx(source, valUnmatched);
            unmatchedEdge->setIntDataAtIdx(e, valUnmatched);
        }
        if (pointAtt->getIntDataAtIdx(target) == 2) {
            unmatchedVertex->setIntDataAtIdx(target, valUnmatched);
            unmatchedEdge->setIntDataAtIdx(e, valUnmatched);
        }
    }
}

void MicrotubuleSpatialGraphAligner::removeIntermediateAndAddLabelForUnmatched(
    const char* matchingLabel) {
    SpatialGraphSelection selectAll(mGraph);
    selectAll.selectAllVerticesAndEdges();
    FindAndConvertVerticesToPointsOperationSet op(mGraph, selectAll, selectAll);
    op.exec();
}
