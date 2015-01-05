#include <hxalignmicrotubules/HxCPDSpatialGraphWarp.h>

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <Inventor/SbLinear.h>
#include <Inventor/nodes/SoTransform.h>

#include <hxcore/HxInterpreter.h>
#include <hxcore/HxMessage.h>
#include <hxcore/HxObjectPool.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxLoc3Uniform.h>
#include <hxfield/HxUniformCoord3.h>
#include <hxfield/HxUniformLabelField3.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformVectorField3.h>
#include <hxlandmark/HxLandmarkSet.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxwarp/HxNearestNeighbourWarp.h>
#include <hxwarp/HxPolygonWarp.h>
#include <mclib/McWatch.h>

#include <hxalignmicrotubules/CoherentPointDriftNLFisherMises.h>
#include <hxalignmicrotubules/CoherentPointDriftRigidFisherMises.h>
#include <hxalignmicrotubules/MovingLeastSquares.h>
#include <hxalignmicrotubules/mtalign.h>

namespace ma = mtalign;

HX_INIT_CLASS(HxCPDSpatialGraphWarp, HxCompModule);

HxCPDSpatialGraphWarp::HxCPDSpatialGraphWarp()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portMethod(this, "method", 4),
      portBeta(this, "betaParam"),
      portLambda(this, "lambdaParam"),
      portW(this, "wParam"),
      portUseDirection(this, "useDirections", 1),
      portUseCoords(this, "useCoords", 1),
      portWithScale(this, "with scaling", 1),
      portCoupleRcRd(this, "coupleRcAndRd", 1),
      portEMParams(this, "emParams", 3),
      portCPDResults(this, "cpdResults", 7),
      portSampleDist(this, "sampleMinDistForRemove"),
      portAlphaForMLS(this, "alphaForMLS"),
      portAction(this, "action") {
    portAction.setAliasName("doIt");
    portMethod.setLabel(0, "NA");
    portMethod.setLabel(1, "NA");
    portMethod.setLabel(2, "elastic");
    portMethod.setLabel(3, "rigid");
    portBeta.setLabel("Beta:");
    portBeta.setMinMax(0, 10);
    portBeta.setValue(10);
    portLambda.setLabel("Lambda:");
    portLambda.setMinMax(0, 10);
    portLambda.setValue(1);
    portW.setLabel("w:");
    portW.setMinMax(0, 10);
    portW.setValue(0.1);
    portUseDirection.setValue(0);
    portUseCoords.setValue(1);
    portUseCoords.setLabel("Use coords (only rigid)");
    portWithScale.setValue(1);
    portWithScale.hide();
    portEMParams.setLabel(0, "maxIter");
    portEMParams.setMinMax(0.0, 1000000.0);
    portEMParams.setValue(0, 100);

    portEMParams.setLabel(1, "sigmaTol");
    portEMParams.setValue(1, 1.e-5);

    portEMParams.setLabel(2, "LLDiffTol");
    portEMParams.setValue(2, 1.e-5);
    portCPDResults.setLabel(0, "sigma^2");
    portCPDResults.setLabel(1, "kappa");
    portCPDResults.setLabel(2, "s");
    portCPDResults.setLabel(3, "rotation");
    portCPDResults.setLabel(4, "tx");
    portCPDResults.setLabel(5, "ty");
    portCPDResults.setLabel(6, "numIter");
    portAlphaForMLS.setMinMax(0, 10);
    portAlphaForMLS.setValue(2);
    portSampleDist.setMinMax(0, 50000000);

    portSampleDist.setValue(0);
    portCoupleRcRd.setValue(1);
    portCoupleRcRd.hide();
    portUseDirection.setLabel("Use orientation");
}

HxCPDSpatialGraphWarp::~HxCPDSpatialGraphWarp() {}

void HxCPDSpatialGraphWarp::update() {}

HxUniformVectorField3* HxCPDSpatialGraphWarp::createOutputVectorDataSet() {
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

HxSpatialGraph* HxCPDSpatialGraphWarp::createOutputDataSet() {
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

    return warpedSG;
}

void HxCPDSpatialGraphWarp::preparePoints(McDArray<McVec3f>& p1,
                                          McDArray<McVec3f>& p2,
                                          SpatialGraphSelection& slice2,
                                          const HxSpatialGraph* spatialGraph) {
    ma::SliceSelector selectionHelper(spatialGraph, "TransformInfo");

    ma::EndPointParams params;
    params.endPointRegion = 30;
    params.projectionPlane = selectionHelper.computeMidPlane(0, 1);
    params.projectionType = McString("Orthogonal");
    params.refSliceNum = 0;
    params.transSliceNum = 1;
    params.useAbsouteValueForEndPointRegion = false;
    params.transformType = 2;
    params.maxDistForAngle = 2000;
    params.angleToPlaneFilter = 0.01;
    SpatialGraphSelection slice1;
    ma::FacingPointSets pr =
        ma::projectEndPoints(spatialGraph, slice1, slice2, params);
    McDArray<McVec3f> refCoords = pr.ref.positions;
    McDArray<McVec3f> transCoords = pr.trans.positions;
    mcassert(refCoords.size() == slice1.getNumSelectedVertices());
    mcassert(transCoords.size() == slice2.getNumSelectedVertices());

    p1.resize(refCoords.size());
    for (int i = 0; i < refCoords.size(); i++) {
        McVec3f coord = refCoords[i];
        p1[i] = McVec3f(coord.x, coord.y, 0);
    }
    p2.resize(transCoords.size());
    for (int i = 0; i < transCoords.size(); i++) {
        McVec3f coord = transCoords[i];
        p2[i] = McVec3f(coord.x, coord.y, 0);
    }
    mcassert(p2.size() == slice2.getNumSelectedVertices());
}

void HxCPDSpatialGraphWarp::preparePointsAndDirections(
    McDArray<McVec3f>& p1, McDArray<McVec3f>& p2, McDArray<McVec3f>& d1,
    McDArray<McVec3f>& d2, SpatialGraphSelection& slice2,
    const HxSpatialGraph* spatialGraph) {

    ma::SliceSelector selectionHelper(spatialGraph, "TransformInfo");

    ma::EndPointParams params;
    params.endPointRegion = 30;
    params.projectionPlane = selectionHelper.computeMidPlane(0, 1);
    params.projectionType = McString("Orthogonal");
    params.refSliceNum = 0;
    params.transSliceNum = 1;
    params.useAbsouteValueForEndPointRegion = false;
    params.transformType = 2;
    params.maxDistForAngle = 2000;
    params.angleToPlaneFilter = 0.01;
    SpatialGraphSelection slice1;
    ma::FacingPointSets pr =
        ma::projectEndPoints(spatialGraph, slice1, slice2, params);
    McDArray<McVec3f> refCoords = pr.ref.positions;
    McDArray<McVec3f> transCoords = pr.trans.positions;
    McDArray<McVec3f> refDirs = pr.ref.directions;
    McDArray<McVec3f> transDirs = pr.trans.directions;
    mcassert(refCoords.size() == slice1.getNumSelectedVertices());
    mcassert(transCoords.size() == slice2.getNumSelectedVertices());

    p1.resize(refCoords.size());
    d1.resize(refCoords.size());
    for (int i = 0; i < refCoords.size(); i++) {
        McVec3f coord = refCoords[i];
        p1[i] = McVec3f(coord.x, coord.y, 0);
        d1[i] = p1[i] + refDirs[i] * 2000;
    }
    p2.resize(transCoords.size());
    d2.resize(transCoords.size());
    for (int i = 0; i < transCoords.size(); i++) {
        McVec3f coord = transCoords[i];
        p2[i] = McVec3f(coord.x, coord.y, 0);
        d2[i] = p2[i] + transDirs[i] * -2000;
    }
    mcassert(p2.size() == slice2.getNumSelectedVertices());
}

void HxCPDSpatialGraphWarp::preparePointsAndDirectionsRigid(
    McDArray<McVec3f>& p1, McDArray<McVec3f>& p2, McDArray<McVec3f>& d1,
    McDArray<McVec3f>& d2, SpatialGraphSelection& slice2,
    const HxSpatialGraph* spatialGraph) {

    ma::SliceSelector selectionHelper(spatialGraph, "TransformInfo");

    ma::EndPointParams params;
    params.endPointRegion = 30;
    params.projectionPlane = selectionHelper.computeMidPlane(0, 1);
    params.projectionType = McString("Orthogonal");
    params.refSliceNum = 0;
    params.transSliceNum = 1;
    params.useAbsouteValueForEndPointRegion = false;
    params.transformType = 2;
    params.maxDistForAngle = 2000;
    params.angleToPlaneFilter = 0.01;
    SpatialGraphSelection slice1;
    ma::FacingPointSets pr =
        ma::projectEndPoints(spatialGraph, slice1, slice2, params);
    McDArray<McVec3f> refCoords = pr.ref.positions;
    McDArray<McVec3f> transCoords = pr.trans.positions;
    McDArray<McVec3f> refDirs = pr.ref.directions;
    McDArray<McVec3f> transDirs = pr.trans.directions;
    mcassert(refCoords.size() == slice1.getNumSelectedVertices());
    mcassert(transCoords.size() == slice2.getNumSelectedVertices());

    p1.resize(refCoords.size());
    d1.resize(refCoords.size());
    for (int i = 0; i < refCoords.size(); i++) {
        McVec3f coord = refCoords[i];
        p1[i] = McVec3f(coord.x, coord.y, 0);
        d1[i] = refDirs[i] * 2000;
    }
    p2.resize(transCoords.size());
    d2.resize(transCoords.size());
    for (int i = 0; i < transCoords.size(); i++) {
        McVec3f coord = transCoords[i];
        p2[i] = McVec3f(coord.x, coord.y, 0);
        d2[i] = transDirs[i] * -2000;
    }
    mcassert(p2.size() == slice2.getNumSelectedVertices());
}

void HxCPDSpatialGraphWarp::compute() {
    if (!portAction.wasHit())
        return;

    const HxSpatialGraph* inputSG = (HxSpatialGraph*)portData.source();

    if (!inputSG)
        return;

    if (portMethod.getValue() == 2)
        computeNL();
    else if (portMethod.getValue() == 3)
        computeRigidVanMises();
    else
        theMsg->printf("Not implemented yet!");
}

void HxCPDSpatialGraphWarp::computeRigidVanMises() {
    const HxSpatialGraph* inputSG = (HxSpatialGraph*)portData.source();

    ma::SliceSelector selectionHelper(inputSG, "TransformInfo");

    HxSpatialGraph* resultSG = createOutputDataSet();
    resultSG->copyFrom(inputSG);
    McDArray<McVec3f> points1, points2;
    McDArray<McVec3f> dirs1, dirs2;
    SpatialGraphSelection slice2;
    preparePointsAndDirectionsRigid(points1, points2, dirs1, dirs2, slice2,
                                    resultSG);
    mcassert(slice2.getNumSelectedVertices() == points2.size());
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.w = portW.getValue();
    cpd.params.withScaling = portWithScale.getValue();
    cpd.params.maxIterations = portEMParams.getValue(0);
    cpd.params.sigmaSquareStop = portEMParams.getValue(1);
    cpd.params.eDiffRelStop = portEMParams.getValue(2);
    cpd.params.useDirections = portUseDirection.getValue();
    cpd.params.usePositions = portUseCoords.getValue();

    if (!portUseDirection.getValue()) {
        dirs1.fill(McVec3f(0, 0, 1));
        dirs2.fill(McVec3f(0, 0, 1));
    }

    cpd.convertMcDArraysToMatrix(points1, dirs1, cpd.Xc, cpd.Xd);
    std::cout << "num X " << points1.size() << ", num Y " << points2.size();
    cpd.convertMcDArraysToMatrix(points2, dirs2, cpd.Yc, cpd.Yd);

    McDMatrix<double> Rc, Rd;
    McDVector<double> t;
    double s;
    McDArray<McVec2i> correspondences;
    double kappa, sigmaSquare;
    const mtalign::AlignInfo info =
        cpd.align(Rc, s, t, Rd, correspondences, sigmaSquare, kappa);

    applyRigidDeformationToSliceVanMises(slice2, resultSG, cpd, Rc, s, t);
    portCPDResults.setValue(0, sigmaSquare);
    portCPDResults.setValue(1, kappa);
    portCPDResults.setValue(2, s);
    double rho = mtalign::rotationAngle2d(Rc);

    portCPDResults.setValue(3, rho);
    portCPDResults.setValue(4, t[0]);
    portCPDResults.setValue(5, t[1]);
    portCPDResults.setValue(6, info.numIterations);

    resultSG->touch();
    resultSG->fire();
}

void HxCPDSpatialGraphWarp::applyRigidDeformationToSliceVanMises(
    SpatialGraphSelection& slice, HxSpatialGraph* spatialGraph,
    const CoherentPointDriftRigidFisherMises& deformation,
    const McDMatrix<double>& R, const double s, const McDVector<double>& t) {

    McWatch watch;
    watch.start();

    ma::SliceSelector ssh(spatialGraph, "TransformInfo");
    SpatialGraphSelection fullSliceSelection;

    ssh.getSlice(ssh.getSliceAttributeValueFromIndex(1), fullSliceSelection);

    for (int i = 0; i < fullSliceSelection.getNumSelectedEdges(); i++) {
        int edge = fullSliceSelection.getSelectedEdge(i);
        McDArray<McVec3f> newEdgePoints;
        for (int j = 0; j < spatialGraph->getNumEdgePoints(edge); j++) {

            McVec3f edgePoint = spatialGraph->getEdgePoint(edge, j);
            McVec3f edgePointWarped;
            deformation.warpPoint(edgePoint, R, s, t, edgePointWarped);
            newEdgePoints.append(edgePointWarped);
        }
        spatialGraph->setEdgePoints(edge, newEdgePoints);
    }

    for (int i = 0; i < fullSliceSelection.getNumSelectedVertices(); i++) {
        int curVertex = fullSliceSelection.getSelectedVertex(i);
        McVec3f curCoord = spatialGraph->getVertexCoords(curVertex);
        McVec3f curCoordWarped;
        deformation.warpPoint(curCoord, R, s, t, curCoordWarped);
        spatialGraph->setVertexCoords(curVertex, curCoordWarped);
    }

    std::cout << "\n Apply deformation took " << watch.stop() << " seconds.";
}

void HxCPDSpatialGraphWarp::computeNL() {
    const HxSpatialGraph* inputSG = (HxSpatialGraph*)portData.source();

    ma::SliceSelector selectionHelper(inputSG, "TransformInfo");

    HxSpatialGraph* resultSG = createOutputDataSet();
    resultSG->copyFrom(inputSG);
    McDArray<McVec3f> points1, points2;
    McDArray<McVec3f> dirs1, dirs2;
    SpatialGraphSelection slice2;
    preparePointsAndDirectionsRigid(points1, points2, dirs1, dirs2, slice2,
                                    resultSG);
    mcassert(slice2.getNumSelectedVertices() == points2.size());
    CoherentPointDriftNLFisherMises cpd;

    cpd.params.beta = portBeta.getValue();
    cpd.params.lambda = portLambda.getValue();
    cpd.params.w = portW.getValue();
    cpd.params.eDiffRelStop = portEMParams.getValue(2);
    cpd.params.sigmaSquareStop = portEMParams.getValue(1);
    cpd.params.maxIterations = portEMParams.getValue(0);
    cpd.params.useDirections = portUseDirection.getValue();
    cpd.convertCoordsToMatrix(points1, cpd.xs);
    cpd.convertCoordsToMatrix(points2, cpd.ys);
    cpd.convertDirectionsToMatrix(dirs1, cpd.xDirs);
    cpd.convertDirectionsToMatrix(dirs2, cpd.yDirs);
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> correspondences;
    const mtalign::AlignInfo info = cpd.align(G, W, correspondences);
    mcassert(slice2.getNumSelectedVertices() == cpd.ys.nRows());
    applyNLDeformationToSlice(slice2, resultSG, cpd, W, G);
    portCPDResults.setValue(0, info.sigmaSquare);
    portCPDResults.setValue(1, info.kappa);
    portCPDResults.setValue(6, info.numIterations);

    resultSG->touch();
    resultSG->fire();
}

void HxCPDSpatialGraphWarp::resamplePairs(McDArray<McVec3f>& p1,
                                          McDArray<McVec3f>& p2) {
    int numPairs = p1.size();
    for (int i = p1.size() - 1; i > -1; i--) {
        bool resetI = false;
        for (int j = i - 1; j > -1; j--) {
            // std::cout<<"\nCompare "<<i<<" and "<<j;
            McVec3f set1Coord = p1[i];
            McVec3f set2Coord = p1[j];
            float dist = (set1Coord - set2Coord).length();
            if (dist < portSampleDist.getValue()) {
                std::cout << "\ndist between " << i << " and " << j << " is "
                          << dist;
            }
            if (dist < portSampleDist.getValue()) {
                p1.remove(j, 1);
                p2.remove(j, 1);
                resetI = true;
            }
        }
        if (resetI) {
            i = p1.size() - 1;
        }
    }
    std::cout << "\n" << p1.size() << " point from " << numPairs
              << " left after resampling.";
}

void HxCPDSpatialGraphWarp::applyNLDeformationToSlice(
    SpatialGraphSelection& slice, HxSpatialGraph* spatialGraph,
    const CoherentPointDriftNLFisherMises& deformation,
    const McDMatrix<double>& W, const McDMatrix<double>& G) {

    McWatch watch;
    watch.start();
    McDMatrix<double> shiftedVertices;
    deformation.shiftYs(deformation.ys, G, W, shiftedVertices);
    deformation.rescaleYs(shiftedVertices);
    McDMatrix<double> origVertices = deformation.ys;
    deformation.rescaleYs(origVertices);
    McDArray<McVec3f> shiftedCoords, origCoords;

    McDArray<McVec3f> origDirections, shiftedDirections;
    deformation.convertMatrixToCoords(shiftedCoords, shiftedVertices);
    deformation.convertMatrixToCoords(origCoords, origVertices);
    resamplePairs(origCoords, shiftedCoords);

    MovingLeastSquares mls;
    mls.setAlpha(portAlphaForMLS.getValue());
    mls.setLandmarks(origCoords, shiftedCoords);
    ma::SliceSelector ssh(spatialGraph, "TransformInfo");
    SpatialGraphSelection fullSliceSelection;
    ssh.getSlice(ssh.getSliceAttributeValueFromIndex(1), fullSliceSelection);

    for (int i = 0; i < fullSliceSelection.getNumSelectedEdges(); i++) {
        int edge = fullSliceSelection.getSelectedEdge(i);
        McDArray<McVec3f> newEdgePoints;
        for (int j = 0; j < spatialGraph->getNumEdgePoints(edge); j++) {

            McVec3f edgePoint = spatialGraph->getEdgePoint(edge, j);
            warpPoint(edgePoint, edgePoint, mls);
            newEdgePoints.append(edgePoint);
        }
        spatialGraph->setEdgePoints(edge, newEdgePoints);
    }

    for (int i = 0; i < fullSliceSelection.getNumSelectedVertices(); i++) {
        int curVertex = fullSliceSelection.getSelectedVertex(i);
        McVec3f curCoord = spatialGraph->getVertexCoords(curVertex);
        McVec3f curCoordWarped;
        warpPoint(curCoord, curCoordWarped, mls);
        spatialGraph->setVertexCoords(curVertex, curCoordWarped);
        if (slice.isSelectedVertex(curVertex)) {
            int newVertex = spatialGraph->addVertex(curCoord);
            McDArray<McVec3f> edgePoints(2);
            edgePoints[0] = curCoord;
            edgePoints[1] = curCoordWarped;
            spatialGraph->addEdge(newVertex, curVertex, edgePoints);
        }
    }

    std::cout << "\n Apply deformation took " << watch.stop() << " seconds.";
}

void HxCPDSpatialGraphWarp::warpPoint(const McVec3f& source, McVec3f& target,
                                      MovingLeastSquares& mlsInterpolator) {
    McVec2d curPoint = McVec2d(source.x, source.y);
    McVec2d warpedPoint = mlsInterpolator.interpolate(curPoint);
    target = McVec3f(warpedPoint.x, warpedPoint.y, source.z);
}
