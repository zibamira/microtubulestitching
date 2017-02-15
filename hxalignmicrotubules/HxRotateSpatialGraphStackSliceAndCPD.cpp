#include <hxalignmicrotubules/HxRotateSpatialGraphStackSliceAndCPD.h>

#include <hxcore/HxObjectPool.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <mclib/McMat3f.h>

#include <hxalignmicrotubules/HxCPDSpatialGraphWarp.h>
#include <hxalignmicrotubules/SpreadSheetWrapper.h>
#include <hxalignmicrotubules/mtalign/SliceSelector.h>

HX_INIT_CLASS(HxRotateSpatialGraphStackSliceAndCDP, HxCompModule);

HxRotateSpatialGraphStackSliceAndCDP::HxRotateSpatialGraphStackSliceAndCDP()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portMethod(this, "method", tr("Method"), 4),
      portBeta(this, "betaParam", tr("Beta Param")),
      portLambda(this, "lambdaParam", tr("Lambda Param")),
      portW(this, "wParam", tr("W Param")),
      portUseDirection(this, "useDirections", tr("Use Direction"), 1),
      portUseCoords(this, "useCoords", tr("Use Coords"), 1),
      portWithScale(this, "with scaling", tr("With Scaling"), 1),
      portCoupleRcRd(this, "coupleRcAndRd", tr("Couple Rc And Rd"), 1),
      portEMParams(this, "emParams", tr("Em Params"), 3),
      portRotationAngles(this, "rotation(from-by-to)", tr("Rotation (from-by-to)"), 3),
      portAction(this, "action", tr("Action")) {
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
    portUseDirection.setLabel("use orientations");
    portUseCoords.setLabel("use coords (only rigid)");
    portWithScale.setValue(1);
    portRotationAngles.setLabel(0, "from");
    portRotationAngles.setLabel(1, "by");
    portRotationAngles.setLabel(2, "to");
    portRotationAngles.setValue(0, -180);
    portRotationAngles.setValue(1, 90);
    portRotationAngles.setValue(2, 180);
    portCoupleRcRd.setValue(true);
    portCoupleRcRd.hide();
    portWithScale.setValue(1);
    portWithScale.hide();
    portEMParams.setLabel(0, "num iterations");
    portEMParams.setLabel(1, "min sigma");
    portEMParams.setLabel(2, "min likelihood diff");
    portEMParams.setValue(0, 500);
    portEMParams.setValue(1, 1.e-6);
    portEMParams.setValue(2, 1.e-6);
}

HxRotateSpatialGraphStackSliceAndCDP::~HxRotateSpatialGraphStackSliceAndCDP() {}

void HxRotateSpatialGraphStackSliceAndCDP::update() {}

void HxRotateSpatialGraphStackSliceAndCDP::compute() {
    if (!portAction.wasHit())
        return;

    const HxSpatialGraph* inputSG = (HxSpatialGraph*)portData.getSource();

    if (!inputSG)
        return;

    HxCPDSpatialGraphWarp* warpModule = HxCPDSpatialGraphWarp::createInstance();
    warpModule->portMethod.setValue(portMethod.getValue());
    warpModule->portBeta.setValue(portBeta.getValue());
    warpModule->portCoupleRcRd.setValue(portCoupleRcRd.getValue());
    warpModule->portEMParams.setValue(0, portEMParams.getValue(0));
    warpModule->portEMParams.setValue(1, portEMParams.getValue(1));
    warpModule->portEMParams.setValue(2, portEMParams.getValue(2));
    warpModule->portLambda.setValue(portLambda.getValue());
    warpModule->portUseCoords.setValue(portUseCoords.getValue());
    warpModule->portUseDirection.setValue(portUseDirection.getValue());
    warpModule->portW.setValue(portW.getValue());
    warpModule->portWithScale.setValue(portWithScale.getValue());
    warpModule->setLabel("SpatialGraphWarp");
    warpModule->fire();
    theObjectPool->addObject(warpModule);

    HxSpatialGraph* rotatedSG = inputSG->duplicate();
    this->setResult(rotatedSG);
    McString rotatedString;
    McHandle<SpreadSheetWrapper> resultSpreadSheet = SpreadSheetWrapper::createInstance();
    int rowID = 0;
    for (double i = portRotationAngles.getValue(0);
         i < portRotationAngles.getValue(2);
         i += portRotationAngles.getValue(1)) {
        HxSpatialGraph* warped = (HxSpatialGraph*)warpModule->getResult();
        if (warped)
            theObjectPool->removeObject(warped);
        rotatedSG->copyFrom(inputSG);
        rotatedString.printf("rotated%fdegree", i);
        rotatedSG->setLabel(rotatedString.getString());
        rotateSlice(rotatedSG, i / 180.0 * M_PI);
        theObjectPool->addObject(rotatedSG);
        warpModule->portData.connect(rotatedSG);
        warpModule->portAction.hit();
        warpModule->fire();
        theMsg->printf("Rotation: %f. Should be: %f.",
                       warpModule->portCPDResults.getValue(3),
                       i / 180.0 * M_PI);
        McString rowIDString;
        rowIDString.printf("%d", rowID);
        rowID++;

        theMsg->printf("Rotation: %f. Should be: %f.",
                       warpModule->portCPDResults.getValue(3),
                       i / 180.0 * M_PI);
        resultSpreadSheet->addEntry("Rotations", "computedAngle", rowIDString,
                                    warpModule->portCPDResults.getValue(3) *
                                        180.0 / M_PI);
        resultSpreadSheet->addEntry("Rotations", "realAngle", rowIDString, i);
        resultSpreadSheet->addEntry("Rotations", "kappa", rowIDString,
                                    warpModule->portCPDResults.getValue(1));
        resultSpreadSheet->addEntry("Rotations", "sigma^2", rowIDString,
                                    warpModule->portCPDResults.getValue(0));
        resultSpreadSheet->addEntry("Rotations", "NumIter", rowIDString,
                                    warpModule->portCPDResults.getValue(6));
    }
    resultSpreadSheet->setLabel("IterationStatistics");
    setResult(resultSpreadSheet);
}

void HxRotateSpatialGraphStackSliceAndCDP::rotateSlice(HxSpatialGraph* graph,
                                                       const double angle) {
    mtalign::SliceSelector ssh(graph, "TransformInfo");
    SpatialGraphSelection fullSliceSelection;
    ssh.getSlice(ssh.getSliceAttributeValueFromIndex(1), fullSliceSelection);
    McMat3f rotMat = McMat3f::IDENTITY;
    rotMat[0][0] = cos(angle);
    rotMat[0][1] = -sin(angle);
    rotMat[1][0] = sin(angle);
    rotMat[1][1] = cos(angle);
    // compute barycenter
    McVec3f barycenter(0, 0, 0);
    for (int i = 0; i < fullSliceSelection.getNumSelectedVertices(); i++) {
        barycenter +=
            graph->getVertexCoords(fullSliceSelection.getSelectedVertex(i));
    }
    barycenter /= fullSliceSelection.getNumSelectedVertices();
    McVec3f rotatedBarycenter;
    rotMat.multMatrixVec(barycenter, rotatedBarycenter);
    McVec3f translation = barycenter - rotatedBarycenter;

    for (int i = 0; i < fullSliceSelection.getNumSelectedEdges(); i++) {
        int edge = fullSliceSelection.getSelectedEdge(i);
        McDArray<McVec3f> newEdgePoints;
        for (int j = 0; j < graph->getNumEdgePoints(edge); j++) {
            McVec3f edgePoint = graph->getEdgePoint(edge, j);
            McVec3f rotatedEdgePoint;
            rotMat.multMatrixVec(edgePoint, rotatedEdgePoint);
            rotatedEdgePoint += translation;
            newEdgePoints.append(rotatedEdgePoint);
        }
        graph->setEdgePoints(edge, newEdgePoints);
    }

    for (int i = 0; i < fullSliceSelection.getNumSelectedVertices(); i++) {
        int curVertex = fullSliceSelection.getSelectedVertex(i);
        McVec3f curCoord = graph->getVertexCoords(curVertex);
        McVec3f rotatedVertex;
        rotMat.multMatrixVec(curCoord, rotatedVertex);
        rotatedVertex += translation;
        graph->setVertexCoords(curVertex, rotatedVertex);
    }
}
