#include <hxalignmicrotubules/MicrotubuleTransformOperation.h>

MicrotubuleTransformOperation::MicrotubuleTransformOperation(
    HxSpatialGraph* sg, const SpatialGraphSelection& selectedElements,
    const SpatialGraphSelection& visibleElements, const SbMatrix mat)
    : Operation(sg, selectedElements, visibleElements), mMat(mat) {
    // include points of non-selected edges, which coincide with a vertex
    mSelection.deselectAllPoints();
    SpatialGraphSelection::Iterator iter(mSelection);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    while (vNum != -1) {
        const IncidenceList& il = graph->getIncidentEdges(vNum);
        for (int i = 0; i < il.size(); ++i) {
            int edgeNum = il[i];
            if (!mSelection.isSelectedEdge(edgeNum)) {
                if (graph->getEdgeSource(edgeNum) == vNum) {
                    mSelection.selectPoint(SpatialGraphPoint(edgeNum, 0));
                }
                if (graph->getEdgeTarget(edgeNum) == vNum) {
                    mSelection.selectPoint(SpatialGraphPoint(
                        edgeNum, graph->getNumEdgePoints(edgeNum) - 1));
                }
            }
        }
        vNum = iter.vertices.nextSelected();
    }
}

MicrotubuleTransformOperation::~MicrotubuleTransformOperation() {}

void MicrotubuleTransformOperation::exec() {
    SpatialGraphSelection::Iterator iter(mSelection);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    while (vNum != -1) {
        McVec3f c = graph->getVertexCoords(vNum);
        SbVec3f t(c.x, c.y, c.z);
        SbVec3f res;
        mMat.multVecMatrix(t, res);
        graph->setVertexCoords(vNum, McVec3f(res[0], res[1], res[2]));
        vNum = iter.vertices.nextSelected();
    }

    iter.edges.reset();
    int eNum = iter.edges.nextSelected();
    while (eNum != -1) {
        McDArray<McVec3f> points = graph->getEdgePoints(eNum);
        for (int p = 0; p < points.size(); p++) {
            SbVec3f t(points[p].x, points[p].y, points[p].z);
            SbVec3f res;
            mMat.multVecMatrix(t, res);
            points[p] = McVec3f(res[0], res[1], res[2]);
        }
        graph->setEdgePoints(eNum, points);
        eNum = iter.edges.nextSelected();
    }
    int numPoints = mSelection.getNumSelectedPoints();
    for (int i = 0; i < numPoints; ++i) {
        SpatialGraphPoint p = mSelection.getSelectedPoint(i);
        McDArray<McVec3f> points = graph->getEdgePoints(p.edgeNum);
        SbVec3f t(points[p.pointNum].x, points[p.pointNum].y,
                  points[p.pointNum].z);
        SbVec3f res;
        mMat.multVecMatrix(t, res);
        points[p.pointNum] = McVec3f(res[0], res[1], res[2]);
        graph->setEdgePoints(p.edgeNum, points);
    }

    // update the transform parameters
    for (int i = 0; i < mTransParams.size(); ++i) {
        appendTransform(mTransParams[i], mMat);
    }
}

void MicrotubuleTransformOperation::undo() {
    SbMatrix invMat = mMat.inverse();
    SpatialGraphSelection::Iterator iter(mSelection);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    while (vNum != -1) {
        McVec3f c = graph->getVertexCoords(vNum);
        SbVec3f t(c.x, c.y, c.z);
        SbVec3f res;
        invMat.multVecMatrix(t, res);
        graph->setVertexCoords(vNum, McVec3f(res[0], res[1], res[2]));
        vNum = iter.vertices.nextSelected();
    }

    iter.edges.reset();
    int eNum = iter.edges.nextSelected();
    while (eNum != -1) {
        McDArray<McVec3f> points = graph->getEdgePoints(eNum);
        for (int p = 0; p < points.size(); p++) {
            SbVec3f t(points[p].x, points[p].y, points[p].z);
            SbVec3f res;
            invMat.multVecMatrix(t, res);
            points[p] = McVec3f(res[0], res[1], res[2]);
        }
        graph->setEdgePoints(eNum, points);
        eNum = iter.edges.nextSelected();
    }
    int numPoints = mSelection.getNumSelectedPoints();
    for (int i = 0; i < numPoints; ++i) {
        SpatialGraphPoint p = mSelection.getSelectedPoint(i);
        McDArray<McVec3f> points = graph->getEdgePoints(p.edgeNum);
        SbVec3f t(points[p.pointNum].x, points[p.pointNum].y,
                  points[p.pointNum].z);
        SbVec3f res;
        invMat.multVecMatrix(t, res);
        points[p.pointNum] = McVec3f(res[0], res[1], res[2]);
        graph->setEdgePoints(p.edgeNum, points);
    }

    // update the transform parameters
    for (int i = 0; i < mTransParams.size(); ++i) {
        appendTransform(mTransParams[i], invMat);
    }
}

SpatialGraphSelection
MicrotubuleTransformOperation::getVisibleSelectionAfterOperation() const {
    return mVisibleSelection;
}

SpatialGraphSelection
MicrotubuleTransformOperation::getHighlightSelectionAfterOperation() const {
    return mSelection;
}

void MicrotubuleTransformOperation::setTransformParameterList(
    const McDArray<HxParameter*> transformParameters) {
    mTransParams = transformParameters;
}

void MicrotubuleTransformOperation::appendTransform(HxParameter* p,
                                                    const SbMatrix& mat) {
    // get existing transform m
    SbMatrix m;
    double res[16];
    if (p->getDimension() == 16) {
        p->getReal(res);
        float* ptr = &m[0][0];
        for (int i = 0; i < 16; ++i) {
            ptr[i] = float(res[i]);
        }
    } else {
        m.makeIdentity();
    }
    // append transform mat and write parameter
    m.multRight(mat);
    p->set(16, &m[0][0]);
}
