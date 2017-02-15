#include <hxalignmicrotubules/HxManualMTAlign.h>

#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/draggers/SoTransformerDragger.h>

#include <hxspatialgraph/internal/HxSpatialGraphInterface.h>

HX_INIT_CLASS(HxManualMTAlign, HxCompModule);

SbRotation HxManualMTAlign::mRot = SbRotation();
SbVec3f HxManualMTAlign::mTrans = SbVec3f();
SbVec3f HxManualMTAlign::mScale = SbVec3f();

HxManualMTAlign::HxManualMTAlign()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      mNodeColorAttribute(0),
      mEdgeColorAttribute(0),
      mNodeScaleAttribute(0),
      mNodeScaleFactor(1.0f),
      mShowNodes(true),
      mShowSegments(true) {}

HxManualMTAlign::~HxManualMTAlign()
{
}

void HxManualMTAlign::compute() {
    if (portData.isNew()) {
        mManualSlice = hxconnection_cast<HxSpatialGraph>(portData);
        mSelectionToRender = SpatialGraphSelection(mManualSlice);

        mManualAlignSep = new SoSeparator;

        mManualTransform = new SoTransform;
        mManualAlignSep->addChild(mManualTransform);

        mTransformedSliceView =
            new SpatialGraphViewNode(mManualSlice->getSpatialGraphInterface());
        mManualAlignSep->addChild(mTransformedSliceView->getSpatialGraphNode());
    }
}

void HxManualMTAlign::startTransform() {
    mcassert(mManualSlice);

    mManip = new SoTransformerManip;
    SoTransformerDragger* dragger =
        dynamic_cast<SoTransformerDragger*>(mManip->getDragger());
    mcassert(dragger);

    // strip all parts dealing with scaling and rotation about x- and y-axes
    dragger->setPart("scale1", new SoSeparator);
    dragger->setPart("scale1active", new SoSeparator);
    dragger->setPart("scale2", new SoSeparator);
    dragger->setPart("scale2active", new SoSeparator);
    dragger->setPart("scale3", new SoSeparator);
    dragger->setPart("scale3active", new SoSeparator);
    dragger->setPart("scale4", new SoSeparator);
    dragger->setPart("scale4active", new SoSeparator);
    dragger->setPart("scale5", new SoSeparator);
    dragger->setPart("scale5active", new SoSeparator);
    dragger->setPart("scale6", new SoSeparator);
    dragger->setPart("scale6active", new SoSeparator);
    dragger->setPart("scale7", new SoSeparator);
    dragger->setPart("scale7active", new SoSeparator);
    dragger->setPart("scale8", new SoSeparator);
    dragger->setPart("scale8active", new SoSeparator);
    dragger->setPart("rotator5", new SoSeparator);
    dragger->setPart("rotator5active", new SoSeparator);
    dragger->setPart("rotator6", new SoSeparator);
    dragger->setPart("rotator6active", new SoSeparator);
    dragger->setPart("zAxisFeedbackActive", new SoSeparator);
    dragger->setPart("zAxisFeedbackSelect", new SoSeparator);
    dragger->setPart("zCrosshairFeedback", new SoSeparator);
    dragger->setPart("xCircleFeedback", new SoSeparator);
    dragger->setPart("yCircleFeedback", new SoSeparator);

    dragger->addMotionCallback(motionCB, this);

    SoSearchAction mySearchAction;
    mySearchAction.setNode(mManualTransform);
    mySearchAction.apply(mManualAlignSep);

    mManip->replaceNode(mySearchAction.getPath());
    showGeom(mManualAlignSep);

    HxManualMTAlign::mRot = dragger->rotation.getValue();
    HxManualMTAlign::mTrans = dragger->translation.getValue();
    HxManualMTAlign::mScale = dragger->scaleFactor.getValue();
}

SbMatrix HxManualMTAlign::stopTransform() {

    hideGeom(mManualAlignSep);

    SoDragger* dragger = mManip->getDragger();
    dragger->removeMotionCallback(motionCB, this);

    SbMatrix mat, invMat;
    mManip->getTranslationSpaceMatrix(mat, invMat);

    SoSearchAction mySearchAction;
    mySearchAction.setNode(mManip);
    mySearchAction.apply(mManualAlignSep);

    mManip->replaceManip(mySearchAction.getPath(), mManualTransform);

    return mat;
}

void HxManualMTAlign::showNodes(const bool show) {
    mShowNodes = show;
    createSlice();
}

void HxManualMTAlign::showSegments(const bool show) {
    mShowSegments = show;
    createSlice();
}

void HxManualMTAlign::createSlice() {
    mSelectionToRender.clear();
    if (mShowNodes) {
        mSelectionToRender.selectAllVertices();
    }
    if (mShowSegments) {
        mSelectionToRender.selectAllEdges();
    }

    mTransformedSliceView->makeVertices(&mSelectionToRender);
    mTransformedSliceView->makeLines(&mSelectionToRender);
    render();
}

void HxManualMTAlign::render() {
    mTransformedSliceView->setVertexColor(&mSelectionToRender,
                                          mNodeColorAttribute);
    mTransformedSliceView->setVertexScale(
        &mSelectionToRender, mNodeScaleAttribute, mNodeScaleFactor);
    mTransformedSliceView->setLineColor(&mSelectionToRender,
                                        mEdgeColorAttribute);
    mTransformedSliceView->setLineWidth(mSegmentWidth);
}

void
HxManualMTAlign::setNodeColorAttribute(const EdgeVertexAttribute* colorAtt) {
    mNodeColorAttribute = colorAtt;
    render();
}

void
HxManualMTAlign::setNodeScaleAttribute(const EdgeVertexAttribute* scaleAtt) {
    mNodeScaleAttribute = scaleAtt;
    render();
}

void HxManualMTAlign::setNodeScaleFactor(const float factor) {
    mNodeScaleFactor = factor;
    render();
}

void
HxManualMTAlign::setEdgeColorAttribute(const EdgeVertexAttribute* colorAtt) {
    mEdgeColorAttribute = colorAtt;
    render();
}

void HxManualMTAlign::setSegmentWidth(const float width) {
    mSegmentWidth = width;
    render();
}

void HxManualMTAlign::motionCB(void* userData, SoDragger* dragger) {
    SoTransformerDragger* d = dynamic_cast<SoTransformerDragger*>(dragger);
    if (d) {
        // constrain translation to xy-plane
        SbVec3f t = d->translation.getValue();
        d->translation.setValue(t[0], t[1], 0.0f);

        // constrain rotation to rotation around z-axis
        SbVec3f axis;
        float angle;
        d->rotation.getValue(axis, angle);
        float tol = 0.0001f;
        if (axis.equals(SbVec3f(0.0f, 0.0f, 1.0f), tol) ||
            axis.equals(SbVec3f(0.0f, 0.0f, -1.0f), tol)) {
            // rotation is valid (around z-axis), store state
            HxManualMTAlign::mRot = d->rotation.getValue();
            HxManualMTAlign::mTrans = d->translation.getValue();
            HxManualMTAlign::mScale = d->scaleFactor.getValue();
        } else {
            // rotation is NOT around z-axis, restore previous valid state
            d->rotation.setValue(HxManualMTAlign::mRot);
            d->translation.setValue(HxManualMTAlign::mTrans);
            d->scaleFactor.setValue(HxManualMTAlign::mScale);
        }
    }
}

