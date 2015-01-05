#pragma once

#include <Inventor/SbLinear.h>
#include <Inventor/manips/SoTransformerManip.h>

#include <hxcore/HxCompModule.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxspatialgraph/SpatialGraphViewNode.h>

#include <hxalignmicrotubules/api.h>

/// `HxManualMTAlign` is used by the filament editor toolbox
/// `QxMicrotubuleAlignSpatialGraphTool` to implement interactive editing of
/// section transformations.
class HXALIGNMICROTUBULES_API HxManualMTAlign : public HxCompModule {
    HX_HEADER(HxManualMTAlign);

  public:
    HxManualMTAlign();
    virtual void compute();

    void startTransform();
    SbMatrix stopTransform();

    void showNodes(const bool show);
    void showSegments(const bool show);

    void setNodeColorAttribute(const EdgeVertexAttribute* colorAtt);
    void setNodeScaleAttribute(const EdgeVertexAttribute* scaleAtt);
    void setNodeScaleFactor(const float factor);
    void setEdgeColorAttribute(const EdgeVertexAttribute* colorAtt);
    void setSegmentWidth(const float width);

  private:
    void createSlice();
    void render();

    McHandle<SoSeparator> mManualAlignSep;
    McHandle<SoTransform> mManualTransform;
    McHandle<SoTransformerManip> mManip;
    McHandle<HxSpatialGraph> mManualSlice;
    McHandle<SpatialGraphViewNode> mTransformedSliceView;
    SpatialGraphSelection mSelectionToRender;

    const EdgeVertexAttribute* mNodeColorAttribute;
    const EdgeVertexAttribute* mEdgeColorAttribute;
    const EdgeVertexAttribute* mNodeScaleAttribute;
    float mNodeScaleFactor;
    float mSegmentWidth;

    bool mShowNodes;
    bool mShowSegments;

    static SbRotation mRot;
    static SbVec3f mTrans;
    static SbVec3f mScale;

    static void motionCB(void* userData, SoDragger* dragger);
};
