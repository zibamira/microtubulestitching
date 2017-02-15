#pragma once

#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/SpatialGraphSelection.h>

#include <hxalignmicrotubules/api.h>

/// Transforms a selection containing only edges and vertices by a
/// transformation
class HXALIGNMICROTUBULES_API MicrotubuleTransformOperation : public Operation {
  public:
    MicrotubuleTransformOperation(HxSpatialGraph* sg,
                                  const SpatialGraphSelection& selectedElements,
                                  const SpatialGraphSelection& visibleElements,
                                  const SbMatrix mat);

    virtual ~MicrotubuleTransformOperation();
    virtual void exec();
    virtual void undo();
    virtual SpatialGraphSelection getVisibleSelectionAfterOperation() const;
    virtual SpatialGraphSelection getHighlightSelectionAfterOperation() const;

    void
    setTransformParameterList(const McDArray<HxParameter*> transformParameters);

  protected:
    SbMatrix mMat;
    McDArray<HxParameter*> mTransParams;
    void appendTransform(HxParameter* p, const SbMatrix& mat);
};
