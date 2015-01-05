#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatSlider.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortToggleList.h>

#include <hxalignmicrotubules/api.h>

class HxSpatialGraph;

///
class HXALIGNMICROTUBULES_API HxRotateSpatialGraphStackSliceAndCDP
    : public HxCompModule {
    HX_HEADER(HxRotateSpatialGraphStackSliceAndCDP);

  public:
    HxRotateSpatialGraphStackSliceAndCDP();
    virtual void update();
    virtual void compute();

    HxPortMultiMenu portMethod;
    HxPortFloatSlider portBeta;
    HxPortFloatSlider portLambda;
    HxPortFloatSlider portW;
    HxPortToggleList portUseDirection;
    HxPortToggleList portUseCoords;
    HxPortToggleList portWithScale;
    HxPortToggleList portCoupleRcRd;
    HxPortFloatTextN portEMParams;
    HxPortFloatTextN portRotationAngles;
    HxPortDoIt portAction;

  protected:
    ~HxRotateSpatialGraphStackSliceAndCDP();

    void rotateSlice(HxSpatialGraph* graph, const double angle);
};
