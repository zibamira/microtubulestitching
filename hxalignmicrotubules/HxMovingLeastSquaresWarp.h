#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatSlider.h>
#include <mclib/McVec2d.h>

#include <hxalignmicrotubules/api.h>

class HxUniformScalarField3;
class HxUniformVectorField3;

/// `HxMovingLeastSquaresWarp` applies a moving least squares warp of the input
/// image onto the given second image using `MovingLeastSquares`.
/// `HxMovingLeastSquaresWarp` is neither used in the supplementary data
/// package nor by one of our protocols in the spindle project.
class HXALIGNMICROTUBULES_API HxMovingLeastSquaresWarp : public HxCompModule {
    HX_HEADER(HxMovingLeastSquaresWarp);

  public:
    HxMovingLeastSquaresWarp();

    HxConnection portFromImage;
    HxConnection portToImage;
    HxPortMultiMenu portMethod;
    HxPortFloatSlider portAlpha;
    HxPortDoIt portAction;

    virtual void update();
    virtual void compute();

  protected:
    ~HxMovingLeastSquaresWarp();

    HxUniformScalarField3* createOutputDataSet();
    HxUniformVectorField3* createOutputVectorDataSet();

    void prepareLandmarks(McDArray<McVec2d>& p1, McDArray<McVec2d>& p2);
};
