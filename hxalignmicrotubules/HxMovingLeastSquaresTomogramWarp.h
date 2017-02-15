#pragma once

#include <map>

#include <QString>

#include <Inventor/SbLinear.h>

#include <hxcore/HxCompModule.h>
#include <hxcore/HxConnection.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortMultiMenu.h>

#include <hxalignmicrotubules/MovingLeastSquares.h>
#include <hxalignmicrotubules/api.h>

class HxSpatialGraph;

class HXALIGNMICROTUBULES_API HxMovingLeastSquaresTomogramWarp
    : public HxCompModule {
    HX_HEADER(HxMovingLeastSquaresTomogramWarp);

  public:

    ///
    HxConnection portTomogram;

    ///
    HxPortMultiMenu portTransform;

    ///
    HxPortDoIt portAction;

    ///
    virtual void update();

    ///
    virtual void compute();

    struct TransformInfos {
        QString name;
        QString filename;
        SbMatrix matrix;
        MovingLeastSquares mls;
    };

  private:
    void extractParameters(const HxSpatialGraph* sg);

    // Slice number -> Transformations.
    std::map<int, TransformInfos> transformations;
};
