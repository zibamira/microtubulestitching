#pragma once

#include <vector>

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatSlider.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxfield/HxUniformVectorField3.h>
#include <mclib/McVec2d.h>

#include <hxalignmicrotubules/api.h>

namespace mtalign {
class SliceSelector;
}

class SbMatrix;
class HxSpatialGraph;
class SpatialGraphSelection;
class MovingLeastSquares;

/// `HxMovingLeastSquaresSpatialGraphWarp` applies a moving least squares warp
/// based on the corresponding landmarks in the spatial graph attribute
/// "WarpPairs".  Pairs can be added via the hotkey `CTRL-W` defined by the
/// microtubule align filament editor plugin (see
/// `QxMicrotubuleAlignSpatialGraphTool`).
/// `HxMovingLeastSquaresSpatialGraphWarp` is neither used in the supplementary
/// data package nor by one of our protocols in the spindle project.
class HXALIGNMICROTUBULES_API HxMovingLeastSquaresSpatialGraphWarp
    : public HxCompModule {
    HX_HEADER(HxMovingLeastSquaresSpatialGraphWarp);

  public:
    struct TransformationOrderInfo {
        int upperSliceIdx;
        int lowerSliceIdx;
        std::vector<int> upperTransformSlices;
        std::vector<int> lowerTransformSlices;
    };

    HxMovingLeastSquaresSpatialGraphWarp();

    HxPortMultiMenu portMethod;
    HxPortFloatSlider portAlpha;
    HxPortDoIt portAction;

    virtual void update();

    virtual void compute();

    static void prepareLandmarks(McDArray<McVec2d>& p1, McDArray<McVec2d>& p2,
                                 const HxSpatialGraph* graph,
                                 const int slice1Num, const int slice2Num);

    static void warpPoint(const McVec3f& source, McVec3f& traget,
                          MovingLeastSquares& mlsInterpolator);
    static void
    getSameValueVertices(std::map<int, std::vector<int> >& sameValueVertices,
                         const HxSpatialGraph* graph,
                         const std::string& attributeName);

    static void generateTransformsForSlicesSimple(
        const HxSpatialGraph* spatialGraph, const int slice1Index,
        const int slice2Index, MovingLeastSquares& lowerTransform,
        MovingLeastSquares& upperTransform, int alpha);

    static void applyDeformationToSlice(SpatialGraphSelection& slice,
                                        HxSpatialGraph* spatialGraph,
                                        MovingLeastSquares& deformation);

  protected:
    ~HxMovingLeastSquaresSpatialGraphWarp();
    HxSpatialGraph* createOutputDataSet();
    HxUniformVectorField3* createOutputVectorDataSet();
    void generateOrderOfDeformationSimple(
        const mtalign::SliceSelector& selectionHelper,
        std::vector<TransformationOrderInfo>& orderOfDeformation);
    static bool isValidPair(const std::vector<int>& vertexIndices,
                            const mtalign::SliceSelector& selectionHelper);
    void applyDeformation(TransformationOrderInfo& transformation,
                          HxSpatialGraph* spatialGraph);
};
