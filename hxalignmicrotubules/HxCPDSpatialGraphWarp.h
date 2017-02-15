#pragma once

#include <vector>

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatSlider.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortToggleList.h>
#include <mclib/internal/McDMatrix.h>
#include <mclib/internal/McDVector.h>
#include <mclib/McVec3.h>

#include <hxalignmicrotubules/api.h>

namespace mtalign {
class CPDLinearAligner;
}
class HxSpatialGraph;
class MovingLeastSquares;
class SpatialGraphSelection;
class HxUniformVectorField3;

/// `HxCPDSpatialGraphWarp` is used in the supplementary data package to apply
/// the registration algorithms.
class HXALIGNMICROTUBULES_API HxCPDSpatialGraphWarp : public HxCompModule {

    HX_HEADER(HxCPDSpatialGraphWarp);

  public:
    struct TransformationOrderInfo {
        int upperSliceIdx;
        int lowerSliceIdx;
        std::vector<int> upperTransformSlices;
        std::vector<int> lowerTransformSlices;
    };

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
    HxPortFloatTextN portCPDResults;

    HxPortFloatSlider portSampleDist;
    HxPortFloatSlider portAlphaForMLS;

    HxPortDoIt portAction;

    static void preparePoints(McDArray<McVec3f>& p1, McDArray<McVec3f>& p2,
                              SpatialGraphSelection& slice2,
                              const HxSpatialGraph* graph);

    void preparePointsAndDirections(McDArray<McVec3f>& p1,
                                    McDArray<McVec3f>& p2,
                                    McDArray<McVec3f>& d1,
                                    McDArray<McVec3f>& d2,
                                    SpatialGraphSelection& slice2,
                                    const HxSpatialGraph* spatialGraph);

    void preparePointsAndDirectionsRigid(McDArray<McVec3f>& p1,
                                         McDArray<McVec3f>& p2,
                                         McDArray<McVec3f>& d1,
                                         McDArray<McVec3f>& d2,
                                         SpatialGraphSelection& slice2,
                                         const HxSpatialGraph* spatialGraph);

    static void warpPoint(const McVec3f& source, McVec3f& traget,
                          MovingLeastSquares& mlsInterpolator);

    void applyNLDeformationToSlice(SpatialGraphSelection& slice,
                                   HxSpatialGraph* spatialGraph,
                                   const McDArray<McVec3f>& origCoords,
                                   const McDArray<McVec3f>& shiftedCoords);

    void applyRigidDeformationToSliceVanMises(
        SpatialGraphSelection& slice, HxSpatialGraph* spatialGraph,
        const mtalign::CPDLinearAligner& deformation,
        const McDMatrix<double>& R, const double s, const McDVector<double>& t);

    void computeRigidVanMises();

  protected:
    HxSpatialGraph* createOutputDataSet();
    HxUniformVectorField3* createOutputVectorDataSet();
    void computeNL();
    void computeRigid();
    void resamplePairs(McDArray<McVec3f>& p1, McDArray<McVec3f>& p2);
};
