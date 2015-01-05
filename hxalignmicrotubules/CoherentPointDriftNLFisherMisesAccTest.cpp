#include <hxalignmicrotubules/CoherentPointDriftNLFisherMises.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <hxcore/TestingData.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxgtest/hxtesting.h>
#include <hxspatialgraph/HxSpatialGraph.h>

#include <hxalignmicrotubules/MicrotubuleSpatialGraphAligner.h>
#include <hxalignmicrotubules/mtalign.h>

#if defined(HX_OS_LINUX)
#define TEST_LINUX(test_case_name, test_name) TEST(test_case_name, test_name)
#else
#define TEST_LINUX(test_case_name, test_name)                                  \
    TEST(DISABLED_##test_case_name, test_name)
#endif

namespace gt = testing;
namespace ht = hxtesting;
namespace ma = mtalign;

static ma::EndPointParams makeEndPointParams(int refSliceNum, int transSliceNum,
                                             float projectionPlane) {
    ma::EndPointParams params;
    params.refSliceNum = refSliceNum;
    params.transSliceNum = transSliceNum;
    params.projectionPlane = projectionPlane;

    // From GUI 'Boundary (%)'.
    params.endPointRegion = 40;

    // From GUI 'Use absolute value'.
    params.useAbsouteValueForEndPointRegion = false;

    // From GUI 'Projection' (String).
    params.projectionType = "Orthogonal";

    // From default value in ctor
    // MicrotubuleSpatialGraphAligner::mNumMaxPointsForInitialTransform.
    params.numMaxPointsForInitTransform = 50;

    // From GUI 'Transform' (String), mapped through
    // MicrotubuleSpatialGraphAligner::TransformTypes to determine index.
    // 0 is 'Rigid'.
    params.transformType = 0;

    // From GUI 'Approx. from dist'.
    params.maxDistForAngle = 2000;

    // Seems to be uninitialized when called from
    // MicrotubuleSpatialGraphAligner.
    params.angleToPlaneFilter = 0;
    return params;
}

static void dontPrint(const char* msg) {
    // Use the following line to temporarily enable printing.
    // printf("%s", msg);
}

// Restrict test to Linux, because sha1 on Windows differs.
TEST_LINUX(CoherentPointDriftNLFisherMises, shouldComputeKnownResult_E5MS) {
    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/fullp0p1.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    ma::SliceSelector selectionHelper(sg, "TransformInfo");
    const int refSliceNum = 0;
    const int transSliceNum = 1;
    const float midPlane =
        selectionHelper.computeMidPlane(refSliceNum, transSliceNum);

    ma::EndPointParams params =
        makeEndPointParams(refSliceNum, transSliceNum, midPlane);
    ma::FacingPointSets opt = ma::projectEndPoints(sg, params);

    // Configure cpd params.
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    cpd.params.beta = 10;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    cpd.params.eDiffRelStop = 1e-5;
    cpd.params.sigmaSquareStop = 1e-7;
    cpd.params.maxIterations = 200;
    cpd.params.useDirections = 1;
    const float alpha = 2;

    // Set fixed points.
    McDArray<McVec3f> directions = opt.ref.directions;
    McDArray<McVec3f> coords = opt.ref.positions;
    if (cpd.params.useDirections) {
        for (int i = 0; i < directions.size(); i++) {
            directions[i].normalize();
        }
    } else {
        directions.fill(McVec3f(0.0, 0.0, 1.0));
    }
    cpd.convertCoordsToMatrix(coords, cpd.xs);
    cpd.convertDirectionsToMatrix(directions, cpd.xDirs);

    // Set moving points.
    directions = opt.trans.directions;
    coords = opt.trans.positions;
    if (cpd.params.useDirections) {
        for (int i = 0; i < directions.size(); i++) {
            directions[i].normalize();
            directions[i] *= -1;
        }
    } else {
        directions.fill(McVec3f(0.0, 0.0, 1.0));
    }
    cpd.convertCoordsToMatrix(coords, cpd.ys);
    cpd.convertDirectionsToMatrix(directions, cpd.yDirs);

    // Solve.
    McDMatrix<double> G;
    McDMatrix<double> W;
    McDArray<McVec2i> correspondences;
    const mtalign::AlignInfo info = cpd.align(G, W, correspondences);
    McDMatrix<double> transCoordsShiftedM;
    cpd.shiftYs(cpd.ys, G, W, transCoordsShiftedM);
    cpd.rescaleYs(transCoordsShiftedM);

    // Get result.
    McDArray<McVec3f> origCoords = opt.trans.positions;
    const McDArray<McVec3f> transCoordsOld = opt.trans.positions;
    McDArray<McVec3f> transCoords;
    cpd.convertMatrixToCoords(transCoords, transCoordsShiftedM);
    for (int i = 0; i < transCoords.size(); i++) {
        transCoords[i].z = transCoordsOld[i].z;
    }

    EXPECT_THAT(info.eDiffRel, gt::Lt(1e-5));
    EXPECT_THAT(info.sigmaSquare, gt::Lt(0.001));
    EXPECT_THAT(info.kappa, gt::Gt(90.902));
    EXPECT_THAT(info.kappa, gt::Lt(90.903));
    EXPECT_THAT(info.e, gt::Gt(1492.199));
    EXPECT_THAT(info.e, gt::Lt(1492.2));
    EXPECT_EQ(45, info.numIterations);

    MovingLeastSquares mls;
    mls.setAlpha(alpha);
    mls.setLandmarks(origCoords, transCoords);
    SpatialGraphSelection applysel(sg);
    selectionHelper.getSlice(
        selectionHelper.getSliceAttributeValueFromIndex(transSliceNum),
        applysel);
    MicrotubuleSpatialGraphAligner::applyMLSToSelection(mls, sg, applysel);

    EXPECT_THAT(sg, ht::EqDataSha1("90eedf2d24408b82b3a3806222f3464f886a1207"));
}
