#include <hxalignmicrotubules/mtalign/CPDElasticAligner.h>

#include <gmock/internal/gmock.h>
#include <gtest/internal/gtest.h>
#include <hxcore/internal/TestingData.h>
#include <hxcore/internal/TestingObjectPoolCleaner.h>
#include <hxgtest/internal/hxtesting.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <hxalignmicrotubules/MicrotubuleSpatialGraphAligner.h>
#include <hxalignmicrotubules/mtalign.h>
#include <hxalignmicrotubules/MovingLeastSquares.h>

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
    params.useAbsoluteValueForEndPointRegion = false;

    // From GUI 'Projection' (String).
    params.projectionType = ma::P_ORTHOGONAL;

    // From default value in ctor
    // MicrotubuleSpatialGraphAligner::mNumMaxPointsForInitialTransform.
    params.numMaxPointsForInitTransform = 50;

    // From GUI 'Approx. from dist'.
    params.maxDistForAngle = 2000;

    // Seems to be uninitialized when called from
    // MicrotubuleSpatialGraphAligner.
    params.angleToPlaneFilter = 0;
    return params;
}

static void testingPrint(QString msg) {
    // Use the following line to temporarily enable printing.
    // puts(qPrintable(msg));
}

static ma::Context makeTestingContext() {
    ma::Context ctx;
    ctx.print = &testingPrint;
    return ctx;
}

static ma::Context testingContext = makeTestingContext();

// Restrict test to Linux, because sha1 on Windows differs.
TEST_LINUX(DISABLED_mtalign__CPDElasticAlignerAccTest, shouldComputeKnownResult_E5MS) {
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
    mtalign::CPDElasticAligner cpd;
    cpd.setContext(&testingContext);;
    cpd.params.beta = 10;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    cpd.params.eDiffRelStop = 1e-5;
    cpd.params.sigmaSquareStop = 1e-7;
    cpd.params.maxIterations = 200;
    cpd.params.useDirections = 1;
    const float alpha = 2;

    for (int i = 0; i < opt.trans.directions.size(); i++) {
        // Keep `normalize()`, because removing it changes sha1.
        opt.trans.directions[i].normalize();
        opt.trans.directions[i] *= -1;
    }
    cpd.setPoints(opt);

    ma::AlignInfo info;
    const McDArray<McVec3f> transCoords = cpd.align(info);

    EXPECT_THAT(info.eDiffRel, gt::Lt(1e-5));
    EXPECT_THAT(info.sigmaSquare, gt::Lt(0.001));
    EXPECT_THAT(info.kappa, gt::Gt(90.902));
    EXPECT_THAT(info.kappa, gt::Lt(90.903));
    EXPECT_THAT(info.e, gt::Gt(1492.199));
    EXPECT_THAT(info.e, gt::Lt(1492.2));
    EXPECT_EQ(45, info.numIterations);

    const McDArray<McVec3f> origCoords = opt.trans.positions;
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
