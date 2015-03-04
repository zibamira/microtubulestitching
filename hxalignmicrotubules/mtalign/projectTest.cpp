#include <hxalignmicrotubules/mtalign/project.h>

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <hxcore/TestingData.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxgtest/hxtesting.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxspatialgraph/SpatialGraphSelection.h>
#include <mclib/TestingDevNullRedirect.h>

#include <hxalignmicrotubules/mtalign.h>

#if defined(HX_OS_LINUX)
#define INSTANTIATE_TEST_CASE_P_LINUX(test_case_name, test_name, test_params)  \
    INSTANTIATE_TEST_CASE_P(test_case_name, test_name, test_params)
#else
#define INSTANTIATE_TEST_CASE_P_LINUX(test_case_name, test_name, test_params)  \
    INSTANTIATE_TEST_CASE_P(DISABLED_##test_case_name, test_name, test_params)
#endif

namespace ht = hxtesting;
namespace ma = mtalign;

static ma::EndPointParams makeParams() {
    ma::EndPointParams params;
    params.refSliceNum = 0;
    params.transSliceNum = 1;
    params.projectionPlane = 0;

    // From GUI 'Boundary (%)'.
    params.endPointRegion = 40;

    // From GUI 'Use absolute value'.
    params.useAbsoluteValueForEndPointRegion = false;

    // From GUI 'Projection' (String) mapped to enum.
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

namespace {

struct Sha1s {
    const char* refCoords;
    const char* refDirections;
    const char* transCoords;
    const char* transDirections;
};

struct Expectation {
    const char* filename;
    ma::EndPointParams params;
    Sha1s sha1s;
};

std::ostream& operator<<(std::ostream& os, const Expectation& e) {
    os << "{ ";
    {
        os << "filename: ";
        os << "'" << e.filename << "'";
        os << ", ";
        os << "projectionType: ";
        os << "'" << e.params.projectionType << "'";
    }
    os << " }";
    return os;
}

class mtalign__projectTestWithTestingData
    : public ::testing::TestWithParam<Expectation> {};

}  // namespace

static void expectSha1(const McDArray<McVec3f>& arr, const char* sha1) {
    int dims[3];
    dims[0] = 3;  // 3 floats per element.
    dims[1] = arr.size();
    dims[2] = 1;
    EXPECT_THAT(&arr[0][0], ht::EqArray3Sha1<float>(dims, sha1));
}

TEST_P(mtalign__projectTestWithTestingData, computesBaseline_E6MS) {
    TestingDevNullRedirect silentout(stdout);
    TestingDevNullRedirect silenterr(stderr);
    Expectation e = GetParam();

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat(e.filename);
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    SpatialGraphSelection refVertexSelection, transVertexSelection;
    ma::SliceSelector selectionHelper(sg, "TransformInfo");
    const int refSliceNum = 0;
    const int transSliceNum = 1;
    selectionHelper.getSlice(refSliceNum + 1, refVertexSelection);
    selectionHelper.getSlice(transSliceNum + 1, transVertexSelection);
    const float projectionPlane =
        selectionHelper.computeMidPlane(refSliceNum, transSliceNum);
    e.params.projectionPlane = projectionPlane;

    ma::FacingPointSets opt = ma::projectEndPoints(sg, e.params);
    expectSha1(opt.ref.positions, e.sha1s.refCoords);
    expectSha1(opt.ref.directions, e.sha1s.refDirections);
    expectSha1(opt.trans.positions, e.sha1s.transCoords);
    expectSha1(opt.trans.directions, e.sha1s.transDirections);
}

std::vector<Expectation> makeExpectations() {
    std::vector<Expectation> es;

    Expectation e;
    e.filename = "spatialgraph/fullp0p1.am";
    e.params = makeParams();
    e.sha1s.refDirections = "8a389ff07b282284a0c19e7965b8f6e213500406";
    e.sha1s.transDirections = "13f6cfb19dbf1fb4349d086ef7b60f93c4f74fed";

    e.params.projectionType = ma::P_ORTHOGONAL;
    e.sha1s.refCoords = "d6cefb5238517d697c9b706e441d0ea57e7845e8";
    e.sha1s.transCoords = "6d86d7bd4242efec27dd54fd14780829ecfb65d1";
    es.push_back(e);

    e.params.projectionType = ma::P_LINEAR;
    e.sha1s.refCoords = "62b6c83c7021c3f8a8be5947ff5cd08eac98c2b8";
    e.sha1s.transCoords = "1574467e28a5f70dec5d62e884315d55c6abbc06";
    es.push_back(e);

    e.params.projectionType = ma::P_TANGENT;
    e.sha1s.refCoords = "d2acb46cf4194afc59d1d922017b5342b053cadb";
    e.sha1s.transCoords = "ab9482357449e87e08a10e26d210b5deb4e43818";
    es.push_back(e);

    e.params.projectionType = ma::P_FIT_0;
    e.sha1s.refCoords = "b24d8efd2cf3164af8728938ed78a38fa11a2755";
    e.sha1s.transCoords = "a5dca2b414c42b18c2b7f0a9606712200dfb56ea";
    es.push_back(e);

    e.params.projectionType = ma::P_FIT_1;
    e.sha1s.refCoords = "ff4b65c6d9379b05600b0396e27405f384e4569f";
    e.sha1s.transCoords = "38112540bb3091920036f3ebc4167fe8a5d50f39";
    es.push_back(e);

    e.params.projectionType = ma::P_FIT_2;
    e.sha1s.refCoords = "26ba469999f374b1f2aecfd2467581a6588a6d14";
    e.sha1s.transCoords = "52294db34ed5bdb1064fee3c2b50b2e3fbc09949";
    es.push_back(e);

    e.params.projectionType = ma::P_FIT_3;
    e.sha1s.refCoords = "6f3257701b9decb0cfea6b52c8117d951cbd24c9";
    e.sha1s.transCoords = "5542eabaabc7b5036c7b2ba9f03d33b837bca13b";
    es.push_back(e);

    e.params.projectionType = ma::P_APPROX_TANGENT;
    e.sha1s.refCoords = "cd1444df24586eb2cd4961bf37ee329d860adaf6";
    e.sha1s.transCoords = "63ca608cdabb8251f338d9175912e8a9715d3432";
    es.push_back(e);

    e.params.projectionType = ma::P_NONE;
    e.sha1s.refCoords = "5dd43e06d3f8a1f9867adbd78a284ca9c4536232";
    e.sha1s.transCoords = "917830282e6220634018236d2f88ed918094dd9f";
    es.push_back(e);

    return es;
}

// Limit tests to Linux, because the sha1 comparison fails on Windows.
INSTANTIATE_TEST_CASE_P_LINUX(WithBaselineForAllProjectionTypes,
                              mtalign__projectTestWithTestingData,
                              testing::ValuesIn(makeExpectations()));
