#include <hxalignmicrotubules/HxCPDSpatialGraphWarp.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <hxcore/TestingData.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxgtest/hxtesting.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/TestingDevNullRedirect.h>

#if defined(HX_OS_LINUX)
#define TEST_F_LINUX(test_case_name, test_name)                                \
    TEST_F(test_case_name, test_name)
#else
#define TEST_F_LINUX(test_case_name, test_name)                                \
    TEST_F(test_case_name, DISABLED_##test_name)
#endif

namespace ht = hxtesting;

namespace {

class HxCPDSpatialGraphWarpTest : public ::testing::Test {
  public:
    TestingObjectPoolCleaner cleaner;
    TestingData sgdat;
    HxSpatialGraph* sg;
    McHandle<HxCPDSpatialGraphWarp> cpd;

    HxCPDSpatialGraphWarpTest() : sg(0) {}

    virtual void SetUp() {
        sgdat = TestingData("spatialgraph/fullp0p1.am");
        ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
        sg = sgdat.get<HxSpatialGraph>();
        cpd = HxCPDSpatialGraphWarp::createInstance();
        cpd->portData.connect(sg);
    }

    void expectComputesSha1(const char* sha1) {
        cpd->portAction.hit();
        {
            TestingDevNullRedirect silentout(stdout);
            cpd->fire();
        }
        HxSpatialGraph* r = mcinterface_cast<HxSpatialGraph>(cpd->getResult());
        ASSERT_TRUE(r);
        EXPECT_THAT(r, ht::EqDataSha1(sha1));
    }
};

}  // namespace

TEST_F_LINUX(HxCPDSpatialGraphWarpTest, computesBaselineElastic_E5MS) {
    cpd->portMethod.setValue(2);
    expectComputesSha1("edff0ccb0796506a1585068bc0b5b9376a50275b");
}

TEST_F_LINUX(HxCPDSpatialGraphWarpTest, computesBaselineRigidFisherMises_E5MS) {
    cpd->portMethod.setValue(3);
    expectComputesSha1("b3590bc0f4ae16f0f574420d602deaaedc0a857f");
}
