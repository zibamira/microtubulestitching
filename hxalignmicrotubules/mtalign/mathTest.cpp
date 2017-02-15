#include <hxalignmicrotubules/mtalign/math.h>

#include <gtest/internal/gtest.h>

namespace ma = mtalign;

TEST(mtalign__trace, computesTrace) {
    McDMatrix<double> idMat(10, 10);
    idMat.makeIdentity();
    EXPECT_EQ(ma::trace(idMat), 10);
}
