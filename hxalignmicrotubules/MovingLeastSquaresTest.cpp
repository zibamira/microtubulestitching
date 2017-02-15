#include <hxalignmicrotubules/MovingLeastSquares.h>

#include <gmock/internal/gmock.h>
#include <gtest/internal/gtest.h>
#include <mclib/internal/McRandom.h>
#include <mclib/McVec3.h>

using ::testing::Lt;
using ::testing::Gt;

namespace {
struct Landmarks {
    McDArray<McVec3f> src;
    McDArray<McVec3f> dst;
};
}

// Create grid with random shifts.
Landmarks randomGrid(const double rbegin, const double rend,
                     const double repsmarks, const double shift,
                     const double sigma) {
    Landmarks marks;
    const int seed = 11832;
    McRandom rand(seed);
    for (float x = rbegin; x < rend; x += repsmarks) {
        for (float y = rbegin; y < rend; y += repsmarks) {
            McVec3f pos(x, y, 0);
            marks.src.append(pos);
            const double delta = shift + sigma * (rand.nextNumber() - 0.5);
            marks.dst.append(pos + McVec3f(delta, delta, 0));
        }
    }
    return marks;
}

TEST(MovingLeastSquares, shouldInterpolateSmoothly) {
    const double rbegin = 0;
    const double rend = 5;
    const double repsmarks = 2.0;  // Spacing for landmarks.
    const double shift = 1;        // Fixed shift.
    const double sigma = 0.2;      // For random shifts.
    const double tolerance = 0.01;

    const Landmarks marks = randomGrid(rbegin, rend, repsmarks, shift, sigma);

    MovingLeastSquares mls;
    mls.setLandmarks(marks.src, marks.dst);

    // Exactly interpolate landmarks.
    for (int i = 0; i < marks.src.size(); i++) {
        const McVec2d pos = McVec2d(marks.src[i][0], marks.src[i][1]);
        const McVec2d tpos = mls.interpolate(pos);
        const McVec2d tdiff = tpos - McVec2d(marks.dst[i][0], marks.dst[i][1]);
        EXPECT_FLOAT_EQ(0, tdiff.length());
    }

    const int seed = 38273;
    McRandom rand(seed);
    // Interpolate points close to landmarks.
    for (int i = 0; i < marks.src.size(); i++) {
        const double d = sigma * (rand.nextNumber() - 0.5);
        const McVec2d delta(d, d);
        const McVec2d pos = McVec2d(marks.src[i][0], marks.src[i][1]) + delta;
        const McVec2d tpos = mls.interpolate(pos);
        const McVec2d tdiff =
            tpos - delta - McVec2d(marks.dst[i][0], marks.dst[i][1]);
        EXPECT_THAT(tdiff.length(), Lt(tolerance));
#if 0  // Useful for debugging.
        printf("tdiff = %f.\n", tdiff.length());
#endif
    }
}

TEST(MovingLeastSquares, isInvertibleByExchangingLandmarks) {
    const double rbegin = 0;
    const double rend = 5;
    const double repsmarks = 2.0;  // Spacing for landmarks.
    const double repstest = 0.1;   // Spacing for test positions.
    const double shift = 1;        // Fixed shift.
    const double sigma = 0.2;      // For random shifts.
    const double tolerance = 0.03;

    const Landmarks marks = randomGrid(rbegin, rend, repsmarks, shift, sigma);

    MovingLeastSquares forward;
    forward.setLandmarks(marks.src, marks.dst);

    MovingLeastSquares backward;
    backward.setLandmarks(marks.dst, marks.src);

    // Test that forward transform does something and forward-backward
    // transform gets back to original position.
    for (double x = rbegin; x < rend; x += repstest) {
        for (double y = rbegin; y < rend; y += repstest) {
            const McVec2d pos(x, y);
            const McVec2d tpos = forward.interpolate(pos);
            const McVec2d ttpos = backward.interpolate(tpos);
            const McVec2d tdiff = tpos - pos;
            const McVec2d ttdiff = ttpos - pos;

            EXPECT_THAT(tdiff.length(), Gt(tolerance));
            EXPECT_THAT(ttdiff.length(), Lt(tolerance));

#if 0  // Useful for debugging.
            printf("tdiff = %f, ttdiff = %f\n", tdiff.length(),
                   ttdiff.length());
#endif
        }
    }
}
