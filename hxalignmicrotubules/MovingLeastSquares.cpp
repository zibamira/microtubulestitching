#include <hxalignmicrotubules/MovingLeastSquares.h>

#include <mclib/McAssert.h>
#include <mclib/McMat2d.h>
#include <mclib/McVec2d.h>

MovingLeastSquares::MovingLeastSquares() : mAlpha(2) {}

void MovingLeastSquares::setLandmarks(const McDArray<McVec2d>& src,
                                      const McDArray<McVec2d>& dst) {
    mcrequire(src.size() == dst.size());
    mPs = src;
    mQs = dst;
}

static McDArray<McVec2d> asMcVec2d(const McDArray<McVec3f>& xsF) {
    McDArray<McVec2d> xs;
    xs.resize(xsF.size());
    for (int i = 0; i < xsF.size(); i++) {
        xs[i] = McVec2d(xsF[i].x, xsF[i].y);
    }
    return xs;
}

void MovingLeastSquares::setLandmarks(const McDArray<McVec3f>& src,
                                      const McDArray<McVec3f>& dst) {
    mcrequire(src.size() == dst.size());
    setLandmarks(asMcVec2d(src), asMcVec2d(dst));
}

void MovingLeastSquares::setAlpha(const double alpha) {
    mcrequire(alpha > 0);
    mAlpha = alpha;
}

static McMat2d transpose(McMat2d& mat) {
    return McMat2d(mat[0][0], mat[1][0], mat[0][1], mat[1][1]);
}

namespace {
struct Weights {
    McDArray<double> weights;
    double sumOfWeights;
    int atIdx;
};
}  // namespace

// Computes the q or p weights as in the equation after equation 1.
static Weights computeWeights(const McDArray<McVec2d>& points,
                              const McVec2d& point, double alpha) {
    Weights cw;
    cw.atIdx = -1;
    cw.weights.resize(points.size());
    cw.sumOfWeights = 0.0;
    float c = 0.0;
    for (int i = 0; i < points.size(); i++) {
        const double dist = (points[i] - point).length();
        if (dist <= 1.e-10) {
            cw.atIdx = i;
            return cw;
        }
        cw.weights[i] = 1.0 / pow(dist, (double)2.0 * alpha);
        // Use Kahan summation to reduce numerical error.
        const float y = cw.weights[i] - c;
        const float t = cw.sumOfWeights + y;
        c = (t - cw.sumOfWeights) - y;
        cw.sumOfWeights = t;
    }
    return cw;
}

static McVec2d computeCentroid(const McDArray<McVec2d>& points,
                               const Weights& weights) {
    McVec2d centroid(0, 0);
    for (int i = 0; i < points.size(); i++) {
        centroid += weights.weights[i] * McVec2d(points[i].x, points[i].y);
    }
    centroid /= weights.sumOfWeights;
    return centroid;
}

// Compute the Ai matrices according to equation (7).
static McDArray<McMat2d> computeAis(const McVec2d& point,
                                    const McDArray<McVec2d>& ps,
                                    const McDArray<double>& weights,
                                    const McVec2d& pCentroid) {
    McDArray<McMat2d> ais;
    ais.resize(ps.size());
    for (int i = 0; i < ais.size(); i++) {
        const McVec2d piHat = ps[i] - pCentroid;
        const McVec2d piHatOrtho = McVec2d(-1.0 * piHat.y, piHat.x);

        McMat2d firstMat(piHat, -1 * piHatOrtho);
        const McVec2d pointMinusPCentroid = point - pCentroid;
        const McVec2d pointMinusPCentroidOrtho =
            McVec2d(-1 * pointMinusPCentroid.y, pointMinusPCentroid.x);
        const McMat2d secondMat(pointMinusPCentroid,
                                -1 * pointMinusPCentroidOrtho);

        firstMat = transpose(firstMat);
        ais[i] = firstMat * secondMat;
        ais[i] *= weights[i];
    }
    return ais;
}

McVec2d MovingLeastSquares::interpolate(const McVec2d& point) const {
    if (mPs.size() == 0)
        return point;

    const Weights weights = computeWeights(mPs, point, mAlpha);
    if (weights.atIdx >= 0) {
        return mQs[weights.atIdx];
    }
    const McVec2d pCentroid = computeCentroid(mPs, weights);
    const McVec2d qCentroid = computeCentroid(mQs, weights);

    // Compute Ai's as in equation (7).
    const McDArray<McMat2d> ais =
        computeAis(point, mPs, weights.weights, pCentroid);

    // Equation before equation 8.
    McVec2d frVec(0, 0);
    for (int i = 0; i < mPs.size(); i++) {
        McVec2d qVec;
        ais[i].multVecMatrix((mQs[i] - qCentroid), qVec);
        frVec += qVec;
    }

    // Equation 8, compute final new position.
    const McVec2d finalPosition =
        (point - pCentroid).length() / frVec.length() * frVec + qCentroid;
    return finalPosition;
}
