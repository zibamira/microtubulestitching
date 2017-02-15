#include <hxalignmicrotubules/MovingLeastSquares.h>

#include <mclib/internal/McAssert.h>
#include <mclib/McMath.h>
#include <mclib/McMat2d.h>
#include <mclib/McVec2.h>
#include <mclib/McVec3.h>

#include <algorithm>


MovingLeastSquares::Grid2D::Cell2D::Cell2D()
    : mElements(0, 10)
{
}


MovingLeastSquares::Grid2D::Grid2D()
{
    initializeEmptyGrid();
}


void MovingLeastSquares::Grid2D::clear()
{
    int n = mCells.size();

    for (int i = 0; i < n; ++i)
    {
        mCells[i].mElements.clear();
    }
}


int MovingLeastSquares::Grid2D::coordToCell(const McVec2i& cellCoord) const
{
    int x = std::min(std::max(0, cellCoord.x), mDimensions.x - 1);
    int y = std::min(std::max(0, cellCoord.y), mDimensions.y - 1);

    return y * mDimensions.x + x;
}


int MovingLeastSquares::Grid2D::coordToCell(const McVec2d& position) const
{
    McVec2d refPoint = position - mGridStart;

    refPoint.x /= mCellSize.x;
    refPoint.y /= mCellSize.y;

    return coordToCell(McVec2i((int) refPoint.x, (int) refPoint.y));
}


const McVec2d& MovingLeastSquares::Grid2D::getCellSize() const
{
    return mCellSize;
}


int MovingLeastSquares::Grid2D::getNeighbors(McDArray< int >& neighbors, const McVec2d& position, double radius, int minNeighbors) const
{
    neighbors.clear();

    McVec2d refPosition = position - mGridStart;

    int sX = (int)((refPosition.x - radius) / mCellSize.x);
    int sY = (int)((refPosition.y - radius) / mCellSize.y);
    int eX = (int)((refPosition.x + radius) / mCellSize.x);
    int eY = (int)((refPosition.y + radius) / mCellSize.y);

    sX = std::min(mDimensions.x - 1, std::max(0, sX));
    eX = std::min(mDimensions.x - 1, std::max(0, eX));
    sY = std::min(mDimensions.y - 1, std::max(0, sY));
    eY = std::min(mDimensions.y - 1, std::max(0, eY));

    for (int i = sX; i <= eX; ++i)
    {
        for (int j = sY; j <= eY; ++j)
        {
            int cell = coordToCell(McVec2i(i, j));

            neighbors.appendArray(mCells[cell].mElements);
        }
    }

    if (neighbors.size() < minNeighbors && radius <= mMaxRadius)
    {
        return getNeighbors(neighbors, position, radius * 2, minNeighbors) + 1;
    }

    return 1;
}


void MovingLeastSquares::Grid2D::initialize(const McDArray< McVec2d >& points)
{
    int n = points.size();

    if (n <= 0)
    {
        initializeEmptyGrid();
        return;
    }

    mGridStart.setValue(DBL_MAX, DBL_MAX);
    mGridEnd.setValue(-DBL_MAX, -DBL_MAX);

    for (int i = 0; i < n; ++i)
    {
        mGridStart.x = std::min(mGridStart.x, points[i].x);
        mGridStart.y = std::min(mGridStart.y, points[i].y);
        mGridEnd.x   = std::max(mGridEnd.x,   points[i].x);
        mGridEnd.y   = std::max(mGridEnd.y,   points[i].y);
    }

    mGridStart -= McVec2d(1.e-5, 1.e-5);
    mGridEnd   += McVec2d(1.e-5, 1.e-5);
    mCellSize   = mGridEnd - mGridStart;

    double x = mCellSize.x / (mCellSize.x + mCellSize.y);
    double y = mCellSize.y / (mCellSize.x + mCellSize.y);
    double c = (double) n * 2.0;

    double m = sqrt( c / (x * y));

    mDimensions.x = std::max((int)(x * m), 1);
    mDimensions.y = std::max((int)(y * m), 1);

    mCells.resize(mDimensions.x * mDimensions.y);

    clear();

    mCellSize.x = mCellSize.x / (double) mDimensions.x;
    mCellSize.y = mCellSize.y / (double) mDimensions.y;

    mMaxRadius = (mGridEnd - mGridStart).length() * 4.0;

    // fill with points

    for (int i = 0; i < n; ++i)
    {
        int cellID = coordToCell(points[i]);

        mCells[cellID].mElements.push(i);
    }
}


void MovingLeastSquares::Grid2D::initializeEmptyGrid()
{
    mCells.resize(1);
    mCellSize.setValue(1.0, 1.0);
    mDimensions.setValue(1, 1);
    mGridStart.setValue(0.0, 0.0);
    mGridEnd.setValue(1.0, 1.0);
    mMaxRadius = 4.0;
}


MovingLeastSquares::MovingLeastSquares()
    : mAlpha(2)
    , mBoxTransform(false)
    , mScale(1.0)
    , mTranslate(0.0, 0.0) {
}

void MovingLeastSquares::setLandmarks(const McDArray<McVec2d>& src,
                                      const McDArray<McVec2d>& dst) {
    mcrequire(src.size() == dst.size());
    mPs = src;
    mQs = dst;
    mBoxTransform = false;
    updateTransformation();
}

void MovingLeastSquares::setLandmarks(const McDArray<McVec2d>& src,
                                      const McDArray<McVec2d>& dst,
                                      const McBox2f& box)
{
    mcrequire(src.size() == dst.size());
    mBox = box;
    mBoxTransform = true;
    mPs = src;
    mQs = dst;
    updateTransformation();
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

void MovingLeastSquares::swapLandmarks() {
    int n = mPs.size();
    for (int i = 0; i < n; ++i) {
        McVec2d tmp = mPs[i];
        mPs[i] = mQs[i];
        mQs[i] = tmp;
    }
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
        cw.weights[i] =  1.0 / std::pow(dist, 2.0 * alpha);

        // Use Kahan summation to reduce numerical error.
        const double y = cw.weights[i] - c;
        const double t = cw.sumOfWeights + y;
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

void MovingLeastSquares::initializeGrid() {
    mGrid.initialize(mPs);
}

McVec2d MovingLeastSquares::interpolate(const McVec2d& point) const {
    if (mPs.size() == 0)
        return point;

    McVec2d pointT = ((double) 1.0 / mScale) * (point - mTranslate);

    const Weights weights = computeWeights(mPs, pointT, mAlpha);
    if (weights.atIdx >= 0) {
        return (mScale * mQs[weights.atIdx]) + mTranslate;
    }
    const McVec2d pCentroid = computeCentroid(mPs, weights);
    const McVec2d qCentroid = computeCentroid(mQs, weights);

    // Compute Ai's as in equation (7).
    const McDArray<McMat2d> ais =
        computeAis(pointT, mPs, weights.weights, pCentroid);

    // Equation before equation 8.
    McVec2d frVec(0, 0);
    for (int i = 0; i < mPs.size(); i++) {
        McVec2d qVec = (mQs[i] - qCentroid) * ais[i];
        frVec += qVec;
    }

    // Equation 8, compute final new position.
    const McVec2d finalPosition =
        (pointT - pCentroid).length() / frVec.length() * frVec + qCentroid;
    return (mScale * finalPosition) + mTranslate;
}

McVec2d MovingLeastSquares::interpolateFast(const McVec2d& point, int minPoints) const {
    if (mPs.size() < 10) return interpolate(point);

    McVec2d pointT = ((double) 1.0 / mScale) * (point - mTranslate);

    // detect close landmarks
    McDArray< int >     landmarkIDs(0, 100);
    McDArray< McVec2d > p(0, 100);
    McDArray< McVec2d > q(0, 100);

    double radius = 1.0 * std::max(mGrid.getCellSize().x, mGrid.getCellSize().y);

    mGrid.getNeighbors(landmarkIDs, pointT, radius, minPoints);

    for (int i = 0; i < landmarkIDs.size(); ++i){
        p.push(mPs[landmarkIDs[i]]);
        q.push(mQs[landmarkIDs[i]]);
    }

    // classical algo with less points

    const Weights weights = computeWeights(p, pointT, mAlpha);
    if (weights.atIdx >= 0) {
        return (mScale * q[weights.atIdx]) + mTranslate;
    }
    const McVec2d pCentroid = computeCentroid(p, weights);
    const McVec2d qCentroid = computeCentroid(q, weights);

    // Compute Ai's as in equation (7).
    const McDArray<McMat2d> ais =
        computeAis(pointT, p, weights.weights, pCentroid);

    // Equation before equation 8.
    McVec2d frVec(0, 0);
    for (int i = 0; i < p.size(); i++) {
        McVec2d qVec = (q[i] - qCentroid) * ais[i];
        frVec += qVec;
    }

    // Equation 8, compute final new position.
    const McVec2d finalPosition =
        (pointT - pCentroid).length() / frVec.length() * frVec + qCentroid;
    return (mScale * finalPosition) + mTranslate;
}

void MovingLeastSquares::updateTransformation()
{
    if (!mBoxTransform)
    {
        mScale = 1.0;
        mTranslate.setValue(0.0, 0.0);
        return;
    }

    McVec2f d = mBox.getMax() - mBox.getMin();
    double  m = std::max(d.x, d.y);

    if (m == 0.0)
    {
        mScale     = 1.0;
        mTranslate = McVec2d(0.0, 0.0);
    }
    else
    {
        mScale     = m;
        mTranslate = McVec2d(mBox.getMin().x, mBox.getMin().y);
    }

    for (int i = 0; i < mPs.size(); ++i)
    {
        mPs[i] = ((double) 1.0 / mScale) * (mPs[i] - mTranslate);
        mQs[i] = ((double) 1.0 / mScale) * (mQs[i] - mTranslate);
    }
}