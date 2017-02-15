#include <hxalignmicrotubules/mtalign/matchingGreedy.h>

#include <mclib/internal/McComparators.h>
#include <mclib/McDArray.h>
#include <mclib/McMath.h>
#include <mclib/McPlane.h>
#include <mclib/internal/McSorter.h>
#include <mclib/McVec3.h>

#include <hxalignmicrotubules/mtalign/NullPointRepresentation.h>
#include <hxalignmicrotubules/mtalign/data.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <pointmatching/GreedyPointMatchingAlgorithm.h>
#include <pointmatching/PointMatchingScoringFunction.h>
#include <pointmatching/Transformation.h>
#include <pointmatching/PointMatchingDataStruct.h>

class CacheObject;

namespace ma = mtalign;

namespace {

class Cache {
  public:
    void setCoordsAndDirections(const McDArray<McVec3f>& coordinates,
                                const McDArray<McVec3f>& directions,
                                const McMat4f& transform) {

        assert(coordinates.size() == directions.size());
        mCoords.resize(coordinates.size());
        mDirections.resize(directions.size());

        // transform coordinates
        for (int i = 0; i < coordinates.size(); i++) {
            McVec3f dir = coordinates[i] + directions[i];
            McVec3f newDir;
            transform.multVecMatrix(dir, newDir);
            transform.multVecMatrix(coordinates[i], mCoords[i]);
            newDir = newDir - mCoords[i];
            newDir.normalize();
            mDirections[i] = newDir;
        }
    }

    const McDArray<McVec3f>& getCoords() const { return mCoords; }

    const McDArray<McVec3f>& getDirections() const { return mDirections; }

    void computeDistancesAndSort(ma::PGMPairWeights weightComputation) {

        pointSet1IsMatched.resize(weightComputation.mCoords1.size());
        pointSet1IsMatched.unsetAll();
        pointSet2IsMatched.resize(weightComputation.mCoords2.size());
        pointSet2IsMatched.unsetAll();
        currListIndex = 0;

        pointCorrespondence.clear();
        pointDist2.clear();

        double dist2;

        for (int i = 0; i < weightComputation.mCoords1.size(); ++i) {
            for (int j = 0; j < weightComputation.mCoords2.size(); ++j) {

                if (weightComputation.canPointsBeAPair(i, j)) {
                    // This clas expects a low value to be a good weight and
                    // large to be bad.
                    dist2 = 1.0 - weightComputation.getWeight(i, j);
                    pointCorrespondence.append(McVec2i(i, j));
                    pointDist2.append(dist2);
                }
            }
        }

        sortedDist2List.resize(pointCorrespondence.size());
        for (int i = 0; i < sortedDist2List.size(); ++i) {
            sortedDist2List[i] = i;
        }

        if (pointDist2.size() > 0) {
            McIndexByValueComparator<double> comparator(&pointDist2[0]);
            sort(&(sortedDist2List[0]), sortedDist2List.size(), comparator);
        }
    }

    bool getNextMatchingPair(int& set1PointIx, int& set2PointIx,
                             double& dist2) {
        if (currListIndex < sortedDist2List.size()) {
            int ix1, ix2;
            do {
                ix1 = pointCorrespondence[sortedDist2List[currListIndex]][0];
                if (pointSet1IsMatched[ix1]) {
                    currListIndex++;
                } else {
                    ix2 =
                        pointCorrespondence[sortedDist2List[currListIndex]][1];
                    if (pointSet2IsMatched[ix2]) {
                        currListIndex++;
                    } else {
                        set1PointIx = ix1;
                        set2PointIx = ix2;
                        dist2 = pointDist2[sortedDist2List[currListIndex]];
                        pointSet1IsMatched.set(ix1);
                        pointSet2IsMatched.set(ix2);
                        return true;
                    }
                }
            } while (currListIndex < sortedDist2List.size());
        }

        return false;
    }

  protected:
    McDArray<McVec3f> mCoords;
    McDArray<McVec3f> mDirections;

  private:
    McDArray<McVec2i> pointCorrespondence;
    McDArray<double> pointDist2;
    McDArray<int> sortedDist2List;
    McBitfield pointSet1IsMatched;
    McBitfield pointSet2IsMatched;
    int currListIndex;
};

class RefRepr : public ma::NullPointRepresentation {
  public:
    RefRepr(const ma::FacingPointSets& pts, const McMat4f& startTransformation,
            ma::PGMPairWeights weightComputation)
        : mNRef(pts.ref.positions.size()), mNTrans(pts.trans.positions.size()) {
        mCache.setCoordsAndDirections(pts.trans.positions, pts.trans.directions,
                                      startTransformation);
        ma::PGMPairWeights newComp(
            weightComputation.mCoords1, mCache.getCoords(),
            weightComputation.mDirections1, mCache.getDirections(),
            weightComputation.mConfig);
        mCache.computeDistancesAndSort(newComp);
    }

    virtual int getNumPoints() const { return mNRef; }

  private:
    virtual void
    setScoreDivisor(const PointRepresentation* pointSet2,
                    PointMatchingScoringFunction* scoringFunction) const {
        scoringFunction->setScoreDivisor(mNRef, mNTrans);
    }

    virtual CacheObject* startGreedyPointMatching(
        const PointRepresentation* pointSet2,
        const Transformation* pointSet2StartTransform) const {
        return 0;
    }

    virtual bool getNextMatchingPair(CacheObject* cacheObject, int& refPointIx,
                                     int& queryPointIx, double& dist2) const {
        return mCache.getNextMatchingPair(refPointIx, queryPointIx, dist2);
    }

    virtual void finishGreedyPointMatching() const {}

  private:
    int mNRef;
    int mNTrans;
    mutable Cache mCache;
};

class TransRepr : public ma::NullPointRepresentation {
  public:
    TransRepr(int n) : mN(n) {}

  private:
    virtual int getNumPoints() const { return mN; }

    int mN;
};
}

ma::Matching ma::matchingGreedy(const ma::FacingPointSets& pts,
                                McMat4f startTransformation,
                                const ma::PGMPairWeightsParams& weightConfig) {

    // carry out point matching for all start transformations
    McHandle<PointMatchingScoringFunction> scoringFunction =
        new PointMatchingScoringFunction();
    scoringFunction->setRMSDScaleFactor(1.e-15);

    GreedyPointMatchingAlgorithm pmAlg;
    pmAlg.setMinPointMatchingSize(3);
    pmAlg.setScoringFunction(scoringFunction.ptr());

    Transformation st;
    st.setTransformation3d(startTransformation);

    const ma::PGMPairWeights weightComputation(
        pts.ref.positions, pts.trans.positions, pts.ref.directions,
        pts.trans.directions, weightConfig);

    const RefRepr refRepresentation(pts, startTransformation,
                                    weightComputation);
    const TransRepr transRepresentation(pts.trans.positions.size());

    PointMatchingDataStruct pointMatching;
    pmAlg.computePointMatching(&refRepresentation, &transRepresentation, &st,
                               &pointMatching);
    Matching matching;
    matching.matchedRefPointIds = pointMatching.getRefPoints();
    matching.matchedTransPointIds = pointMatching.getQueryPoints();
    return matching;
}
