#include <hxalignmicrotubules/mtalign/matchingExact.h>

#include <hxcore/HxMessage.h>
#include <mclib/internal/McAssert.h>
#include <mclib/McDArray.h>
#include <mclib/McMath.h>
#include <mclib/McVec3.h>

#include <hxalignmicrotubules/mtalign/NullPointRepresentation.h>
#include <hxalignmicrotubules/mtalign/data.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <pointmatching/ExactPointMatchingAlgorithm.h>
#include <pointmatching/Transformation.h>
#include <pointmatching/PointMatchingDataStruct.h>

namespace ma = mtalign;

namespace {

class RefRepr : public ma::NullPointRepresentation {
  public:
    RefRepr(int nRef, int nTrans, const ma::PGMPairWeights& weightComputation)
        : mNRef(nRef), mNTrans(nTrans), mWeightComputation(weightComputation) {}

  private:
    virtual int getNumPoints() const { return mNRef; }

    virtual void
    setScoreDivisor(const PointRepresentation* pointSet2,
                    PointMatchingScoringFunction* scoringFunction) const {
        scoringFunction->setScoreDivisor(mNRef, mNTrans);
    }

    virtual CacheObject*
    startExactPointMatching(const PointRepresentation* pointRep2,
                            const Transformation* pointRep2Transform) const {
        return 0;
    }

    virtual void finishExactPointMatching() const {}

    virtual bool canPointsBeMatched(const int pointRep1PointIdx,
                                    const PointRepresentation* pointRep2,
                                    const int pointRep2PointIdx,
                                    CacheObject* cacheObject,
                                    double& squaredEdgeWeight) const {
        return mWeightComputation.canPointsBeAPair(pointRep1PointIdx,
                                                   pointRep2PointIdx);
    }

    virtual double getSquaredEdgeWeight(const int pointRep1PointIdx,
                                        const PointRepresentation* pointRep2,
                                        const int pointRep2PointIdx,
                                        CacheObject* cacheObject) const {
        return 1.0 -
               mWeightComputation.getWeight(pointRep1PointIdx,
                                            pointRep2PointIdx);
    }

  private:
    int mNRef;
    int mNTrans;
    const ma::PGMPairWeights& mWeightComputation;
};

class TransRepr : public ma::NullPointRepresentation {
  public:
    TransRepr(int n) : mN(n) {}

  private:
    virtual int getNumPoints() const { return mN; }

    int mN;
};

}  // namespace

static McDArray<McDArray<double> >
makeEdgeWeight(const ma::FacingPointSets& pts,
               const ma::PGMPairWeights& weightComputer) {
    McDArray<McDArray<double> > edgeWeights;
    for (int i = 0; i < pts.ref.positions.size(); i++) {
        McDArray<double> weights;
        for (int j = 0; j < pts.trans.positions.size(); j++) {
            weights.append(1.0 - weightComputer.getWeight(i, j));
        }
        edgeWeights.append(weights);
    }
    return edgeWeights;
}

ma::Matching ma::matchingExact(const ma::FacingPointSets& pts,
                               McMat4f startTransformation,
                               const ma::PGMPairWeightsParams& weightConfig) {

    Transformation st;
    st.setTransformation3d(startTransformation);

    const ma::PGMPairWeights weightComputation(
        pts.ref.positions, pts.trans.positions, pts.ref.directions,
        pts.trans.directions, weightConfig);

    const RefRepr refRepresentation(pts.ref.positions.size(),
                                    pts.trans.positions.size(),
                                    weightComputation);

    const TransRepr transRepresentation(pts.trans.positions.size());

    ExactPointMatchingAlgorithm epma;
    epma.setUserDefinedEdgeWeights(makeEdgeWeight(pts, weightComputation));
    PointMatchingDataStruct pointMatching;
    epma.computeMaxCardinalityPointMatching(
        &refRepresentation, &transRepresentation, &st, &pointMatching);
    Matching matching;
    matching.matchedRefPointIds = pointMatching.getRefPoints();
    matching.matchedTransPointIds = pointMatching.getQueryPoints();
    return matching;
}
