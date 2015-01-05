#pragma once

#include <pointmatching/PointRepresentation.h>

class PointMatchingDataStruct;
class PointMatchingScoringFunction;

namespace mtalign {

/// `NullPointRepresentation` is used internally by the `matching*()`
/// functions.
class NullPointRepresentation : public PointRepresentation {
  protected:
    virtual ~NullPointRepresentation();

  private:
    virtual int getNumPoints() const;

    virtual void
    computeTransformation(const PointRepresentation* pointRep2,
                          const PointMatchingDataStruct* pointMatching,
                          Transformation* pointRep2Transform) const;

    virtual void
    setScoreDivisor(const PointRepresentation* pointRep2,
                    PointMatchingScoringFunction* scoringFunction) const;

    virtual CacheObject* startGreedyPointMatching(
        const PointRepresentation* pointRep2,
        const Transformation* pointRep2StartTransform) const;

    virtual CacheObject* startExactPointMatching(
        const PointRepresentation* pointRep2,
        const Transformation* pointRep2StartTransform) const;

    virtual bool getNextMatchingPair(CacheObject* cacheObject, int& refPointIx,
                                     int& queryPointIx, double& dist2) const;

    virtual const McDArray<McVec3f>& getCoords() const;

    virtual void
    computeCorrespondenceGraphVertices(const PointRepresentation* pointSet2,
                                       McDArray<McVec2i>& corrGraph) const;

    virtual void
    computeCorrespondenceGraphEdges(const PointRepresentation* pointSet2,
                                    const McDArray<McVec2i>& corrGraph,
                                    McDArray<McBitfield>& connected) const;

    virtual void finishGreedyPointMatching() const;

    virtual void finishExactPointMatching() const;

    virtual const McDArray<McDArray<double> >& getDistMatrix() const;

    virtual bool canPointsBeMatched(const int pointRep1PointIdx,
                                    const PointRepresentation* pointRep2,
                                    const int pointRep2PointIdx,
                                    CacheObject* cacheObject,
                                    double& squaredEdgeWeight) const;

    virtual double getSquaredEdgeWeight(const int pointRep1PointIdx,
                                        const PointRepresentation* pointRep2,
                                        const int pointRep2PointIdx,
                                        CacheObject* cacheObject) const;
};

}  // namespace mtalign
