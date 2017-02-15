#include <hxalignmicrotubules/mtalign/NullPointRepresentation.h>

#include <mclib/internal/McAssert.h>

namespace ma = mtalign;

ma::NullPointRepresentation::~NullPointRepresentation() {}

int ma::NullPointRepresentation::getNumPoints() const {
    mcerror("Missing implementation");
    return 0;
}

void ma::NullPointRepresentation::computeTransformation(
    const PointRepresentation* pointRep2,
    const PointMatchingDataStruct* pointMatching,
    Transformation* pointRep2Transform) const {
    mcerror("Missing implementation");
}

void ma::NullPointRepresentation::setScoreDivisor(
    const PointRepresentation* pointRep2,
    PointMatchingScoringFunction* scoringFunction) const {
    mcerror("Missing implementation");
}

CacheObject* ma::NullPointRepresentation::startGreedyPointMatching(
    const PointRepresentation* pointRep2,
    const Transformation* pointRep2StartTransform) const {
    mcerror("Missing implementation");
    return 0;
}

CacheObject* ma::NullPointRepresentation::startExactPointMatching(
    const PointRepresentation* pointRep2,
    const Transformation* pointRep2StartTransform) const {
    mcerror("Missing implementation");
    return 0;
}

bool ma::NullPointRepresentation::getNextMatchingPair(CacheObject* cacheObject,
                                                      int& refPointIx,
                                                      int& queryPointIx,
                                                      double& dist2) const {
    mcerror("Missing implementation");
    return false;
}

const McDArray<McVec3f>& ma::NullPointRepresentation::getCoords() const {
    mcerror("Missing implementation");
    static McDArray<McVec3f> dummy;
    return dummy;
}

void ma::NullPointRepresentation::computeCorrespondenceGraphVertices(
    const PointRepresentation* pointSet2, McDArray<McVec2i>& corrGraph) const {
    mcerror("Missing implementation");
}

void ma::NullPointRepresentation::computeCorrespondenceGraphEdges(
    const PointRepresentation* pointSet2, const McDArray<McVec2i>& corrGraph,
    McDArray<McBitfield>& connected) const {
    mcerror("Missing implementation");
}

void ma::NullPointRepresentation::finishGreedyPointMatching() const {
    mcerror("Missing implementation");
}

void ma::NullPointRepresentation::finishExactPointMatching() const {
    mcerror("Missing implementation");
}

const McDArray<McDArray<double> >&
ma::NullPointRepresentation::getDistMatrix() const {
    mcerror("Missing implementation");
    static McDArray<McDArray<double> > dummy;
    return dummy;
}

bool ma::NullPointRepresentation::canPointsBeMatched(
    const int pointRep1PointIdx, const PointRepresentation* pointRep2,
    const int pointRep2PointIdx, CacheObject* cacheObject,
    double& squaredEdgeWeight) const {
    mcerror("Missing implementation");
    return false;
}

double ma::NullPointRepresentation::getSquaredEdgeWeight(
    const int pointRep1PointIdx, const PointRepresentation* pointRep2,
    const int pointRep2PointIdx, CacheObject* cacheObject) const {
    mcerror("Missing implementation");
    return 0;
}
