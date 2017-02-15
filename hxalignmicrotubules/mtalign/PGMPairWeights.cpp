#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

#include <mclib/McPlane.h>
#include <mclib/McLine.h>
#include <mclib/McMath.h>

namespace ma = mtalign;

ma::PGMPairWeights::PGMPairWeights(const McDArray<McVec3f>& coords1,
                                   const McDArray<McVec3f>& coords2,
                                   const McDArray<McVec3f>& directions1,
                                   const McDArray<McVec3f>& directions2,
                                   ma::PGMPairWeightsParams config)
    : mCoords1(coords1),
      mCoords2(coords2),
      mDirections1(directions1),
      mDirections2(directions2) {
    mConfig = config;
}

double ma::PGMPairWeights::getWeight(const int point1, const int point2) const {
    if (!canPointsBeAPair(point1, point2))
        return 0.0;
    double weight = 1.0;
    if (mConfig.useDist3dWeight)
        weight *= get3dDistanceWeight(point1, point2);
    if (mConfig.useProjectedDistWeight)
        weight *= getProjectedDistanceWeight(point1, point2);
    if (mConfig.useAngleWeight)
        weight *= getAngleWeight(point1, point2);
    return weight;
}

double ma::PGMPairWeights::getDummyWeight() const {
    double weight = 1.0;
    const double dist3dDummy =
        -1.0 * mConfig.dist3dParam * log(mConfig.dummySignificance);
    const double distProjectedDummy =
        -1.0 * mConfig.distProjectedParam * log(mConfig.dummySignificance);
    const double angleDummy =
        -1.0 * mConfig.angleWeightParam * log(mConfig.dummySignificance);
    if (mConfig.useDist3dWeight)
        weight *= computeWeightFromRawValue(dist3dDummy, mConfig.dist3dParam);
    if (mConfig.useProjectedDistWeight)
        weight *= computeWeightFromRawValue(distProjectedDummy,
                                            mConfig.distProjectedParam);
    if (mConfig.useAngleWeight)
        weight *=
            computeWeightFromRawValue(angleDummy, mConfig.angleWeightParam);
    return weight;
}

double ma::PGMPairWeights::get3dDistanceWeight(const int point1,
                                               const int point2) const {
    const double distance = get3dDistance(point1, point2);
    return computeWeightFromRawValue(distance, mConfig.dist3dParam);
}

double ma::PGMPairWeights::getProjectedDistanceWeight(const int point1,
                                                      const int point2) const {
    const double projectedDistance = getProjectedDistance(point1, point2);
    return computeWeightFromRawValue(projectedDistance,
                                     mConfig.distProjectedParam);
}

double ma::PGMPairWeights::getAngleWeight(const int point1,
                                          const int point2) const {
    const double angle = getAngle(point1, point2);
    return computeWeightFromRawValue(angle, mConfig.angleWeightParam);
}

double ma::PGMPairWeights::get3dDistance(const int point1,
                                         const int point2) const {
    return (mCoords1[point1] - mCoords2[point2]).length();
}

double ma::PGMPairWeights::getProjectedDistance(const int point1,
                                                const int point2) const {
    return getMinNormalProjectedDistance(mCoords1[point1], mCoords2[point2],
                                         mDirections1[point1],
                                         mDirections2[point2]);
}

double ma::PGMPairWeights::getAngle(const int point1, const int point2) const {
    double angle = mDirections1[point1].angle(mDirections2[point2]);
    angle = angle / M_PI * 180;
    angle = 180 - angle;
    return angle;
}

bool ma::PGMPairWeights::canPointsBeAPair(const int point1,
                                          const int point2) const {
    const double angle = getAngle(point1, point2);
    if (mConfig.useAngleThreshold && (angle > mConfig.angleThreshold))
        return false;

    const float dist3d = get3dDistance(point1, point2);
    if (mConfig.useDistanceThreshold3d && (dist3d > mConfig.distanceThreshold3d))
        return false;

    const float distProjected = getProjectedDistance(point1, point2);
    if (mConfig.useDistanceThresholdProjected &&
        (distProjected > mConfig.distanceThresholdProjected))
        return false;

    return true;
}

double ma::PGMPairWeights::computeWeightFromRawValue(double rawWeight,
                                                     double param) const {
    double weight = -1;
    if (mConfig.weightType == ma::PGMPairWeightsParams::EXPONENTIAL)
        // We do not multiply the parameter, although in an exponential
        // distribution this should be so.  Exact and greedy matching expect a
        // value <= 1.0.
        weight = exp((double)(-1.0 * rawWeight / param));
    else {
        weight = param - rawWeight;
        if (weight < 0)
            weight = 0;
        weight /= param;
    }
    return weight;
}

float ma::PGMPairWeights::getMinNormalProjectedDistance(const McVec3f coord1,
                                                        const McVec3f coord2,
                                                        const McVec3f dir1,
                                                        const McVec3f dir2) {
    // project coord1 on plane perpendicular to dir2 positioned at coord2
    const float projectedDist12 =
        getNormalProjectedDistance(coord1, coord2, dir1, dir2);
    const float projectedDist21 =
        getNormalProjectedDistance(coord2, coord1, dir2, dir1);
    if (projectedDist12 > projectedDist21)
        return projectedDist21;
    else
        return projectedDist12;
}

float ma::PGMPairWeights::getNormalProjectedDistance(const McVec3f coord1,
                                                     const McVec3f coord2,
                                                     const McVec3f dir1,
                                                     const McVec3f dir2) {
    const McPlane plane(dir2, coord2);
    const McLine line(coord1, coord1 + dir1);
    McVec3f intersectionPoint;
    if (plane.intersect(line, intersectionPoint)) {
        return (coord2 - intersectionPoint).length();
    } else {
        return FLT_MAX;
    }
}
