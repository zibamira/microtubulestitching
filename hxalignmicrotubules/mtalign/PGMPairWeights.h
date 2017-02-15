#pragma once

#include <mclib/McDArray.h>
#include <mclib/McVec3.h>

#include <hxalignmicrotubules/api.h>

namespace mtalign {

/// `PGMPairWeightsParams` controls the factor values for matching microtubule
/// endpoints.  See paper at <http://dx.doi.org/10.1371/journal.pone.0113222>.
/// A parameter is ignored if the corresponding `use*` is `false`.
struct HXALIGNMICROTUBULES_API PGMPairWeightsParams {
    /// `distanceThreshold3d` is `t_c` from the paper.
    double distanceThreshold3d;
    bool useDistanceThreshold3d;

    /// `distanceThresholdProjected` is `t_p` from the paper.
    double distanceThresholdProjected;
    bool useDistanceThresholdProjected;

    /// `angleThreshold` is `t_alpha` from the paper.
    double angleThreshold;
    bool useAngleThreshold;

    /// `angleWeightParam` is the inverse weight `lambda_alpha^-1` from the
    /// paper when `weightType` is `EXPONENTIAL`.
    double angleWeightParam;
    bool useAngleWeight;

    /// `dist3dParam` is the inverse weight `lambda_c^-1` from the paper when
    /// `weightType` is `EXPONENTIAL`.
    double dist3dParam;
    bool useDist3dWeight;

    /// `distProjectedParam` is the inverse weight `lambda_p^-1` from the paper
    /// when `weightType` is `EXPONENTIAL`.
    double distProjectedParam;
    bool useProjectedDistWeight;

    /// `dummySignificance` is the placeholder significance parameter `r` from
    /// the paper.
    double dummySignificance;

    /// Use `weightType = EXPONENTIAL` for the algorithms from the paper.
    enum WeightType { LINEAR = 0, EXPONENTIAL };
    WeightType weightType;
};

/// `PGMPairWeights` implements the factor values for matching microtubules
/// endpoints (see paper).
class HXALIGNMICROTUBULES_API PGMPairWeights {
  public:
    PGMPairWeights(const McDArray<McVec3f>& coords1,
                   const McDArray<McVec3f>& coords2,
                   const McDArray<McVec3f>& directions1,
                   const McDArray<McVec3f>& directions2,
                   mtalign::PGMPairWeightsParams config);

    double get3dDistance(const int point1, const int point2) const;
    double get3dDistanceWeight(const int point1, const int point2) const;
    double getAngle(const int point1, const int point2) const;
    double getAngleWeight(const int point1, const int point2) const;
    double getProjectedDistance(const int point1, const int point2) const;
    double getProjectedDistanceWeight(const int point1, const int point2) const;
    bool canPointsBeAPair(const int point1, const int point2) const;
    double computeWeightFromRawValue(const double rawWeight,
                                     const double parameter) const;
    double getWeight(const int point1, const int point2) const;

    McDArray<McVec3f> mCoords1;
    McDArray<McVec3f> mCoords2;
    McDArray<McVec3f> mDirections1;
    McDArray<McVec3f> mDirections2;
    mtalign::PGMPairWeightsParams mConfig;

    static float getMinNormalProjectedDistance(const McVec3f coord1,
                                               const McVec3f coord2,
                                               const McVec3f dir1,
                                               const McVec3f dir2);
    static float getNormalProjectedDistance(const McVec3f coord1,
                                            const McVec3f coord2,
                                            const McVec3f dir1,
                                            const McVec3f dir2);

    double getDummyWeight() const;
};

}  // namespace mtalign
