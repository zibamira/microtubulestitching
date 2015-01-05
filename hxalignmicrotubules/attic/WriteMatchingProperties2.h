#pragma once

#include <map>
#include <vector>

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatTextN.h>

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

#include <hxalignmicrotubules/SpreadSheetWrapper.h>
#include <hxalignmicrotubules/mtalign/SliceSelector.h>

namespace mtalign {
struct FacingPointSets;
}

/// `WriteMatchingProperties2` might be useful to understand the factor values
/// for matching microtubule endpoints as implemented by
/// `mtalign::PGMPairWeights`
/// Use it for example on `data/spatialgraph/small3SlicesWithWarp.am`.
class WriteMatchingProperties2 : public HxCompModule {

    HX_HEADER(WriteMatchingProperties2);

  public:
    WriteMatchingProperties2();
    ~WriteMatchingProperties2();

    HxPortFloatTextN portDistanceThreshold;
    HxPortFloatTextN portProjectedDistanceThreshold;
    HxPortFloatTextN portAngleThreshold;
    HxPortDoIt portAction;

    virtual void compute();

  private:
    mtalign::SliceSelector* mSelectionHelper;

    int getMatchedVertex(const int vertex);

    SpreadSheetWrapper* getResultSpreadSheet(const std::string& outputName);

    void writeParameters();

    void writeParametersForTwoConsecutiveSlices(const int slice1Idx,
                                                const int slice2Idx);

    void createWeightConfig(mtalign::PGMPairWeightsParams& config);

    void simulateAlignPairFunction(const mtalign::FacingPointSets&,
                                   SpatialGraphSelection& refSelection,
                                   SpatialGraphSelection& transSelection);

    std::vector<int> mPairMap;
    std::map<std::string, int> mOutputNamesMap;

    void initPairMap();

    void writeSingleProbParameters(const int variableIdx,
                                   const int assignmentIdx,
                                   const mtalign::PGMPairWeights& weighter);

    void addResultEntry(const std::string& whichTable, const float value);

    McVec3f getNextPointAtDist(const McDArray<McVec3f> edgePoints,
                               const float maxDist);

    void getDirections(const SpatialGraphSelection& selectedVertices,
                       McDArray<McVec3f>& directions);
};
