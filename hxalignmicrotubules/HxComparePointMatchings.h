#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortMultiMenu.h>

#include <hxalignmicrotubules/api.h>

class HxSpatialGraph;
class EdgeVertexAttribute;

/// `HxComparePointMatchings` is used in the supplementary data package.
class HXALIGNMICROTUBULES_API HxComparePointMatchings : public HxCompModule {

    HX_HEADER(HxComparePointMatchings);

  public:
    HxComparePointMatchings(void);

    ~HxComparePointMatchings(void);

    void compute();

    int parse(Tcl_Interp* t, int argc, char** argv);

    void update();

    HxPortMultiMenu portPairLabel;
    HxPortIntTextN portFalsePositivePairs;
    HxPortIntTextN portFalseNegativePairs;
    HxPortIntTextN portDisagreementPairs;
    HxPortIntTextN portTruePositivePairs;
    HxPortIntTextN portTreeNodes;
    HxPortIntTextN portCriticalNodes;
    HxPortIntTextN portNumberOfManualPairs;
    HxPortIntTextN portNumberOfAutoPairs;

    HxPortDoIt mDoIt;

  private:
    int createUserDefinedEdges(HxSpatialGraph* sg,
                               McDArray<int>& userDefinedEdges);
    int removeUserDefinedMatchingsThatAreNotInAdjacentSlices(
        const HxSpatialGraph* sg, McDArray<int>& allEdgesWithNoTransformInfo);
    void getAllUnmatchedPairs(const HxSpatialGraph* sg,
                              const EdgeVertexAttribute* matchingAttribute,
                              const McDArray<int>& userDefinedMatches,
                              McDArray<McVec2i>& falseNegativePairs);
    void getAllWronglyMatchedPairs(const HxSpatialGraph* sg,
                                   const EdgeVertexAttribute* matchingAttribute,
                                   McDArray<McVec2i>& falsePositivePairs);
    bool hasNodeConnectionToOtherSlice(const HxSpatialGraph* sg,
                                       const int nodeIdx);
    void getAllAutoPairs(const HxSpatialGraph* sg,
                         const EdgeVertexAttribute* matchingAttribute,
                         McDArray<McVec2i>& autoPairs);
    void
    getAllCorrectlyMatchedPairs(const HxSpatialGraph* sg,
                                const EdgeVertexAttribute* matchingAttribute,
                                const McDArray<int>& userDefinedMatches,
                                McDArray<McVec2i>& correctPairs);
    void getAllDisagreeingMatches(const HxSpatialGraph* sg,
                                  const EdgeVertexAttribute* matchingAttribute,
                                  const McDArray<int>& userDefinedMatches,
                                  McDArray<McVec2i>& disagreeingMatches);
    void updatePairLabelPort();
    void rewritePairs(HxSpatialGraph* graph, const McDArray<McVec2i>& newPairs,
                      const McString& pairsName, const bool writeLabels);
    int countLabels(HxSpatialGraph* graph, McString labelName);
};
