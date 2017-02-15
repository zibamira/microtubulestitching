#pragma once

#include <vector>

#include <mclib/McDArray.h>
#include <mclib/internal/McDMatrix.h>
#include <mclib/McVec2i.h>
#include <mclib/McVec3.h>

#include <dai/alldai.h>
#include <dai/bipgraph.h>
#include <dai/bp.h>
#include <dai/daialg.h>
#include <dai/factor.h>
#include <dai/factorgraph.h>

#include <hxalignmicrotubules/api.h>

class HXALIGNMICROTUBULES_API BruteForceOptMatching {
  public:
    struct ConnectedFactorGraph {

      public:
        std::vector<dai::TFactor<dai::Real> > factors;
        McDArray<int> variables;
        McDMatrix<int> adjacencyMatrix;
    };

    BruteForceOptMatching(
        const McDArray<McVec3f>& coords1, const McDArray<McVec3f>& coords2,
        const McDArray<McVec3f>& directions1,
        const McDArray<McVec3f>& directions2, const float distanceThreshold3d,
        const float distanceThresholdProjected, const float angleThreshold,
        const McDArray<McVec2i>& evidence, bool projectPoints = false);
    BruteForceOptMatching(const McDArray<McVec3f>& coords1,
                          const McDArray<McVec3f>& coords2,
                          const McDArray<McVec3f>& directions1,
                          const McDArray<McVec3f>& directions2,
                          const float distanceThreshold3d,
                          const float distanceThresholdProjected,
                          const float angleThreshold,
                          bool projectPoints = false);

    ~BruteForceOptMatching(void);

    void createProjectedDistanceMatrix(McDMatrix<float>& projMatrix);

    void create3dDistanceMatrix(McDMatrix<float>& distance3dMatrix);
    void
    createDistanceThreshold3dMatrix(const McDMatrix<float>& distance3dMatrix,
                                    McDMatrix<int>& threshMatrix);
    void createDistanceThresholdProjectedMatrix(
        McDMatrix<int>& threshMatrix,
        const McDMatrix<float>& projectionDistMatrix);
    void createAngleMatrix(McDMatrix<float>& angleMatrix);
    void createAngleThresholdMatrix(McDMatrix<int>& threshMatrix,
                                    const McDMatrix<float>& angleMatrix);

    void createConnectionFactors(const McDMatrix<int>& adjMat,
                                 const McDArray<int> connectedComps,
                                 std::vector<dai::Factor>& pairFactorList);
    void createSingletonFactors(const McDArray<int>& connectedComp,
                                std::vector<dai::Factor>& singletonFactorList);

    static void
    createAdjacenceMatrix(const McDMatrix<int>& variableAssignmentMat,
                          McDMatrix<int>& adjacenceMatrix);
    static bool getConnectedComponent(
        const McDMatrix<int>& adjacenceMatrix,
        McDMatrix<int>& adjacenceMatrixWithoutConnctedComponent,
        McDArray<int>& connComp);
    void getAllConnectedComponents(std::vector<ConnectedFactorGraph>& graphs);

    void fillSingletonVals(const McDMatrix<float>& angleMatrix,
                           const McDMatrix<float>& projDistanceMatrix,
                           const McDMatrix<float>& distanceMatrix3d,
                           const McDMatrix<int>& variableAssignmentMat,
                           std::vector<dai::Factor>& singletonFactors,
                           ConnectedFactorGraph& graph);
    void getSingletonProbs(const McDMatrix<float>& angleMatrix,
                           const McDMatrix<float>& projDistanceMatrix,
                           const McDMatrix<float>& distanceMatrix3d,
                           const McDMatrix<int>& variableAssignmentMat,
                           const McDArray<int>& assignments, const int varLabel,
                           McDArray<float>& probs);
    void getSingletonProbs(const McDMatrix<float>& angleMatrix,
                           const McDMatrix<float>& projDistanceMatrix,
                           const McDMatrix<float>& distanceMatrix3d,
                           const McDMatrix<int>& variableAssignmentMat,
                           const McDArray<int>& assignments, const int varLabel,
                           McDArray<double>& probs);
    void fillPairVals(const McDMatrix<int>& variableAssignmentMat,
                      const McDMatrix<float>& projDistanceMatrix,
                      std::vector<dai::Factor>& pairFactors,
                      ConnectedFactorGraph& graph);

    void getAssignedValuesForVar(const McDMatrix<float>& allValues,
                                 const McDMatrix<int>& variableAssignmentMat,
                                 const McDArray<int>& possibleAssignemnts,
                                 const int label, const float zeroVal,
                                 McDArray<float>& values);
    void computeAngleProbs(McDArray<float>& angles);
    void computeProjDistProbs(McDArray<float>& dists);
    void compute3dDistProbs(McDArray<float>& dists);
    void createSameShiftMatrix(const ConnectedFactorGraph& graph,
                               const dai::Factor& curFac,
                               const McDMatrix<float>& projecedDistMatrix,
                               McDMatrix<float>& shiftEntries);

    bool matchPoints(ConnectedFactorGraph& graph, const int numIter,
                     const float maxDiff, const float damping,
                     McDArray<McVec2i>& matchedPointPairs,
                     McDArray<int>& ambiguousAssignments);

    void getMaxProbAssignments(const dai::BP& ia, const dai::FactorGraph& fg,
                               const ConnectedFactorGraph& graph,
                               McDArray<McVec2i>& pairs);
    void checkAmbiguities(const dai::BP& ia, const dai::FactorGraph& fg,
                          const ConnectedFactorGraph& graph,
                          McDArray<int>& ambiguities);
    void
    checkAmbiguitiesInAssignments(const ConnectedFactorGraph& graph,
                                  const McDArray<McVec2i>& matchedPointPairs,
                                  McDArray<int>& ambiguities);
    void createFactorGraphFactors(const McDMatrix<int>& variableAssignmentMat,
                                  const McDMatrix<int>& variableAdjacenceMatrix,
                                  const McDMatrix<float>& angleMatrix,
                                  const McDMatrix<float>& projDistanceMatrix,
                                  const McDMatrix<float>& distance3dMatrix,
                                  ConnectedFactorGraph& factorGraph);

    const McDMatrix<int>& getVariableAssignmentMat() {
        return mVariableAssignmentMat;
    };

    McDArray<McVec3f> mCoords1;
    McDArray<McVec3f> mCoords2;
    McDArray<McVec3f> mDirections1;
    McDArray<McVec3f> mDirections2;
    McDArray<McVec2i> mEvidence;
    McDArray<McVec2i> mQueerEvidence;

    McDMatrix<int> mVariableAssignmentMat;

    float mDistanceThreshold3d;
    float mDistanceThresholdProjected;
    float mAngleThreshold;

    float getMinNormalProjectedDistance(const McVec3f coord1,
                                        const McVec3f coord2,
                                        const McVec3f dir1, const McVec3f dir2);
    float getNormalProjectedDistance(const McVec3f coord1, const McVec3f coord2,
                                     const McVec3f dir1, const McVec3f dir2);

    int getEvidenceAssignment(const ConnectedFactorGraph& graph,
                              const int varLabel, int& evidenceIndexInGraph);
    void projectToPlaneApproxDirection(const McDArray<McVec3f>& vertices,
                                       const McDArray<McVec3f>& directions,
                                       const float planeZ,
                                       McDArray<McVec3f>& result);
    float getMedianZ(const McDArray<McVec3f>& vertices);
    void preparePoints(void);
    void get3dShiftProbs(const McVec3f pointCoord1, const McVec3f pointCoord2,
                         const McDArray<McVec3f>& assignmentCoordsVar1,
                         const McDArray<McVec3f>& assignmentCoordsVar2,
                         McDMatrix<float>& probs);

    void outputSingleFactorValues(const ConnectedFactorGraph& graph);
    void outputDoubleFactorValues(const ConnectedFactorGraph& graph);
    void handleEvidenceAssignment(const int varLabel,
                                  const int assignmentInEvidence,
                                  const int assignemtnIndexInWholeModel);

    // computes the index of the entry in coords2 (the possible assignments) for
    // the assigment i for this variable
    int mapVariableAssignmentToIndexInVertexList(const int variableLabel,
                                                 const int assigment);
    int mapIndexInVertexListToVariableAssignment(const int variableLabel,
                                                 const int indexInCoords2Array);
    void getAssignmentsForVariable(const int variableLabel,
                                   McDArray<int>& possibleAssignments);
};
