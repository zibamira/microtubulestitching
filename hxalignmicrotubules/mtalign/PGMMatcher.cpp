#include <hxalignmicrotubules/mtalign/PGMMatcher.h>

#include <QString>

#include <dai/bipgraph.h>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/alldai.h>
#include <dai/bp.h>
#include <dai/factor.h>

using namespace std;
using namespace dai;

namespace ma = mtalign;

ma::PGMMatcher::PGMMatcher(const McDArray<McVec3f>& coords1,
                           const McDArray<McVec3f>& coords2,
                           const McDArray<McVec3f>& directions1,
                           const McDArray<McVec3f>& directions2,
                           const McDArray<McVec2i>& evidence,
                           PGMPairWeights& weightComputationFunction,
                           double shiftProbParam)
    : mCoords1(coords1),
      mCoords2(coords2),
      mDirections1(directions1),
      mDirections2(directions2),
      mEvidence(evidence),
      mShiftProbParam(shiftProbParam),
      mWeightComputationFunction(weightComputationFunction) {
    mContext = &defaultContext();
}

void ma::PGMMatcher::print(QString msg) {
    mContext->print(msg);
}

void ma::PGMMatcher::setContext(ma::Context* ctx) {
    mContext = ctx;
}

void ma::PGMMatcher::createVariableAssignmentMat() {
    mVariableAssignmentMat.resize(mCoords1.size(), mCoords2.size());
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {
            mVariableAssignmentMat[i][j] =
                mWeightComputationFunction.canPointsBeAPair(i, j) ? 1 : 0;
        }
    }
}

void ma::PGMMatcher::initConnectedComponentIteration() {
    mQueerEvidence.resize(0);

    createVariableAssignmentMat();
    createAdjacenceMatrix(mVariableAssignmentMat,
                          mVariableAdjacenceMatrixForConnCompIteration);
}
bool ma::PGMMatcher::getNextConnectedComponent(ConnectedFactorGraph& graph) {
    McDMatrix<int> adjacenceMatrixWithoutConnctedComponent;
    McDArray<int> connComp;

    bool gotConnectedComponent = getConnectedComponent(
        mVariableAdjacenceMatrixForConnCompIteration,
        adjacenceMatrixWithoutConnctedComponent, connComp);

    graph.variables = connComp;
    graph.adjacencyMatrix = mVariableAdjacenceMatrixForConnCompIteration -
                            adjacenceMatrixWithoutConnctedComponent;

    createFactorGraphFactors(mVariableAssignmentMat,
                             mVariableAdjacenceMatrixForConnCompIteration,
                             graph);

    mVariableAdjacenceMatrixForConnCompIteration =
        adjacenceMatrixWithoutConnctedComponent;
    return gotConnectedComponent;
}

// Finds all connected components in adjacenceMatrix, that is, all variable
// labels that are connected in the pgm are after processinf an array in the
// connComps array.

void ma::PGMMatcher::getAllConnectedComponents(
    vector<ConnectedFactorGraph>& graphs) {

    mQueerEvidence.resize(0);

    createVariableAssignmentMat();
    McDMatrix<int> variableAdjacenceMatrix;
    createAdjacenceMatrix(mVariableAssignmentMat, variableAdjacenceMatrix);

    McDMatrix<int> adjacenceMatrixWithoutConnctedComponent;
    McDArray<int> connComp;

    while (getConnectedComponent(variableAdjacenceMatrix,
                                 adjacenceMatrixWithoutConnctedComponent,
                                 connComp)) {
        ConnectedFactorGraph graph;
        graph.variables = connComp;
        graph.adjacencyMatrix =
            variableAdjacenceMatrix - adjacenceMatrixWithoutConnctedComponent;

        createFactorGraphFactors(mVariableAssignmentMat,
                                 variableAdjacenceMatrix, graph);
        graphs.push_back(graph);

        variableAdjacenceMatrix = adjacenceMatrixWithoutConnctedComponent;
    }
}

// Computes all variables that form a connected component in the
// adjacenceMatrix.
// The connected component chosen is arbitrary if variableInConnectedComponent
// is set to -1.
// If variableInConnectedComponent is set, it finds the connected component this
// variable belongs to.
bool ma::PGMMatcher::getConnectedComponent(
    const McDMatrix<int>& adjacenceMatrix,
    McDMatrix<int>& adjacenceMatrixWithoutConnctedComponent,
    McDArray<int>& connComp, int variableInConnectedComponent) {

    adjacenceMatrixWithoutConnctedComponent.resize(adjacenceMatrix.nRows(),
                                                   adjacenceMatrix.nCols());
    memcpy(adjacenceMatrixWithoutConnctedComponent.dataPtr(),
           adjacenceMatrix.dataPtr(),
           sizeof(int) * adjacenceMatrix.nRows() * adjacenceMatrix.nCols());

    // find first a startpoint

    connComp.resize(0);
    if (variableInConnectedComponent == -1) {
        for (int i = 0; i < adjacenceMatrix.nRows(); i++) {
            for (int j = i; j < adjacenceMatrix.nCols(); j++) {
                if (adjacenceMatrix[i][j] == 1) {
                    variableInConnectedComponent = i;
                    break;
                }
            }
        }
    }
    if (variableInConnectedComponent == -1)
        return false;

    McDArray<int> queue;
    queue.append(variableInConnectedComponent);
    connComp.clear();
    while (queue.size() > 0) {
        int cur = queue.last();
        connComp.append(cur);
        queue.pop_back();
        for (int i = 0; i < adjacenceMatrixWithoutConnctedComponent.nCols();
             i++) {
            if (adjacenceMatrixWithoutConnctedComponent[cur][i] == 1) {
                queue.push(i);
                adjacenceMatrixWithoutConnctedComponent[cur][i] = 0;
            }
        }
    }
    // remove duplicates
    connComp.sort(&mcStandardCompare);
    int cur = connComp.last();
    for (int i = connComp.size() - 2; i >= 0; i--) {
        if (cur == connComp[i])
            connComp.remove(i, 1);
        else
            cur = connComp[i];
    }
    return true;
}

// Finds variables that could compete for the same assignment.

void ma::PGMMatcher::createAdjacenceMatrix(
    const McDMatrix<int>& variableAssignmentMat,
    McDMatrix<int>& adjacenceMatrix) {
    // create adjacence matrix for variables
    adjacenceMatrix.resize(variableAssignmentMat.nRows(),
                           variableAssignmentMat.nRows());
    adjacenceMatrix.fill(0);

    // rows represent variables and columns the values that the variables take
    // so we iterate over all rows/variables
    for (int i = 0; i < variableAssignmentMat.nRows(); i++) {
        // and all other variables including itself
        for (int j = i; j < variableAssignmentMat.nRows(); j++) {
            // and see, if they can potentially take the same value at all.
            for (int k = 0; k < variableAssignmentMat.nCols(); k++) {
                if (variableAssignmentMat[i][k] > 1.e-6 &&
                    variableAssignmentMat[j][k] > 1.e-6) {
                    adjacenceMatrix[i][j] = 1;
                    adjacenceMatrix[j][i] = 1;
                }
            }
        }
    }
}

void ma::PGMMatcher::createFactorGraphFactors(
    const McDMatrix<int>& variableAssignmentMat,
    const McDMatrix<int>& variableAdjacenceMatrix,
    ConnectedFactorGraph& factorGraph) {

    vector<Factor> pairFactors;
    createConnectionFactors(variableAdjacenceMatrix, factorGraph.variables,
                            pairFactors);
    vector<Factor> singletonFactors;
    createSingletonFactors(factorGraph.variables, singletonFactors);

    fillSingletonVals(variableAssignmentMat, singletonFactors, factorGraph);
    fillPairVals(variableAssignmentMat, pairFactors, factorGraph);
    // outputSingleFactorValues(factorGraph);
}

int ma::PGMMatcher::mapVariableAssignmentToIndexInVertexList(
    const int variableLabel, const int assignment) {
    McDArray<int> possibleAssignments;
    getAssignmentsForVariable(variableLabel, possibleAssignments);
    // see if the assignment was the dummy variable
    if (assignment == possibleAssignments.size())
        return -1;

    return possibleAssignments[assignment];
}

int ma::PGMMatcher::mapIndexInVertexListToVariableAssignment(
    const int variableLabel, const int indexInCoords2Array) {
    McDArray<int> possibleAssignments;
    getAssignmentsForVariable(variableLabel, possibleAssignments);
    return possibleAssignments.findSorted(indexInCoords2Array,
                                          &mcStandardCompare);
}

void ma::PGMMatcher::getAssignmentsForVariable(
    const int variableLabel,
    McDArray<int>& indicesThisVariableCanBeAssignedTo) {
    indicesThisVariableCanBeAssignedTo.resize(0);
    for (int i = 0; i < mVariableAssignmentMat.nCols(); i++) {
        if (mVariableAssignmentMat[variableLabel][i] > 0)
            indicesThisVariableCanBeAssignedTo.append(i);
    }
}

// creates a pairwise factor for each adjacent variables in one connected
// component
// does not set the values yet
void ma::PGMMatcher::createConnectionFactors(const McDMatrix<int>& adjMat,
                                             const McDArray<int> connectedComp,
                                             vector<Factor>& pairFactorList) {
    for (int i = 0; i < connectedComp.size(); i++) {
        int curVar = connectedComp[i];
        McDArray<int> assigmentsForIthVariable;
        getAssignmentsForVariable(curVar, assigmentsForIthVariable);
        Var pi(curVar, assigmentsForIthVariable.size() + 1);
        for (int j = curVar + 1; j < adjMat.nCols(); j++) {
            if (adjMat[curVar][j] == 1) {
                McDArray<int> assigmentsForJthVariable;
                getAssignmentsForVariable(j, assigmentsForJthVariable);

                Var pj(j, assigmentsForJthVariable.size() + 1);
                Factor facij(VarSet(pi, pj));
                pairFactorList.push_back(facij);
            }
        }
    }
}

// creates all singleton factors for one connected component
// does not set the values
void
ma::PGMMatcher::createSingletonFactors(const McDArray<int>& connectedComp,
                                       vector<Factor>& singletonFactorList) {

    for (int i = 0; i < connectedComp.size(); i++) {
        McDArray<int> assigmentsForIthVariable;
        int curVar = connectedComp[i];
        getAssignmentsForVariable(curVar, assigmentsForIthVariable);

        Var pi(curVar, assigmentsForIthVariable.size() + 1);
        Factor faci(pi);
        singletonFactorList.push_back(faci);
    }
}

void ma::PGMMatcher::handleEvidenceAssignment(const int varLabel,
                                              const int errorCode,
                                              int assignemtnIndexInWholeModel) {
    if (errorCode > -2) {
        return;
    } else if (errorCode == -2) {
        print(QString(
                  "Variable with label %1 was assigned to assignment %2 "
                  "but it these two are not in the same connected component!")
                  .arg(varLabel)
                  .arg(assignemtnIndexInWholeModel));
    } else if (errorCode == -3) {
        print(QString(
                  "Variable with label %1 was assigned to assignment %2 "
                  "but this conflicts with another assignment - pair factors "
                  "not OK!")
                  .arg(varLabel)
                  .arg(assignemtnIndexInWholeModel));
    }

    mQueerEvidence.append(McVec2i(varLabel, assignemtnIndexInWholeModel));
}

// fills the singleton values with a propability computed from the angles the
// adjacent nodes have

void
ma::PGMMatcher::fillSingletonVals(const McDMatrix<int>& variableAssignmentMat,
                                  vector<Factor>& singletonFactors,
                                  ConnectedFactorGraph& graph) {

    while (!singletonFactors.empty()) {

        Factor curFac = singletonFactors.back();
        singletonFactors.pop_back();

        int varLabel = curFac.vars().front().label();
        int assignmentInEvidence;
        int evidenceAssgnmentInWholeModel =
            getEvidenceAssignment(graph, varLabel, assignmentInEvidence);
        handleEvidenceAssignment(varLabel, assignmentInEvidence,
                                 evidenceAssgnmentInWholeModel);
        // Get all possible assignments for variable
        McDArray<int> possibleAssignmentsForVariable;
        getAssignmentsForVariable(varLabel, possibleAssignmentsForVariable);
        if (assignmentInEvidence < 0) {
            McDArray<float> singletonProbs;
            getSingletonProbs(possibleAssignmentsForVariable, varLabel,
                              singletonProbs);

            // set values of factors: Multiply angle and dist threshold
            for (int j = 0; j < curFac.vars().front().states(); j++) {
                curFac.set(j, singletonProbs[j]);
            }
        } else {
            for (int j = 0; j < curFac.vars().front().states(); j++) {
                if (j == assignmentInEvidence) {
                    curFac.set(j, 1);
                } else {
                    curFac.set(j, 0);
                }
            }
        }
        graph.factors.push_back(curFac);
    }
}

void ma::PGMMatcher::getSingletonProbs(const McDArray<int>& assignments,
                                       const int varLabel,
                                       McDArray<double>& probs) {
    McDArray<float> floatProbs;
    getSingletonProbs(assignments, varLabel, floatProbs);

    probs.resize(floatProbs.size());
    for (int i = 0; i < floatProbs.size(); ++i)
        probs[i] = floatProbs[i];
}

void ma::PGMMatcher::normalize(McDArray<float>& probs) {
    double sum = 0.0;
    for (int i = 0; i < probs.size(); i++)
        sum += probs[i];
    if (sum > 0.0) {
        for (int i = 0; i < probs.size(); i++)
            probs[i] /= sum;
    }
}
void ma::PGMMatcher::getSingletonProbs(const McDArray<int>& assignments,
                                       const int varLabel,
                                       McDArray<float>& probs) {

    getWeightsForVariable(assignments, varLabel, probs);
    // append dummy
    probs.append(mWeightComputationFunction.getDummyWeight());
    normalize(probs);
}
// evidenceAssignmentIndexInGraph will be assigned to: >= 0 if assignment was
// found - in that case the value is the index of the assignemnt
//                                          -1 if no assigned evidence was found
//                                          -2 is an assigment was found - but
//                                          the assignment is not in the graph!

int ma::PGMMatcher::getEvidenceAssignment(const ConnectedFactorGraph& graph,
                                          const int varLabel,
                                          int& evidenceAssignmentIndexInGraph) {
    int evidenceAssignmentIndexInWholeModel = -1;
    McVec2i evidence(-1, -1);
    for (int i = 0; i < mEvidence.size(); i++) {
        if (mEvidence[i].x == varLabel) {
            evidence = mEvidence[i];
        }
    }
    // get assignment
    if (evidence.x == varLabel) {
        evidenceAssignmentIndexInWholeModel = evidence.y;
    } else {
        // if there was no assignment, we return -1
        evidenceAssignmentIndexInGraph = -1;
        return evidenceAssignmentIndexInWholeModel;
    }

    // The varaible was assigned, let's see to which one.
    McDArray<int> assignementsForVariable;
    getAssignmentsForVariable(varLabel, assignementsForVariable);
    if (evidenceAssignmentIndexInWholeModel == -1) {

        // if the assigned vaue is the dummy value,
        // assign the dummy value, but return -1 anyway, because there is no
        // label for the assignment

        evidenceAssignmentIndexInGraph = assignementsForVariable.size();
        return evidenceAssignmentIndexInWholeModel;
    }

    for (int i = 0; i < assignementsForVariable.size(); i++) {
        if (assignementsForVariable[i] == evidenceAssignmentIndexInWholeModel) {
            // assign the index of assignment in the current graph but return
            // the index in the full model
            evidenceAssignmentIndexInGraph = i;
            return evidenceAssignmentIndexInWholeModel;
        }
    }
    // not OK
    evidenceAssignmentIndexInGraph = -2;
    return evidenceAssignmentIndexInWholeModel;
}

// This method takes the entries in @allValues at row @label at entries listed
// in @possibleAssignment.
// If @variableAssignmentMat in row @label is 0, meaning that the assignment is
// impossible anyways,
// the entry is replaced by the value @zeroVal.

void
ma::PGMMatcher::getWeightsForVariable(const McDArray<int>& possibleAssignments,
                                      const int label,
                                      McDArray<float>& values) {
    values.resize(0);
    for (int i = 0; i < possibleAssignments.size(); i++) {
        values.append(mWeightComputationFunction.getWeight(
            label, possibleAssignments[i]));
    }
}

// This method takes the entries in @allValues at row @label at entries listed
// in @possibleAssignment.
// If @variableAssignmentMat in row @label is 0, meaning that the assignment is
// impossible anyways,
// the entry is replaced by the value @zeroVal.

void ma::PGMMatcher::getAssignedValuesForVar(
    const McDMatrix<float>& allValues,
    const McDMatrix<int>& variableAssignmentMat,
    const McDArray<int>& possibleAssignments, const int label,
    const float zeroVal, McDArray<float>& values) {
    values.resize(0);
    for (int i = 0; i < possibleAssignments.size(); i++) {
        if (variableAssignmentMat[label][possibleAssignments[i]] > 1.e-6) {
            values.append(allValues[label][possibleAssignments[i]]);
        } else
            values.append(zeroVal);
    }
}

// This method sets the mutual exclusive constraint for each pair factor
void ma::PGMMatcher::fillPairVals(const McDMatrix<int>& variableAssignmentMat,
                                  vector<Factor>& pairFactors,
                                  ConnectedFactorGraph& graph) {

    while (!pairFactors.empty()) {
        Factor curFac = pairFactors.back();
        pairFactors.pop_back();
        int numStatesVar1 = curFac.vars().front().states();
        int numStatesVar2 = curFac.vars().back().states();

        // get the evidence assignment for the two vars, to make sure that the
        // evidence was assigned correctly!
        int var1 = curFac.vars().elements()[0].label();
        int var2 = curFac.vars().elements()[1].label();
        int assignmentForVar1, assignmentForVar2;
        int assignmentIndexInModelForVar1 =
            getEvidenceAssignment(graph, var1, assignmentForVar1);
        int assignmentIndexInModelForVar2 =
            getEvidenceAssignment(graph, var2, assignmentForVar2);

        // create shift matrix, ensuring same shift for vertices are weighted
        // higher
        McDMatrix<float> sameShift(numStatesVar1, numStatesVar2);
        createSameShiftMatrix(graph, curFac, sameShift);
        // we must check, if the hard coded assignemnts in the evidence do not
        // conflict with the pair factors.
        // this can happen, if evidence was given, that would result in 0
        // probability according to parameters
        if (assignmentForVar1 > -1 && assignmentForVar2 > -1) {
            // both were assigned. Now, check if pair entry is 0
            if (sameShift[assignmentForVar1][assignmentForVar2] == 0.0) {
                // add the assignments to the queer evidence
                handleEvidenceAssignment(var1, -3,
                                         assignmentIndexInModelForVar1);
                handleEvidenceAssignment(var2, -3,
                                         assignmentIndexInModelForVar2);
                // set the shift probability to some value >0, the value does
                // not matter
                sameShift[assignmentForVar1][assignmentForVar2] = 1.e-5;
            }
        }

        // create the actual probability matrix

        McDMatrix<float> finalProb = sameShift;  // here woe could multiply
                                                 // other probability factors as
                                                 // well....

        for (int i = 0; i < numStatesVar1 - 1; i++) {
            for (int j = 0; j < numStatesVar2 - 1; j++) {
                int assignmentIndexInArrayForVar1 =
                    mapVariableAssignmentToIndexInVertexList(
                        curFac.vars().front().label(), i);
                int assignmentIndexInArrayForVar2 =
                    mapVariableAssignmentToIndexInVertexList(
                        curFac.vars().back().label(), j);
                if (assignmentIndexInArrayForVar1 ==
                    assignmentIndexInArrayForVar2)
                    finalProb[i][j] = 0.0;
            }
        }

        // set values for factor

        finalProb = finalProb.transpose();

        for (int i = 0; i < numStatesVar1 * numStatesVar2; i++) {
            curFac.set(i, finalProb.dataPtr()[i]);
        }
        graph.factors.push_back(curFac);
    }
}

float ma::PGMMatcher::computeShiftDistance(const McVec3f variableCoord1,
                                           const McVec3f variableCoord2,
                                           const McVec3f assignmentCoord1,
                                           const McVec3f assignmentCoord2) {
    McVec3f dir1 = variableCoord1 - assignmentCoord1;
    McVec3f dir2 = variableCoord2 - assignmentCoord2;

    return (dir1 - dir2).length();
}
void
ma::PGMMatcher::getShiftDistances(const McVec3f pointCoord1,
                                  const McVec3f pointCoord2,
                                  const McDArray<McVec3f>& assignmentCoordsVar1,
                                  const McDArray<McVec3f>& assignmentCoordsVar2,
                                  McDMatrix<float>& distances) {
    distances.resize(assignmentCoordsVar1.size() + 1,
                     assignmentCoordsVar2.size() + 1);

    double shiftProbDummy =
        -1.0 * mShiftProbParam *
        std::log(
            (double)(mWeightComputationFunction.mConfig.dummySignificance));
    for (int i = 0; i < distances.nRows(); i++) {
        for (int j = 0; j < distances.nCols(); j++) {
            if (i == distances.nRows() - 1 || j == distances.nCols() - 1)
                distances[i][j] = shiftProbDummy;
            else {
                distances[i][j] = computeShiftDistance(pointCoord1, pointCoord2,
                                                       assignmentCoordsVar1[i],
                                                       assignmentCoordsVar2[j]);
            }
        }
    }
}

void ma::PGMMatcher::get3dShiftProbs(
    const McVec3f pointCoord1, const McVec3f pointCoord2,
    const McDArray<McVec3f>& assignmentCoordsVar1,
    const McDArray<McVec3f>& assignmentCoordsVar2, McDMatrix<float>& probs) {

    getShiftDistances(pointCoord1, pointCoord2, assignmentCoordsVar1,
                      assignmentCoordsVar2, probs);

    for (int i = 0; i < probs.nRows(); i++) {
        for (int j = 0; j < probs.nCols(); j++) {
            probs[i][j] = mWeightComputationFunction.computeWeightFromRawValue(
                probs[i][j], mShiftProbParam);
        }
    }
    // normalize
    float sumDists = 0.0;
    for (int i = 0; i < probs.nCols() * probs.nRows(); i++)
        sumDists += probs.dataPtr()[i];

    for (int i = 0; i < probs.nCols() * probs.nRows(); i++)
        probs.dataPtr()[i] /= sumDists;
}

void ma::PGMMatcher::createSameShiftMatrix(const ConnectedFactorGraph& graph,
                                           const Factor& curFac,
                                           McDMatrix<float>& shiftEntries) {
    const VarSet& vars = curFac.vars();
    Var var1 = vars.elements()[0];
    Var var2 = vars.elements()[1];
    McVec3f pointCoord1 = mCoords1[var1.label()];
    McVec3f pointCoord2 = mCoords1[var2.label()];
    McVec3f ori1 = mDirections1[var1.label()];
    McVec3f ori2 = mDirections1[var2.label()];
    McDArray<McVec3f> assignmentCoordsVar1;
    McDArray<McVec3f> assignmentCoordsVar2;
    McDArray<McVec3f> assignmentDirections1;
    McDArray<McVec3f> assignmentDirections2;
    for (int i = 0; i < var1.states() - 1; i++) {
        int assignmentCoordIndex =
            mapVariableAssignmentToIndexInVertexList(var1.label(), i);
        assignmentCoordsVar1.append(mCoords2[assignmentCoordIndex]);
        assignmentDirections1.append(mDirections2[assignmentCoordIndex]);
    }
    for (int i = 0; i < var2.states() - 1; i++) {
        int assignmentCoordIndex =
            mapVariableAssignmentToIndexInVertexList(var2.label(), i);
        assignmentCoordsVar2.append(mCoords2[assignmentCoordIndex]);
        assignmentDirections2.append(mDirections2[assignmentCoordIndex]);
    }

    McDMatrix<float> shift3dEntries = shiftEntries;
    get3dShiftProbs(pointCoord1, pointCoord2, assignmentCoordsVar1,
                    assignmentCoordsVar2, shift3dEntries);
    for (int i = 0; i < var1.states(); i++) {
        for (int j = 0; j < var2.states(); j++) {
            shiftEntries[i][j] = shift3dEntries[i][j];
        }
    }
}

bool ma::PGMMatcher::matchPoints(
    ConnectedFactorGraph& graph, const int numIter, const float maxDiff,
    const float damping, McDArray<McVec2i>& matchedPointPairs,
    McDArray<int>& ambiguousAssignments,
    McDArray<int>& ambiguousAssignmentsMarginalDifference,
    McDArray<float>& maxMessageDiffsForNodes,
    McDArray<float>& beliefEntropies) {
    FactorGraph fg(graph.factors);

    mculong maxNumberOfAssignmentsForExactInference = 1.e6;

    print(QString("Graph has %1 factors and %2 variables")
              .arg(graph.factors.size())
              .arg(graph.variables.size()));
    mculong numAssignments = 1;
    McDArray<int> possibleAssignments;
    for (int i = 0;
         (i < graph.variables.size()) &&
             (numAssignments < maxNumberOfAssignmentsForExactInference);
         i++) {
        getAssignmentsForVariable(graph.variables[i], possibleAssignments);
        numAssignments *= ((mculong)(possibleAssignments.size() + 1));
        possibleAssignments.clear();
    }
    bool useExactInf = numAssignments < maxNumberOfAssignmentsForExactInference;
    if (useExactInf) {
        ;
    } else {
        print(QString("PD has more than %1 assignments")
                  .arg(maxNumberOfAssignmentsForExactInference));
    }

    PropertySet opts;
    opts.set("tol", (Real)maxDiff);
    opts.set("maxiter", (size_t)numIter);
    opts.set("maxtime", (Real)300);
    if (useExactInf)
        opts.set("verbose", (size_t)0);
    else
        opts.set("verbose", (size_t)3);
    opts.set("updates", string("SEQMAX"));
    opts.set("logdomain", (bool)true);
    opts.set("inference", string("MAXPROD"));
    opts.set("damping", (Real)damping);

    DAIAlgFG* infAlgorithm;

    if (useExactInf) {
        print("Running exact inference...");
        infAlgorithm = new ExactInf(fg, opts);
    } else {
        print("Running belief propagation...");
        infAlgorithm = new BP(fg, opts);
    }
    infAlgorithm->init();
    try {
        infAlgorithm->run();
    } catch (Exception e) {
        outputSingleFactorValues(graph);
        outputDoubleFactorValues(graph);
        throw e;
    }

    getMaxProbAssignments(infAlgorithm, fg, graph, matchedPointPairs);
    checkAmbiguities(infAlgorithm, fg, graph,
                     ambiguousAssignmentsMarginalDifference);
    checkAmbiguitiesInAssignments(graph, matchedPointPairs,
                                  ambiguousAssignments);
    bool success;

    if (useExactInf)
        success = true;
    else {
        success = (infAlgorithm->maxDiff() < maxDiff);
        if (!success) {
            getMaxMessageDiffs(infAlgorithm, maxMessageDiffsForNodes);
            getBeliefEntropies(infAlgorithm, beliefEntropies);
        }
    }

    delete infAlgorithm;
    return success;
}

void ma::PGMMatcher::getBeliefEntropies(DAIAlgFG* infAlgorithm,
                                        McDArray<float>& beliefEntropies) {

    const std::vector<Factor>& factors = infAlgorithm->factors();
    for (int i = 0; i < factors.size(); i++) {

        Factor curFactor = infAlgorithm->beliefF(i);
        const std::vector<Var>& curVars = curFactor.vars().elements();
        for (int j = 0; j < curVars.size(); j++) {
            Factor curMaxMarginal = curFactor.maxMarginal(curVars[j]);
            float curEntropie = curMaxMarginal.normalized().entropy();
            if (beliefEntropies[curVars[j].label()] < curEntropie)
                beliefEntropies[curVars[j].label()] = curEntropie;
        }
    }
}

void
ma::PGMMatcher::getMaxMessageDiffs(DAIAlgFG* infAlgorithm,
                                   McDArray<float>& maxMessageDiffsForNodes) {
    std::map<int, std::vector<Factor> > maxMarginals;
    std::map<int, double> maxDiffs;

    const std::vector<Factor>& factors = infAlgorithm->factors();
    for (int i = 0; i < factors.size(); i++) {

        Factor curFactor = infAlgorithm->beliefF(i);
        const std::vector<Var>& curVars = curFactor.vars().elements();
        for (int j = 0; j < curVars.size(); j++) {
            Factor curMaxMarginal = curFactor.maxMarginal(curVars[j]);
            std::map<int, std::vector<Factor> >::iterator oldMaxMarginal =
                maxMarginals.find(curVars[j].label());
            std::map<int, std::vector<Factor> >::iterator end =
                maxMarginals.end();
            if (oldMaxMarginal == end) {
                std::vector<Factor> newVector;
                newVector.push_back(curMaxMarginal);
                maxMarginals.insert(pair<int, std::vector<Factor> >(
                    curVars[j].label(), newVector));
            } else {
                std::vector<Factor>& oldMarginals = (*oldMaxMarginal).second;
                for (int k = 0; k < oldMarginals.size(); k++) {

                    double distance =
                        dist(oldMarginals[k].normalized(),
                             curMaxMarginal.normalized(), DISTLINF);

                    if (distance > maxMessageDiffsForNodes[curVars[j].label()])
                        maxMessageDiffsForNodes[curVars[j].label()] = distance;
                }
                oldMarginals.push_back(curMaxMarginal);
            }
        }
    }
}

void ma::PGMMatcher::getMaxProbAssignments(const DAIAlgFG* ia,
                                           const FactorGraph& fg,
                                           const ConnectedFactorGraph& graph,
                                           McDArray<McVec2i>& pairs) {
    try {
        std::vector<std::size_t> maxes = ia->findMaximum();
        McDArray<int> possibleAssignments;

        for (int i = 0; i < maxes.size(); i++) {
            int indexOfAssignmentInVertexList =
                mapVariableAssignmentToIndexInVertexList(graph.variables[i],
                                                         maxes[i]);
            McVec2i pair =
                McVec2i(graph.variables[i], indexOfAssignmentInVertexList);

            pairs.append(pair);
        }

    } catch (Exception e) {
        print(e.getDetailedMsg().c_str());

        for (int i = 0; i < graph.variables.size(); i++) {

            McDArray<int> possibleAssignments;
            getAssignmentsForVariable(graph.variables[i], possibleAssignments);
            Factor belief = ia->belief(
                Var(graph.variables[i], possibleAssignments.size() + 1));
            float maxVal = -1 * FLT_MAX;
            int maxIdx = -1;

            for (int j = 0; j < possibleAssignments.size() + 1; j++) {
                if (belief.get(j) > maxVal) {
                    maxVal = belief.get(j);
                    maxIdx = j;
                }
            }
            int indexOfAssignmentInVertexList =
                mapVariableAssignmentToIndexInVertexList(graph.variables[i],
                                                         maxIdx);

            McVec2i pair =
                McVec2i(graph.variables[i], indexOfAssignmentInVertexList);

            pairs.append(pair);
        }
    }
    // outputSingleFactorValues(graph);

    //    std::vector<std::size_t> maxes= ia.findMaximum();
    //    vector<std::size_t>::iterator it=maxes.begin();
}

void
ma::PGMMatcher::outputSingleFactorValues(const ConnectedFactorGraph& graph) {
    // output factor values
    for (int j = 0; j < graph.factors.size(); j++) {

        Factor fac = graph.factors[j];
        if (fac.vars().size() != 1)
            continue;
        print(QString("singvals for var %1:").arg(fac.vars().front().label()));
        for (int k = 0; k < fac.nrStates(); k++) {
            print(QString("%1").arg(fac.get(k)));
        }
    }
}

void
ma::PGMMatcher::outputDoubleFactorValues(const ConnectedFactorGraph& graph) {
    // output factor values
    for (int j = 0; j < graph.factors.size(); j++) {

        Factor fac = graph.factors[j];
        if (fac.vars().size() != 2)
            continue;
        print(QString("Potentials for vars %1 - %2")
                  .arg(fac.vars().front().label())
                  .arg(fac.vars().back().label()));
        for (int k = 0; k < fac.vars().front().states(); k++) {

            for (int l = 0; l < fac.vars().front().states(); l++) {
                print(QString("%1")
                          .arg(fac.get(k * fac.vars().front().states() + l)));
            }
        }
    }
}

void ma::PGMMatcher::checkAmbiguitiesInAssignments(
    const ConnectedFactorGraph& graph,
    const McDArray<McVec2i>& matchedPointPairs, McDArray<int>& ambiguities) {

    McBitfield assignedAlready(mCoords2.size());
    assignedAlready.unsetAll();

    for (int i = 0; i < matchedPointPairs.size(); i++) {
        if (matchedPointPairs[i].y < 0)
            continue;
        if (assignedAlready[matchedPointPairs[i].y])
            ambiguities.append(matchedPointPairs[i].x);
        assignedAlready.set(matchedPointPairs[i].y);
    }
}

void ma::PGMMatcher::checkAmbiguities(const DAIAlgFG* ia, const FactorGraph& fg,
                                      const ConnectedFactorGraph& graph,
                                      McDArray<int>& ambiguities) {
    for (int h = 0; h < graph.variables.size(); h++) {
        McDArray<int> possibleAssignments;
        getAssignmentsForVariable(graph.variables[h], possibleAssignments);
        Factor belief =
            ia->belief(Var(graph.variables[h], possibleAssignments.size() + 1));

        float maxProb = belief.max();
        int countSame = 0;

        for (int k = 0; k < possibleAssignments.size() + 1; k++) {
            float curProb = belief.get(k);
            if (fabs(curProb - maxProb) < 0.2)
                countSame++;
        }

        if (countSame > 1) {
            // oh no! We found an ambiguous assignment!

            ambiguities.append(graph.variables[h]);

            // print it out:
            print(QString("Found an ambiguous assignment to variable %1")
                      .arg(graph.variables[h]));
            for (int k = 0; k < possibleAssignments.size() + 1; k++) {
                float curProb = belief.get(k);
                print(QString("%1").arg(curProb));
            }
        }
    }
}
