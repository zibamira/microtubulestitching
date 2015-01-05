#include <hxalignmicrotubules/attic/BruteForceOptMatching.h>

#include <mclib/McLine.h>
#include <mclib/McPlane.h>

#include <dai/alldai.h>
#include <dai/bipgraph.h>
#include <dai/bp.h>
#include <dai/daialg.h>
#include <dai/factor.h>
#include <dai/factorgraph.h>

using namespace std;
using namespace dai;

BruteForceOptMatching::BruteForceOptMatching(
    const McDArray<McVec3f>& coords1, const McDArray<McVec3f>& coords2,
    const McDArray<McVec3f>& directions1, const McDArray<McVec3f>& directions2,
    const float distanceThreshold3d, const float distanceThresholdProjected,
    const float angleThreshold, const McDArray<McVec2i>& evidence,
    bool projectPoints)
    : mCoords1(coords1),
      mCoords2(coords2),
      mDirections1(directions1),
      mDirections2(directions2),
      mEvidence(evidence),
      mDistanceThreshold3d(distanceThreshold3d),
      mDistanceThresholdProjected(distanceThresholdProjected),
      mAngleThreshold(angleThreshold) {
    if (projectPoints)
        preparePoints();
}

BruteForceOptMatching::BruteForceOptMatching(
    const McDArray<McVec3f>& coords1, const McDArray<McVec3f>& coords2,
    const McDArray<McVec3f>& directions1, const McDArray<McVec3f>& directions2,
    const float distanceThreshold3d, const float distanceThresholdProjected,
    const float angleThreshold, bool projectPoints)
    : mCoords1(coords1),
      mCoords2(coords2),
      mDirections1(directions1),
      mDirections2(directions2),
      mDistanceThreshold3d(distanceThreshold3d),
      mDistanceThresholdProjected(distanceThresholdProjected),
      mAngleThreshold(angleThreshold) {
    if (projectPoints)
        preparePoints();
}

BruteForceOptMatching::~BruteForceOptMatching(void) {}

void BruteForceOptMatching::preparePoints(void) {
    float medianZ1 = getMedianZ(mCoords1);
    McDArray<McVec3f> oldCoords = mCoords1;
    mCoords1.clear();
    projectToPlaneApproxDirection(oldCoords, mDirections1, medianZ1, mCoords1);

    float medianZ2 = getMedianZ(mCoords2);
    oldCoords = mCoords2;
    mCoords2.clear();
    projectToPlaneApproxDirection(oldCoords, mDirections2, medianZ2, mCoords2);
}

void BruteForceOptMatching::createDistanceThreshold3dMatrix(
    const McDMatrix<float>& distance3dMatrix, McDMatrix<int>& threshMatrix) {
    threshMatrix.resize(mCoords1.size(), mCoords2.size());
    threshMatrix.fill(0);
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {
            float dist = distance3dMatrix(i, j);
            if (dist > mDistanceThreshold3d)
                threshMatrix(i, j) = 0;
            else
                threshMatrix(i, j) = 1;
        }
    }
}

void BruteForceOptMatching::create3dDistanceMatrix(
    McDMatrix<float>& distance3dMatrix) {
    distance3dMatrix.resize(mCoords1.size(), mCoords2.size());
    distance3dMatrix.fill(0);
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {
            if ((mCoords1[i] - mCoords2[j]).length() >
                2.0 * mDistanceThreshold3d) {
                distance3dMatrix(i, j) = (mCoords1[i] - mCoords2[j]).length();
            } else {
                // remove gap component
                // first: gap distance
                float z = fabs(mCoords1[i].z - mCoords2[j].z);

                // second: Figure out angle and compute how long the gap
                // bridging would be
                McVec3f dir1 = mDirections1[i];
                McVec3f dir2 = mDirections2[j];
                // compute the mean angle
                dir1 *= -1;
                dir1.normalize();
                dir2.normalize();
                McVec3f dir = dir1 + dir2;
                float angle = McVec3f(0, 0, 1).angle(dir);
                if (angle > M_PI / 2.0)
                    angle = M_PI - angle;
                angle = M_PI / 2.0 - angle;
                float projDist = z / sin(angle);

                distance3dMatrix(i, j) =
                    (mCoords1[i] - mCoords2[j]).length() - projDist;
                if (distance3dMatrix(i, j) < 0)
                    distance3dMatrix(i, j) = 0;
            }
        }
    }
}

void BruteForceOptMatching::createDistanceThresholdProjectedMatrix(
    McDMatrix<int>& threshMatrix,
    const McDMatrix<float>& projectionDistMatrix) {

    threshMatrix.resize(mCoords1.size(), mCoords2.size());
    threshMatrix.fill(0);
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {

            if (projectionDistMatrix(i, j) > mDistanceThresholdProjected)
                threshMatrix(i, j) = 0;
            else
                threshMatrix(i, j) = 1;
        }
    }
}

void BruteForceOptMatching::createAngleThresholdMatrix(
    McDMatrix<int>& threshMatrix, const McDMatrix<float>& angleMatrix) {

    threshMatrix.resize(angleMatrix.nRows(), angleMatrix.nCols());
    threshMatrix.fill(0);
    for (int i = 0; i < angleMatrix.nRows(); i++) {
        for (int j = 0; j < angleMatrix.nCols(); j++) {

            if (angleMatrix(i, j) > mAngleThreshold)
                threshMatrix(i, j) = 1;
            else
                threshMatrix(i, j) = 0;
        }
    }
}

void BruteForceOptMatching::createProjectedDistanceMatrix(
    McDMatrix<float>& projMatrix) {
    projMatrix.resize(mCoords1.size(), mCoords2.size());
    projMatrix.fill(0);
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {
            projMatrix(i, j) = getMinNormalProjectedDistance(
                mCoords1[i], mCoords2[j], mDirections1[i], mDirections2[j]);
        }
    }
}

void BruteForceOptMatching::createAngleMatrix(McDMatrix<float>& angleMatrix) {
    angleMatrix.resize(mCoords1.size(), mCoords2.size());
    angleMatrix.fill(0);
    for (int i = 0; i < mCoords1.size(); i++) {
        for (int j = 0; j < mCoords2.size(); j++) {
            float angle = mDirections1[i].angle(mDirections2[j]);
            // angle = acos(angle);
            // angle = angle/(2*M_PI)*360.0;
            angleMatrix(i, j) = angle / M_PI * 180;
        }
    }
}

float BruteForceOptMatching::getMinNormalProjectedDistance(const McVec3f coord1,
                                                           const McVec3f coord2,
                                                           const McVec3f dir1,
                                                           const McVec3f dir2) {
    // project coord1 on plane perpendicular to dir2 positioned at coord2

    // check, if the two lines overlapp:
    // if  dir1 points into positive z, then origPoint2 must be below
    // origPoint1,

    float todoParam = 200;
    if (dir1.z >= 0) {
        if (coord1.z + todoParam < coord2.z)
            return FLT_MAX;
    } else {
        if (coord1.z - todoParam > coord2.z)
            return FLT_MAX;
    }

    float projectedDist12 =
        getNormalProjectedDistance(coord1, coord2, dir1, dir2);
    float projectedDist21 =
        getNormalProjectedDistance(coord2, coord1, dir2, dir1);
    if (projectedDist12 > projectedDist21)
        return projectedDist21;
    else
        return projectedDist12;
}

float BruteForceOptMatching::getNormalProjectedDistance(const McVec3f coord1,
                                                        const McVec3f coord2,
                                                        const McVec3f dir1,
                                                        const McVec3f dir2) {
    McPlane plane(dir2, coord2);
    McLine line(coord1, coord1 + dir1);
    McVec3f intersectionPoint;
    if (plane.intersect(line, intersectionPoint)) {
        return (coord2 - intersectionPoint).length();
    } else {
        return FLT_MAX;
    }
}

// Finds all connected components in adjacenceMatrix, that is, all variable
// labels that are connected in the pgm are after processinf an array in the
// connComps array.

void BruteForceOptMatching::getAllConnectedComponents(
    vector<ConnectedFactorGraph>& graphs) {

    mQueerEvidence.resize(0);
    McDMatrix<int> treshMat3d;
    McDMatrix<float> distMat3d;
    create3dDistanceMatrix(distMat3d);
    createDistanceThreshold3dMatrix(distMat3d, treshMat3d);
    McDMatrix<float> projDistMatrix;
    createProjectedDistanceMatrix(projDistMatrix);
    McDMatrix<int> projDistThreshMat;
    createDistanceThresholdProjectedMatrix(projDistThreshMat, projDistMatrix);

    McDMatrix<float> angleMatrix;
    createAngleMatrix(angleMatrix);

    McDMatrix<int> angleThreshMatrix;
    createAngleThresholdMatrix(angleThreshMatrix, angleMatrix);
    mVariableAssignmentMat = projDistThreshMat.multElementWise(treshMat3d);
    mVariableAssignmentMat =
        mVariableAssignmentMat.multElementWise(angleThreshMatrix);

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
                                 variableAdjacenceMatrix, angleMatrix,
                                 projDistMatrix, distMat3d, graph);
        graphs.push_back(graph);

        variableAdjacenceMatrix = adjacenceMatrixWithoutConnctedComponent;
    }
}

// Computes all variables that form a connected component in the
// adjacenceMatrix.
// The connected component chosen is arbitrary - it takes the first it finds.
bool BruteForceOptMatching::getConnectedComponent(
    const McDMatrix<int>& adjacenceMatrix,
    McDMatrix<int>& adjacenceMatrixWithoutConnctedComponent,
    McDArray<int>& connComp) {
    adjacenceMatrixWithoutConnctedComponent.resize(adjacenceMatrix.nRows(),
                                                   adjacenceMatrix.nCols());
    memcpy(adjacenceMatrixWithoutConnctedComponent.dataPtr(),
           adjacenceMatrix.dataPtr(),
           sizeof(int) * adjacenceMatrix.nRows() * adjacenceMatrix.nCols());

    // find first a startpoint
    int start = -1;
    connComp.resize(0);
    for (int i = 0; i < adjacenceMatrix.nRows(); i++) {
        for (int j = i; j < adjacenceMatrix.nCols(); j++) {
            if (adjacenceMatrix[i][j] == 1) {
                start = i;
                break;
            }
        }
    }
    if (start == -1)
        return false;

    McDArray<int> queue;
    queue.append(start);
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

void BruteForceOptMatching::createAdjacenceMatrix(
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

void BruteForceOptMatching::createFactorGraphFactors(
    const McDMatrix<int>& variableAssignmentMat,
    const McDMatrix<int>& variableAdjacenceMatrix,
    const McDMatrix<float>& angleMatrix,
    const McDMatrix<float>& projDistanceMatrix,
    const McDMatrix<float>& distance3dMatrix,
    ConnectedFactorGraph& factorGraph) {

    vector<Factor> pairFactors;
    createConnectionFactors(variableAdjacenceMatrix, factorGraph.variables,
                            pairFactors);
    vector<Factor> singletonFactors;
    createSingletonFactors(factorGraph.variables, singletonFactors);

    fillSingletonVals(angleMatrix, projDistanceMatrix, distance3dMatrix,
                      variableAssignmentMat, singletonFactors, factorGraph);
    fillPairVals(variableAssignmentMat, projDistanceMatrix, pairFactors,
                 factorGraph);
}

int BruteForceOptMatching::mapVariableAssignmentToIndexInVertexList(
    const int variableLabel, const int assignment) {
    McDArray<int> possibleAssignments;
    getAssignmentsForVariable(variableLabel, possibleAssignments);
    // see if the assignment was the dummy variable
    if (assignment == possibleAssignments.size())
        return -1;

    return possibleAssignments[assignment];
}

int BruteForceOptMatching::mapIndexInVertexListToVariableAssignment(
    const int variableLabel, const int indexInCoords2Array) {
    McDArray<int> possibleAssignments;
    getAssignmentsForVariable(variableLabel, possibleAssignments);
    return possibleAssignments.findSorted(indexInCoords2Array,
                                          &mcStandardCompare);
}

void BruteForceOptMatching::getAssignmentsForVariable(
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
void BruteForceOptMatching::createConnectionFactors(
    const McDMatrix<int>& adjMat, const McDArray<int> connectedComp,
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
void BruteForceOptMatching::createSingletonFactors(
    const McDArray<int>& connectedComp, vector<Factor>& singletonFactorList) {

    for (int i = 0; i < connectedComp.size(); i++) {
        McDArray<int> assigmentsForIthVariable;
        int curVar = connectedComp[i];
        getAssignmentsForVariable(curVar, assigmentsForIthVariable);

        Var pi(curVar, assigmentsForIthVariable.size() + 1);
        Factor faci(pi);
        singletonFactorList.push_back(faci);
    }
}

void BruteForceOptMatching::handleEvidenceAssignment(
    const int varLabel, const int errorCode, int assignemtnIndexInWholeModel) {
    if (errorCode > -2)
        return;

    if (errorCode == -2)
        cout << "\nVariable with label " << varLabel
             << " was assigned to assignment " << assignemtnIndexInWholeModel
             << " but it these two are not in the same connected component!";
    if (errorCode == -3)
        cout << "\nVariable with label " << varLabel
             << " was assigned to assignment " << assignemtnIndexInWholeModel
             << " but this conflicts with another assignment - pair factors "
                "not OK!";

    mQueerEvidence.append(McVec2i(varLabel, assignemtnIndexInWholeModel));
}

// fills the singleton values with a propability computed from the angles the
// adjacent nodes have

void BruteForceOptMatching::fillSingletonVals(
    const McDMatrix<float>& angleMatrix,
    const McDMatrix<float>& projDistanceMatrix,
    const McDMatrix<float>& distanceMatrix3d,
    const McDMatrix<int>& variableAssignmentMat,
    vector<Factor>& singletonFactors, ConnectedFactorGraph& graph) {

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
            getSingletonProbs(angleMatrix, projDistanceMatrix, distanceMatrix3d,
                              variableAssignmentMat,
                              possibleAssignmentsForVariable, varLabel,
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

void BruteForceOptMatching::getSingletonProbs(
    const McDMatrix<float>& angleMatrix,
    const McDMatrix<float>& projDistanceMatrix,
    const McDMatrix<float>& distanceMatrix3d,
    const McDMatrix<int>& variableAssignmentMat,
    const McDArray<int>& assignments, const int varLabel,
    McDArray<double>& probs) {
    McDArray<float> floatProbs;
    getSingletonProbs(angleMatrix, projDistanceMatrix, distanceMatrix3d,
                      variableAssignmentMat, assignments, varLabel, floatProbs);

    probs.resize(floatProbs.size());
    for (int i = 0; i < floatProbs.size(); ++i)
        probs[i] = floatProbs[i];
}

void BruteForceOptMatching::getSingletonProbs(
    const McDMatrix<float>& angleMatrix,
    const McDMatrix<float>& projDistanceMatrix,
    const McDMatrix<float>& distanceMatrix3d,
    const McDMatrix<int>& variableAssignmentMat,
    const McDArray<int>& assignments, const int varLabel,
    McDArray<float>& probs) {

    McDArray<float> angleValues;
    getAssignedValuesForVar(angleMatrix, variableAssignmentMat, assignments,
                            varLabel, 0, angleValues);
    mcassert(angleValues.size() == assignments.size());
    // add dummy
    angleValues.append(mAngleThreshold / 2.0);
    // compute actual prob representation
    computeAngleProbs(angleValues);

    McDArray<float> projDistValues;
    getAssignedValuesForVar(projDistanceMatrix, variableAssignmentMat,
                            assignments, varLabel, FLT_MAX, projDistValues);
    // add dummy
    projDistValues.append(mDistanceThresholdProjected / 2.0);
    // compute actual prob representation
    computeProjDistProbs(projDistValues);

    McDArray<float> distValues3d;
    getAssignedValuesForVar(distanceMatrix3d, variableAssignmentMat,
                            assignments, varLabel, FLT_MAX, distValues3d);
    // add dummy
    distValues3d.append(mDistanceThreshold3d / 2.0);
    // compute actual prob representation
    compute3dDistProbs(distValues3d);

    probs.resize(assignments.size() + 1);
    // set values of factors: Multiply angle and dist threshold
    for (int j = 0; j < probs.size(); j++) {
        probs[j] = projDistValues[j] * angleValues[j];
    }
}
// evidenceAssignmentIndexInGraph will be assigned to: >= 0 if assignment was
// found - in that case the value is the index of the assignemnt
//                                          -1 if no assigned evidence was found
//                                          -2 is an assigment was found - but
//                                          the assignment is not in the graph!

int BruteForceOptMatching::getEvidenceAssignment(
    const ConnectedFactorGraph& graph, const int varLabel,
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

// Normalizes an array representing distances

void BruteForceOptMatching::compute3dDistProbs(McDArray<float>& dists) {
    float sumDists = 0.0;
    // This should be a user parameter...
    float probParam = 500.0;
    for (int i = 0; i < dists.size(); i++) {
        if (dists[i] != FLT_MAX) {
            dists[i] = exp(-dists[i] / probParam);
        } else
            dists[i] = 0;
        sumDists += dists[i];
    }
    for (int i = 0; i < dists.size(); i++)
        dists[i] /= sumDists;
}

// Normalizes an array representing distances

void BruteForceOptMatching::computeProjDistProbs(McDArray<float>& dists) {
    float sumDists = 0.0;
    for (int i = 0; i < dists.size(); i++) {
        if (dists[i] != FLT_MAX) {
            dists[i] = exp((float)(-dists[i] / 100.0));
        } else
            dists[i] = 0;
        sumDists += dists[i];
    }
    for (int i = 0; i < dists.size(); i++)
        dists[i] /= sumDists;
}

// This method normalizes the array @angles. This assumes that a high angle
// means high probability.

void BruteForceOptMatching::computeAngleProbs(McDArray<float>& angles) {
    float sumAngles = 0;
    for (int i = 0; i < angles.size(); i++) {
        angles[i] = 180.0 - angles[i];
        angles[i] = exp((float)(-angles[i] / 10.0));
        sumAngles += angles[i];
    }
    for (int i = 0; i < angles.size(); i++)
        angles[i] /= sumAngles;
}
// This method takes the entries in @allValues at row @label at entries listed
// in @possibleAssignment.
// If @variableAssignmentMat in row @label is 0, meaning that the assignment is
// impossible anyways,
// the entry is replaced by the value @zeroVal.

void BruteForceOptMatching::getAssignedValuesForVar(
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
void
BruteForceOptMatching::fillPairVals(const McDMatrix<int>& variableAssignmentMat,
                                    const McDMatrix<float>& projDistanceMatrix,
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
        createSameShiftMatrix(graph, curFac, projDistanceMatrix, sameShift);
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
                    finalProb[i][j] = 0;
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

void BruteForceOptMatching::get3dShiftProbs(
    const McVec3f pointCoord1, const McVec3f pointCoord2,
    const McDArray<McVec3f>& assignmentCoordsVar1,
    const McDArray<McVec3f>& assignmentCoordsVar2, McDMatrix<float>& probs) {

    float probParam = 50.0;
    float dummyDistance = mDistanceThresholdProjected / 2.0;

    for (int i = 0; i < assignmentCoordsVar1.size(); i++) {
        for (int j = 0; j < assignmentCoordsVar2.size(); j++) {
            McVec3f dir1 = pointCoord1 - assignmentCoordsVar1[i];
            McVec3f dir2 = pointCoord2 - assignmentCoordsVar2[j];

            float dist1 = dir1.length();
            float dist2 = dir2.length();

            dir1.normalize();
            dir2.normalize();
            dir1 *= -1;
            dir2 *= -1;
            McVec3f p2Projected = dir1 * dist2 + pointCoord2;
            McVec3f p1Projected = dir2 * dist1 + pointCoord1;

            float distP1Projected =
                (assignmentCoordsVar1[i] - p1Projected).length();
            float distP2Projected =
                (assignmentCoordsVar2[j] - p2Projected).length();
            float maxDist = distP1Projected > distP2Projected ? distP1Projected
                                                              : distP2Projected;
            probs[i][j] = exp(-maxDist / probParam);
        }
    }
    for (int i = 0; i <= assignmentCoordsVar1.size(); i++) {
        probs[i][assignmentCoordsVar2.size()] = exp(-dummyDistance / probParam);
    }
    for (int i = 0; i <= assignmentCoordsVar2.size(); i++) {
        probs[(int)(assignmentCoordsVar1.size())][i] =
            exp(-dummyDistance / probParam);
    }

    // normalize
    float sumDists = 0.0;
    for (int i = 0; i < probs.nCols() * probs.nRows(); i++)
        sumDists += probs.dataPtr()[i];

    for (int i = 0; i < probs.nCols() * probs.nRows(); i++)
        probs.dataPtr()[i] /= sumDists;
}

void BruteForceOptMatching::createSameShiftMatrix(
    const ConnectedFactorGraph& graph, const Factor& curFac,
    const McDMatrix<float>& projecedDistMatrix,
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

bool BruteForceOptMatching::matchPoints(ConnectedFactorGraph& graph,
                                        const int numIter, const float maxDiff,
                                        const float damping,
                                        McDArray<McVec2i>& matchedPointPairs,
                                        McDArray<int>& ambiguousAssignments) {
    FactorGraph fg(graph.factors);

    cout << "\n run BP for network with " << graph.factors.size()
         << " factors and " << graph.factors[0].vars().front().states()
         << "states";
    PropertySet opts;
    opts.set("tol", (Real)maxDiff);
    opts.set("maxiter", (size_t)numIter);
    opts.set("maxtime", (Real)180);
    opts.set("verbose", (size_t)1);
    opts.set("updates", string("SEQMAX"));
    opts.set("logdomain", (bool)true);
    opts.set("inference", string("MAXPROD"));
    opts.set("damping", (Real)damping);
    BP ia(fg, opts);
    ia.init();
    try {
        ia.run();
    } catch (Exception e) {
        outputSingleFactorValues(graph);
        outputDoubleFactorValues(graph);
        throw e;
    }

    getMaxProbAssignments(ia, fg, graph, matchedPointPairs);
    checkAmbiguities(ia, fg, graph, ambiguousAssignments);
    checkAmbiguitiesInAssignments(graph, matchedPointPairs,
                                  ambiguousAssignments);
    return (ia.maxDiff() < maxDiff);
}

void BruteForceOptMatching::getMaxProbAssignments(
    const BP& ia, const FactorGraph& fg, const ConnectedFactorGraph& graph,
    McDArray<McVec2i>& pairs) {

    for (int i = 0; i < graph.variables.size(); i++) {

        McDArray<int> possibleAssignments;
        getAssignmentsForVariable(graph.variables[i], possibleAssignments);
        Factor belief =
            ia.belief(Var(graph.variables[i], possibleAssignments.size() + 1));
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
    outputSingleFactorValues(graph);

    //    std::vector<std::size_t> maxes= ia.findMaximum();
    //    vector<std::size_t>::iterator it=maxes.begin();
}

void BruteForceOptMatching::outputSingleFactorValues(
    const ConnectedFactorGraph& graph) {
    // output factor values
    for (int j = 0; j < graph.factors.size(); j++) {

        Factor fac = graph.factors[j];
        if (fac.vars().size() != 1)
            continue;
        cout << "singvals for var " << fac.vars().front().label() << " :\n";
        for (int k = 0; k < fac.nrStates(); k++) {
            cout << fac.get(k) << " ";
        }
        cout << "\n";
    }
}

void BruteForceOptMatching::outputDoubleFactorValues(
    const ConnectedFactorGraph& graph) {
    // output factor values
    for (int j = 0; j < graph.factors.size(); j++) {

        Factor fac = graph.factors[j];
        if (fac.vars().size() != 2)
            continue;
        cout << "\nPotentials for vars " << fac.vars().front().label() << " - "
             << fac.vars().back().label() << "\n";
        for (int k = 0; k < fac.vars().front().states(); k++) {

            for (int l = 0; l < fac.vars().front().states(); l++) {
                cout << fac.get(k * fac.vars().front().states() + l) << " ";
            }
            cout << "\n";
        }
    }
}

void BruteForceOptMatching::checkAmbiguitiesInAssignments(
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

void BruteForceOptMatching::checkAmbiguities(const BP& ia,
                                             const FactorGraph& fg,
                                             const ConnectedFactorGraph& graph,
                                             McDArray<int>& ambiguities) {
    for (int h = 0; h < graph.variables.size(); h++) {
        McDArray<int> possibleAssignments;
        getAssignmentsForVariable(graph.variables[h], possibleAssignments);
        Factor belief =
            ia.belief(Var(graph.variables[h], possibleAssignments.size() + 1));

        float maxProb = belief.max();
        int countSame = 0;

        for (int k = 0; k < possibleAssignments.size() + 1; k++) {
            float curProb = belief.get(k);
            if (fabs(curProb - maxProb) < 0.1)
                countSame++;
        }

        /////
        cout << "\n Belief for var " << graph.variables[h] << "\n";
        for (int k = 0; k < possibleAssignments.size() + 1; k++) {
            float curProb = belief.get(k);
            cout << curProb << " ";
        }
        cout << "\n";

        ////

        if (countSame > 1) {
            // oh no! We found an ambiguos assignment!

            ambiguities.append(graph.variables[h]);

            // print it out:
            cout << "Found an ambiguous assignemnt to variable "
                 << graph.variables[h] << "\n";
            for (int k = 0; k < possibleAssignments.size() + 1; k++) {
                float curProb = belief.get(k);
                cout << curProb << " ";
            }
            cout << "\n";
        }
    }
}

void BruteForceOptMatching::projectToPlaneApproxDirection(
    const McDArray<McVec3f>& vertices, const McDArray<McVec3f>& directions,
    const float planeZ, McDArray<McVec3f>& result) {

    for (int i = 0; i < vertices.size(); i++) {
        McPlane theZPlane(McVec3f(0, 0, 1), planeZ);
        McVec3f vertexCoord = vertices[i];
        McVec3f dir = directions[i] * -1;
        dir.normalize();
        float angle = dir.angle(McVec3f(0, 0, 1));
        if (angle > M_PI / 2.0)
            angle = M_PI - angle;

        McLine theLine(vertexCoord, vertexCoord + directions[i]);

        McVec3f intersectionPoint;
        bool intersected = theZPlane.intersect(theLine, intersectionPoint);
        // if(fabs(angle)<0.1)
        //   cout<<"\n Angle for vertex "<<i <<" too low: "<< angle;
        if (intersected && (fabs(angle) < (M_PI / 2.0 - M_PI / 8.0))) {
            result.append(intersectionPoint);
        } else {
            result.append(McVec3f(vertexCoord.x, vertexCoord.y, planeZ));
        }
    }
}

float BruteForceOptMatching::getMedianZ(const McDArray<McVec3f>& vertices) {

    if (!vertices.size())
        return -1 * FLT_MAX;
    McDArray<float> zs;
    float mean = 0.0;

    for (int i = 0; i < vertices.size(); i++) {
        zs.append(vertices[i].z);
        mean += vertices[i].z;
    }
    zs.sort(&mcStandardCompare);
    int medianIdx = zs.size() / 2.0;
    // cout <<"MeanZ: "<<mean/vertices.size();
    // return mean/vertices.size();
    return zs[medianIdx];
}
