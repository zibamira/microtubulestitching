#include <hxalignmicrotubules/mtalign/matchingPGM.h>

#include <QString>

#include <mclib/McVec2i.h>
#include <mclib/McBitfield.h>

#include <dai/exceptions.h>
#include <hxalignmicrotubules/mtalign/data.h>
#include <hxalignmicrotubules/mtalign/PGMMatcher.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

namespace ma = mtalign;

static bool graphIsEmpty(const McDMatrix<int>& mat) {
    for (int i = 0; i < mat.nRows(); i++) {
        for (int j = 0; j < mat.nCols(); j++) {
            if (mat[i][j] > 0)
                return false;
        }
    }
    return true;
}

static void pruneGraph(McDMatrix<int>& mat) {
    for (int i = 0; i < mat.nRows(); i++) {
        int numEdges = 0;
        for (int j = 0; j < mat.nCols(); j++) {
            numEdges += mat[i][j];
        }
        if (numEdges == 2 || numEdges == 1) {
            // We found a branch!  Prune it.
            for (int j = 0; j < mat.nCols(); j++) {
                mat[i][j] = 0;
                mat[j][i] = 0;
            }
            i = 0;
        }
    }
}

static int pruneMostConnectedVar(McDMatrix<int>& mat) {
    int maxConn = -1;
    int maxIdx = -1;
    for (int i = 0; i < mat.nRows(); i++) {
        int numConn = 0;
        for (int j = 0; j < mat.nCols(); j++) {
            numConn += mat[i][j];
        }
        if (numConn > maxConn) {
            maxConn = numConn;
            maxIdx = i;
        }
    }
    for (int i = 0; i < mat.nRows(); i++) {
        mat[i][maxIdx] = 0;
        mat[maxIdx][i] = 0;
    }
    return maxIdx;
}

static int getCriticalVariable(McDArray<float>& maxMessageDiffsForNodes,
                               const McDArray<int>& vars,
                               const McBitfield& isEvidence) {
    int maximumIndex = -1;
    float maxDist = 0.0;
    for (int i = 0; i < vars.size(); i++) {
        const double distance = maxMessageDiffsForNodes[vars[i]];
        if (distance > maxDist && !isEvidence[vars[i]]) {
            maxDist = distance;
            maximumIndex = vars[i];
        }
    }
    mcassert(maximumIndex > -1);
    return maximumIndex;
}

static McDArray<int>
getPointsToRemoveToMakeATree(const mtalign::ConnectedFactorGraph& graph,
                             const McDMatrix<int>& variableAssignmentMat,
                             const McBitfield& isEvidence) {

    McDMatrix<int> prunedGraph = graph.adjacencyMatrix;
    pruneGraph(prunedGraph);

    McDArray<int> assignToMakeTree;
    while (!graphIsEmpty(prunedGraph)) {
        const int mostConnectedVar = pruneMostConnectedVar(prunedGraph);
        if (!isEvidence[mostConnectedVar]) {
            assignToMakeTree.append(mostConnectedVar);
        }
        pruneGraph(prunedGraph);
    }

    return assignToMakeTree;
}

static void findPointsToAssignWithAdaptiveThreshold(
    const mtalign::ConnectedFactorGraph& graph, McDArray<float>& maxDiff,
    McDArray<int>& resultNodes) {
    for (int i = 0; i < graph.variables.size(); i++) {
        int curVar = graph.variables[i];
        float curVal = maxDiff[curVar];
        if (curVal < 1.e-3)
            continue;
        bool assign = true;
        // find the neighbours for the current node
        for (int j = 0; j < graph.adjacencyMatrix.nCols(); j++) {
            if (graph.adjacencyMatrix[curVar][j] ||
                graph.adjacencyMatrix[j][curVar]) {
                if (maxDiff[j] > curVal)
                    assign = false;
            }
        }
        if (assign)
            resultNodes.append(curVar);
    }
}

ma::MatchingPGM ma::matchingPGM(const ma::FacingPointSets& pts,
                                const MatchingPGMParams& params,
                                ma::Context* ctx) {
    if (!ctx) {
        ctx = &defaultContext();
    }
    ma::MatchingPGM res;
    const mclong nRefPoints = pts.ref.positions.size();

    McBitfield isEvidence(nRefPoints, McBitfield::INIT_UNSET);
    for (int i = 0; i < params.evidence.size(); i++) {
        isEvidence.set(params.evidence[i][0]);
    }

    ma::PGMPairWeights weighter(pts.ref.positions, pts.trans.positions,
                                pts.ref.directions, pts.trans.directions,
                                params.weightConfig);
    mtalign::PGMMatcher bfom(pts.ref.positions, pts.trans.positions,
                             pts.ref.directions, pts.trans.directions,
                             params.evidence, weighter, params.pairFactorParam);
    bfom.setContext(ctx);

    McDArray<int> notConverged;
    McDArray<int> ambiguousAssignments;
    McDArray<int> criticalVariables;
    McDArray<int> groupsWithAmbiguities;
    McDArray<int> pointsToRemove;
    McDArray<int> errorGroups;
    McDArray<int> nodesToAssignWithAdaptiveThreshold;
    McDArray<int> ambiguousAssignmentsMarginalDifference;
    McDArray<int> maxEntropieVariables;
    McDArray<dai::Exception> exceptions;

    McDArray<float> maxMessageDiffOfLastIteration(nRefPoints);
    maxMessageDiffOfLastIteration.fill(0);

    McDArray<float> beliefEntropies(nRefPoints);
    beliefEntropies.fill(0);

    int groupCounter = 0;
    bfom.initConnectedComponentIteration();
    mtalign::ConnectedFactorGraph graph;
    int numVarsWithNeighbour = 0;
    while (bfom.getNextConnectedComponent(graph)) {
        if (graph.variables.size() > 1)
            numVarsWithNeighbour += graph.variables.size();

        dai::FactorGraph fg(graph.factors);
        if (graph.factors.size() > 0) {
            McDArray<McVec2i> matchedPointPairs;
            int ambSize = ambiguousAssignments.size();
            bool converged = false;

            try {
                converged = bfom.matchPoints(
                    graph, 5, 1.e-7, 0.0, matchedPointPairs,
                    ambiguousAssignments,
                    ambiguousAssignmentsMarginalDifference,
                    maxMessageDiffOfLastIteration, beliefEntropies);
            } catch (dai::Exception e) {
                ctx->print("InfAlgorithm threw an exception!");
                ctx->print(e.getMsg().c_str());
                exceptions.append(e);
                converged = false;
                errorGroups.append(groupCounter);
            }

            if (!converged) {
                notConverged.append(groupCounter);
                pointsToRemove.appendArray(getPointsToRemoveToMakeATree(
                    graph, bfom.getVariableAssignmentMat(), isEvidence));
                criticalVariables.append(
                    getCriticalVariable(maxMessageDiffOfLastIteration,
                                        graph.variables, isEvidence));
                maxEntropieVariables.append(getCriticalVariable(
                    beliefEntropies, graph.variables, isEvidence));
                findPointsToAssignWithAdaptiveThreshold(
                    graph, maxMessageDiffOfLastIteration,
                    nodesToAssignWithAdaptiveThreshold);
            }

            if (ambSize != ambiguousAssignments.size())
                groupsWithAmbiguities.append(groupCounter);

            for (int h = 0; h < matchedPointPairs.size(); h++) {
                const McVec2i p = matchedPointPairs[h];
                if (p.y >= 0) {
                    res.matchedRefPointIds.append(p[0]);
                    res.matchedTransPointIds.append(p[1]);
                }
            }
            groupCounter++;
        }
        graph.adjacencyMatrix.fill(0);
        graph.factors.clear();
        graph.variables.clear();
    }

    for (int h = 0; h < bfom.mQueerEvidence.size(); h++) {
        const McVec2i e = bfom.mQueerEvidence[h];
        if (e.y >= 0) {
            res.queerEvidence.append(e);
        }
    }
    res.ambiguities = ambiguousAssignments;
    res.criticalNodes = criticalVariables;
    res.assignToMakeTree = pointsToRemove;
    res.criticalNodesWithAdaptiveThreshold = nodesToAssignWithAdaptiveThreshold;
    res.ambiguitiesFromMarginalDiff = ambiguousAssignmentsMarginalDifference;
    res.maxEntropieNodes = maxEntropieVariables;
    res.maxDiff = maxMessageDiffOfLastIteration;
    res.entropie = beliefEntropies;

    ctx->print("");
    ctx->print(
        QString("%1 networks did not converge.").arg(notConverged.size()));
    for (int i = 0; i < notConverged.size(); i++) {
        ctx->print(QString("Group %1").arg(notConverged[i]));
    }
    for (int i = 0; i < exceptions.size(); i++) {
        ctx->print(
            QString("Exception occured in group %1").arg(errorGroups[i]));
        ctx->print(exceptions[i].getMsg().c_str());
    }
    ctx->print(
        QString("%1 variables had a neighbour.").arg(numVarsWithNeighbour));
    ctx->print(
        QString("%1 ambiguous assignments").arg(ambiguousAssignments.size()));
    for (int i = 0; i < groupsWithAmbiguities.size(); i++) {
        ctx->print(QString("In group %1").arg(groupsWithAmbiguities[i]));
    }
    if (res.queerEvidence.size() > 0) {
        ctx->print(
            QString("WARNING! Found %1 evidence assignment that did not match "
                    "the parameters.").arg(res.queerEvidence.size()));
    }

    return res;
}

