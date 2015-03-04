#include <hxalignmicrotubules/mtalign/PGMMatcher.h>

#include <fstream>

#include <gtest/gtest.h>
#include <hxcore/HxObjectPool.h>
#include <hxcore/TestingData.h>
#include <hxcore/HxResource.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/TestingDevNullRedirect.h>

#include <dai/bipgraph.h>
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/alldai.h>
#include <dai/bp.h>
#include <dai/factor.h>

using namespace std;
using namespace dai;

// get a point on a sg edge that is maxDist from first point in edgePoints list
static McVec3f getNextPointAtDist(const McDArray<McVec3f> edgePoints,
                                  const float maxDist) {
    float dist = 0;
    int actIndex = 1;

    while ((dist < maxDist) && (actIndex < edgePoints.size())) {
        dist += (edgePoints[actIndex - 1] - edgePoints[actIndex]).length();
        actIndex++;
    }
    // if the max Dist was not reached at all
    if (dist <= maxDist) {
        return edgePoints[actIndex - 1];
    }
    // we interpolate the last point
    else {
        // get difference
        float diffDist = dist - maxDist;
        McVec3f x2 = edgePoints[actIndex - 1];
        McVec3f x1 = edgePoints[actIndex - 2];
        float distLastPoints = (x2 - x1).length();
        float frac = diffDist / distLastPoints;
        mcassert(frac < 1);

        return (1 - frac) * x2 + frac * x1;
    }
}

// computed the direction vectors for the edges of the selected vertices
static void getDirections(const SpatialGraphSelection& selectedVertices,
                          HxSpatialGraph* graph,
                          McDArray<McVec3f>& directions) {

    directions.resize(selectedVertices.getNumSelectedVertices());

    for (int i = 0; i < selectedVertices.getNumSelectedVertices(); i++) {
        int selectedVertex = selectedVertices.getSelectedVertex(i);
        IncidenceList connectedEdges = graph->getIncidentEdges(selectedVertex);

        McVec3f point1 = graph->getVertexCoords(selectedVertex);
        McVec3f point2;

        int connectedEdge = connectedEdges[0];
        int source = graph->getEdgeSource(connectedEdge);
        McDArray<McVec3f> edgePoints = graph->getEdgePoints(connectedEdge);
        if (source != selectedVertex) {
            edgePoints.reverse();
        }

        point2 = getNextPointAtDist(edgePoints, 1000);
        McVec3f direction = point2 - point1;
        direction.normalize();
        directions[i] = direction;
    }
}

// selects the points that should be matched from upper and lower slice

static void selectMatchingVertices(SpatialGraphSelection& selectedVertices,
                                   HxSpatialGraph* graph, int sliceNum) {
    EdgeVertexAttribute* transformAtt = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, "TransformInfo"));
    EdgeVertexAttribute* matchingAtt = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, "Matching"));
    selectedVertices.clear();
    for (int i = 0; i < graph->getNumVertices(); i++) {
        int* trafoDat = transformAtt->intDataAtIdx(i);
        if (*trafoDat == sliceNum) {
            int* pDat = matchingAtt->intDataAtIdx(i);
            if (*pDat == 1) {
                selectedVertices.selectVertex(i);
            }
        }
    }
}

static void makeDefaultConfig(mtalign::PGMPairWeightsParams& config) {
    config.distanceThreshold3d = 1000;
    config.useDistanceThreshold3d = true;
    config.distanceThresholdProjected = 600;
    config.useDistanceThresholdProjected = true;
    config.useAngleThreshold = true;
    config.angleThreshold = 30;

    config.angleWeightParam = 10;
    config.useAngleWeight = true;
    config.dist3dParam = 250;
    config.useDist3dWeight = true;
    config.distProjectedParam = 150;
    config.useProjectedDistWeight = true;
    config.dummySignificance = 0.01;

    // linear or exponential?
    config.weightType = mtalign::PGMPairWeightsParams::EXPONENTIAL;
}

static void getMatchings(McString filename, McDArray<McVec2i>& pairs,
                         McDArray<McVec2i>& evidence) {
    TestingObjectPoolCleaner cleaner;

    TestingData data(filename.dataPtr());
    ASSERT_TRUE(data.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = data.get<HxSpatialGraph>();

    SpatialGraphSelection matchingSelection1(sg);
    selectMatchingVertices(matchingSelection1, sg, 1);

    SpatialGraphSelection matchingSelection2(sg);
    selectMatchingVertices(matchingSelection2, sg, 2);

    McDArray<McVec3f> pointSet2;
    McDArray<McVec3f> pointSet1;

    SpatialGraphSelection::Iterator iter1(matchingSelection1);
    iter1.vertices.reset();
    int vertex = iter1.vertices.nextSelected();
    McDArray<int> points1Indices, points2Indices;
    while (vertex >= 0) {
        pointSet1.append(sg->getVertexCoords(vertex));
        points1Indices.append(vertex);
        cout << "\nlabel " << points1Indices.size() - 1 << " maps to "
             << vertex;
        vertex = iter1.vertices.nextSelected();
    }
    SpatialGraphSelection::Iterator iter2(matchingSelection2);
    iter2.vertices.reset();
    vertex = iter2.vertices.nextSelected();
    while (vertex >= 0) {
        pointSet2.append(sg->getVertexCoords(vertex));
        points2Indices.append(vertex);
        cout << "\nstate " << points2Indices.size() - 1 << " maps to "
             << vertex;
        vertex = iter2.vertices.nextSelected();
    }

    McDArray<McVec3f> directions1;
    McDArray<McVec3f> directions2;
    McDArray<McVec2i> evidenceDummy;

    cout << "getting dirs\n";
    getDirections(matchingSelection1, sg, directions1);
    getDirections(matchingSelection2, sg, directions2);

    cout << "got dirs\n";
    mtalign::PGMPairWeightsParams config;
    makeDefaultConfig(config);
    mtalign::PGMPairWeights ppwc(pointSet1, pointSet2, directions1, directions2,
                              config);
    mtalign::PGMMatcher bfom(pointSet1, pointSet2, directions1, directions2,
                             evidenceDummy, ppwc, 150);

    vector<mtalign::ConnectedFactorGraph> graphs;

    bfom.getAllConnectedComponents(graphs);
    cout << "got " << graphs.size() << " conn comp\n";

    McDArray<int> notConverged;
    McDArray<int> ambiguousAssignments;
    McDArray<int> groupsWithAmbiguities;
    McDArray<int> ambiguousIndexes;
    for (unsigned int l = 0; l < graphs.size(); l++) {
        FactorGraph fg(graphs[l].factors);

        ofstream dotfile;
        McString x;
        char* v = x.printf("%s-group%d.dot", filename.dataPtr(), l);
        cout << "\n*********************************************\n\n" << v;
        dotfile.open(v);
        fg.printDot(dotfile);
        dotfile.close();

        if (graphs[l].factors.size() > 0) {
            McDArray<McVec2i> matchedPointPairs;

            int ambSize = ambiguousAssignments.size();
            McDArray<int> marginalDifferenceAmbiguities;
            McDArray<float> messageDiffs(bfom.mCoords1.size()),
                beliefEntropies(bfom.mCoords1.size());

            bool converged = bfom.matchPoints(
                graphs[l], 10, 1.e-7, 0.0, matchedPointPairs,
                ambiguousAssignments, marginalDifferenceAmbiguities,
                messageDiffs, beliefEntropies);

            cout << "\nconverged:" << converged;
            if (!converged) {
                notConverged.append(l);
            }

            if (ambSize != ambiguousAssignments.size())
                groupsWithAmbiguities.append(l);

            // map pairs to indices:
            for (int h = 0; h < matchedPointPairs.size(); h++) {
                cout << "\n iter matches";
                if (matchedPointPairs[h].y >= 0) {
                    pairs.append(
                        McVec2i(points1Indices[matchedPointPairs[h].x],
                                points2Indices[matchedPointPairs[h].y]));
                    cout << "\nAdd pair "
                         << points1Indices[matchedPointPairs[h].x] << " - "
                         << points2Indices[matchedPointPairs[h].y];
                }
            }
        }
    }

    HxSpatialGraph* dbgSG = new HxSpatialGraph();
    for (int i = 0; i < bfom.mCoords1.size(); i++) {
        int v1 = dbgSG->addVertex(bfom.mCoords1[i]);
        int v2 =
            dbgSG->addVertex(bfom.mCoords1[i] + bfom.mDirections1[i] * 2000);
        dbgSG->addEdge(v1, v2);
    }
    //     HxSpatialGraph* dbgSG=new HxSpatialGraph();
    for (int i = 0; i < bfom.mCoords2.size(); i++) {
        int v1 = dbgSG->addVertex(bfom.mCoords2[i]);
        int v2 =
            dbgSG->addVertex(bfom.mCoords2[i] + bfom.mDirections2[i] * 2000);
        dbgSG->addEdge(v1, v2);
    }
    dbgSG->saveAmiraMesh("Projected.am");
}

static void getMatchings(McString filename, McDArray<McVec2i>& pairs) {
    McDArray<McVec2i> evidence;
    getMatchings(filename, pairs, evidence);
}

TEST(mtalign__PGMPairMatcherTest, testLBPResultReallyKlp7Group207_GUI_E6MS) {
    TestingDevNullRedirect silentout(stdout);

    McDArray<McVec2i> pairs;
    getMatchings("spatialgraph/klp7-207.am", pairs);
    EXPECT_EQ(pairs.size(), 3);
    if (pairs.size() != 3)
        return;
    EXPECT_EQ(pairs[2].x, 1);
    EXPECT_EQ(pairs[2].y, 6);
    EXPECT_EQ(pairs[1].x, 3);
    EXPECT_EQ(pairs[1].y, 10);
    EXPECT_EQ(pairs[0].x, 4);
    EXPECT_EQ(pairs[0].y, 9);
}

TEST(mtalign__PGMPairMatcherTest, testObviousMatchTwoOnlyReally_GUI_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    McDArray<McVec2i> pairs;
    McDArray<McVec2i> evidence;
    getMatchings("spatialgraph/obviousMatchTwoOnly.am", pairs, evidence);

    EXPECT_EQ(pairs.size(), 2);
    if (pairs.size() != 2)
        return;
    EXPECT_EQ(pairs[0].x, 0);
    EXPECT_EQ(pairs[0].y, 7);
    EXPECT_EQ(pairs[1].x, 2);
    EXPECT_EQ(pairs[1].y, 5);
}

TEST(mtalign__PGMPairMatcherTest, testThreeShouldMatchReally_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    McDArray<McVec2i> pairs;
    McDArray<McVec2i> evidence;
    getMatchings("spatialgraph/threeShouldMatch.am", pairs, evidence);

    EXPECT_EQ(pairs.size(), 3);
    if (pairs.size() != 3)
        return;
    EXPECT_EQ(pairs[0].x, 0);
    EXPECT_EQ(pairs[0].y, 7);
    EXPECT_EQ(pairs[1].x, 2);
    EXPECT_EQ(pairs[1].y, 9);
    EXPECT_EQ(pairs[2].x, 4);
    EXPECT_EQ(pairs[2].y, 11);
}

TEST(mtalign__PGMPairMatcherTest, testLBPResultReallyInaNoConnIII_GUI_E4MS) {
    TestingDevNullRedirect silentout(stdout);
    TestingDevNullRedirect silenterr(stderr);

    McDArray<McVec2i> pairs;
    McDArray<McVec2i> evidence;

    getMatchings("spatialgraph/inaNoConnIII.am", pairs, evidence);

    EXPECT_EQ(pairs.size(), 11);
}

TEST(mtalign__PGMPairMatcherTest,
     testLBPResultReallyInaNoConnectionWhy_GUI_E6MS) {
    TestingDevNullRedirect silentout(stdout);

    McDArray<McVec2i> pairs;
    getMatchings("spatialgraph/inaNoConnectionWhy.am", pairs);

    EXPECT_EQ(pairs.size(), 4);
    if (pairs.size() != 4)
        return;

    EXPECT_EQ(pairs[0].x, 0);
    EXPECT_EQ(pairs[0].y, 11);
    EXPECT_EQ(pairs[3].x, 2);
    EXPECT_EQ(pairs[3].y, 9);
    EXPECT_EQ(pairs[1].x, 4);
    EXPECT_EQ(pairs[1].y, 13);
    EXPECT_EQ(pairs[2].x, 6);
    EXPECT_EQ(pairs[2].y, 15);
}

TEST(mtalign__PGMPairMatcherTest, checkConnCompReally_GUI_E4MS) {
    int adj[25] = { 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,
                    0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1 };
    McDMatrix<int> adjMat(5, 5, adj);
    McDMatrix<int> adjMatWo;
    McDArray<int> connComp;
    mtalign::PGMMatcher::getConnectedComponent(adjMat, adjMatWo, connComp);
    int expected[25] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
    for (int i = 0; i < 25; i++) {
        EXPECT_EQ(adjMatWo.dataPtr()[i], expected[i]);
    }
    EXPECT_EQ(connComp.size(), 4);
    int expConnComp[4] = { 0, 1, 2, 4 };
    for (int i = 0; i < 4; i++) {
        EXPECT_EQ(expConnComp[i], connComp[i]);
    }
}

TEST(mtalign__PGMPairMatcherTest, checkAdjMatrixReally_GUI_E4MS) {

    McDMatrix<int> adjMat;
    int vals[25] = { 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1,
                     0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0 };
    McDMatrix<int> variableAssignmentMat(5, 5, vals);
    mtalign::PGMMatcher::createAdjacenceMatrix(variableAssignmentMat, adjMat);
    int expected[25] = { 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1,
                         0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1 };
    for (int i = 0; i < 25; i++) {
        EXPECT_EQ(adjMat.dataPtr()[i], expected[i]);
    }
}

static void setDiagOff(Factor& fac) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            fac.set(i * 5 + j, 1);
            if (i == j)
                fac.set(i * 5 + j, 0);
        }
    }
}

static void setSameMove(Factor& fac, int varCard) {
    for (int i = 0; i < varCard; i++) {
        for (int j = 0; j < varCard; j++) {
            fac.set(i * varCard + j, 1);
            if (abs(i - j) != 1)
                fac.set(i * varCard + j, 0);
        }
    }
}

static int checkDiagOff(Factor& fac) {
    TProb<double>& vals = fac.p();
    float expVals[25] = { 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0,
                          1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0 };
    for (int i = 0; i < 25; i++) {
        if (expVals[i] != vals[i])
            return 0;
    }
    return 1;
}

static int checkSameMove(Factor& fac) {
    TProb<double>& vals = fac.p();
    float expVals[25] = { 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0,
                          1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0 };
    for (int i = 0; i < 25; i++) {
        if (expVals[i] != vals[i])
            return 0;
    }
    return 1;
}

TEST(mtalign__PGMPairMatcherTest, simpleFactorGraph_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);
    TestingDevNullRedirect silenterr(stderr);

    Var p1(0, 5);
    Var p2(1, 5);
    Var p3(2, 5);
    Factor p1Fac(p1);
    Factor p3Fac(p3);
    Factor p2Fac(p2);

    p1Fac.set(0, 0.1);
    p1Fac.set(1, 0.4);
    p1Fac.set(2, 0.3);
    p1Fac.set(3, 0.2);
    p1Fac.set(4, 0.1);
    p2Fac.set(0, 0.1);
    p2Fac.set(1, 0.5);
    p2Fac.set(2, 0.4);
    p2Fac.set(3, 0.3);
    p2Fac.set(4, 0.1);
    p3Fac.set(0, 0.1);
    p3Fac.set(1, 0.4);
    p3Fac.set(2, 0.5);
    p3Fac.set(3, 0.4);
    p3Fac.set(4, 0.1);

    Factor p1p2Fac(VarSet(p1, p2));
    Factor p2p3Fac(VarSet(p2, p3));
    Factor p1p3Fac(VarSet(p1, p3));
    setDiagOff(p1p2Fac);
    setDiagOff(p1p3Fac);
    setDiagOff(p2p3Fac);
    Factor p1p2FacMove(VarSet(p1, p2));
    Factor p2p3FacMove(VarSet(p2, p3));
    Factor p1p3FacMove(VarSet(p1, p3));
    setSameMove(p1p2FacMove, 5);
    setSameMove(p1p3FacMove, 5);
    setSameMove(p2p3FacMove, 5);
    EXPECT_EQ(checkSameMove(p2p3FacMove), 1);

    // create factor graph
    vector<Factor> pFactors;
    pFactors.push_back(p1Fac);
    pFactors.push_back(p2Fac);
    pFactors.push_back(p3Fac);
    pFactors.push_back(p1p2Fac);
    pFactors.push_back(p1p3Fac);
    pFactors.push_back(p2p3Fac);

    FactorGraph pNetwork(pFactors);
    cout << "\n is tree " << pNetwork.isTree() << "\n";

    PropertySet opts;
    opts.set("maxiter", string("10000"));  // Maximum number of iterations
    opts.set("tol", string("1.e-8"));      // Tolerance for convergence
    opts.set("verbose", string("1"));  // Verbosity (amount of output generated)
    BP ia(pNetwork,
          opts("updates", string("SEQMAX"))("logdomain", false)(
              "inference", string("MAXPROD"))("damping", string("0.0")));
    ia.init();
    ia.run();

    Factor x = p1Fac;
    x *= p2Fac;
    x *= p3Fac;
    x *= p1p2Fac;
    x *= p1p3Fac;
    x *= p2p3Fac;
    cout << "\nmax:" << x.max();
    float maxVal = x.max();
    EXPECT_FLOAT_EQ(0.064, maxVal);

    Factor beliefP1 = ia.belief(p1);
    Factor beliefP2 = ia.belief(p2);
    Factor beliefP3 = ia.belief(p3);
    cout << "\n beliefP1" << beliefP1;
    cout << "\n beliefP2" << beliefP2;
    cout << "\n beliefP3" << beliefP3;

    EXPECT_FLOAT_EQ(1, checkDiagOff(p1p2Fac));
    EXPECT_FLOAT_EQ(1, checkDiagOff(p1p3Fac));
    EXPECT_FLOAT_EQ(1, checkDiagOff(p2p3Fac));

    // now: clamp one value, see if it is a tree!

    FactorGraph clampedGraph = pNetwork.clamped(0, 1);
    cout << "\nafter clamping  is tree " << clampedGraph.isTree() << "\n";

    BP iaClamped(
        clampedGraph,
        opts("updates", string("SEQMAX"))("logdomain", false)(
            "inference", string("MAXPROD"))("damping", string("0.001")));

    iaClamped.init();
    iaClamped.run();
    std::vector<std::size_t> maxesClamped = iaClamped.findMaximum();
    vector<std::size_t>::iterator itClamped = maxesClamped.begin();
    for (; itClamped < maxesClamped.end(); itClamped++) {
        cout << "\nclamped:" << *itClamped;
    }
    Factor beliefP1Clamped = iaClamped.belief(p1);
    Factor beliefP2Clamped = iaClamped.belief(p2);
    Factor beliefP3Clamped = iaClamped.belief(p3);
    cout << "\n beliefP1" << beliefP1Clamped;
    cout << "\n beliefP2" << beliefP2Clamped;
    cout << "\n beliefP3" << beliefP3Clamped;
}
