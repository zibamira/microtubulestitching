#include <set>

#include <QToolBar>
#include <QLineEdit>
#include <QCheckBox>
#include <QAction>
#include <QComboBox>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <hxcore/HxMain.h>
#include <hxcore/HxObjectPool.h>
#include <hxcore/TestingData.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxgtest/hxtesting.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/TestingDevNullRedirect.h>

#include <hxalignmicrotubules/HxMovingLeastSquaresSpatialGraphWarp.h>
#include <hxalignmicrotubules/attic/WriteMatchingProperties2.h>
#include <hxalignmicrotubules/mtalign/SliceSelector.h>

#if defined(HX_OS_WIN)
#define TEST_UNIX(test_case_name, test_name)                                   \
    TEST(DISABLED_##test_case_name, test_name)
#define TEST_P_UNIX(test_case_name, test_name)                                 \
    TEST_P(test_case_name, DISABLED_##test_name)
#else
#define TEST_UNIX(test_case_name, test_name) TEST(test_case_name, test_name)
#define TEST_P_UNIX(test_case_name, test_name) TEST_P(test_case_name, test_name)
#endif

namespace ht = hxtesting;

static void qxtSelectFilamentEditor() {
    QList<QToolBar*> subAppToolBar =
        theMain->findChildren<QToolBar*>("subApplicationToolBar");
    ASSERT_EQ(1, subAppToolBar.size());

    QList<QAction*> filamentEditorAction =
        subAppToolBar[0]->findChildren<QAction*>("Filament Editor");
    ASSERT_EQ(1, filamentEditorAction.size());

    filamentEditorAction[0]->trigger();
}

static QAction* findNamedAction(QToolBar* toolbar) {
    foreach (QAction* a, toolbar->actions()) {
        if (a->text() == "MicrotubuleAlignTool") {
            return a;
        }
    }
    return 0;
}

static void qxtSelectAlignTool() {
    QList<QToolBar*> toolbar =
        theMain->findChildren<QToolBar*>("FilamentEditorToolboxToolbar");
    ASSERT_EQ(1, toolbar.size());

    QAction* alignToolAction = findNamedAction(toolbar[0]);
    ASSERT_TRUE(alignToolAction);
    alignToolAction->trigger();
}

static void qxtActivateAlignTool() {
    qxtSelectFilamentEditor();
    qxtSelectAlignTool();
}

namespace {

class QxtButton {
  public:
    QAbstractButton* const q;

    explicit QxtButton(QAbstractButton* w) : q(w) {}

    void check() {
        ASSERT_TRUE(q->isCheckable());
        q->setChecked(true);
    }

    void uncheck() {
        ASSERT_TRUE(q->isCheckable());
        q->setChecked(false);
    }

    void click() { q->click(); }
};

class QxtCheckBox {
  public:
    QCheckBox* const q;

    explicit QxtCheckBox(QCheckBox* w) : q(w) {}

    void check() {
        ASSERT_TRUE(q->isCheckable());
        q->setChecked(true);
    }

    void uncheck() {
        ASSERT_TRUE(q->isCheckable());
        q->setChecked(false);
    }

    void checkPartially() {
        ASSERT_TRUE(q->isTristate());
        q->setCheckState(Qt::PartiallyChecked);
    }

    void click() { q->click(); }
};

class QxtLineEdit {
  public:
    QLineEdit* const q;

    explicit QxtLineEdit(QLineEdit* w) : q(w) {}

    void setText(const char* text) { q->setText(text); }
};

class QxtComboBox {
  public:
    QComboBox* const q;

    explicit QxtComboBox(QComboBox* w) : q(w) {}

    void setCurrentIndex(const int theIndex) { q->setCurrentIndex(theIndex); }
};

class QxtWidget {
  public:
    QWidget* const q;

    explicit QxtWidget(QWidget* w) : q(w) {}

    QxtWidget find(const char* name) {
        QList<QWidget*> widgets = q->findChildren<QWidget*>(name);
        EXPECT_EQ(1, widgets.size());
        return QxtWidget(widgets[0]);
    }

    QxtButton findButton(const char* name) {
        QList<QAbstractButton*> button =
            q->findChildren<QAbstractButton*>(name);
        EXPECT_EQ(1, button.size());
        return QxtButton(button[0]);
    }

    QxtLineEdit findLineEdit(const char* name) {
        QList<QLineEdit*> te = q->findChildren<QLineEdit*>(name);
        EXPECT_EQ(1, te.size());
        return QxtLineEdit(te[0]);
    }
    QxtComboBox findComboBox(const char* name) {
        QList<QComboBox*> te = q->findChildren<QComboBox*>(name);
        EXPECT_EQ(1, te.size());
        return QxtComboBox(te[0]);
    }
    QxtCheckBox findCheckBox(const char* name) {
        QList<QCheckBox*> te = q->findChildren<QCheckBox*>(name);
        EXPECT_EQ(1, te.size());
        return QxtCheckBox(te[0]);
    }
};

class QxtAlignMicrotubuleTool {
  public:
    QxtWidget tool;

    QxtAlignMicrotubuleTool(QxtWidget w) : tool(w) {}

    QxtAlignMicrotubuleTool defaults() {
        tool.findButton("allRadioButton").check();
        tool.findComboBox("projectionComboBox").setCurrentIndex(8);
        tool.findLineEdit("angularMatchingSmoothing").setText("2000");
        tool.findButton("radioButtonOptAlign").check();
        tool.findComboBox("pmAlgorithmComboBox").setCurrentIndex(2);
        tool.findLineEdit("angleThres").setText("30");
        tool.findLineEdit("includeLineEdit").setText("25");
        tool.findComboBox("transformComboBox").setCurrentIndex(3);
        tool.findLineEdit("scaleTestMin").setText("1.0");
        tool.findLineEdit("scaleTestMax").setText("1.0");
        tool.findLineEdit("dMaxLineEdit").setText("0");
        tool.findLineEdit("dMinLineEdit").setText("0");
        tool.findLineEdit("deltaHLineEdit").setText("1");
        tool.findButton("GapSizeAutoCheckBox").uncheck();
        tool.findLineEdit("threeDDistanceThres").setText("1500");
        tool.findLineEdit("projectedDistanceThres").setText("600");
        tool.findButton("use3dDistanceWeight").uncheck();
        tool.findLineEdit("threeDDistanceWeight").setText("1000");
        tool.findLineEdit("angleWeight").setText("10");
        tool.findLineEdit("projectedDistanceWeight").setText("300");
        tool.findLineEdit("pairFactorParam").setText("200");
        return *this;
    }

    QxtAlignMicrotubuleTool withScaling() {
        tool.findLineEdit("scaleTestMin").setText("0.95");
        tool.findLineEdit("scaleTestMax").setText("1.05");
        tool.findLineEdit("scaleTestIncrement").setText("0.01");
        return *this;
    }

    QxtAlignMicrotubuleTool noScaling() {
        tool.findLineEdit("scaleTestMin").setText("1");
        tool.findLineEdit("scaleTestMax").setText("1");
        return *this;
    }

    QxtAlignMicrotubuleTool withLabels() {
        tool.findButton("createMatchingLabelsCheckBox").check();
        return *this;
    }

    QxtAlignMicrotubuleTool noLabels() {
        tool.findButton("createMatchingLabelsCheckBox").uncheck();
        return *this;
    }

    QxtAlignMicrotubuleTool varyGap(const char* minGap, const char* maxGap,
                                    const char* stepSize) {
        tool.findCheckBox("GapSizeAutoCheckBox").checkPartially();
        tool.findLineEdit("dMinLineEdit").setText(minGap);
        tool.findLineEdit("dMaxLineEdit").setText(maxGap);
        tool.findLineEdit("deltaHLineEdit").setText(maxGap);
        return *this;
    }

    QxtAlignMicrotubuleTool align() {
        tool.findButton("alignButton").click();
        return *this;
    }
};

}  // namespace

static QxtWidget qxtRoot() { return QxtWidget(theMain); }

static QxtAlignMicrotubuleTool qxtAlignSpatialGraphTool() {
    return QxtAlignMicrotubuleTool(
        qxtRoot().find("MicrotubuleAlignSpatialGraphTool"));
}

static int countMatches(HxSpatialGraph* sg, McString attName) {
    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, attName));
    std::set<int> labels;
    for (int i = 0; i < att->size(); i++) {
        labels.insert(att->getIntDataAtIdx(i));
    }
    return labels.size() - 2;
}

static int getNumUsed(HxSpatialGraph* sg, McString attName) {
    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, attName));
    int num = 0;
    for (int i = 0; i < sg->getNumVertices(); i++) {
        int curAtt = att->getIntDataAtIdx(i);
        if (curAtt > 1)
            num++;
    }
    return num;
}

static void checkAllMatchedPairsInAdjacentSlices(HxSpatialGraph* sg,
                                                 McString labelName) {
    EdgeVertexAttribute* pmAtt = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, labelName));
    EdgeVertexAttribute* tiAtt = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::EDGE, "TransformInfo"));
    assert(tiAtt);
    std::set<int> labels;
    for (int i = 0; i < pmAtt->size(); i++) {
        int matchingP1 = pmAtt->getIntDataAtIdx(i);
        if (matchingP1 < 3)
            continue;
        int numMatches = 0;
        for (int j = 0; j < pmAtt->size(); j++) {
            if (i != j) {
                int matchingP2 = pmAtt->getIntDataAtIdx(j);
                int edgeP1 = sg->getIncidentEdges(i)[0];
                int edgeP2 = sg->getIncidentEdges(j)[0];
                int sliceNumP2 = tiAtt->getIntDataAtIdx(edgeP2);
                int sliceNumP1 = tiAtt->getIntDataAtIdx(edgeP1);
                if (matchingP2 == matchingP1) {
                    if (abs(sliceNumP2 - sliceNumP1) != 1)
                        printf("match label: %d \n points: %d %d \n slices: %d "
                               "%d \n",
                               matchingP2, i, j, sliceNumP1, sliceNumP2);
                    EXPECT_TRUE(abs(sliceNumP2 - sliceNumP1) == 1);
                    numMatches++;
                }
            }
        }
        EXPECT_EQ(numMatches, 1);
    }
}

static void checkMatchingsCorrect(McDArray<McVec2i>& expectedPairs,
                                  HxSpatialGraph* sg) {
    EdgeVertexAttribute* pmAtt = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, "PGMAssignedPairs"));
    assert(pmAtt);
    McDArray<int> allMatchedVertices;
    for (int i = 0; i < expectedPairs.size(); i++) {
        int match1 = pmAtt->getIntDataAtIdx(expectedPairs[i].x);
        int match2 = pmAtt->getIntDataAtIdx(expectedPairs[i].y);
        EXPECT_EQ(match1, match2);
        allMatchedVertices.insertSorted(expectedPairs[i].x, mcStandardCompare);
        allMatchedVertices.insertSorted(expectedPairs[i].y, mcStandardCompare);
    }
    for (int i = 0; i < sg->getNumVertices(); i++) {
        int matchAtt = pmAtt->getIntDataAtIdx(i);
        if (allMatchedVertices.findSorted(i, mcStandardCompare) > -1)
            continue;
        EXPECT_LE(matchAtt, 2);
    }
}

static void checkSameCoord(const McVec2d coord1, const McVec2d coord2) {
    float eps = 2.e-2;
    ASSERT_NEAR(coord1.x, coord2.x, eps);
    ASSERT_NEAR(coord1.y, coord2.y, eps);
}

static void checkCoordinatesTheSame(HxSpatialGraph* graph,
                                    HxSpatialGraph* otherGraph) {
    float eps = 2.e-2;
    ASSERT_EQ(graph->getNumVertices(), otherGraph->getNumVertices());
    for (int i = 0; i < graph->getNumVertices(); i++) {
        McVec3f coord1 = graph->getVertexCoords(i);
        McVec3f coord2 = otherGraph->getVertexCoords(i);
        ASSERT_NEAR(coord1.x, coord2.x, eps);
        ASSERT_NEAR(coord1.y, coord2.y, eps);
        ASSERT_NEAR(coord1.z, coord2.z, eps);
    }
    for (int i = 0; i < graph->getNumEdges(); i++) {
        for (int j = 0; j < graph->getNumEdgePoints(i); j++) {
            McVec3f coord1 = graph->getEdgePoint(i, j);
            McVec3f coord2 = otherGraph->getEdgePoint(i, j);
            ASSERT_NEAR(coord1.x, coord2.x, eps);
            ASSERT_NEAR(coord1.y, coord2.y, eps);
            ASSERT_NEAR(coord1.z, coord2.z, eps);
        }
    }
}

static void checkEvidenceNotAssignedTwice(HxSpatialGraph* graph) {
    EdgeVertexAttribute* evidenceAttrib = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, "Evidence"));
    EdgeVertexAttribute* criticalAttrib = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, "CriticalNodes"));
    EdgeVertexAttribute* treeAttrib = dynamic_cast<EdgeVertexAttribute*>(
        graph->findAttribute(HxSpatialGraph::VERTEX, "AssignToMakeTree"));
    for (int i = 0; i < graph->getNumVertices(); i++) {
        if (evidenceAttrib->getIntDataAtIdx(i) > 0) {
            EXPECT_EQ(criticalAttrib->getIntDataAtIdx(i), 0);
            EXPECT_EQ(treeAttrib->getIntDataAtIdx(i), 0);
        }
    }
}

// writes for each node either -1 if no matching pair found or the vertex index
// of the pair opposite vertex
static void getPairs(HxSpatialGraph* graph,
                     const std::string& pairMatchingAttributeName,
                     std::vector<int>& pairs) {
    std::map<int, int> pairLabelToVertexIdxMap;
    EdgeVertexAttribute* pairAttrib =
        dynamic_cast<EdgeVertexAttribute*>(graph->findAttribute(
            HxSpatialGraph::VERTEX, pairMatchingAttributeName.c_str()));
    ASSERT_TRUE(pairAttrib != NULL);
    pairs.resize(graph->getNumVertices());
    // init pairs
    for (int i = 0; i < graph->getNumVertices(); i++) {
        pairs[i] = -1;
    }
    for (int i = 0; i < graph->getNumVertices(); i++) {
        int currentLabel = pairAttrib->getIntDataAtIdx(i);
        if (currentLabel <= 2)
            continue;
        // is there already a node with the same label?
        std::map<int, int>::iterator otherVertex =
            pairLabelToVertexIdxMap.find(currentLabel);
        if (otherVertex == pairLabelToVertexIdxMap.end())
            pairLabelToVertexIdxMap.insert(
                std::pair<int, int>(currentLabel, i));
        else {
            int otherVertexIdx = otherVertex->second;
            EXPECT_EQ(pairs[otherVertexIdx], -1);
            if (pairs[otherVertexIdx] == -1) {
                pairs[otherVertexIdx] = i;
                pairs[i] = otherVertexIdx;
            }
        }
    }
}

static void clearAllAttributes(HxSpatialGraph* graph) {
    int pAttrib = graph->numAttributes(HxSpatialGraph::POINT);
    for (int i = pAttrib - 1; i >= 0; i--) {
        GraphAttribute* attrib = graph->attribute(HxSpatialGraph::POINT, i);
        graph->deleteLabelGroup(attrib->getName());
        graph->deleteAttribute(HxSpatialGraph::POINT, i);
    }
    pAttrib = graph->numAttributes(HxSpatialGraph::EDGE);
    for (int i = pAttrib - 1; i >= 0; i--) {
        GraphAttribute* attrib = graph->attribute(HxSpatialGraph::EDGE, i);
        graph->deleteLabelGroup(attrib->getName());
        graph->deleteAttribute(HxSpatialGraph::EDGE, i);
    }
    pAttrib = graph->numAttributes(HxSpatialGraph::VERTEX);
    for (int i = pAttrib - 1; i >= 0; i--) {
        GraphAttribute* attrib = graph->attribute(HxSpatialGraph::VERTEX, i);
        graph->deleteLabelGroup(attrib->getName());
        graph->deleteAttribute(HxSpatialGraph::VERTEX, i);
    }
    graph->touch();
    graph->fire();
}

static void checkPairsAreTheSame(HxSpatialGraph* graph1, HxSpatialGraph* graph2,
                                 const std::string& pairMatchingAttributeName) {
    std::vector<int> pairs1(graph1->getNumVertices());
    std::vector<int> pairs2(graph2->getNumVertices());
    getPairs(graph1, pairMatchingAttributeName, pairs1);
    getPairs(graph2, pairMatchingAttributeName, pairs2);
    for (unsigned int i = 0; i < pairs1.size(); i++) {
        EXPECT_EQ(pairs1[i], pairs2[i]);
    }
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSGWithoutTransformInfo_TestFEDoesNotCrash_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    clearAllAttributes(sg);
    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults();
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_TestGreedyAndGapVariationOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    TestingData sg_expectedResult_dat(
        "spatialgraph/"
        "spatialgraphstack_3slices_2kvertices_matchedGreedy_withGap.am");
    ASSERT_TRUE(sg_expectedResult_dat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg_expectedResult =
        sg_expectedResult_dat.get<HxSpatialGraph>();

    HxSpatialGraph* origSg = sg->duplicate();
    theObjectPool->addObject(origSg);
    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool()
        .tool.findButton("radioButtonInitialAlign")
        .check();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("transformComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool().align();

    qxtAlignSpatialGraphTool().tool.findButton("radioButtonOptAlign").check();
    qxtAlignSpatialGraphTool().tool.findLineEdit("dMaxLineEdit").setText("500");
    qxtAlignSpatialGraphTool().tool.findLineEdit("dMinLineEdit").setText(
        "-500");
    qxtAlignSpatialGraphTool().tool.findLineEdit("deltaHLineEdit").setText(
        "500");

    qxtAlignSpatialGraphTool().align();
    checkCoordinatesTheSame(sg, sg_expectedResult);
    checkPairsAreTheSame(sg, sg_expectedResult, "GreedyAssignedPairs");
}

TEST_UNIX(MicrotubuleSpatialGraphAligner, givenSG_TestExactAndGapVariationOK_GUI_E3MS) {
    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();
    HxSpatialGraph* origSg = sg->duplicate();
    theObjectPool->addObject(origSg);

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(1);
    qxtAlignSpatialGraphTool()
        .tool.findButton("radioButtonInitialAlign")
        .check();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("transformComboBox")
        .setCurrentIndex(0);

    {
        TestingDevNullRedirect silentout(stdout);
        qxtAlignSpatialGraphTool().align();
    }

    qxtAlignSpatialGraphTool().tool.findButton("radioButtonOptAlign").check();
    qxtAlignSpatialGraphTool().tool.findLineEdit("dMaxLineEdit").setText("500");
    qxtAlignSpatialGraphTool().tool.findLineEdit("dMinLineEdit").setText(
        "-500");
    qxtAlignSpatialGraphTool().tool.findLineEdit("deltaHLineEdit").setText(
        "500");

    {
        TestingDevNullRedirect silentout(stdout);
        qxtAlignSpatialGraphTool().align();
    }

    EXPECT_THAT(sg, ht::EqDataSha1("723736a0d1c7c99648dd5970ab863417488f4904"));
}

TEST(
    MicrotubuleSpatialGraphAligner,
    givenSG_noMatchingStack_TestNoTransformationAppliedWhenTestingScalingForInitTransform_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/notMatchingStack.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();
    HxSpatialGraph* origSg = sg->duplicate();

    theObjectPool->addObject(origSg);
    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findButton("radioButtonInitialAlign")
        .check();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool().tool.findLineEdit("scaleTestMin").setText(
        "0.95");
    qxtAlignSpatialGraphTool().align();
    checkCoordinatesTheSame(sg, origSg);
}

// This test is crashing from time to time on cherry.  It's disabled to avoid
// spurious build failures.  It still fails.
TEST(
    MicrotubuleSpatialGraphAligner,
    DISABLED_givenSG_verySmallGroundTruth_TestWritingMatchingPropertiesOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    float eps = 1.e-5;
    TestingData sgdat("spatialgraph/verySmallGroundTruth.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    WriteMatchingProperties2* writer = new WriteMatchingProperties2();
    theObjectPool->addObject(writer);
    writer->portData.connect(sg);
    writer->portAction.hit();
    writer->fire();
    SpreadSheetWrapper* wrapper =
        dynamic_cast<SpreadSheetWrapper*>(writer->getResult(2));

    McString result;
    wrapper->getEntry("MutualAngles", "MutualAngles", "0", result);
    float expected = 5.431047;
    bool dummy;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);

    wrapper->getEntry("MutualAngles", "MutualAngles", "1", result);
    expected = 9.967590;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);

    wrapper = dynamic_cast<SpreadSheetWrapper*>(writer->getResult(0));
    wrapper->getEntry("3dDistance", "3dDistance", "0", result);
    expected = 381.41741943359375;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);

    wrapper = dynamic_cast<SpreadSheetWrapper*>(writer->getResult(1));
    wrapper->getEntry("ProjectedDistance", "ProjectedDistance", "0", result);
    expected = 377.99539184570313;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);

    wrapper = dynamic_cast<SpreadSheetWrapper*>(writer->getResult(3));
    wrapper->getEntry("PairwiseShiftDiff", "PairwiseShiftDiff", "0", result);
    expected = 56.932323455810547;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);
    wrapper->getEntry("PairwiseShiftDiff", "PairwiseShiftDiff", "1", result);
    expected = 203.24520874023437;
    EXPECT_NEAR(result.toFloat(dummy), expected, eps);
}

TEST(
    MicrotubuleSpatialGraphAligner,
    givenSG_example5WithQueerEvidence_TestNoAssignmentToEvidenceTwice_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);
    TestingDevNullRedirect silenterr(stderr);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample5WithQueerEvidence.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    EdgeVertexAttribute* evidenceAttrib = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, "Evidence"));
    ASSERT_TRUE(evidenceAttrib);

    // check the evidence is there
    EXPECT_GE(evidenceAttrib->getIntDataAtIdx(2), 1);

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    checkEvidenceNotAssignedTwice(sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_example6_TestNoTransformationApplied_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample6.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();
    HxSpatialGraph* origSg = sg->duplicate();
    theObjectPool->addObject(origSg);

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    checkCoordinatesTheSame(sg, origSg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_example6_TestPGMMatchingOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample6.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(0, 7));
    pairs.append(McVec2i(2, 9));
    pairs.append(McVec2i(4, 11));

    checkMatchingsCorrect(pairs, sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_example3_TestPGMMatchingOK_GUI_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample3.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(0, 7));
    pairs.append(McVec2i(2, 5));

    checkMatchingsCorrect(pairs, sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_example2_TestPGMMatchingOK_GUI_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample2.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(1, 6));
    pairs.append(McVec2i(3, 10));
    pairs.append(McVec2i(4, 9));

    checkMatchingsCorrect(pairs, sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_example6_checkAttributeNotExists_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample6.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    EdgeVertexAttribute* pgmAtt = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, "PGMAssignedPairs"));
    EXPECT_TRUE(pgmAtt == NULL);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_randomSmallExample4_TestPGMMatchingOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample4.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(1, 12));
    pairs.append(McVec2i(3, 16));
    pairs.append(McVec2i(5, 10));
    pairs.append(McVec2i(6, 14));
    pairs.append(McVec2i(8, 18));

    checkMatchingsCorrect(pairs, sg);
}
TEST(MicrotubuleSpatialGraphAligner,
     givenSG_randomSmallExample4_TestCleaningPMAttributeOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample6.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    qxtAlignSpatialGraphTool().tool.findLineEdit("includeLineEdit").setText(
        "-1");
    qxtAlignSpatialGraphTool().align();

    EdgeVertexAttribute* pmAtt = dynamic_cast<EdgeVertexAttribute*>(
        sg->findAttribute(HxSpatialGraph::VERTEX, "PGMAssignedPairs"));
    for (int i = 0; i < sg->getNumVertices(); i++) {
        EXPECT_EQ(pmAtt->getIntDataAtIdx(i), 1);
    }
}

TEST(MicrotubuleSpatialGraphAligner,
     DISABLED_givenSG_randomSmallExample5_TestPGMMatchingOK_GUI_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample5.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(0, 25));
    pairs.append(McVec2i(2, 27));
    pairs.append(McVec2i(4, 33));
    pairs.append(McVec2i(8, 23));
    pairs.append(McVec2i(10, 29));

    pairs.append(McVec2i(12, 31));
    pairs.append(McVec2i(14, 39));
    pairs.append(McVec2i(16, 35));
    pairs.append(McVec2i(18, 37));

    checkMatchingsCorrect(pairs, sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSG_smallExample1_TestPGMMatchingOK_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample1.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    McDArray<McVec2i> pairs;
    pairs.append(McVec2i(0, 11));
    pairs.append(McVec2i(2, 9));
    pairs.append(McVec2i(4, 13));
    pairs.append(McVec2i(6, 15));

    checkMatchingsCorrect(pairs, sg);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSpatialGraphStackCheckOptInitOptMatchingOK_greedy_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    float origBB[6];
    float eps = 1.e-7;
    sg->getBoundingBox(origBB);

    qxtActivateAlignTool();

    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("transformComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "GreedyAssignedPairs");

    EXPECT_EQ(getNumUsed(sg, "GreedyAssignedPairs"), 983);

    EXPECT_EQ(countMatches(sg, "GreedyAssignedPairs"), 163);

    float afterFirstOptBB[6];
    sg->getBoundingBox(afterFirstOptBB);
    EXPECT_NEAR(origBB[4], afterFirstOptBB[4], eps);
    EXPECT_NEAR(origBB[5], afterFirstOptBB[5], eps);

    qxtAlignSpatialGraphTool()
        .tool.findButton("radioButtonInitialAlign")
        .check();
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "GreedyAssignedPairs");

    EXPECT_EQ(getNumUsed(sg, "GreedyAssignedPairs"), 188);

    EXPECT_EQ(countMatches(sg, "GreedyAssignedPairs"), 58);

    float afterFirstInitBB[6];
    sg->getBoundingBox(afterFirstInitBB);
    EXPECT_NEAR(origBB[4], afterFirstInitBB[4], eps);
    EXPECT_NEAR(origBB[5], afterFirstInitBB[5], eps);

    qxtAlignSpatialGraphTool().tool.findButton("radioButtonOptAlign").check();
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "GreedyAssignedPairs");

    EXPECT_EQ(getNumUsed(sg, "GreedyAssignedPairs"), 983);

    EXPECT_EQ(countMatches(sg, "GreedyAssignedPairs"), 329);

    float afterSecondOptBB[6];
    sg->getBoundingBox(afterSecondOptBB);
    EXPECT_NEAR(origBB[4], afterSecondOptBB[4], eps);
    EXPECT_NEAR(origBB[5], afterSecondOptBB[5], eps);
}

TEST(
    MicrotubuleSpatialGraphAligner,
    givenSpatialGraphStackCheckMatchingsInAdjacentSlicesAligning1To2_pgm_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();

    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(2);
    qxtAlignSpatialGraphTool().tool.findButton("pairRadioButton").check();
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "PGMAssignedPairs");
}
TEST(
    MicrotubuleSpatialGraphAligner,
    givenSpatialGraphStackCheckMatchingsInAdjacentSlicesAligning1To2_greedy_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();

    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool().tool.findButton("pairRadioButton").check();
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "GreedyAssignedPairs");
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSpatialGraphStackCheckMatchingsInAdjacentSlices_greedy_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();

    qxtAlignSpatialGraphTool().defaults();
    qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(0);
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "GreedyAssignedPairs");
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSpatialGraphStackCheckMatchingsInAdjacentSlices_pgm_GUI_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    qxtActivateAlignTool();

    qxtAlignSpatialGraphTool().defaults();
    // qxtAlignSpatialGraphTool().tool.findComboBox("pmAlgorithmComboBox").setCurrentIndex(1);
    qxtAlignSpatialGraphTool().align();

    checkAllMatchedPairsInAdjacentSlices(sg, "PGMAssignedPairs");
}

TEST(MicrotubuleSpatialGraphAligner,
     sliceSelectorReturnsCorrectNumberOfSlices_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/spatialgraphstack_3slices_2kvertices.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();

    mtalign::SliceSelector selectionHelper(sg, "TransformInfo");
    EXPECT_EQ(selectionHelper.getNumSlices(), 3);
}

TEST(MicrotubuleSpatialGraphAligner,
     givenSpatialWithWarpPairsGetCorrectLandmarks_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    TestingObjectPoolCleaner cleaner;
    TestingData sgdat("spatialgraph/randomSmallExample5.am");
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();
    McDArray<McVec2d> p1;
    McDArray<McVec2d> p2;
    HxMovingLeastSquaresSpatialGraphWarp::prepareLandmarks(p1, p2, sg, 0, 1);
    EXPECT_EQ(p1.size(), p2.size());
    EXPECT_EQ(p1.size(), 3);
    checkSameCoord(p1[0], McVec2d(10802.271484375, 18836.2109375));
    checkSameCoord(p1[1], McVec2d(10095.4697265625, 17387.271484375));
    checkSameCoord(p1[2], McVec2d(8882.12890625, 19154.271484375));
    checkSameCoord(p2[0], McVec2d(10535.7744140625, 19023.763671875));
    checkSameCoord(p2[1], McVec2d(9737.3486328125, 17568.169921875));
    checkSameCoord(p2[2], McVec2d(8674.9599609375, 19236.84765625));
}

namespace {

struct Expectation {
    const char* package;
    const char* path;
    const char* sha1;
};

std::ostream& operator<<(std::ostream& os, const Expectation& e) {
    os << "{ ";
    {
        os << "package: ";
        os << "'" << e.package << "'";
        os << ", ";
        os << "path: ";
        os << "'" << e.path << "'";
        os << ", ";
        os << "sha1: ";
        os << "'" << e.sha1 << "'";
    }
    os << " }";
    return os;
}

class MicrotubuleSpatialGraphAlignerWithRegressionTestingData
    : public ::testing::TestWithParam<Expectation> {};
}

TEST_P_UNIX(MicrotubuleSpatialGraphAlignerWithRegressionTestingData,
            PGMMatchingComputesExpectedBaseline_GUI_E3MS) {
    TestingDevNullRedirect silentout(stdout);
    TestingDevNullRedirect silenterr(stderr);
    TestingObjectPoolCleaner cleaner;

    Expectation e = GetParam();
    TestingData sgdat(e.package, e.path);
    ASSERT_TRUE(sgdat.dataOk<HxSpatialGraph>());

    qxtActivateAlignTool();
    qxtAlignSpatialGraphTool().defaults().align();
    HxSpatialGraph* sg = sgdat.get<HxSpatialGraph>();
    EXPECT_THAT(sg, ht::EqDataSha1(e.sha1));
}

static std::vector<Expectation> makeExpectations() {
    std::vector<Expectation> es;

    Expectation e;
    e.package = "hxalignmicrotubules";

    e.path = "test/data/randomSmallExample1.am";
    e.sha1 = "b781989e77a11541f0ab92048744a159a70a2391";
    es.push_back(e);

    e.path = "test/data/randomSmallExample2.am";
    e.sha1 = "a04203f369ae2e38fcd5d3823d57f5393f31cfff";
    es.push_back(e);

    e.path = "test/data/randomSmallExample3.am";
    e.sha1 = "bc74ac19b69f8c274b6572f9637b6968c068eab3";
    es.push_back(e);

    e.path = "test/data/randomSmallExample4.am";
    e.sha1 = "4ebca0231d0113647d7001a0f45990835506af8f";
    es.push_back(e);

    e.path = "test/data/randomSmallExample5.am";
    e.sha1 = "33459a2bc95851cb76954f65bb6f2547323c21f1";
    es.push_back(e);

    // DISABLED, because the result is unstable on Linux.  The test worked on
    // MacX.  The status on Windows is unknown.
    e.path = "test/data/randomSmallExample5WithQueerEvidence.am";
    e.sha1 = "2aad523db667a1320cd2bb90cb869b21954a2e50";
    // es.push_back(e);

    e.path = "test/data/randomSmallExample5WithUnassignedAndAssignedEvidence.am";
    e.sha1 = "da918504f5ef9c7ad1029d8830436109c507bcc0";
    es.push_back(e);

    e.path = "test/data/randomSmallExample6.am";
    e.sha1 = "21576ec36d1dfeb09b17a14852eabc0bf7c4e5ef";
    es.push_back(e);

    return es;
}

INSTANTIATE_TEST_CASE_P(WithRadomSmallExamples,
                        MicrotubuleSpatialGraphAlignerWithRegressionTestingData,
                        testing::ValuesIn(makeExpectations()));
