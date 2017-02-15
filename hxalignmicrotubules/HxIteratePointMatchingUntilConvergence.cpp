#include <hxalignmicrotubules/HxIteratePointMatchingUntilConvergence.h>

#include <QString>

#include <hxcore/HxMessage.h>
#include <hxcore/HxObjectPool.h>
#include <hxspatialgraph/internal/HierarchicalLabels.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <hxalignmicrotubules/HxComparePointMatchings.h>
#include <hxalignmicrotubules/HxTestPointMatching.h>
#include <hxalignmicrotubules/SpreadSheetWrapper.h>

HX_INIT_CLASS(HxIteratePointMatchingUntilConvergence, HxCompModule);

HxIteratePointMatchingUntilConvergence::HxIteratePointMatchingUntilConvergence()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portEvidenceHeuristic(this, "EvidenceHeuristic", tr("Evidence Heuristic"), 1),
      portNumIterations(this, "numIterations", tr("Num Iterations"), 1),
      connectionToFEModule(this, "FESteeringModule", tr("FE Steering Module"), HxCompModule::getTypeId()),
      connectionToPMEvalModule(this, "PMComparisonModule", tr("PM Comparaison Module"),
                               HxCompModule::getTypeId()),
      mDoIt(this, "apply", tr("Apply")) {
    HxTestPointMatching* testModule = HxTestPointMatching::createInstance();
    theObjectPool->addObject(testModule);
    connectionToFEModule.connect(testModule);
    testModule->setLabel("FEModule");
    connectionToFEModule.setTightness(1);

    HxComparePointMatchings* compareModule = HxComparePointMatchings::createInstance();
    theObjectPool->addObject(compareModule);
    connectionToPMEvalModule.connect(compareModule);
    compareModule->setLabel("ComparisonModule");
    connectionToPMEvalModule.setTightness(1);
    portNumIterations.setValue(10);
}

HxIteratePointMatchingUntilConvergence::
    ~HxIteratePointMatchingUntilConvergence() {}

void HxIteratePointMatchingUntilConvergence::compute() {
    if (!mDoIt.wasHit()) {
        return;
    }
    if (portData.getSource() == NULL) {
        return;
    }

    int numberOfNeededEvidenceAssignments = 1;

    McHandle<SpreadSheetWrapper> resultSpreadSheet = SpreadSheetWrapper::createInstance();

    int numIterations = 0;

    while (numberOfNeededEvidenceAssignments != 0 &&
           numIterations < portNumIterations.getValue()) {
        // run the pm
        HxTestPointMatching* feModule = dynamic_cast<HxTestPointMatching*>(
            hxconnection_cast<HxCompModule>(connectionToFEModule));
        feModule->mDoIt.hit();
        feModule->fire();
        // do the comparison
        HxComparePointMatchings* compareModule =
            dynamic_cast<HxComparePointMatchings*>(
                hxconnection_cast<HxCompModule>(connectionToPMEvalModule));
        compareModule->mDoIt.hit();
        compareModule->fire();
        // assign the evidence

        assignNeededEvidence();

        numberOfNeededEvidenceAssignments = getNumberOfPointsToAssign();
        std::cout << "\n Assigning " << numberOfNeededEvidenceAssignments
                  << " nodes.";
        McString rowId;
        rowId.printf("Iteration %d", numIterations);
        McString itEntry;
        itEntry.printf("%d", numIterations);
        resultSpreadSheet->addEntry(
            qPrintable(QString(portEvidenceHeuristic.getLabel(
                portEvidenceHeuristic.getValue()))),
            "Iteration", rowId, itEntry);
        McString numAssigned;
        numAssigned.printf("%d", numberOfNeededEvidenceAssignments);
        resultSpreadSheet->addEntry(
            qPrintable(QString(portEvidenceHeuristic.getLabel(
                portEvidenceHeuristic.getValue()))),
            "Assigned", rowId, numAssigned);
        QString alignResultsName(portData.getSource()->getLabel());
        alignResultsName += "-alignResults";
        theMsg->printf("Remove %s", qPrintable(alignResultsName));
        HxObject* obj = theObjectPool->findObject(alignResultsName);
        if (obj)
            theObjectPool->removeObject(obj);
        theObjectPool->removeObject(compareModule->getResult());
        numIterations++;
    }

    resultSpreadSheet->setLabel("IterationStatistics");
    setResult(resultSpreadSheet);
    HxComparePointMatchings* compareModule =
        dynamic_cast<HxComparePointMatchings*>(
            hxconnection_cast<HxCompModule>(connectionToPMEvalModule));
    compareModule->select();
}

void HxIteratePointMatchingUntilConvergence::assignNeededEvidence() {
    McDArray<int> nodesToAssign;
    getListOfNeededEvidenceNodes(nodesToAssign);
    addEvidenceForNodes(nodesToAssign);
}

int HxIteratePointMatchingUntilConvergence::getNumberOfPointsToAssign() {
    McDArray<int> nodesToAssign;
    getListOfNeededEvidenceNodes(nodesToAssign);
    return nodesToAssign.size();
}

void HxIteratePointMatchingUntilConvergence::addEvidenceForNodes(
    const McDArray<int>& nodesToAssign) {
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    const EdgeVertexAttribute* evidenceAttrib =
        dynamic_cast<const EdgeVertexAttribute*>(graph->findAttribute(
            HxSpatialGraph::VERTEX, "UserDefinedMatchings"));
    McDArray<int> pairNodes(nodesToAssign.size());
    pairNodes.fill(-1);

    for (int j = 0; j < nodesToAssign.size(); j++) {
        int nodeJLabel = evidenceAttrib->getIntDataAtIdx(nodesToAssign[j]);
        if (nodeJLabel == 0) {
            pairNodes[j] = -1;
            continue;
        }
        for (int i = 0; i < graph->getNumVertices(); i++) {
            if (i != nodesToAssign[j]) {
                int nodeILabel = evidenceAttrib->getIntDataAtIdx(i);
                if (nodeILabel == nodeJLabel) {
                    pairNodes[j] = i;
                }
            }
        }
    }
    for (int j = 0; j < nodesToAssign.size(); j++) {
        addEvidence(nodesToAssign[j], pairNodes[j]);
    }
}

void HxIteratePointMatchingUntilConvergence::getListOfNeededEvidenceNodes(
    McDArray<int>& nodesToAssign) {
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    McString evidenceLabel = qPrintable(QString(
        portEvidenceHeuristic.getLabel(portEvidenceHeuristic.getValue())));
    const EdgeVertexAttribute* evidenceAttrib =
        dynamic_cast<const EdgeVertexAttribute*>(
            graph->findAttribute(HxSpatialGraph::VERTEX, evidenceLabel));
    for (int i = 0; i < graph->getNumVertices(); i++) {
        float test = evidenceAttrib->getIntDataAtIdx(i);
        if (test > 0) {
            nodesToAssign.append(i);
        }
    }
}

int HxIteratePointMatchingUntilConvergence::parse(Tcl_Interp* t, int argc,
                                                  char** argv) {
    return HxCompModule::parse(t, argc, argv);
}

void HxIteratePointMatchingUntilConvergence::update() {
    if (portData.isNew()) {
        HxCompModule* cm =
            hxconnection_cast<HxCompModule>(connectionToPMEvalModule);
        cm->portData.connect(portData.getSource());
        cm = hxconnection_cast<HxCompModule>(connectionToFEModule);
        cm->portData.connect(portData.getSource());
        updateEvidenceToAssignPort();
    }
}

void HxIteratePointMatchingUntilConvergence::updateEvidenceToAssignPort() {
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    if (graph == NULL) {
        return;
    }
    int numVertexAttributes = graph->numAttributes(HxSpatialGraph::VERTEX);
    int vertexItemIdx = 0;
    int curState = portEvidenceHeuristic.getValue();
    for (int i = 0; i < numVertexAttributes; ++i) {
        const EdgeVertexAttribute* attrib =
            dynamic_cast<const EdgeVertexAttribute*>(
                graph->attribute(HxSpatialGraph::VERTEX, i));
        if (attrib->nDataVar() == 1) {
            portEvidenceHeuristic.setNum(portEvidenceHeuristic.getNum() + 1);
            portEvidenceHeuristic.setLabel(vertexItemIdx++, attrib->getName());
        }
    }
    if (portEvidenceHeuristic.getNum() > curState)
        portEvidenceHeuristic.setValue(curState);
}

void HxIteratePointMatchingUntilConvergence::addEvidence(const int node1,
                                                         const int node2) {
    // create label group or get if if it exists already
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);
    graph->addNewLabelGroup("Evidence", false, true);
    HierarchicalLabels* labelGroup = graph->getLabelGroup("Evidence");
    int labelId = labelGroup->getMaxLabelId();
    McString s;
    s.printf("Assignment%d", labelId + 1);
    SbColor color;
    color[0] = float(rand()) / float(RAND_MAX);
    color[1] = float(rand()) / float(RAND_MAX);
    color[2] = float(rand()) / float(RAND_MAX);
    int val = graph->addLabel("Evidence", 0, s.getString(), color);

    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(graph->addAttribute(
            "Evidence", HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));
    att->setIntDataAtIdx(node1, val);
    if (node2 > -1)
        att->setIntDataAtIdx(node2, val);

    graph->touch(HxData::NEW_PARAMETERS);
}
