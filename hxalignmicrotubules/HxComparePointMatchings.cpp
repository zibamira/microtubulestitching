#include <hxalignmicrotubules/HxComparePointMatchings.h>

#include <cstdlib>

#include <hxcore/HxMessage.h>
#include <hxspatialgraph/HierarchicalLabels.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/McMath.h>

HX_INIT_CLASS(HxComparePointMatchings, HxCompModule);

HxComparePointMatchings::HxComparePointMatchings()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portPairLabel(this, "MatchingLabelToCompareTo", 1),
      portFalsePositivePairs(this, "falsePositivePairs", 1),
      portFalseNegativePairs(this, "falseNegativePairs", 1),
      portDisagreementPairs(this, "disagreementPairs", 1),
      portTruePositivePairs(this, "correctPairs", 1),
      portTreeNodes(this, "treeNodes", 1),
      portCriticalNodes(this, "criticalNodes", 1),
      portNumberOfManualPairs(this, "numberOfManualPairs"),
      portNumberOfAutoPairs(this, "portNumberOfAutoPairs"),
      mDoIt(this, "apply") {
    portTreeNodes.hide();
    portCriticalNodes.hide();
}

HxComparePointMatchings::~HxComparePointMatchings() {}

void HxComparePointMatchings::compute() {
    int pairLabelIdx = portPairLabel.getValue(0);
    if (pairLabelIdx < 0)
        return;
    McString matchingLabel =
        qPrintable(QString(portPairLabel.getLabel(pairLabelIdx)));

    if (!mDoIt.wasHit()) {
        return;
    }

    if (portData.source() == NULL) {
        return;
    }
    // Recompute all

    HxSpatialGraph* origGraph =
        dynamic_cast<HxSpatialGraph*>(portData.source());
    HxSpatialGraph* sg = origGraph->duplicate();

    if (sg == NULL) {
        return;
    }

    McDArray<int> allEdgesWithNoTransformInfo;
    int numUserDefinedMatchings =
        createUserDefinedEdges(sg, allEdgesWithNoTransformInfo);
    std::cout << "\nFound " << numUserDefinedMatchings
              << " user defined edges.";
    int numRemoved = removeUserDefinedMatchingsThatAreNotInAdjacentSlices(
        sg, allEdgesWithNoTransformInfo);
    std::cout << "\nRemoved " << numRemoved
              << " user defined edges, because they are in same slice.";

    McDArray<int> userDefinedMatches = allEdgesWithNoTransformInfo;

    EdgeVertexAttribute* matchingAttribute =
        (EdgeVertexAttribute*)(sg->findVertexAttribute(
            matchingLabel.dataPtr()));
    if (!matchingAttribute) {
        theMsg->printf("Could not find matchingLabel.");
        return;
    }

    McDArray<McVec2i> autoPairs;
    getAllAutoPairs(sg, matchingAttribute, autoPairs);
    portNumberOfAutoPairs.setValue(autoPairs.size());
    portNumberOfManualPairs.setValue(userDefinedMatches.size());

    McDArray<McVec2i> falseNegativePairs;
    getAllUnmatchedPairs(sg, matchingAttribute, userDefinedMatches,
                         falseNegativePairs);
    std::cout << "\nFound " << falseNegativePairs.size()
              << " pairs that were not matched in label "
              << matchingLabel.dataPtr();
    rewritePairs(sg, falseNegativePairs, McString("falseNegativePairs"), false);
    portFalseNegativePairs.setValue(falseNegativePairs.size());

    McDArray<McVec2i> falsePositivePairs;
    getAllWronglyMatchedPairs(sg, matchingAttribute, falsePositivePairs);
    std::cout << "\nFound " << falsePositivePairs.size()
              << " pairs that were matched in label " << matchingLabel.dataPtr()
              << " but should not be matched";
    rewritePairs(sg, falsePositivePairs, McString("falsePositivePairs"), false);
    portFalsePositivePairs.setValue(falsePositivePairs.size());

    McDArray<McVec2i> truePositivePairs;
    getAllCorrectlyMatchedPairs(sg, matchingAttribute, userDefinedMatches,
                                truePositivePairs);
    std::cout << "\nFound " << truePositivePairs.size()
              << " pairs that are correct in label " << matchingLabel.dataPtr();
    rewritePairs(sg, truePositivePairs, McString("truePositivePairs"), false);
    portTruePositivePairs.setValue(truePositivePairs.size());

    McDArray<McVec2i> disagreementPairs;
    getAllDisagreeingMatches(sg, matchingAttribute, userDefinedMatches,
                             disagreementPairs);
    std::cout << "\nFound " << disagreementPairs.size()
              << " pairs that were differently macthed by user and label "
              << matchingLabel.dataPtr();
    rewritePairs(sg, disagreementPairs, McString("disagreementPairs"), true);
    portDisagreementPairs.setValue(disagreementPairs.size());

    portTreeNodes.setValue(countLabels(sg, McString("AssignToMakeTree")));
    portCriticalNodes.setValue(countLabels(sg, McString("CriticalNodes")));
    sg->setLabel(origGraph->getLabel() + McString("-with-GT-connections"));
    setResult(sg);
}

int HxComparePointMatchings::countLabels(HxSpatialGraph* graph,
                                         McString attribName) {
    EdgeVertexAttribute* attribute =
        (EdgeVertexAttribute*)(graph->findVertexAttribute(
            attribName.dataPtr()));
    if (!attribute) {
        theMsg->printf("Could not find ambiguity attribute.");
        return 0;
    }

    int numLabels = 0;
    for (int i = 0; i < graph->getNumVertices(); ++i) {
        int val = attribute->getIntDataAtIdx(i);
        if (val > 0)
            numLabels++;
    }
    return numLabels;
}

void HxComparePointMatchings::rewritePairs(HxSpatialGraph* graph,
                                           const McDArray<McVec2i>& newPairs,
                                           const McString& pairsName,
                                           const bool writeLabels) {

    EdgeVertexAttribute* att = dynamic_cast<EdgeVertexAttribute*>(
        graph->addAttribute(pairsName.dataPtr(), HxSpatialGraph::VERTEX,
                            McPrimType::mc_int32, 1));
    HierarchicalLabels* labelGroup = NULL;
    if (writeLabels) {
        graph->addNewLabelGroup(pairsName.dataPtr(), false, true);

        labelGroup = graph->getLabelGroup(pairsName.dataPtr());
    }

    // find all existing pairs

    int maxAttVal = 0;
    for (int i = 0; i < graph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (pairNum > maxAttVal)
            maxAttVal = pairNum;
    }
    McDArray<McVec2i> existingPairs(maxAttVal + 1);
    for (int i = 0; i < existingPairs.size(); i++) {
        existingPairs[i] = McVec2i(-1, -1);
    }

    // find mapping from unsorted labels to new labels
    for (int i = 0; i < graph->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (pairNum > 0) {
            if (existingPairs[pairNum].x == -1) {
                existingPairs[pairNum].x = i;
            } else if (existingPairs[pairNum].y == -1) {
                existingPairs[pairNum].y = i;
            }
        }
    }

    existingPairs.appendArray(newPairs);

    // rewrite labels
    if (writeLabels)
        labelGroup->removeChildLabels();

    int labelCounter = 0;
    for (int i = 0; i < existingPairs.size(); ++i) {

        if ((existingPairs[i].x > -1) && (existingPairs[i].y > -1)) {
            mcassert(existingPairs[i].x != existingPairs[i].y);
            McString s;
            s.printf("Pair%d", labelCounter);
            SbColor color;
            color[0] = float(rand()) / float(RAND_MAX);
            color[1] = float(rand()) / float(RAND_MAX);
            color[2] = float(rand()) / float(RAND_MAX);
            int val;
            if (writeLabels)
                val = graph->addLabel(pairsName.dataPtr(), 0, s.getString(),
                                      color);
            else
                val = i;
            att->setIntDataAtIdx(existingPairs[i].x, val);
            att->setIntDataAtIdx(existingPairs[i].y, val);

            labelCounter++;
        }
    }
}

void HxComparePointMatchings::getAllDisagreeingMatches(
    const HxSpatialGraph* sg, const EdgeVertexAttribute* matchingAttribute,
    const McDArray<int>& userDefinedMatches,
    McDArray<McVec2i>& disagreeingMatches) {
    for (int i = 0; i < userDefinedMatches.size(); i++) {
        int edgeSource = sg->getEdgeSource(userDefinedMatches[i]);
        int edgeTarget = sg->getEdgeTarget(userDefinedMatches[i]);
        int pairLabelSource = matchingAttribute->getIntDataAtIdx(edgeSource);
        int pairLabelTarget = matchingAttribute->getIntDataAtIdx(edgeTarget);
        if (((pairLabelSource > 2) || (pairLabelTarget > 2)) &&
            (pairLabelTarget != pairLabelSource)) {
            disagreeingMatches.append(McVec2i(edgeSource, edgeTarget));
        }
    }
}

void HxComparePointMatchings::updatePairLabelPort() {
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portData);

    if (graph == NULL) {
        return;
    }

    int numVertexAttributes = graph->numAttributes(HxSpatialGraph::VERTEX);
    int vertexItemIdx = 0;
    for (int i = 0; i < numVertexAttributes; ++i) {
        const EdgeVertexAttribute* attrib =
            dynamic_cast<const EdgeVertexAttribute*>(
                graph->attribute(HxSpatialGraph::VERTEX, i));
        if (attrib->nDataVar() == 1) {
            portPairLabel.setNum(portPairLabel.getNum() + 1);
            portPairLabel.setLabel(vertexItemIdx++, attrib->getName());
        }
    }
}

void HxComparePointMatchings::getAllCorrectlyMatchedPairs(
    const HxSpatialGraph* sg, const EdgeVertexAttribute* matchingAttribute,
    const McDArray<int>& userDefinedMatches, McDArray<McVec2i>& correctPairs) {
    for (int i = 0; i < userDefinedMatches.size(); i++) {
        int edgeSource = sg->getEdgeSource(userDefinedMatches[i]);
        int edgeTarget = sg->getEdgeTarget(userDefinedMatches[i]);
        int pairLabelSource = matchingAttribute->getIntDataAtIdx(edgeSource);
        int pairLabelTarget = matchingAttribute->getIntDataAtIdx(edgeTarget);
        if ((pairLabelSource > 2) && (pairLabelTarget > 2) &&
            (pairLabelTarget == pairLabelSource)) {
            correctPairs.append(McVec2i(edgeSource, edgeTarget));
        }
    }
}
void HxComparePointMatchings::getAllWronglyMatchedPairs(
    const HxSpatialGraph* sg, const EdgeVertexAttribute* matchingAttribute,
    McDArray<McVec2i>& falsePositivePairs) {
    McDArray<McVec2i> autoPairs;
    getAllAutoPairs(sg, matchingAttribute, autoPairs);
    std::cout << "\nFound " << autoPairs.size() << " automatedPairs";
    for (int i = 0; i < autoPairs.size(); i++) {
        if (!hasNodeConnectionToOtherSlice(sg, autoPairs[i].x) &&
            !hasNodeConnectionToOtherSlice(sg, autoPairs[i].y)) {
            falsePositivePairs.append(autoPairs[i]);
        }
    }
}
bool
HxComparePointMatchings::hasNodeConnectionToOtherSlice(const HxSpatialGraph* sg,
                                                       const int nodeIdx) {
    EdgeVertexAttribute* transInf =
        (EdgeVertexAttribute*)(sg->findEdgeAttribute("TransformInfo"));
    if (!transInf) {
        theMsg->printf("Could not find TransformInfo.");
        return false;
    }
    McSmallArray<int, 2> incidentEdges = sg->getIncidentEdges(nodeIdx);
    if (incidentEdges.size() == 0)
        return false;

    int curSlice = transInf->getIntDataAtIdx(incidentEdges[0]);
    for (int i = 1; i < incidentEdges.size(); i++) {
        int otherSlice = transInf->getIntDataAtIdx(incidentEdges[i]);
        if (otherSlice != curSlice)
            return true;
    }
    return false;
}
void HxComparePointMatchings::getAllAutoPairs(
    const HxSpatialGraph* sg, const EdgeVertexAttribute* matchingAttribute,
    McDArray<McVec2i>& autoPairs) {
    for (int i = 0; i < sg->getNumVertices(); i++) {
        int pairLabel = matchingAttribute->getIntDataAtIdx(i);
        if (pairLabel > 2) {
            for (int j = i + 1; j < sg->getNumVertices(); j++) {
                int otherPairLabel = matchingAttribute->getIntDataAtIdx(j);
                if (pairLabel == otherPairLabel)
                    autoPairs.append(McVec2i(i, j));
            }
        }
    }
}

void HxComparePointMatchings::getAllUnmatchedPairs(
    const HxSpatialGraph* sg, const EdgeVertexAttribute* matchingAttribute,
    const McDArray<int>& userDefinedMatches,
    McDArray<McVec2i>& falseNegativePairs) {
    for (int i = 0; i < userDefinedMatches.size(); i++) {
        int edgeSource = sg->getEdgeSource(userDefinedMatches[i]);
        int edgeTarget = sg->getEdgeTarget(userDefinedMatches[i]);
        int pairLabelSource = matchingAttribute->getIntDataAtIdx(edgeSource);
        int pairLabelTarget = matchingAttribute->getIntDataAtIdx(edgeTarget);
        if ((pairLabelSource < 3) && (pairLabelTarget < 3)) {
            falseNegativePairs.append(McVec2i(edgeSource, edgeTarget));
        }
    }
}
int
HxComparePointMatchings::removeUserDefinedMatchingsThatAreNotInAdjacentSlices(
    const HxSpatialGraph* sg, McDArray<int>& allEdgesWithNoTransformInfo) {
    int inputNum = allEdgesWithNoTransformInfo.size();
    EdgeVertexAttribute* transInf =
        (EdgeVertexAttribute*)(sg->findVertexAttribute("TransformInfo"));
    if (!transInf) {
        theMsg->printf("Could not find TransformInfo.");
        return -1;
    }
    for (int i = 0; i < allEdgesWithNoTransformInfo.size(); i++) {
        int edgeSource = sg->getEdgeSource(allEdgesWithNoTransformInfo[i]);
        int edgeTarget = sg->getEdgeTarget(allEdgesWithNoTransformInfo[i]);
        int sliceNumSource = transInf->getIntDataAtIdx(edgeSource);
        int sliceNumTarget = transInf->getIntDataAtIdx(edgeTarget);
        // this is not perfect because the label need not be 1,2,3,4.... but
        // could also be 1,4,7,2,...
        // but usually it is ok.
        if (std::abs(sliceNumSource - sliceNumTarget) != 1) {
            allEdgesWithNoTransformInfo.remove(i);
            i--;
        }
    }
    return inputNum - allEdgesWithNoTransformInfo.size();
}

int HxComparePointMatchings::createUserDefinedEdges(
    HxSpatialGraph* sg, McDArray<int>& userDefinedEdges) {
    EdgeVertexAttribute* userMatchingsAttrib =
        (EdgeVertexAttribute*)(sg->findVertexAttribute("UserDefinedMatchings"));
    if (!userMatchingsAttrib) {
        theMsg->printf("Could not find ground truth labels.");
        return -1;
    }
    for (int i = 0; i < sg->getNumVertices(); i++) {
        int attrib1 = userMatchingsAttrib->getIntDataAtIdx(i);
        if (attrib1 > 0) {
            for (int j = i + 1; j < sg->getNumVertices(); j++) {
                int attrib2 = userMatchingsAttrib->getIntDataAtIdx(j);
                if (attrib1 == attrib2) {
                    userDefinedEdges.append(sg->addEdge(i, j));
                    break;
                }
            }
        }
    }
    return userDefinedEdges.size();
}

int HxComparePointMatchings::parse(Tcl_Interp* t, int argc, char** argv) {
    return HxCompModule::parse(t, argc, argv);
}

void HxComparePointMatchings::update() {
    if (portData.isNew()) {
        updatePairLabelPort();
    }
}
