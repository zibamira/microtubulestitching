#include <hxalignmicrotubules/attic/HxCopyLabelsFromSpatialGraph.h>

#include <hxspatialgraph/HxSpatialGraph.h>

#include <QString>

HX_INIT_CLASS(HxCopyLabelsFromSpatialGraph, HxCompModule);

HxCopyLabelsFromSpatialGraph::HxCopyLabelsFromSpatialGraph(void)
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      mDoIt(this, "apply"),
      portPairLabel(this, "LabelToCopy", 1),
      portOtherSG(this, "otherSG", HxSpatialGraph::getClassTypeId()) {}

HxCopyLabelsFromSpatialGraph::~HxCopyLabelsFromSpatialGraph(void) {}

void HxCopyLabelsFromSpatialGraph::updatePairLabelPort() {
    HxSpatialGraph* graph = hxconnection_cast<HxSpatialGraph>(portOtherSG);

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

int HxCopyLabelsFromSpatialGraph::parse(Tcl_Interp* t, int argc, char** argv) {
    return HxCompModule::parse(t, argc, argv);
}

void HxCopyLabelsFromSpatialGraph::update() {
    if (portOtherSG.isNew()) {
        updatePairLabelPort();
    }
}

void HxCopyLabelsFromSpatialGraph::compute(void) {
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

    HxSpatialGraph* graph = dynamic_cast<HxSpatialGraph*>(portData.source());
    HxSpatialGraph* otherGraph =
        dynamic_cast<HxSpatialGraph*>(portOtherSG.source());

    if (otherGraph == NULL || graph == NULL) {
        return;
    }

    EdgeVertexAttribute* attribute =
        (EdgeVertexAttribute*)(otherGraph->findVertexAttribute(
            matchingLabel.dataPtr()));
    if (!attribute) {
        theMsg->printf("Could not find label.");
        return;
    }
    if ((attribute->primType() != McPrimType::mc_int32)) {
        theMsg->printf("Can only copy int labels on vertices.");
    }
    if (graph->getNumVertices() != otherGraph->getNumVertices()) {
        theMsg->printf("Number of vertices differ. Cannot copy attribute.");
    }
    EdgeVertexAttribute* copiedAtt = (EdgeVertexAttribute*)(graph->addAttribute(
        "CopiedLabels", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    for (int i = 0; i < graph->getNumVertices(); i++)
        copiedAtt->setIntDataAtIdx(i, attribute->getIntDataAtIdx(i));
}
