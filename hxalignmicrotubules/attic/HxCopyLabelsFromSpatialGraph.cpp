#include <hxalignmicrotubules/attic/HxCopyLabelsFromSpatialGraph.h>

#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <QString>

HX_INIT_CLASS(HxCopyLabelsFromSpatialGraph, HxCompModule);

HxCopyLabelsFromSpatialGraph::HxCopyLabelsFromSpatialGraph(void)
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      mDoIt(this, "apply", tr("Apply")),
      portPairLabel(this, "LabelToCopy", tr("Label To Copy"), 1),
      portOtherSG(this, "otherSG", tr("Other SG"), HxSpatialGraph::getClassTypeId()) {}

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

    if (portData.getSource() == NULL) {
        return;
    }
    // Recompute all

    HxSpatialGraph* graph = dynamic_cast<HxSpatialGraph*>(portData.getSource());
    HxSpatialGraph* otherGraph =
        dynamic_cast<HxSpatialGraph*>(portOtherSG.getSource());

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
    if ((attribute->primType() != McPrimType::MC_INT32)) {
        theMsg->printf("Can only copy int labels on vertices.");
    }
    if (graph->getNumVertices() != otherGraph->getNumVertices()) {
        theMsg->printf("Number of vertices differ. Cannot copy attribute.");
    }
    EdgeVertexAttribute* copiedAtt = (EdgeVertexAttribute*)(graph->addAttribute(
        "CopiedLabels", HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));
    for (int i = 0; i < graph->getNumVertices(); i++)
        copiedAtt->setIntDataAtIdx(i, attribute->getIntDataAtIdx(i));
}
