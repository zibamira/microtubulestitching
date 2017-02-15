#include <hxalignmicrotubules/hxtesting.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/HierarchicalLabels.h>

static void addSegment(HxSpatialGraph* sg, McDArray<McVec3f>& points) {
    const int node0 = sg->addVertex(points[0]);
    const int node1 = sg->addVertex(points.last());
    sg->addEdge(node0, node1, points);
}

static McHandle<HxSpatialGraph> makeSpatialGraph() {
    McHandle<HxSpatialGraph> sg(HxSpatialGraph::createInstance());

    for (int z = 0; z < 3; z++) {
        McDArray<McVec3f> points;
        points.append(McVec3f(0, 0, float(z) + 0));
        points.append(McVec3f(0.5, 0.5, float(z) + 0.5));
        points.append(McVec3f(0.9, 0.9, float(z) + 0.9));
        addSegment(sg, points);
    }

    return sg;
}

static void addOrderedNodeIntAttribute(HxSpatialGraph* sg) {
    const char* attrName = "section";
    const int nDataVar = 1;
    EdgeVertexAttribute* attr =
        static_cast<EdgeVertexAttribute*>(sg->addAttribute(
            attrName, HxSpatialGraph::VERTEX, McPrimType::MC_INT32, nDataVar));
    mcassert(attr);
    for (int i = 0; i < sg->getNumVertices(); i++) {
        // Each section contains two nodes.
        attr->setIntDataAtIdx(i, 1 + i / 2);
    }
}

static void addUnorderedNodeIntAttribute(HxSpatialGraph* sg) {
    const char* attrName = "unordered";
    const int nDataVar = 1;
    EdgeVertexAttribute* attr =
        static_cast<EdgeVertexAttribute*>(sg->addAttribute(
            attrName, HxSpatialGraph::VERTEX, McPrimType::MC_INT32, nDataVar));
    mcassert(attr);
    for (int i = 0; i < sg->getNumVertices(); i++) {
        attr->setIntDataAtIdx(i, i % 3);
    }
}

McHandle<HxSpatialGraph> hxtesting::makeSpatialGraphThreeStackedSections() {
    McHandle<HxSpatialGraph> sg(makeSpatialGraph());
    addOrderedNodeIntAttribute(sg);
    addUnorderedNodeIntAttribute(sg);
    return sg;
}
