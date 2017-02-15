#include <hxalignmicrotubules/mtalign/matchingPGM_sg.h>

#include <limits>

#include <QString>

#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/SpatialGraphSelection.h>
#include <mclib/McDArray.h>
#include <mclib/internal/McString.h>

#include <hxalignmicrotubules/mtalign.h>
#include <hxspatialgraph/internal/HierarchicalLabels.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

namespace ma = mtalign;

static int findAttrVal(const EdgeVertexAttribute* attr,
                       const McDArray<int>& vertexIds, const int val) {
    for (int j = 0; j < vertexIds.size(); j++) {
        if (val == attr->getIntDataAtIdx(vertexIds[j])) {
            return j;
        }
    }
    return -1;
}

// `findUnassigned` returns an array of increasing indices into `first` of
// vertices that have an attribute value defined without corresponding
// attribute value for vertices in `second`.
static McDArray<int> findUnassigned(const EdgeVertexAttribute* attr,
                                    const McDArray<int>& first,
                                    const McDArray<int>& second) {
    McDArray<int> unassigned;
    for (int i = 0; i < first.size(); i++) {
        const int val = attr->getIntDataAtIdx(first[i]);
        if (val <= 0) {  // 0 indicates undefined; ignore it.
            continue;
        }
        if (findAttrVal(attr, second, val) == -1) {
            unassigned.append(i);
        }
    }
    return unassigned;
}

static void removeUnassigned(const EdgeVertexAttribute* attr,
                             McDArray<int>& refSlicePointsToMatch,
                             McDArray<McVec3f>& refSlicePointsToMatchCoords,
                             McDArray<McVec3f>& refDirections,
                             SpatialGraphSelection& refSelection,
                             McDArray<int>& transSlicePointsToMatch,
                             McDArray<McVec3f>& transSlicePointsToMatchCoords,
                             McDArray<McVec3f>& transDirections,
                             SpatialGraphSelection& transSelection) {

    const McDArray<int> refUnassigned =
        findUnassigned(attr, refSlicePointsToMatch, transSlicePointsToMatch);
    for (int i = refUnassigned.size() - 1; i >= 0; i--) {
        const int idx = refUnassigned[i];
        refSelection.deselectVertex(refSlicePointsToMatch[idx]);
        refSlicePointsToMatch.remove(idx);
        refSlicePointsToMatchCoords.remove(idx);
        refDirections.remove(idx);
    }

    const McDArray<int> transUnassigned =
        findUnassigned(attr, transSlicePointsToMatch, refSlicePointsToMatch);
    for (int i = transUnassigned.size() - 1; i >= 0; i--) {
        const int idx = transUnassigned[i];
        transSelection.deselectVertex(transSlicePointsToMatch[idx]);
        transSlicePointsToMatch.remove(idx);
        transSlicePointsToMatchCoords.remove(idx);
        transDirections.remove(idx);
    }
}

static McDArray<McVec2i>
getAssignments(const EdgeVertexAttribute* attr,
               const McDArray<int>& refSlicePointsToMatch,
               const McDArray<int>& transSlicePointsToMatch) {
    McDArray<McVec2i> evidence;
    for (int i = 0; i < transSlicePointsToMatch.size(); i++) {
        const int val = attr->getIntDataAtIdx(transSlicePointsToMatch[i]);
        if (val <= 0) {
            continue;
        }
        const int assignment = findAttrVal(attr, refSlicePointsToMatch, val);
        if (assignment != -1) {
            evidence.append(McVec2i(assignment, i));
        }
    }
    return evidence;
}

static void clearAttribute(HxSpatialGraph* sg, const SpatialGraphSelection& sel,
                           const char* attName, const McPrimType primType) {
    mcrequire(primType == McPrimType::MC_INT32 ||
              primType == McPrimType::MC_FLOAT);

    EdgeVertexAttribute* attrib = dynamic_cast<EdgeVertexAttribute*>(
        sg->addAttribute(attName, HxSpatialGraph::VERTEX, primType, 1));
    if (!attrib)
        return;

    for (int i = 0; i < sel.getNumSelectedVertices(); i++) {
        const int curVertex = sel.getSelectedVertex(i);
        if (primType == McPrimType::MC_INT32) {
            attrib->setIntDataAtIdx(curVertex, 0);
        } else {
            attrib->setFloatDataAtIdx(curVertex, 0);
        }
    }
}

static void clearAttributes(HxSpatialGraph* sg,
                            const SpatialGraphSelection& sel) {
    struct Descr {
        const char* name;
        char ty;
    };
    const Descr attrs[] = { { "Ambiguities", 'i' },
                            { "AmbiguitiesFromMarginalDiff", 'i' },
                            { "AssignToMakeTree", 'i' },
                            { "AssignThesePointsStep", 'i' },
                            { "QueerEvidence", 'i' },
                            { "CriticalNodes", 'i' },
                            { "CriticalNodesWithAdaptiveThreshold", 'i' },
                            { "MaxEntropieNodes", 'i' },
                            { "entropie", 'f' },
                            { "maxDiff", 'f' },
                            { "entropieIsZero", 'i' },
                            { "maxDiffIsZero", 'i' } };
    for (unsigned int i = 0; i < sizeof(attrs) / sizeof(attrs[0]); ++i) {
        const Descr a = attrs[i];
        const McPrimType primType =
            (a.ty == 'i') ? McPrimType::MC_INT32 : McPrimType::MC_FLOAT;
        clearAttribute(sg, sel, a.name, primType);
    }
}

static void rewritePairs(HxSpatialGraph* sg, const McDArray<McVec2i>& newPairs,
                         const char* pairsName) {

    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(sg->addAttribute(
            pairsName, HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));

    HierarchicalLabels* labelGroup = NULL;
    sg->addNewLabelGroup(pairsName, false, true);
    labelGroup = sg->getLabelGroup(pairsName);

    // find all existing pairs

    int maxAttVal = 0;
    for (int i = 0; i < sg->getNumVertices(); ++i) {
        int pairNum = att->getIntDataAtIdx(i);
        if (pairNum > maxAttVal)
            maxAttVal = pairNum;
    }
    McDArray<McVec2i> existingPairs(maxAttVal + 1);
    for (int i = 0; i < existingPairs.size(); i++) {
        existingPairs[i] = McVec2i(-1, -1);
    }

    // find mapping from unsorted labels to new labels
    for (int i = 0; i < sg->getNumVertices(); ++i) {
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
    labelGroup->removeChildLabels();

    int labelCounter = 3;
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
            val = sg->addLabel(pairsName, 0, s.getString(), color);
            att->setIntDataAtIdx(existingPairs[i].x, val);
            att->setIntDataAtIdx(existingPairs[i].y, val);

            labelCounter++;
        }
    }
}

// Find nodes for which attribute `labelName` is `true` and also set attribute
// for `newPositives` to true.
static void rewriteBinaryLabel(HxSpatialGraph* sg,
                               const McDArray<int>& newPositives,
                               const char* labelName) {
    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(sg->addAttribute(
            labelName, HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));
    mcassert(att);
    sg->addNewLabelGroup(labelName, false, true);
    HierarchicalLabels* labelGroup = sg->getLabelGroup(labelName);
    if (!labelGroup) {
        return;
    }

    McDArray<int> existingPositives;
    for (int i = 0; i < sg->getNumVertices(); ++i) {
        const int ambNum = att->getIntDataAtIdx(i);
        if (ambNum > 0) {
            existingPositives.append(i);
        }
    }
    existingPositives.appendArray(newPositives);

    labelGroup->removeChildLabels();
    const SbColor color(1, 0, 0);
    const int val = sg->addLabel(labelName, 0, "P", color);
    for (int i = 0; i < existingPositives.size(); ++i) {
        att->setIntDataAtIdx(existingPositives[i], val);
    }
}

static void writeFloatAttribute(HxSpatialGraph* sg,
                                const McDArray<int>& onNodeIds,
                                const McDArray<float>& values,
                                const char* labelName,
                                const char* labelNameZero) {
    EdgeVertexAttribute* attrib =
        dynamic_cast<EdgeVertexAttribute*>(sg->addAttribute(
            labelName, HxSpatialGraph::VERTEX, McPrimType::MC_FLOAT, 1));
    mcassert(attrib);

    EdgeVertexAttribute* attribZero =
        dynamic_cast<EdgeVertexAttribute*>(sg->addAttribute(
            labelNameZero, HxSpatialGraph::VERTEX, McPrimType::MC_INT32, 1));
    mcassert(attribZero);

    for (int i = 0; i < onNodeIds.size(); i++) {
        attrib->setFloatDataAtIdx(onNodeIds[i], values[i]);
        if (values[i] < 1.e-10)
            attribZero->setIntDataAtIdx(onNodeIds[i], 1);
        else
            attribZero->setIntDataAtIdx(onNodeIds[i], 0);
    }
}

static McDArray<int> getSelectedNodeIds(const SpatialGraphSelection& sel) {
    McDArray<int> ids;
    for (int i = 0; i < sel.getNumSelectedVertices(); i++)
        ids.append(sel.getSelectedVertex(i));
    return ids;
}

// Check user-provided assignments: remove facing points that the user
// explicitly marked to be unmatched and return assigned pairs as indices into
// the `FacingPointSets` from which explicitly unassigned points have been
// removed.
static McDArray<McVec2i>
handleEvidence(const EdgeVertexAttribute* evidenceAttrib,
               ma::FacingPointSets& pts, SpatialGraphSelection& refSelection,
               SpatialGraphSelection& transSelection) {
    mcrequire(pts.ref.positions.size() ==
              refSelection.getNumSelectedVertices());
    mcrequire(pts.trans.positions.size() ==
              transSelection.getNumSelectedVertices());

    McDArray<int> refNodeIds = getSelectedNodeIds(refSelection);
    McDArray<int> transSlicePointsToMatch = getSelectedNodeIds(transSelection);

    removeUnassigned(evidenceAttrib, refNodeIds, pts.ref.positions,
                     pts.ref.directions, refSelection, transSlicePointsToMatch,
                     pts.trans.positions, pts.trans.directions, transSelection);
    return getAssignments(evidenceAttrib, refNodeIds, transSlicePointsToMatch);
}

McDArray<int> remapIndices(const McDArray<int>& ids,
                           const McDArray<int>& refIds) {
    McDArray<int> r(ids.size());
    for (mclong i = 0; i < ids.size(); i++) {
        r[i] = refIds[ids[i]];
    }
    return r;
}

McDArray<McVec2i> remapIndices(const McDArray<McVec2i>& pairs,
                               const McDArray<int>& refIds,
                               const McDArray<int>& transIds) {
    McDArray<McVec2i> r(pairs.size());
    for (mclong h = 0; h < pairs.size(); h++) {
        const McVec2i e = pairs[h];
        r[h] = McVec2i(refIds[e.x], transIds[e.y]);
    }
    return r;
}

static void reportEvidenceDetails(const ma::MatchingPGM& res,
                                  const McDArray<McVec2i>& evidence,
                                  const McDArray<int>& refSlicePointsToMatch,
                                  const McDArray<int>& transSlicePointsToMatch,
                                  ma::Context* ctx) {
    mcrequire(ctx);
    for (int h = 0; h < res.queerEvidence.size(); h++) {
        const McVec2i e = res.queerEvidence[h];
        ctx->print(QString("WARNING! This evidence does not fit the "
                           "parameters: %1 %2 (%3 - %4)")
                       .arg(refSlicePointsToMatch[e.x])
                       .arg(transSlicePointsToMatch[e.y])
                       .arg(e.x)
                       .arg(e.y));
        ctx->print("Check the queer evidence label and remove evidence that is "
                   "queer or adjust your parameters.");
    }
    for (int i = 0; i < evidence.size(); i++) {
        ctx->print(QString("Evidence pair %1: %2 - %3")
                       .arg(i)
                       .arg(refSlicePointsToMatch[evidence[i].x])
                       .arg(transSlicePointsToMatch[evidence[i].y]));
    }
}

static McDArray<int> concat(McDArray<int> first, const McDArray<int>& second) {
    first.appendArray(second);
    return first;
}

ma::MatchingPGM
ma::matchingPGM(const FacingPointSets& pts, MatchingPGMParams mparams,
                HxSpatialGraph* sg, SpatialGraphSelection& refVertexSelection,
                SpatialGraphSelection& transVertexSelection,
                const SpatialGraphSelection& aroundMidPlaneSelection,
                Context* ctx) {
    if (!ctx) {
        ctx = &defaultContext();
    }

    clearAttributes(sg, aroundMidPlaneSelection);

    ma::FacingPointSets ptsWithoutUnassigned = pts;
    const EdgeVertexAttribute* evidenceAttrib =
        dynamic_cast<EdgeVertexAttribute*>(
            sg->findAttribute(HxSpatialGraph::VERTEX, "Evidence"));
    if (evidenceAttrib) {
        mparams.evidence =
            handleEvidence(evidenceAttrib, ptsWithoutUnassigned,
                           refVertexSelection, transVertexSelection);
    }
    mcassert(ptsWithoutUnassigned.ref.positions.size() ==
             refVertexSelection.getNumSelectedVertices());
    mcassert(ptsWithoutUnassigned.trans.positions.size() ==
             transVertexSelection.getNumSelectedVertices());

    const ma::MatchingPGM matching =
        matchingPGM(ptsWithoutUnassigned, mparams, ctx);

    const McDArray<int> refNodeIds = getSelectedNodeIds(refVertexSelection);
    const McDArray<int> transNodeIds = getSelectedNodeIds(transVertexSelection);
    reportEvidenceDetails(matching, mparams.evidence, refNodeIds, transNodeIds,
                          ctx);

    struct RemapWriter {
        HxSpatialGraph* sg;
        const McDArray<int>& refNodeIds;
        const McDArray<int>& transNodeIds;

        void writeFloats(const McDArray<float>& vals, const char* attr,
                         const char* attrZero) {
            writeFloatAttribute(sg, refNodeIds, vals, attr, attrZero);
        }

        void writePairs(const McDArray<McVec2i>& vals, const char* attr) {
            rewritePairs(sg, remapIndices(vals, refNodeIds, transNodeIds),
                         attr);
        }

        void writeBool(const McDArray<int>& vals, const char* attr) {
            rewriteBinaryLabel(sg, remapIndices(vals, refNodeIds), attr);
        }
    };

    RemapWriter w = { sg, refNodeIds, transNodeIds };
    w.writeFloats(matching.maxDiff, "maxDiff", "maxDiffIsZero");
    w.writeFloats(matching.entropie, "entropie", "entropieIsZero");
    w.writePairs(matching.queerEvidence, "QueerEvidence");
    w.writeBool(matching.ambiguities, "Ambiguities");
    w.writeBool(matching.ambiguitiesFromMarginalDiff,
                "AmbiguitiesFromMarginalDiff");
    w.writeBool(matching.assignToMakeTree, "AssignToMakeTree");
    w.writeBool(matching.criticalNodes, "CriticalNodes");
    w.writeBool(matching.criticalNodesWithAdaptiveThreshold,
                "CriticalNodesWithAdaptiveThreshold");
    w.writeBool(concat(matching.ambiguities, matching.assignToMakeTree),
                "AssignThesePointsStep");
    w.writeBool(matching.maxEntropieNodes, "MaxEntropieNodes");

    sg->touch(HxData::NEW_PARAMETERS);

    return matching;
}
