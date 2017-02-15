#include <hxalignmicrotubules/mtalign/project.h>

#include <limits>

#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <hxspatialgraph/internal/SpatialGraphSelection.h>
#include <mclib/internal/McDMatrix.h>
#include <mclib/McPlane.h>
#include <mclib/McLine.h>

#include <hxalignmicrotubules/mtalign/SliceSelector.h>
#include <hxalignmicrotubules/mtalign/data.h>

namespace ma = mtalign;

static void selectSliceByIdx(const ma::SliceSelector& sections,
                             const int sliceNum, SpatialGraphSelection& slice) {
    sections.getSlice(sections.getSliceAttributeValueFromIndex(sliceNum),
                      slice);
}

static void getZRange(const ma::SliceSelector& sections, const bool isReference,
                      float& minVal, float& maxVal,
                      const ma::EndPointParams& params) {

    SpatialGraphSelection referenceSelection;
    selectSliceByIdx(sections, params.refSliceNum, referenceSelection);
    const ma::MinMax refmm = sections.getZRange(referenceSelection);

    SpatialGraphSelection transSelection;
    selectSliceByIdx(sections, params.transSliceNum, transSelection);
    const ma::MinMax transmm = sections.getZRange(transSelection);

    float const refThickness =
        (!params.useAbsoluteValueForEndPointRegion)
            ? (params.endPointRegion / 100.0) * (refmm.max - refmm.min)
            : params.endPointRegion;
    float const transThickness =
        (!params.useAbsoluteValueForEndPointRegion)
            ? (params.endPointRegion / 100.0) * (transmm.max - transmm.min)
            : params.endPointRegion;

    if (isReference) {
        if (refmm.min < transmm.min) {
            minVal = refmm.max - refThickness;
            maxVal = refmm.max;
        } else {
            minVal = refmm.min;
            maxVal = refmm.min + refThickness;
        }
    } else {
        if (refmm.min < transmm.min) {
            minVal = transmm.min;
            maxVal = transmm.min + transThickness;
        } else {
            minVal = transmm.max - transThickness;
            maxVal = transmm.max;
        }
    }
}

namespace {

struct VertexWithAngle {
    int vertexId;
    float angle;
    VertexWithAngle(int vI, float a) : vertexId(vI), angle(a) {}
};

}  // namespace

static int leVertexAngleStruct(const VertexWithAngle& one,
                               const VertexWithAngle& other) {
    return one.angle > other.angle ? 1 : -1;
}

static void selectVertices(const HxSpatialGraph* sg,
                           const ma::SliceSelector& sections,
                           const bool isReference, const int maxNumberOfPoints,
                           SpatialGraphSelection& selectedVertices,
                           const ma::EndPointParams& params) {

    SpatialGraphSelection sliceSelection;
    selectSliceByIdx(sections,
                     isReference ? params.refSliceNum : params.transSliceNum,
                     sliceSelection);
    // A previous version contained code at this place whose name suggested
    // that it removed potential loops.  But it had no effect and has been
    // removed.
    float zMin, zMax;
    getZRange(sections, isReference, zMin, zMax, params);

    selectedVertices.resize(sg->getNumVertices(), sg->getNumEdges());
    selectedVertices.clear();

    SpatialGraphSelection::Iterator iter(sliceSelection);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    McDArray<VertexWithAngle> vertexAngleList;
    while (vNum != -1) {
        McVec3f v = sg->getVertexCoords(vNum);
        if (v.z >= zMin && v.z <= zMax) {
            selectedVertices.selectVertex(vNum);
        }
        vNum = iter.vertices.nextSelected();
    }

    // Remove non-endpoints and compute angle and so on.
    const SpatialGraphSelection currSelection = selectedVertices;
    SpatialGraphSelection::Iterator iter3(currSelection);
    iter3.vertices.reset();
    vNum = iter3.vertices.nextSelected();
    while (vNum != -1) {
        if (sg->getIncidentEdges(vNum).size() != 1) {
            selectedVertices.deselectVertex(vNum);
        } else {
            // Compute angle.
            int edgeNum = sg->getIncidentEdges(vNum)[0];
            McVec3f p1 = sg->getVertexCoords(sg->getEdgeSource(edgeNum));
            McVec3f p2 = sg->getVertexCoords(sg->getEdgeTarget(edgeNum));
            McVec3f zDir(0, 0, 1);
            McVec3f edgeDir = p1 - p2;
            float angle = zDir.angle(edgeDir);
            if (angle > M_PI * 0.5)
                angle = M_PI - angle;
            angle = angle * 180.0f / M_PI;
            vertexAngleList.append(VertexWithAngle(vNum, angle));
        }
        vNum = iter3.vertices.nextSelected();
    }

    vertexAngleList.sort(&leVertexAngleStruct);

    selectedVertices.clear();
    for (int i = 0; (i < maxNumberOfPoints) && (i < vertexAngleList.size());
         i++) {
        selectedVertices.selectVertex(vertexAngleList[i].vertexId);
    }
    theMsg->printf("Selected %d vertices",
                   selectedVertices.getNumSelectedVertices());
}

static McVec3f getNextPointAtDist(const McDArray<McVec3f>& edgePoints,
                                  const float maxDist) {
    float dist = 0;
    int actIndex = 1;

    while ((dist < maxDist) && (actIndex < edgePoints.size())) {
        dist += (edgePoints[actIndex - 1] - edgePoints[actIndex]).length();
        actIndex++;
    }
    if (dist <= maxDist) {
        // The max Dist was not reached at all.
        return edgePoints[actIndex - 1];
    } else {
        // Interpolate the last point.
        const McVec3f x2 = edgePoints[actIndex - 1];
        const McVec3f x1 = edgePoints[actIndex - 2];
        const float diffDist = dist - maxDist;
        const float distLastPoints = (x2 - x1).length();
        const float frac = diffDist / distLastPoints;
        mcassert(frac < 1);

        return (1 - frac) * x2 + frac * x1;
    }
}

static McDArray<McVec3f>
getDirections(const HxSpatialGraph* sg,
              const SpatialGraphSelection& selectedVertices,
              const float maxDist) {
    McDArray<McVec3f> directions;
    directions.resize(selectedVertices.getNumSelectedVertices());
    for (int i = 0; i < selectedVertices.getNumSelectedVertices(); i++) {
        int const selectedVertex = selectedVertices.getSelectedVertex(i);
        IncidenceList connectedEdges = sg->getIncidentEdges(selectedVertex);

        McVec3f const point1 = sg->getVertexCoords(selectedVertex);

        int const connectedEdge = connectedEdges[0];
        int const source = sg->getEdgeSource(connectedEdge);
        McDArray<McVec3f> edgePoints = sg->getEdgePoints(connectedEdge);
        if (source != selectedVertex) {
            edgePoints.reverse();
        }

        McVec3f const point2 = getNextPointAtDist(edgePoints, maxDist);
        McVec3f direction = point2 - point1;
        direction.normalize();
        directions[i] = direction;
    }
    return directions;
}

// `deselectEndpointsWithWrongDirection()` deselects lines that violate the
// model of the section boundary.  `directions` point from the boundary along
// the line; good directions point towards the section center.  Endpoints are
// rejected if the angle with the vector from the midplane to the section
// center gets too large.  An greater `angleToPlaneFilter` reduces the
// threshold so that fewer endpoints are accepted.
static void deselectEndpointsWithWrongDirection(
    const HxSpatialGraph* sg, SpatialGraphSelection& vertices,
    bool const isAboveMidPlane, const ma::EndPointParams& params) {
    const McDArray<McVec3f> directions =
        getDirections(sg, vertices, params.maxDistForAngle);
    const McVec3f midPlaneToSection =
        McVec3f(0.0, 0.0, isAboveMidPlane ? 1.0 : -1.0);
    const float threshold = (M_PI * 0.5 - M_PI * params.angleToPlaneFilter);
    for (int i = vertices.getNumSelectedVertices() - 1; i >= 0; i--) {
        if (directions[i].angle(midPlaneToSection) > threshold) {
            vertices.deselectVertex(vertices.getSelectedVertex(i));
        }
    }
}

static McDArray<McVec3f>
projectToPlaneFit(const HxSpatialGraph* sg,
                  const SpatialGraphSelection& vertices, const float planeZ,
                  const int degree) {

    theMsg->printf("using degree %d", degree);

    McDArray<McVec3f> result(0);

    SpatialGraphSelection::Iterator iter(vertices);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    while (vNum != -1) {
        int const edgeNum = sg->getIncidentEdges(vNum)[0];
        McDArray<McVec3f> edgePoints = sg->getEdgePoints(edgeNum);
        // Use only the last numPoints-fraction points for fitting.
        int fraction = int(0.5f * float(edgePoints.size()));
        if (edgePoints.size() - fraction < 2) {
            fraction = edgePoints.size() - 2;
        }
        if (sg->getEdgeSource(edgeNum) == vNum) {
            edgePoints.removeLast(fraction);
        } else {
            edgePoints.remove(0, fraction);
        }

        int const m = edgePoints.size();
        int matDim = degree + 1;
        if (m < matDim) {
            matDim = m;
            theMsg->printf("Warning: too few points (%d). Using polynomial of "
                           "degree %d instead.",
                           m, matDim);
        }
        McDMatrix<double> mat(2 * matDim, 2 * matDim);
        mat.fill(0.0);
        for (int r = 0; r < matDim; ++r) {
            for (int c = 0; c < matDim; ++c) {
                double sum = 0.0;
                for (int j = 0; j < m; ++j) {
                    sum += std::pow(edgePoints[j][2], c + r);
                }
                mat[r][c] = sum;
                mat[r + matDim][c + matDim] = sum;
            }
        }

        printf("Determinant: %f\n", mat.det());
        if (abs(mat.det()) < 0.001)
            mat.print();

        McDVector<double> xy(2 * matDim);
        for (int i = 0; i < matDim; ++i) {
            xy[i] = 0.0f;
            xy[i + matDim] = 0.0f;
            for (int j = 0; j < m; ++j) {
                xy[i] += std::pow(edgePoints[j][2], i) * edgePoints[j][0];
                xy[i + matDim] +=
                    std::pow(edgePoints[j][2], i) * edgePoints[j][1];
            }
        }

        McDMatrix<double> const invMat = mat.inverse();

        McDVector<double> ab(2 * matDim);
        invMat.multVec(xy.dataPtr(), ab.dataPtr());
        printf("ab\n");
        ab.print();

        McVec3f proj(0.0f, 0.0f, planeZ);
        for (int i = 0; i < matDim; ++i) {
            proj[0] += ab[i] * std::pow(planeZ, i);
            proj[1] += ab[i + matDim] * std::pow(planeZ, i);
        }

        bool ortho = false;
        if ((fabs(proj[0] - sg->getVertexCoords(vNum)[0]) > 1000.f) ||
            (fabs(proj[1] - sg->getVertexCoords(vNum)[1]) > 1000.f)) {
            theMsg->printf(
                "Warning: projection out of range. using ortho projection");
            proj[0] = sg->getVertexCoords(vNum)[0];
            proj[1] = sg->getVertexCoords(vNum)[1];
            ortho = true;
        }
        result.append(proj);

        // Add sampled line to debugging lineset.
        McDArray<int> pointIds;
        if (ortho) {
            McVec3f q = sg->getVertexCoords(vNum);

            q[2] = planeZ;

        } else {
            float minZ = std::numeric_limits<float>::infinity();
            float maxZ = -std::numeric_limits<float>::infinity();
            for (int i = 0; i < m; ++i) {
                if (edgePoints[i][2] < minZ)
                    minZ = edgePoints[i][2];
                if (edgePoints[i][2] > maxZ)
                    maxZ = edgePoints[i][2];
            }
            minZ = std::min(planeZ, minZ);
            maxZ = std::max(planeZ, maxZ);
            int const numSteps = 50;
            float hStep = (maxZ - minZ) / float(numSteps - 1);
            if (fabs(hStep) < 0.01) {
                theMsg->printf("Warning: small hstep");
                hStep = 0.01;
            }
            if (maxZ < minZ) {
                theMsg->printf("Warning: switch minZ and maxZ");
                float const tmp = maxZ;
                maxZ = minZ;
                minZ = tmp;
            }
            for (float h = minZ; h <= maxZ; h += hStep) {
                McVec3f proj(0.0f, 0.0f, h);
                for (int i = 0; i < matDim; ++i) {
                    proj[0] += ab[i] * std::pow(h, i);
                    proj[1] += ab[i + matDim] * std::pow(h, i);
                }
            }
        }

        vNum = iter.vertices.nextSelected();
    }
    theMsg->printf("Projected %ld points", result.size());
    return result;
}

static McDArray<McVec3f>
projectToPlaneTangent(const HxSpatialGraph* sg,
                      const SpatialGraphSelection& vertices,
                      const float planeZ) {
    McDArray<McVec3f> result(0);

    SpatialGraphSelection::Iterator iter(vertices);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    while (vNum != -1) {
        // Get the edge points of the edge connected to vertex vNum.
        int const edgeNum = sg->getIncidentEdges(vNum)[0];
        McDArray<McVec3f> edgePoints = sg->getEdgePoints(edgeNum);
        int const smoothNum = 5;
        McDArray<McVec3f> smoothed = edgePoints;
        for (int i = 0; i < smoothNum; ++i) {
            for (int j = 1; j < edgePoints.size() - 1; ++j) {
                McVec3f newPos(0.0f, 0.0f, 0.0f);
                for (int k = -1; k <= 1; ++k) {
                    newPos += edgePoints[j + k];
                }
                newPos /= 3.0f;
                smoothed[j] = newPos;
            }
            edgePoints = smoothed;
        }
        McVec3f const p1 = sg->getVertexCoords(vNum);
        McVec3f p2;
        if (vNum == sg->getEdgeSource(edgeNum)) {
            if (edgePoints.size() >= 50) {
                p2 = edgePoints[49];
            } else {
                p2 = edgePoints[1];
            }
        } else {
            if (edgePoints.size() >= 50) {
                p2 = edgePoints[edgePoints.size() - 48];
            } else {
                p2 = edgePoints[edgePoints.size() - 2];
            }
        }
        McVec3f dir = p2 - p1;

        // Line parallel to plane, use orthogonal projection.
        if (fabs(dir.z) < 0.001) {
            McVec3f p = sg->getVertexCoords(vNum);
            p.z = planeZ;
            result.append(p);

        }
        // Intersection of line through edge source and target with planeZ.
        else {
            float const t = (planeZ - p1.z) / dir.z;
            McVec3f p = dir * t + p1;
            result.append(p);
        }

        vNum = iter.vertices.nextSelected();
    }
    return result;
}

static McDArray<McVec3f>
projectToPlaneApproxDirection(const HxSpatialGraph* sg,
                              const SpatialGraphSelection& vertices,
                              const float planeZ, const float maxDist) {

    McDArray<McVec3f> result(0);

    const McDArray<McVec3f> directions = getDirections(sg, vertices, maxDist);

    SpatialGraphSelection::Iterator iter(vertices);
    iter.vertices.reset();
    int vNum = iter.vertices.nextSelected();
    int counter = 0;
    while (vNum != -1) {
        McPlane const theZPlane(McVec3f(0, 0, 1), planeZ);
        McVec3f const vertexCoord = sg->getVertexCoords(vNum);
        McVec3f dir = directions[counter];
        dir.normalize();
        float const angle = dir.dot(McVec3f(0, 0, 1));

        McLine theLine(vertexCoord, vertexCoord + directions[counter]);

        McVec3f intersectionPoint;
        bool const intersected =
            theZPlane.intersect(theLine, intersectionPoint);
        if (intersected && (fabs(angle) > 0.1)) {
            result.append(intersectionPoint);
        } else {
            result.append(McVec3f(vertexCoord.x, vertexCoord.y, planeZ));
        }

        vNum = iter.vertices.nextSelected();
        counter++;
    }
    return result;
}

static McDArray<McVec3f> projectToPlane(const HxSpatialGraph* sg,
                                        const SpatialGraphSelection& vertices,
                                        const float planeZ,
                                        const ma::EndPointParams& params) {
    McDArray<McVec3f> result(0);
    // Orthogonal projection: simply set the z value to planeZ.
    if (params.projectionType == ma::P_ORTHOGONAL) {
        SpatialGraphSelection::Iterator iter(vertices);
        iter.vertices.reset();
        int vNum = iter.vertices.nextSelected();
        while (vNum != -1) {
            result.append(sg->getVertexCoords(vNum));
            result.last().z = planeZ;
            vNum = iter.vertices.nextSelected();
        }
    }

    // Linear projection.
    else if (params.projectionType == ma::P_LINEAR) {
        SpatialGraphSelection::Iterator iter(vertices);
        iter.vertices.reset();
        int vNum = iter.vertices.nextSelected();
        while (vNum != -1) {
            // Get the edge points of the edge connected to vertex vNum.
            int const edgeNum = sg->getIncidentEdges(vNum)[0];
            McVec3f const p1 = sg->getVertexCoords(sg->getEdgeSource(edgeNum));
            McVec3f const p2 = sg->getVertexCoords(sg->getEdgeTarget(edgeNum));
            McVec3f const dir = p2 - p1;

            // Line parallel to plane, use orthogonal projection.
            if (fabs(dir.z) < 0.001) {
                McVec3f p = sg->getVertexCoords(vNum);
                p.z = planeZ;
                result.append(p);

            }
            // Intersection of line through edge source and target with planeZ.
            else {
                float const t = (planeZ - p1.z) / dir.z;
                McVec3f const p = dir * t + p1;
                result.append(p);
            }

            vNum = iter.vertices.nextSelected();
        }
    } else if (params.projectionType == ma::P_TANGENT) {
        result = projectToPlaneTangent(sg, vertices, planeZ);
    } else if (params.projectionType == ma::P_FIT_0) {
        result = projectToPlaneFit(sg, vertices, planeZ, 0);
    } else if (params.projectionType == ma::P_FIT_1) {
        result = projectToPlaneFit(sg, vertices, planeZ, 1);
    } else if (params.projectionType == ma::P_FIT_2) {
        result = projectToPlaneFit(sg, vertices, planeZ, 2);
    } else if (params.projectionType == ma::P_FIT_3) {
        result = projectToPlaneFit(sg, vertices, planeZ, 3);
    } else if (params.projectionType == ma::P_APPROX_TANGENT) {
        result = projectToPlaneApproxDirection(sg, vertices, planeZ,
                                               params.maxDistForAngle);
    }
    // No projection.
    else if (params.projectionType == ma::P_NONE) {
        SpatialGraphSelection::Iterator iter(vertices);
        iter.vertices.reset();
        int vNum = iter.vertices.nextSelected();
        while (vNum != -1) {
            result.append(sg->getVertexCoords(vNum));
            vNum = iter.vertices.nextSelected();
        }
    } else {
        theMsg->printf("Error: not existing ProjectionType");
    }
    return result;
}

static bool isSliceAboveMidPlane(const ma::SliceSelector& sections,
                                 const SpatialGraphSelection& sliceSelection1,
                                 const SpatialGraphSelection& sliceSelection2,
                                 const int thisSlice) {
    const ma::MinMax mm1 = sections.getZRange(sliceSelection1);
    const ma::MinMax mm2 = sections.getZRange(sliceSelection2);
    if (thisSlice == 0)
        return mm1.min > mm2.min;
    else
        return mm2.min > mm1.min;
}

static ma::FacingPointSets
projectNumEndPoints(const HxSpatialGraph* sg, const ma::EndPointParams& params,
                    const int nPoints, SpatialGraphSelection& refSel,
                    SpatialGraphSelection& transSel) {
    ma::SliceSelector sections(sg, "TransformInfo");

    mcrequire(params.refSliceNum >= 0);
    mcrequire(params.refSliceNum < sections.getNumSlices());
    mcrequire(params.transSliceNum >= 0);
    mcrequire(params.transSliceNum < sections.getNumSlices());
    // endPointRegion must be an absolute value or %.
    mcrequire(params.useAbsoluteValueForEndPointRegion ||
              params.endPointRegion >= 0);
    mcrequire(params.useAbsoluteValueForEndPointRegion ||
              params.endPointRegion <= 100);
    // angleToPlaneFilter is proportion of 180 degrees.  90 degrees is the
    // largest reasonable filter value.
    mcrequire(params.angleToPlaneFilter <= 0.5);

    const bool isReference = true;
    selectVertices(sg, sections, isReference, nPoints, refSel, params);
    selectVertices(sg, sections, !isReference, nPoints, transSel, params);

    bool const isAbove = isSliceAboveMidPlane(sections, refSel, transSel, 0);
    deselectEndpointsWithWrongDirection(sg, refSel, isAbove, params);
    deselectEndpointsWithWrongDirection(sg, transSel, !isAbove, params);

    ma::FacingPointSets f;
    f.ref.positions =
        projectToPlane(sg, refSel, params.projectionPlane, params);
    f.ref.directions = getDirections(sg, refSel, params.maxDistForAngle);
    f.trans.positions =
        projectToPlane(sg, transSel, params.projectionPlane, params);
    f.trans.directions = getDirections(sg, transSel, params.maxDistForAngle);
    return f;
}

ma::FacingPointSets ma::projectEndPoints(const HxSpatialGraph* sg,
                                         const ma::EndPointParams& params) {
    SpatialGraphSelection refSel;
    SpatialGraphSelection transSel;
    return projectEndPoints(sg, refSel, transSel, params);
}

ma::FacingPointSets ma::projectEndPoints(const HxSpatialGraph* sg,
                                         SpatialGraphSelection& refSel,
                                         SpatialGraphSelection& transSel,
                                         const ma::EndPointParams& params) {
    const int nPoints = sg->getNumVertices();
    return projectNumEndPoints(sg, params, nPoints, refSel, transSel);
}

ma::FacingPointSets ma::projectEndPointsSubset(
    const HxSpatialGraph* sg, SpatialGraphSelection& refSel,
    SpatialGraphSelection& transSel, const ma::EndPointParams& params) {
    const int nPoints = params.numMaxPointsForInitTransform;
    return projectNumEndPoints(sg, params, nPoints, refSel, transSel);
}
