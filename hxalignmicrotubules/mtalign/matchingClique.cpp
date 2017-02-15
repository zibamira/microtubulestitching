#include <hxalignmicrotubules/mtalign/matchingClique.h>

#include <QTime>
#include <limits>
#include <algorithm>

#include <hxcore/HxMessage.h>
#include <mclib/internal/McAssert.h>
#include <mclib/McBitfield.h>
#include <mclib/McDArray.h>
#include <mclib/McMath.h>
#include <mclib/McRot.h>
#include <mclib/McVec2i.h>
#include <mclib/McVec3.h>

#include <pointmatching/StartTransformationGenerator3d.h>
#include <pointmatching/Transformation.h>
#include <hxalignmicrotubules/mtalign/NullPointRepresentation.h>

namespace ma = mtalign;

McDArray<McDArray<double> > computeDistMatrix(const McDArray<McVec3f>& points) {
    McDArray<McDArray<double> > dist;

    const int n = points.size();
    dist.resize(n);
    for (int i = 0; i < n; i++) {
        dist[i].resize(n);
    }

    for (int i = 0; i < n; i++) {
        dist[i][i] = 0.0f;
        for (int j = i + 1; j < n; j++) {
            const double d = (points[i] - points[j]).length();
            dist[i][j] = d;
            dist[j][i] = d;
        }
    }

    return dist;
}

namespace {

class Repr : public ma::NullPointRepresentation {
  public:
    Repr(const McDArray<McVec3f>& positions,
         const McDArray<McVec3f>& directions, const double cliqueDistThreshold,
         const double maxAngleDistForInitMatching)
        : mCoords(positions), mDirections(directions) {
        mCliqueDistThreshold = cliqueDistThreshold;
        mMaxAngleDistForInitMatching = maxAngleDistForInitMatching;
        distMatrix = computeDistMatrix(positions);
    }

  private:
    virtual int getNumPoints() const;

    const McDArray<McVec3f>& getDirections() const;

    virtual const McDArray<McVec3f>& getCoords() const;

    virtual void
    computeCorrespondenceGraphVertices(const PointRepresentation* pointSet2,
                                       McDArray<McVec2i>& corrGraph) const;

    virtual void
    computeCorrespondenceGraphEdges(const PointRepresentation* pointSet2,
                                    const McDArray<McVec2i>& corrGraph,
                                    McDArray<McBitfield>& connected) const;

  private:
    const McDArray<McDArray<double> >& getDistMatrix() const;

    bool shouldConnectByEdge(const int vertex1, const int vertex2,
                             const PointRepresentation* otherPoints,
                             const McDArray<McVec2i>& corrGraph) const;

    McDArray<McVec3f> mCoords;
    McDArray<McVec3f> mDirections;
    McDArray<McDArray<double> > distMatrix;
    double mCliqueDistThreshold;
    double mMaxAngleDistForInitMatching;
};
}

int Repr::getNumPoints() const { return mCoords.size(); }

const McDArray<McVec3f>& Repr::getCoords() const { return mCoords; }

const McDArray<McVec3f>& Repr::getDirections() const { return mDirections; }

const McDArray<McDArray<double> >& Repr::getDistMatrix() const {
    return distMatrix;
}

// All points can be matched.
void
Repr::computeCorrespondenceGraphVertices(const PointRepresentation* pointSet2,
                                         McDArray<McVec2i>& corrGraph) const {
    const int numPoints1 = getNumPoints();
    const int numPoints2 = pointSet2->getNumPoints();
    corrGraph.remax(numPoints1 * numPoints2);
    for (int i = 0; i < numPoints1; i++) {
        for (int j = 0; j < numPoints2; j++) {
            corrGraph.append(McVec2i(i, j));
        }
    }
}

void
Repr::computeCorrespondenceGraphEdges(const PointRepresentation* pointSet2,
                                      const McDArray<McVec2i>& corrGraph,
                                      McDArray<McBitfield>& connected) const {
    connected.resize(corrGraph.size());

    for (int i = 0; i < connected.size(); i++) {
        connected[i].resize(corrGraph.size());
        connected[i].setAll();
    }

    for (int i = 0; i < corrGraph.size(); i++) {
        for (int j = i + 1; j < corrGraph.size(); j++) {
            if (!shouldConnectByEdge(i, j, pointSet2, corrGraph)) {
                connected[i].unset(j);
                connected[j].unset(i);
            }
        }
    }
}

bool Repr::shouldConnectByEdge(const int vertex1, const int vertex2,
                               const PointRepresentation* otherPoints,
                               const McDArray<McVec2i>& corrGraph) const {
    const McDArray<McDArray<double> >& distMatrix1 = getDistMatrix();
    const McDArray<McDArray<double> >& distMatrix2 =
        otherPoints->getDistMatrix();

    const float d1 = distMatrix1[corrGraph[vertex1][0]][corrGraph[vertex2][0]];
    const float d2 = distMatrix2[corrGraph[vertex1][1]][corrGraph[vertex2][1]];
    const bool distanceThresholdOK = !(fabs(d1 - d2) > mCliqueDistThreshold);
    const Repr* const otherPointsA = dynamic_cast<const Repr*>(otherPoints);
    const McDArray<McVec3f>& directions1 = getDirections();
    const McDArray<McVec3f>& directions2 = otherPointsA->getDirections();

    double angleDiff1 = directions1[corrGraph[vertex1][0]].dot(
        directions1[corrGraph[vertex2][0]]);
    double angleDiff2 = directions2[corrGraph[vertex1][1]].dot(
        directions2[corrGraph[vertex2][1]]);
    angleDiff1 = acos(angleDiff1);
    angleDiff2 = acos(angleDiff2);

    double angleDiffDiff = fabs(angleDiff1 - angleDiff2);
    angleDiffDiff = angleDiffDiff / (2.0 * M_PI) * 360.0;

    const bool directionDifferenceOK =
        angleDiffDiff < mMaxAngleDistForInitMatching;

    return directionDifferenceOK && distanceThresholdOK;
}

// Remove any z shift and return transforms as McMat4f.
McDArray<McMat4f> asMats(const McDArray<Transformation>& ts) {
    McDArray<McMat4f> mats;
    for (int i = 0; i < ts.size(); i++) {
        McRotation rotation, so;
        McVec3f translation, scale;
        ts[i].getTransformation3d().getTransform(translation, rotation, scale,
                                                 so);
        translation.z = 0.0;
        McMat4f mat = McMat4f::IDENTITY;
        mat.setTransform(translation, rotation, scale);
        mats.append(mat);
    }
    return mats;
}

McDArray<McMat4f>
ma::matchingCliqueTransforms(const ma::FacingPointSets& pts,
                             const ma::MatchingParams& params) {
    mcassert(params.transformType == TF_RIGID ||
             params.transformType == TF_RIGID_ISO_SCALE);

    McDArray<int> cliqueSizes;
    StartTransformationGenerator3d startTransfGen;
    startTransfGen.setMaxNumStartTransformations(params.maxNumCliques);

    // `StartTransformationGenerator3d` expects transform type as integer in
    // the same way as `mcAlignPointSets()`:
    const int rigid = 0;
    const int rigidIsoScale = 1;
    startTransfGen.setTransformType(
        (params.transformType == TF_RIGID) ? rigid : rigidIsoScale);

    // The minimum clique size is determined by the smaller set.
    int mcs =
        MC_MIN2(int(pts.ref.positions.size() * params.minCliqueSizeFraction),
                int(pts.trans.positions.size() * params.minCliqueSizeFraction));

    mcs = std::max(mcs, 3);
    printf("Minimum clique size: %d\n", mcs);
    startTransfGen.setMinimumCliqueSize(mcs);

    Repr refRepresentation(pts.ref.positions, pts.ref.directions,
                           params.maxDistanceForGraphConstruction,
                           params.maxAngleDiffForInitMatching);

    Repr transRepresentation(pts.trans.positions, pts.trans.directions,
                             params.maxDistanceForGraphConstruction,
                             params.maxAngleDiffForInitMatching);

    McDArray<Transformation> transforms;
    if (!startTransfGen.compute(&refRepresentation, &transRepresentation,
                                transforms, cliqueSizes)) {
        printf("MaxNumStartTransformations (%d) reached\n",
               (int)transforms.size());
    }

    printf("Number of start transformations: %ld\n", transforms.size());

    return asMats(transforms);
}
