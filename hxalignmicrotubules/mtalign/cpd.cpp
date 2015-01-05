#include <hxalignmicrotubules/mtalign/cpd.h>

#include <iostream>

#include <mclib/McException.h>

#include <hxalignmicrotubules/CoherentPointDriftNLFisherMises.h>
#include <hxalignmicrotubules/CoherentPointDriftRigidFisherMises.h>
#include <hxalignmicrotubules/mtalign/data.h>

namespace ma = mtalign;

static void resamplePairs(McDArray<McVec3f>& p1, McDArray<McVec3f>& p2,
                          float sampleDist) {
    int numPairs = p1.size();
    for (int i = p1.size() - 1; i > -1; i--) {
        bool resetI = false;
        for (int j = i - 1; j > -1; j--) {
            McVec3f set1Coord = p1[i];
            McVec3f set2Coord = p1[j];
            float dist = (set1Coord - set2Coord).length();
            if (dist < sampleDist) {
                p1.remove(j, 1);
                p2.remove(j, 1);
                resetI = true;
            }
        }
        if (resetI)
            i = p1.size() - 1;
    }
    std::cout << "\n" << p1.size() << " point from " << numPairs
              << " left after resampling.";
}

void ma::cpdElastic(const ma::FacingPointSets& points,
                    const ma::CPDParams& params, ma::WarpResult& warpResult) {

    // Configure cpd params.
    CoherentPointDriftNLFisherMises cpd;
    cpd.params = params.elastic;

    // Set fixed points.
    McDArray<McVec3f> coords = points.ref.positions;
    McDArray<McVec3f> directions = points.ref.directions;
    if (params.elastic.useDirections) {
        for (int i = 0; i < directions.size(); i++) {
            directions[i].normalize();
        }
    } else {
        directions.fill(McVec3f(0.0, 0.0, 1.0));
    }
    cpd.convertCoordsToMatrix(coords, cpd.xs);
    cpd.convertDirectionsToMatrix(directions, cpd.xDirs);

    // Set moving points.
    coords = points.trans.positions;
    directions = points.trans.directions;
    if (params.elastic.useDirections) {
        for (int i = 0; i < directions.size(); i++) {
            directions[i].normalize();
            directions[i] *= -1;
        }
    } else {
        directions.fill(McVec3f(0.0, 0.0, 1.0));
    }
    cpd.convertCoordsToMatrix(coords, cpd.ys);
    cpd.convertDirectionsToMatrix(directions, cpd.yDirs);

    // Solve.
    McDMatrix<double> G;
    McDMatrix<double> W;
    McDArray<McVec2i> correspondences;
    const AlignInfo info = cpd.align(G, W, correspondences);
    McDMatrix<double> transCoordsShiftedM;
    cpd.shiftYs(cpd.ys, G, W, transCoordsShiftedM);
    cpd.rescaleYs(transCoordsShiftedM);

    // Get result.
    const McDArray<McVec3f> transCoordsOld = points.trans.positions;
    McDArray<McVec3f> transCoords;
    cpd.convertMatrixToCoords(transCoords, transCoordsShiftedM);
    for (int i = 0; i < transCoords.size(); i++) {
        transCoords[i].z = transCoordsOld[i].z;
    }
    McDArray<McVec3f> origCoords = points.trans.positions;
    resamplePairs(origCoords, transCoords,
                  params.elastic.sampleDistForWarpingLandmarks);
    warpResult.mls.setAlpha(params.alphaForMLS);
    warpResult.mls.setLandmarks(origCoords, transCoords);
    warpResult.alignInfo = info;
}

void ma::cpdRigidFisherMises(const ma::FacingPointSets& points,
                             const ma::CPDParams& params,
                             ma::WarpResult& warpResult) {
    McDMatrix<double> refCoordsM, transCoordsM;

    CoherentPointDriftRigidFisherMises cpd;
    cpd.params = params.rigid;

    McDArray<McVec3f> directions = points.ref.directions;
    for (int i = 0; i < directions.size(); i++) {
        directions[i].normalize();
    }
    cpd.convertMcDArraysToMatrix(points.ref.positions, directions, cpd.Xc,
                                 cpd.Xd);

    directions = points.trans.directions;
    for (int i = 0; i < directions.size(); i++) {
        directions[i] *= -1.0;
        directions[i].normalize();
    }
    cpd.convertMcDArraysToMatrix(points.trans.positions, directions, cpd.Yc,
                                 cpd.Yd);

    McDMatrix<double> R, Rd;
    McDVector<double> t;
    double s;

    McDArray<McVec2i> correspondences;
    const AlignInfo info = cpd.align(R, s, t, Rd, correspondences);

    // construct 4d matrix
    warpResult.transformMatrix.makeIdentity();
    t *= cpd.mMeansAndStds.stdC;
    t[0] += cpd.mMeansAndStds.meanC[0];
    t[1] += cpd.mMeansAndStds.meanC[1];
    McDVector<double> muRotScaled(2), mu(2);
    mu[0] = cpd.mMeansAndStds.meanC[0];
    mu[1] = cpd.mMeansAndStds.meanC[1];
    R.multVec(mu, muRotScaled);
    muRotScaled *= s;
    t -= muRotScaled;
    R.transpose();
    warpResult.transformMatrix[0][0] = R[0][0] * s;
    warpResult.transformMatrix[0][1] = R[0][1] * s;
    warpResult.transformMatrix[1][0] = R[1][0] * s;
    warpResult.transformMatrix[1][1] = R[1][1] * s;
    warpResult.transformMatrix[3][0] = t[0];
    warpResult.transformMatrix[3][1] = t[1];
    warpResult.alignInfo = info;
}

void ma::cpd(const ma::FacingPointSets& points, const ma::CPDParams& params,
             ma::WarpResult& warpResult) {
    if (params.type == CPD_RIGID_FISHER_MISES) {
        cpdRigidFisherMises(points, params, warpResult);
    } else if (params.type == CPD_ELASTIC) {
        cpdElastic(points, params, warpResult);
    } else {
        mcthrow("CPDType not implemented!");
    }
}
