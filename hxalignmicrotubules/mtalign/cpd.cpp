#include <hxalignmicrotubules/mtalign/cpd.h>

#include <iostream>

#include <mclib/McException.h>

#include <hxalignmicrotubules/mtalign/CPDElasticAligner.h>
#include <hxalignmicrotubules/mtalign/CPDLinearAligner.h>
#include <hxalignmicrotubules/mtalign/data.h>

namespace ma = mtalign;

// `resamplePairs()` removes points that are closer than `sampleDist` in `p1`.
// The corresponding points are removed from `p2`.  The two arrays are required
// to have the same size.
static void resamplePairs(McDArray<McVec3f>& p1, McDArray<McVec3f>& p2,
                          const float sampleDist) {
    mcrequire(p1.size() == p2.size());
    const int numPairs = p1.size();
    for (int i = p1.size() - 1; i > -1; i--) {
        bool resetI = false;
        for (int j = i - 1; j > -1; j--) {
            const McVec3f set1Coord = p1[i];
            const McVec3f set2Coord = p1[j];
            const float dist = (set1Coord - set2Coord).length();
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

static McDArray<McVec2d> asVec2dArray(const McDArray<McVec3f>& a) {
    McDArray<McVec2d> b;
    b.resize(a.size());
    for (long i = 0; i < a.size(); i++) {
        b[i].x = a[i].x;
        b[i].y = a[i].y;
    }
    return b;
}

void ma::cpdElastic(const ma::FacingPointSets& points,
                    const ma::CPDParams& params, ma::WarpResult& warpResult) {

    // Configure cpd.
    CPDElasticAligner cpd;
    cpd.params = params.elastic;

    if (params.elastic.useDirections) {
        FacingPointSets copy = points;
        for (int i = 0; i < copy.trans.directions.size(); i++) {
            copy.trans.directions[i] *= -1;
        }
        cpd.setPoints(copy);
    } else {
        cpd.setPoints(points);
    }

    // Solve.
    AlignInfo info;
    McDArray<McVec3f> transCoords = cpd.align(info);

    McDArray<McVec3f> origCoords = points.trans.positions;
    resamplePairs(origCoords, transCoords,
                  params.elastic.sampleDistForWarpingLandmarks);
    warpResult.mlsParams.alpha = params.alphaForMLS;
    warpResult.mlsParams.ps = asVec2dArray(origCoords);
    warpResult.mlsParams.qs = asVec2dArray(transCoords);
    warpResult.alignInfo = info;
}

void ma::cpdLinear(const ma::FacingPointSets& points,
                   const ma::CPDParams& params, ma::WarpResult& warpResult) {
    McDMatrix<double> refCoordsM, transCoordsM;

    CPDLinearAligner cpd;
    cpd.params = params.linear;

    if (params.linear.useDirections) {
        FacingPointSets copy = points;
        for (int i = 0; i < copy.trans.directions.size(); i++) {
            copy.trans.directions[i] *= -1;
        }
        cpd.setPoints(copy);
    } else {
        cpd.setPoints(points);
    }

    McDMatrix<double> R, Rd;
    McDVector<double> t;
    double s;
    const AlignInfo info = cpd.align(R, s, t, Rd);
    warpResult.transformMatrix = cpd.getTransformMat4f(R, s, t);
    warpResult.alignInfo = info;
}

void ma::cpd(const ma::FacingPointSets& points, const ma::CPDParams& params,
             ma::WarpResult& warpResult) {
    if (params.type == CPD_LINEAR) {
        cpdLinear(points, params, warpResult);
    } else if (params.type == CPD_ELASTIC) {
        cpdElastic(points, params, warpResult);
    } else {
        mcthrow("CPDType not implemented!");
    }
}
