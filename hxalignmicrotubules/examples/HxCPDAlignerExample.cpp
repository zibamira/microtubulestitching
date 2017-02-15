// The module is not declared in an rc file.  You can create an instance from
// the TCL console with the following commands:
//
//     dso open libhxalignmicrotubules.so
//     create HxCPDAlignerExample

#include <hxcore/HxCompModule.h>
#include <hxcore/HxMessage.h>
#include <hxcore/HxObjectPool.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>

#include <hxalignmicrotubules/MovingLeastSquares.h>
#include <hxalignmicrotubules/mtalign.h>

namespace ma = mtalign;

class HxCPDAlignerExample : public HxCompModule {
    HX_HEADER(HxCPDAlignerExample);
  public:

    HxPortDoIt portAction;

    void compute();
};

HX_INIT_CLASS(HxCPDAlignerExample, HxCompModule);

HxCPDAlignerExample::HxCPDAlignerExample()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portAction(this, "action", tr("Action"), 1) {}

// `printMsg` prints to stdout and to the Amira console.
static void printMsg(QString msg) {
    printf("%s\n", qPrintable(msg));
    theMsg->printf(msg);
    theApp->handlePendingEvents();
}

HxCPDAlignerExample::~HxCPDAlignerExample()
{
}

void HxCPDAlignerExample::compute() {
    if (!portAction.wasHit()) {
        return;
    }

    HxSpatialGraph* in = mcinterface_cast<HxSpatialGraph>(portData.getSource());
    if (!in) {
        return;
    }

    // Duplicate input, because it is modified below.
    McHandle<HxSpatialGraph> stack = in->duplicate();

    // Configure an optional Context object to redirect printing to the Amira
    // console.  It will be passed to functions below to configure their
    // behavior.
    ma::Context ctx;
    ctx.print = &printMsg;

    // Determine slices based on attribute `TransformInfo`, and compute midplane
    // between reference slice and transformed slice.
    const int refSliceNum = 0;
    const int transSliceNum = 1;
    ma::SliceSelector selector(stack, "TransformInfo");
    const float midPlane = selector.computeMidPlane(refSliceNum, transSliceNum);

    // Configure projection params, and compute projected points.
    ma::EndPointParams prjParams;
    prjParams.refSliceNum = refSliceNum;
    prjParams.transSliceNum = transSliceNum;
    prjParams.projectionPlane = midPlane;
    prjParams.projectionType = ma::P_ORTHOGONAL;
    // Specify boundary region as 40% of slice thickness.
    prjParams.useAbsoluteValueForEndPointRegion = false;
    prjParams.endPointRegion = 40;
    // Second line point for direction must be closer than 2000 Angstrom.
    prjParams.maxDistForAngle = 2000;
    // Reject lines that are close to parallel to boundary.
    prjParams.angleToPlaneFilter = 0.01;
    ma::FacingPointSets points = ma::projectEndPoints(stack, prjParams);
    printMsg("============================================================");
    printMsg(QString("Projected %1 ref points and %2 trans points.")
                 .arg(points.ref.positions.size())
                 .arg(points.trans.positions.size()));
    printMsg("============================================================");

    // Compute linear CPD.  Set the params explicitly for illustration,
    // although the defaults could be used.
    ma::CPDParams cpdParams;
    ma::WarpResult deformation;
    cpdParams = ma::CPDParams::defaultsLinear();
    cpdParams.type = ma::CPD_LINEAR;
    cpdParams.linear.withScaling = true;
    cpdParams.linear.usePositions = true;
    cpdParams.linear.useDirections = true;
    cpdParams.linear.w = 0.1;
    cpdParams.linear.maxIterations = 200;
    cpdParams.linear.eDiffRelStop = 1.e-5;
    cpdParams.linear.sigmaSquareStop = 1.e-7;
    ma::cpd(points, deformation, cpdParams, &ctx);
    printMsg("============================================================");
    printMsg(
        QString("Linear CPD converged after %1 iterations (took %2 seconds).")
            .arg(deformation.alignInfo.numIterations)
            .arg(deformation.alignInfo.timeInSec));
    printMsg(QString("Relative expectation difference upon convergence: %1.")
                 .arg(deformation.alignInfo.eDiffRel));
    printMsg("============================================================");

    // Apply linear transform.  Apply the matrix to all selected vertices and
    // edges.  `getSlice()` does not select individual points.  The code below
    // assumes that the stack consists of a single pair of sections.  For a
    // larger stack, transformations would need to be consecutively applied
    // across all section pairs.
    {
        SpatialGraphSelection sel(stack);
        selector.getSlice(selector.getSliceAttributeValueFromIndex(1), sel);
        mcassert(sel.getNumSelectedPoints() == 0);
        SpatialGraphSelection::Iterator iter(sel);
        const McMat4f& mat = deformation.transformMatrix;
        int idx;
        while ((idx = iter.vertices.nextSelected()) != -1) {
            const McVec3f c = stack->getVertexCoords(idx);
            McVec3f res;
            mat.multVecMatrix(c, res);
            stack->setVertexCoords(idx, res);
        }
        while ((idx = iter.edges.nextSelected()) != -1) {
            McDArray<McVec3f> points = stack->getEdgePoints(idx);
            for (int p = 0; p < points.size(); p++) {
                McVec3f res;
                mat.multVecMatrix(points[p], res);
                points[p] = res;
            }
            stack->setEdgePoints(idx, points);
        }
        theObjectPool->addObject(stack->duplicate());
    }

    // Re-project points.
    points = ma::projectEndPoints(stack, prjParams);
    printMsg("============================================================");
    printMsg(QString("Projected %1 ref points and %2 trans points.")
                 .arg(points.ref.positions.size())
                 .arg(points.trans.positions.size()));
    printMsg("============================================================");

    // Compute elastic CPD.  Set the params explicitly for illustration,
    // although the defaults could be used.
    cpdParams = ma::CPDParams::defaultsElastic();
    cpdParams.type = ma::CPD_ELASTIC;
    cpdParams.elastic.beta = 10;
    cpdParams.elastic.lambda = 1;
    cpdParams.elastic.useDirections = true;
    cpdParams.elastic.w = 0.1;
    cpdParams.elastic.maxIterations = 200;
    cpdParams.elastic.eDiffRelStop = 1.e-5;
    cpdParams.elastic.sigmaSquareStop = 1.e-7;
    // Disable landmark resampling.
    cpdParams.elastic.sampleDistForWarpingLandmarks = 0;
    cpdParams.alphaForMLS = 2;
    ma::cpd(points, deformation, cpdParams, &ctx);
    printMsg("============================================================");
    printMsg(
        QString("Elastic CPD converged after %1 iterations (took %2 seconds).")
            .arg(deformation.alignInfo.numIterations)
            .arg(deformation.alignInfo.timeInSec));
    printMsg(QString("Relative expectation difference upon convergence: %1.")
                 .arg(deformation.alignInfo.eDiffRel));
    printMsg("============================================================");

    // Apply MLS transformation.
    MovingLeastSquares mls;
    mls.setAlpha(deformation.mlsParams.alpha);
    mls.setLandmarks(deformation.mlsParams.ps, deformation.mlsParams.qs);
    {
        SpatialGraphSelection sel(stack);
        selector.getSlice(selector.getSliceAttributeValueFromIndex(1), sel);
        mcassert(sel.getNumSelectedPoints() == 0);
        SpatialGraphSelection::Iterator iter(sel);
        int idx;
        while ((idx = iter.vertices.nextSelected()) != -1) {
            const McVec3f c = stack->getVertexCoords(idx);
            const McVec2d res = mls.interpolate(McVec2d(c.x, c.y));
            stack->setVertexCoords(idx, McVec3f(res.x, res.y, c.z));
        }
        while ((idx = iter.edges.nextSelected()) != -1) {
            McDArray<McVec3f> points = stack->getEdgePoints(idx);
            for (int p = 0; p < points.size(); p++) {
                const McVec3f c = points[p];
                const McVec2d res = mls.interpolate(McVec2d(c.x, c.y));
                points[p] = McVec3f(res.x, res.y, c.z);
            }
            stack->setEdgePoints(idx, points);
        }
        theObjectPool->addObject(stack->duplicate());
    }

    // Re-project points.
    points = ma::projectEndPoints(stack, prjParams);
    printMsg("============================================================");
    printMsg(QString("Projected %1 ref points and %2 trans points.")
                 .arg(points.ref.positions.size())
                 .arg(points.trans.positions.size()));
    printMsg("============================================================");

    // Compute PGM matching.
    ma::MatchingPGMParams pgmparams;

    // No user assignments.
    pgmparams.evidence.clear();

    // The params are inverted compared to the paper and specified in physical
    // distances.  The values are the ones that are used in the [Weber 2014].
    // They are explicitly set here for illustration, although the defaults
    // could be used.
    pgmparams = ma::MatchingPGMParams::defaultsWeber2014();
    pgmparams.weightConfig.weightType = ma::PGMPairWeightsParams::EXPONENTIAL;
    pgmparams.weightConfig.dist3dParam = 243;  // `lambda_c^-1`.
    pgmparams.weightConfig.distProjectedParam = 170;  // `lambda_p^-1`.
    pgmparams.weightConfig.angleWeightParam = 5.8;  // `lambda_alpha^-1`.
    pgmparams.pairFactorParam = 293;  // `lambda_s^1`.
    pgmparams.weightConfig.dummySignificance = 0.01;  // `r`.

    // `-log(r)` times params from above; values from GUI.
    pgmparams.weightConfig.distanceThreshold3d = 1119;
    pgmparams.weightConfig.distanceThresholdProjected = 783;
    pgmparams.weightConfig.angleThreshold = 26.7;

    // Activate all factors.
    pgmparams.weightConfig.useDist3dWeight = true;
    pgmparams.weightConfig.useProjectedDistWeight = true;
    pgmparams.weightConfig.useAngleWeight = true;
    pgmparams.weightConfig.useDistanceThreshold3d = true;
    pgmparams.weightConfig.useDistanceThresholdProjected = true;
    pgmparams.weightConfig.useAngleThreshold = true;

    // Use the low-level matching API.  The alternative is to use the function
    // `matchingPGM()` that derives evidence from a spatial graph and writes
    // back the matching result as spatial graph attributes.  See source for
    // details.
    const ma::MatchingPGM matching = matchingPGM(points, pgmparams, &ctx);
    printMsg("============================================================");
    printMsg(QString("%1 matched pairs, %2 ambiguities, %3 critical nodes.")
                 .arg(matching.matchedRefPointIds.size())
                 .arg(matching.ambiguities.size())
                 .arg(matching.criticalNodes.size()));
    printMsg("============================================================");

    printMsg("Done.");
}
