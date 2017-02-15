#include <hxalignmicrotubules/HxMovingLeastSquaresTomogramWarp.h>

#include <map>

#include <QFileInfo>
#include <QRegExp>

#include <hxcore/internal/HxWorkArea.h>
#include <hxfield/HxFieldEvaluator.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <mclib/McHandle.h>

#include <hxtemplatematchingutil/orientation3intcode.h>

#ifdef _OPENMP
#include <omp.h>
#endif

HX_INIT_CLASS(HxMovingLeastSquaresTomogramWarp, HxCompModule);

HxMovingLeastSquaresTomogramWarp::HxMovingLeastSquaresTomogramWarp()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portTomogram(this, "tomogram", tr("Tomogram"), HxUniformScalarField3::getClassTypeId()),
      portTransform(this, "transform", tr("Transform"), 1),
      portAction(this, "action", tr("Action")) {}

HxMovingLeastSquaresTomogramWarp::~HxMovingLeastSquaresTomogramWarp()
{
}

void HxMovingLeastSquaresTomogramWarp::update() {
    HxSpatialGraph* sg = hxconnection_cast<HxSpatialGraph>(portData);
    if (portData.isNew()) {
        portTransform.setNum(0, 0);
        extractParameters(sg);
        portTransform.setNum(0, transformations.size());
        for (int i = 0; i < int(transformations.size()); i++) {
            portTransform.setLabel(
                0, i, qPrintable(transformations[i].name + " (" +
                                 transformations[i].filename + ")"));
        }
    }
}

// Create result that covers spatial graph in xy and shifted input tomogram
// in z.  Use input tomogram's voxel size.
static McHandle<HxUniformScalarField3>
makeFieldCoveringXYAndTranslatedInZ(const HxUniformScalarField3* tomogram,
                                    const HxSpatialData* sg, float zshift) {
    const McDim3l& dimsTomo = tomogram->lattice().getDims();
    const McVec3f vs = tomogram->getVoxelSize();
    const McBox3f& bboxTomo = tomogram->getBoundingBox();

    const McBox3f& bboxSG = sg->getBoundingBox();

    int dims[3];
    for (int d = 0; d < 2; d++) {
        dims[d] = ceil((bboxSG[2 * d + 1] - bboxSG[2 * d]) / vs[d] + 1.0);
    }
    dims[2] = dimsTomo[2];

    float bbox[6];
    bbox[0] = bboxSG[0];
    bbox[2] = bboxSG[2];
    bbox[4] = bboxTomo[4] + zshift;
    for (int d = 0; d < 3; d++) {
        bbox[2 * d + 1] = bbox[2 * d] + vs[d] * (dims[d] - 1);
    }

    McHandle<HxUniformScalarField3> result(
        new HxUniformScalarField3(dims, tomogram->primType()));
    result->lattice().setBoundingBox(bbox);

    return result;
}

// Apply transformations for sliceNum; sample from tomogram and store in
// result.  Uses standard trilinear interpolation.
static void
applyMLS(HxUniformScalarField3* tomogram,
         std::map<int, HxMovingLeastSquaresTomogramWarp::TransformInfos>&
             transformations,
         int sliceNum, HxUniformScalarField3* result) {

    bool aborted = false;
    int numIteration = 0;
    theWorkArea->startWorking("Applying MLS...");

    // There are 3 coordinate systems:
    //
    //  - the original tomogram.
    //  - the transformed tomogram.
    //  - the result.
    //
    //  From original to transformed tomogram by infos.matrix.
    //  From transformed to original tomogram by infos.matrix.inverse().
    //  From result to transformed tomogram by applying mls as described below.

    const HxMovingLeastSquaresTomogramWarp::TransformInfos& infos =
        transformations[sliceNum];
    const SbMatrix transformedTomogramToTomogram = infos.matrix.inverse();
    const HxCoord3* coords = result->lattice().coords();
    const McDim3l& dims = result->lattice().getDims();

#pragma omp parallel
    {
        HxLocation3* locationIn = tomogram->createLocation();

#pragma omp for
        for (int j = 0; j < dims[1]; j++) {
            if (aborted)
                continue;

            for (int i = 0; i < dims[0]; i++) {
                McVec3f pos = coords->pos(McVec3i(i, j, 0));

                // MicrotubuleSpatialGraphAligner::applyMLS() applies
                // transforms T0, ..., Tn to create the transformed point t
                // from the original point p as follows:
                //
                //   t = T0(...(Tn(p))...)
                //
                // Here we need to do the inverse:
                //
                //   p = Tn^-1(...(T0^-1(t))...)
                //
                // Plug in t from above to see that p == p.
                //
                // The following applies the inverse transforms, which are
                // stored in mls, starting with the lowest slice number.
                McVec2d warpedCoord(pos[0], pos[1]);
                for (int s = 1; s <= sliceNum; s++) {
                    const HxMovingLeastSquaresTomogramWarp::TransformInfos&
                        info = transformations[s];
                    warpedCoord = info.mls.interpolate(warpedCoord);
                }

                for (int k = 0; k < dims[2]; k++) {
                    pos = coords->pos(McVec3i(i, j, k));
                    const SbVec3f posTransformedTomogram(warpedCoord.x,
                                                         warpedCoord.y, pos[2]);
                    SbVec3f posTomogram;
                    transformedTomogramToTomogram.multVecMatrix(
                        posTransformedTomogram, posTomogram);
                    // Fall back to 0 if position outside tomogram.
                    float resultValue = 0;
                    if (locationIn->set(posTomogram[0], posTomogram[1],
                                        posTomogram[2])) {
                        tomogram->eval(*locationIn, &resultValue);
                    }
                    result->set(i, j, k, resultValue);
                }

#pragma omp atomic
                numIteration++;
            }

#ifdef _OPENMP
            if (omp_get_thread_num() == 0)  // Main thread.
#endif
            {
                theWorkArea->setProgressValue((float)numIteration /
                                              (float)(dims[0] * dims[1]));
                if (theWorkArea->wasInterrupted()) {
                    aborted = true;
                }
            }
        }
        delete locationIn;
    }
    theWorkArea->stopWorking();
}

static double defaultDegRange(double phi) {
    while (phi < 0.0) {
        phi += 360.0;
    }
    while (phi > 360.0) {
        phi -= 360.0;
    }
    return phi;
}

double phiDeg(McVec2d p) {
    p.normalize();
    double phi = acos(p.x) * 180.0 / M_PI;
    if (p.y < 0.0) {
        phi = 360.0 - phi;
    }
    mcassert(defaultDegRange(phi) == phi);
    return phi;
}

using orientation3intcode::Orientation;

static Orientation rotatedByPhiDeg(Orientation o, float dPhiDeg) {
    o.phi = defaultDegRange(o.phi + dPhiDeg);
    return o;
}

// `applyMLSIntVecs()` correctly transforms int-encoded vectors.  See
// applyMLS() above for general idea.  Differences:
//
//  - `applyMLSIntVecs()` uses nearest neighbor interpolation and avoids cast
//    to float, because the encoded value must not be interpolated.
//
//  - `applyMLSIntVecs()` corrects the orientations by the rotation of the
//    transform.
//
static void
applyMLSIntVecs(HxUniformScalarField3* tomogram,
                std::map<int, HxMovingLeastSquaresTomogramWarp::TransformInfos>&
                    transformations,
                int sliceNum, HxUniformScalarField3* result) {
    mcrequire(tomogram->primType() == McPrimType::MC_UINT32);
    mcrequire(result->primType() == McPrimType::MC_UINT32);

    bool aborted = false;
    int numIteration = 0;
    theWorkArea->startWorking("Applying MLS...");

    const HxMovingLeastSquaresTomogramWarp::TransformInfos& infos =
        transformations[sliceNum];
    const SbMatrix transformedTomogramToTomogram = infos.matrix.inverse();

    // Angel by which orientations are rotated from the original tomogram to
    // the transformed tomogram.
    float dPhiTransform;
    {
        // origex is the result of transforming ex from the transformed to the
        // original tomogram.  The orientations need to be rotated in the
        // reverse direction, from the original to the transformed tomogram.
        // This is achieved by the negative of origex's phi.
        SbVec3f o(0.0, 0.0, 0.0);
        SbVec3f ex(1.0, 0.0, 0.0);
        SbVec3f origo;
        SbVec3f origstepex;
        transformedTomogramToTomogram.multVecMatrix(o, origo);
        transformedTomogramToTomogram.multVecMatrix(o + ex, origstepex);
        SbVec3f origex = origstepex - origo;
        dPhiTransform = defaultDegRange(-phiDeg(McVec2d(origex[0], origex[1])));
    }

    const HxCoord3* coords = result->lattice().coords();
    const McDim3l& dimsInt = result->lattice().getDims();
    const McDim3l& dims = result->lattice().getDims();

    mcuint32* resdat = static_cast<mcuint32*>(result->lattice().dataPtr());

#pragma omp parallel
    {
        HxFieldEvaluator* eval = HxEvalNN::createEval(tomogram);
        HxLocation3* loc = eval->createLocation();

#pragma omp for
        for (int j = 0; j < dimsInt[1]; j++) {
            if (aborted)
                continue;

            for (int i = 0; i < dimsInt[0]; i++) {
                McVec3f pos = coords->pos(McVec3i(i, j, 0));

                McVec2d warpedCoord(pos[0], pos[1]);
                const McVec2d ex(1.0, 0);
                McVec2d stepx = warpedCoord + ex;
                for (int s = 1; s <= sliceNum; s++) {
                    const HxMovingLeastSquaresTomogramWarp::TransformInfos&
                        info = transformations[s];
                    warpedCoord = info.mls.interpolate(warpedCoord);
                    stepx = info.mls.interpolate(stepx);
                }
                // See above for explanation.  Here, origex is the result of
                // transforming ex from the result to the transformed tomogram.
                // The orientations need to be rotated in the reverse direction.
                const McVec2d origex = stepx - warpedCoord;
                const float dPhiWarp = defaultDegRange(-phiDeg(origex));

                for (int k = 0; k < dimsInt[2]; k++) {
                    pos = coords->pos(McVec3i(i, j, k));
                    const SbVec3f posTransformedTomogram(warpedCoord.x,
                                                         warpedCoord.y, pos[2]);
                    SbVec3f posTomogram;
                    transformedTomogramToTomogram.multVecMatrix(
                        posTransformedTomogram, posTomogram);
                    // Fall back to 0 if position outside tomogram.
                    mcuint32 origOInt = 0;
                    if (loc->set(posTomogram[0], posTomogram[1],
                                 posTomogram[2])) {
                        eval->evalNative(loc, &origOInt);
                    }
                    const Orientation o = orientation3intcode::decode(origOInt);
                    const mculong idx = (k * dims[1] + j) * dims[0] + i;
                    resdat[idx] =
                        encode(rotatedByPhiDeg(o, dPhiTransform + dPhiWarp));
                }

#pragma omp atomic
                numIteration++;
            }

#ifdef _OPENMP
            if (omp_get_thread_num() == 0)  // Main thread.
#endif
            {
                theWorkArea->setProgressValue((float)numIteration /
                                              (float)(dims[0] * dims[1]));
                if (theWorkArea->wasInterrupted()) {
                    aborted = true;
                }
            }
        }
        delete loc;
        delete eval;
    }
    theWorkArea->stopWorking();
}

void HxMovingLeastSquaresTomogramWarp::compute() {
    if (!portAction.wasHit())
        return;

    HxSpatialGraph* sg = hxconnection_cast<HxSpatialGraph>(portData);
    HxUniformScalarField3* tomogram =
        hxconnection_cast<HxUniformScalarField3>(portTomogram);

    if (!sg || !tomogram || !transformations.size())
        return;

    const int sliceNum = portTransform.getValue(0);
    const TransformInfos& infos = transformations[sliceNum];
    SbVec3f translation;
    {
        SbRotation ignore1;
        SbVec3f ignore2;
        SbRotation ignore3;
        infos.matrix.getTransform(translation, ignore1, ignore2, ignore3);
    }

    McHandle<HxUniformScalarField3> result(
        makeFieldCoveringXYAndTranslatedInZ(tomogram, sg, translation[2]));
    result->composeLabel(tomogram->getLabel(), "warped");

    if (tomogram->primType() == McPrimType::MC_UINT32) {
        theMsg->printf("Assuming int-encoded vectors.");
        applyMLSIntVecs(tomogram, transformations, sliceNum, result);
    } else {
        applyMLS(tomogram, transformations, sliceNum, result);
    }

    setResult(result);
}

void
HxMovingLeastSquaresTomogramWarp::extractParameters(const HxSpatialGraph* sg) {
    transformations.clear();

    if (!sg)
        return;

    // Extract global slices transformations.
    const HxParamBundle* transformInfoPB =
        sg->parameters.getBundle("TransformInfo");
    if (!transformInfoPB) {
        theMsg->printf("Error: could not find TransformInfo parameter bundle");
        return;
    }

    for (int i = 0; i < transformInfoPB->getNumberOfBundles(); i++) {
        HxParamBundle* sliceInfo = transformInfoPB->getBundle(i);
        QString name = sliceInfo->getName();
        QRegExp rx("Slice[0-9]{4}");
        if (rx.exactMatch(name)) {
            const mclong index = name.mid(5).toInt() - 1;

            TransformInfos& info = transformations[index];
            const HxParameter* transfo = sliceInfo->find("Transform");
            if (transfo && transfo->getDimension() == 16) {
                info.name = name;

                double res[16];
                transfo->getReal(res);
                float* ptr = &info.matrix[0][0];
                for (int i = 0; i < 16; ++i) {
                    ptr[i] = float(res[i]);
                }
            }
            const HxParameter* filename = sliceInfo->find("Filename");
            if (filename) {
                info.filename = QFileInfo(filename->getString()).fileName();
            }
        }
    }

    // Sanity check - slice numbers should be consecutive
    for (int i = 0; i < int(transformations.size()); i++) {
        if (transformations.find(i) == transformations.end()) {
            theMsg->printf("Error: missing global transform info for slice %d",
                           i);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Extract MLS infos
    ////////////////////////////////////////////////////////////////////////////////////////////////
    const HxParamBundle* CPDTransformLandmarks =
        sg->parameters.getBundle("CPDTransformLandmarks");
    if (!CPDTransformLandmarks) {
        theMsg->printf(
            "Error: could not find CPDTransformLandmarks parameter bundle");
        return;
    }

    float alpha = 2;
    if (!CPDTransformLandmarks->findReal("alpha", alpha)) {
        theMsg->printf("Warning: missing alpha in CPDTransformLandmarks, using "
                       "default 2.");
    }

    for (int i = 0; i < CPDTransformLandmarks->getNumberOfBundles(); i++) {
        const HxParamBundle* sliceInfo = CPDTransformLandmarks->getBundle(i);
        const QString name = sliceInfo->getName();

        const QRegExp rx("Slices([0-9]+)-([0-9]+)");
        if (rx.exactMatch(name)) {
            const int fromSlice = rx.cap(1).toInt();
            const int toSlice = rx.cap(2).toInt();

            if (toSlice != fromSlice + 1) {
                theMsg->printf("Error: invalid '%s' in CPDTransformLandmarks, "
                               "slices need to be numbered consecutively",
                               qPrintable(name));
            }

            // Sanity check
            if (transformations.find(toSlice - 1) == transformations.end()) {
                theMsg->printf(
                    "Error: missing global transform info for slice %d",
                    toSlice - 1);
                continue;
            }

            TransformInfos info = transformations[toSlice - 1];

            info.mls.setAlpha(alpha);
            McDArray<McVec2d> ps;
            ps.resize(sliceInfo->getSize());
            McDArray<McVec2d> qs;
            qs.resize(sliceInfo->getSize());
            bool ok = true;
            for (int j = 0; j < sliceInfo->getSize(); j++) {
                HxParameter* PQ = static_cast<HxParameter*>((*sliceInfo)[j]);
                if (PQ->getDimension() != 4) {
                    theMsg->printf("Error: invalid number of entries in "
                                   "CPDTransformLandmarks %s, %s (expected 4).",
                                   qPrintable(name), qPrintable(PQ->getName()));
                    ok = false;
                    break;
                }
                double res[4];
                PQ->getReal(res);

                // !! We read PQ but store QP since we need the invert
                // transform. !!
                int index = QString(PQ->getName()).mid(2).toInt();
                qs[index].setValue(res[0], res[1]);
                ps[index].setValue(res[2], res[3]);
            }
            info.mls.setLandmarks(ps, qs);
            if (ok) {
                transformations[toSlice - 1] = info;
            }
        }
    }
}
