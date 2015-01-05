#include <hxalignmicrotubules/HxMovingLeastSquaresWarp.h>

#include <hxcore/HxMessage.h>
#include <hxfield/HxUniformLabelField3.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformVectorField3.h>
#include <hxlandmark/HxLandmarkSet.h>

#ifdef _OPENMP
#include <omp.h>
#include <hxcore/HxSettingsMgr.h>
#endif

#include <hxalignmicrotubules/MovingLeastSquares.h>

static mculong latticePos(int i, int j, int k, const int* dims) {
    return (mculong)dims[0] * (mculong)dims[1] * (mculong)k +
           (mculong)j * (mculong)dims[0] + (mculong)i;
}

HX_INIT_CLASS(HxMovingLeastSquaresWarp, HxCompModule);

HxMovingLeastSquaresWarp::HxMovingLeastSquaresWarp()
    : HxCompModule(HxLandmarkSet::getClassTypeId()),
      portFromImage(this, "from", HxUniformScalarField3::getClassTypeId()),
      portToImage(this, "to", HxUniformScalarField3::getClassTypeId()),
      portMethod(this, "method", 1),
      portAlpha(this, "parameter"),
      portAction(this, "action") {
    portAction.setAliasName("doIt");
    portMethod.setLabel(0, "Rigid");
    portAlpha.setLabel("Alpha:");
    portAlpha.setMinMax(0, 10);
    portAlpha.setValue(4);
}

HxMovingLeastSquaresWarp::~HxMovingLeastSquaresWarp() {}

void HxMovingLeastSquaresWarp::update() {}

HxUniformVectorField3* HxMovingLeastSquaresWarp::createOutputVectorDataSet() {
    HxUniformScalarField3* fromImage =
        dynamic_cast<HxUniformScalarField3*>(portFromImage.source());
    HxUniformVectorField3* output =
        dynamic_cast<HxUniformVectorField3*>(getResult(1));

    if (!fromImage)
        return (0);

    int dims[3];
    float bbox[6];

    if (fromImage) {
        memcpy(dims, fromImage->lattice.dimsInt(), 3 * sizeof(int));
        fromImage->getBoundingBox(bbox);
    }
    if (!output || output->lattice.dimsInt()[0] != dims[0] ||
        output->lattice.dimsInt()[1] != dims[1] ||
        output->lattice.dimsInt()[2] != dims[2])
        output = 0;

    if (!output) {
        output = new HxUniformVectorField3(dims, McPrimType::mc_float);
    }
    output->lattice.setBoundingBox(bbox);
    output->composeLabel(fromImage->getLabel().getString(), "Displacement");
    setResult(1, output);
    return (output);
}

HxUniformScalarField3* HxMovingLeastSquaresWarp::createOutputDataSet() {
    HxUniformScalarField3* fromImage =
        dynamic_cast<HxUniformScalarField3*>(portFromImage.source());
    HxUniformScalarField3* toImage =
        dynamic_cast<HxUniformScalarField3*>(portToImage.source());
    HxUniformScalarField3* warpedImage =
        dynamic_cast<HxUniformScalarField3*>(getResult(0));

    if (!fromImage)
        return (0);

    int dims[3];
    float bbox[6];

    if (toImage) {
        memcpy(dims, toImage->lattice.dimsInt(), 3 * sizeof(int));
        toImage->getBoundingBox(bbox);
    }
    if (!warpedImage || warpedImage->lattice.dimsInt()[0] != dims[0] ||
        warpedImage->lattice.dimsInt()[1] != dims[1] ||
        warpedImage->lattice.dimsInt()[2] != dims[2])
        warpedImage = 0;

    if (!warpedImage) {
        if (toImage->isOfType(HxUniformLabelField3::getClassTypeId())) {
            warpedImage = new HxUniformLabelField3(dims);
            ((HxUniformLabelField3*)warpedImage)->parameters =
                ((HxUniformLabelField3*)fromImage)->parameters;
        } else
            warpedImage =
                new HxUniformScalarField3(dims, fromImage->primType());
    }

    memcpy(warpedImage->bbox(), bbox, 6 * sizeof(float));
    warpedImage->composeLabel(fromImage->getLabel().getString(), "Warped");
    setResult(0, warpedImage);
    return (warpedImage);
}

void HxMovingLeastSquaresWarp::prepareLandmarks(McDArray<McVec2d>& p1,
                                               McDArray<McVec2d>& p2) {
    int set1 = 0;
    int set2 = 1;

    HxLandmarkSet* pointSet = hxconnection_cast<HxLandmarkSet>(portData);

    if (!pointSet)
        return;

    p1.resize(0);
    p2.resize(0);
    int nPoints = pointSet->getNumMarkers();
    for (int i = 0; i < nPoints; i++) {
        p1.append(McVec2d(pointSet->getCoords(set1)[i].x,
                          pointSet->getCoords(set1)[i].y));
        p2.append(McVec2d(pointSet->getCoords(set2)[i].x,
                          pointSet->getCoords(set2)[i].y));
    }
}

void HxMovingLeastSquaresWarp::compute() {
    if (!portAction.wasHit())
        return;

    HxUniformScalarField3* fromImage =
        (HxUniformScalarField3*)portFromImage.source();

    HxLandmarkSet* pointSet = hxconnection_cast<HxLandmarkSet>(portData);

    if (!fromImage)
        return;

    /* It is ok to warp without landmarks, if method is rigid and
       input has a transformation */
    if (!pointSet) {
        return;
    } else if (pointSet->getNumSets() < 2) {
        theMsg->printf(
            "Error: LandmarkWarp data has to contain at least 2 sets.");
        return;
    }

    HxUniformScalarField3* outputImage;

    HxUniformVectorField3* outputVectorField = 0;

    outputImage = createOutputDataSet();
    outputVectorField = createOutputVectorDataSet();

    float* outputImageData = (float*)outputImage->lattice.dataPtr();

    McDArray<McVec2d> landmarks1, landmarks2;

    prepareLandmarks(landmarks1, landmarks2);

    MovingLeastSquares mls;
    mls.setAlpha(portAlpha.getValue());
    mls.setLandmarks(landmarks2, landmarks1);

    const int* dims = outputImage->lattice.dimsInt();

    McVec3f voxelSizeInOutputImage = outputImage->getVoxelSize();
    float bboxOfOutputImage[6];
    outputImage->getBoundingBox(bboxOfOutputImage);
    float bboxOfFromImage[6];
    fromImage->getBoundingBox(bboxOfFromImage);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < dims[0]; i++) {
        HxLocation3* locationInFromImage = fromImage->createLocation();
        std::cout << "\n" << i << " of " << dims[0];
        for (int j = 0; j < dims[1]; j++) {
            McVec2d currentPositionInOutput = McVec2d(
                bboxOfOutputImage[0] + (float)(i)*voxelSizeInOutputImage.x,
                bboxOfOutputImage[2] + (float)(j)*voxelSizeInOutputImage.y);
            McVec2d warpedCoordInFromImage =
                mls.interpolate(currentPositionInOutput);
            McVec3f displacement = McVec3f(
                warpedCoordInFromImage.x - currentPositionInOutput.x,
                warpedCoordInFromImage.y - currentPositionInOutput.y, 0);

            displacement = displacement * -1;
            locationInFromImage->move(McVec3f(warpedCoordInFromImage.x,
                                              warpedCoordInFromImage.y,
                                              bboxOfFromImage[4]));
            float resultValue[1];
            fromImage->eval(locationInFromImage, resultValue);
            unsigned long pos = latticePos(i, j, 0, dims);
            outputImageData[pos] = resultValue[0];
            outputVectorField->lattice.set(i, j, 0, displacement.getValue());
        }
        delete locationInFromImage;
    }

    outputImage->touch();
    outputImage->fire();
}
