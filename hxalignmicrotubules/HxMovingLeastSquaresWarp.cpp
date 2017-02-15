#include <hxalignmicrotubules/HxMovingLeastSquaresWarp.h>

#include <hxcore/HxMessage.h>
#include <hxfield/HxUniformLabelField3.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxUniformVectorField3.h>
#include <hxlandmark/internal/HxLandmarkSet.h>

#ifdef _OPENMP
#include <omp.h>
#include <hxcore/HxSettingsMgr.h>
#endif

#include <hxalignmicrotubules/MovingLeastSquares.h>

static mculong latticePos(int i, int j, int k, const McDim3l& dims) {
    return (mculong)dims[0] * (mculong)dims[1] * (mculong)k +
           (mculong)j * (mculong)dims[0] + (mculong)i;
}

HX_INIT_CLASS(HxMovingLeastSquaresWarp, HxCompModule);

HxMovingLeastSquaresWarp::HxMovingLeastSquaresWarp()
    : HxCompModule(HxLandmarkSet::getClassTypeId()),
      portFromImage(this, "from", tr("From"), HxUniformScalarField3::getClassTypeId()),
      portToImage(this, "to", tr("To"), HxUniformScalarField3::getClassTypeId()),
      portMethod(this, "method", tr("Method"), 1),
      portAlpha(this, "parameter", tr("Parameter")),
      portAction(this, "action", tr("Action")) {
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
        dynamic_cast<HxUniformScalarField3*>(portFromImage.getSource());
    HxUniformVectorField3* output =
        dynamic_cast<HxUniformVectorField3*>(getResult(1));

    if (!fromImage)
        return (0);

    McDim3l dims;
    McBox3f bbox;

    if (fromImage) {
        dims = fromImage->lattice().getDims();
        bbox = fromImage->getBoundingBox();
    }
    if (!output || output->lattice().getDims() != dims)
        output = 0;

    if (!output) {
        output = new HxUniformVectorField3(dims, McPrimType::MC_FLOAT);
    }
    output->lattice().setBoundingBox(bbox);
    output->composeLabel(fromImage->getLabel(), "Displacement");
    setResult(1, output);
    return (output);
}

HxUniformScalarField3* HxMovingLeastSquaresWarp::createOutputDataSet() {
    HxUniformScalarField3* fromImage =
        dynamic_cast<HxUniformScalarField3*>(portFromImage.getSource());
    HxUniformScalarField3* toImage =
        dynamic_cast<HxUniformScalarField3*>(portToImage.getSource());
    HxUniformScalarField3* warpedImage =
        dynamic_cast<HxUniformScalarField3*>(getResult(0));

    if (!fromImage)
        return (0);

    McDim3l dims;
    McBox3f bbox;

    if (toImage) {
        dims = toImage->lattice().getDims();
        bbox = toImage->getBoundingBox();
    }
    if (!warpedImage || warpedImage->lattice().getDims()[0] != dims[0] ||
        warpedImage->lattice().getDims()[1] != dims[1] ||
        warpedImage->lattice().getDims()[2] != dims[2])
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

    warpedImage->setBoundingBox(bbox);
    warpedImage->composeLabel(fromImage->getLabel(), "Warped");
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
        (HxUniformScalarField3*)portFromImage.getSource();

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

    float* outputImageData = (float*)outputImage->lattice().dataPtr();

    McDArray<McVec2d> landmarks1, landmarks2;

    prepareLandmarks(landmarks1, landmarks2);

    MovingLeastSquares mls;
    mls.setAlpha(portAlpha.getValue());
    mls.setLandmarks(landmarks2, landmarks1);

    const McDim3l& dims = outputImage->lattice().getDims();

    McVec3f voxelSizeInOutputImage = outputImage->getVoxelSize();
    const McBox3f& bboxOfOutputImage = outputImage->getBoundingBox();
    const McBox3f& bboxOfFromImage = fromImage->getBoundingBox();
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
            fromImage->eval(*locationInFromImage, resultValue);
            unsigned long pos = latticePos(i, j, 0, dims);
            outputImageData[pos] = resultValue[0];
            outputVectorField->lattice().set(i, j, 0, displacement.getValue());
        }
        delete locationInFromImage;
    }

    outputImage->touch();
    outputImage->fire();
}
