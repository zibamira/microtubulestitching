#include <hxalignmicrotubules/CoherentPointDriftRigidFisherMises.h>

#include <stdio.h>
#include <iostream>

#include <gtest/gtest.h>
#include <hxcore/TestingData.h>
#include <hxcore/TestingObjectPoolCleaner.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/McMat3f.h>
#include <mclib/TestingDevNullRedirect.h>

#include <hxalignmicrotubules/mtalign/rotation.h>

static void printMatrix(const McDMatrix<double>& mat, McString name) {
    std::cout << "\n" << name << " =";
    std::cout << "[";
    for (int i = 0; i < mat.nRows(); i++) {
        std::cout << "\n";
        for (int j = 0; j < mat.nCols(); j++) {
            if (j < mat.nCols() - 1)
                std::cout << mat[i][j] << ", ";
            else
                std::cout << mat[i][j];
        }
        if (i < mat.nRows() - 1)
            std::cout << ";";
    }
    std::cout << "]\n";
}

static CoherentPointDriftRigidFisherMises initCPDdata1NoDir() {
    CoherentPointDriftRigidFisherMises cpd;

    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = false;
    // create data

    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(10, 10, 0));
    tempD.append(McVec3f(0, 0, 0));

    for (int i = 0; i < tempD.size(); i++) {
        xd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yc.append(xc[i] + McVec3f(1.0, 1.0, 0.0));
    }
    yc.append(yc[2] + McVec3f(1.0, 1.0, 0.0));
    yd.append(McVec3f(0, 0, 0));
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdata1() {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(1, 1, 1));
    tempD.append(McVec3f(1, 4, 0));
    tempD.append(McVec3f(1, 4, 8));

    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        xd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yc.append(xc[i] + McVec3f(1.0, 1.0, 0.0));
    }
    yc.append(yc[2] + McVec3f(1.0, 1.1, 0.0));
    yd.append(McVec3f(tempD[2].x, tempD[2].y, tempD[2].z));
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdata2() {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(1, 1, 1));
    tempD.append(McVec3f(1, 4, 0));
    tempD.append(McVec3f(1, 4, 8));
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        xd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yd.append(McVec3f(tempD[i].x, tempD[i].y, tempD[i].z));
        yc.append(xc[i] + McVec3f(1.0, 1.0, 0.0));
    }
    xc.append(yc[2] + McVec3f(1.0, 1.1, 0.0));
    xd.append(McVec3f(tempD[2].x, tempD[2].y, tempD[2].z));
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdata2Rotated() {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    xc.append(McVec3f(1, 1, 0));
    tempD.append(McVec3f(1, 1, 1));

    tempD.append(McVec3f(1, 1, 16));
    tempD.append(McVec3f(1, 4, 0));
    tempD.append(McVec3f(1, 4, 8));
    McMat3f R;
    R.setValue(cos(0.3), -1 * sin(0.3), 0, sin(0.3), cos(0.3), 0, 0, 0, 1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        McVec3f transCoord, rotCoord, rotDir;
        transCoord = McVec3f(1.0, 1.0, 0.0) + xc[i];
        std::cout << "\n transCoord: " << transCoord.x << " " << transCoord.y;
        R.multMatrixVec(transCoord, rotCoord);
        std::cout << "\n rotCoord:   " << rotCoord.x << " " << rotCoord.y;
        R.multMatrixVec(tempD[i], rotDir);
        xd.append(tempD[i]);
        yd.append(rotDir);
        yc.append(rotCoord);
    }
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdata2NoDirRotated() {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = false;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    tempD.append(McVec3f(0, 0, 0));
    McMat3f R;
    R.setValue(cos(0.3), -sin(0.3), 0, sin(0.3), cos(0.3), 0, 0, 0, 1);
    for (int i = 0; i < tempD.size(); i++) {
        McVec3f transCoord, rotCoord, rotDir;
        transCoord = McVec3f(1.0, 1.0, 0.0) + xc[i];
        R.multMatrixVec(transCoord, rotCoord);
        R.multMatrixVec(tempD[i], rotDir);
        xd.append(tempD[i]);
        yd.append(rotDir);
        yc.append(rotCoord);
    }
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdataVMOnlyScaled(const double scale) {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(1, 2, 3));
    tempD.append(McVec3f(0, 10, 0));
    tempD.append(McVec3f(4, 0, 0));
    McMat3f R;
    double ang = 0.0;
    R.setValue(cos(ang), -1 * sin(ang), 0, sin(ang), cos(ang), 0, 0, 0, 1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        McVec3f transCoordRotated, transDirRotated;
        R.multMatrixVec(xc[i], transCoordRotated);
        R.multMatrixVec(tempD[i], transDirRotated);
        transCoordRotated *= scale;
        xd.append(tempD[i]);
        yd.append(transDirRotated);
        yc.append(transCoordRotated);
    }

    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises
initCPDdataVMOnlyRotated(const double ang) {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(1, 2, 3));
    tempD.append(McVec3f(0, 10, 0));
    tempD.append(McVec3f(4, 0, 0));
    McMat3f R;
    R.setValue(cos(ang), -1 * sin(ang), 0, sin(ang), cos(ang), 0, 0, 0, 1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        McVec3f transCoordRotated, transDirRotated;

        R.multMatrixVec(xc[i], transCoordRotated);
        R.multMatrixVec(tempD[i], transDirRotated);
        xd.append(tempD[i]);
        yd.append(transDirRotated);
        yc.append(transCoordRotated);
    }
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises
initCPDdataVMScaledAndRotated(const double scale, const double angle) {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(1, 1, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(1, 1, 0));
    tempD.append(McVec3f(1, 0, 0));
    tempD.append(McVec3f(0, 1, 0));

    McMat3f R;

    R.setValue(cos(angle), -1 * sin(angle), 0, sin(angle), cos(angle), 0, 0, 0,
               1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();

        // here longer?
        McVec3f transCoordRotated, transDirRotated;

        R.multMatrixVec(xc[i], transCoordRotated);
        R.multMatrixVec(tempD[i], transDirRotated);
        xd.append(tempD[i]);
        yd.append(transDirRotated);
        transCoordRotated *= scale;
        yc.append(transCoordRotated);
    }
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises
initCPDdataVMScaledAndRotatedMoreComplex(const double scale,
                                         const double angle) {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(1, 10, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    xc.append(McVec3f(10, 10, 0));
    xc.append(McVec3f(10, 2, 0));
    xc.append(McVec3f(1, 6, 0));
    tempD.append(McVec3f(1, 1, 0));
    tempD.append(McVec3f(1, 90, 2));
    tempD.append(McVec3f(0, 1, 0));
    tempD.append(McVec3f(1, 3, 0));
    tempD.append(McVec3f(5, 0, 5));
    tempD.append(McVec3f(0, 1, 0));

    McMat3f R;

    R.setValue(cos(angle), -1 * sin(angle), 0, sin(angle), cos(angle), 0, 0, 0,
               1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        McVec3f transCoordRotated, transDirRotated;

        R.multMatrixVec(xc[i], transCoordRotated);
        R.multMatrixVec(tempD[i], transDirRotated);
        xd.append(tempD[i]);
        yd.append(transDirRotated);
        transCoordRotated *= scale;
        yc.append(transCoordRotated);
    }
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

static CoherentPointDriftRigidFisherMises initCPDdataVMRotated() {
    CoherentPointDriftRigidFisherMises cpd;
    cpd.params.withScaling = true;
    cpd.params.usePositions = true;
    cpd.params.useDirections = true;
    // create data
    McDArray<McVec3f> xc, xd, yc, yd;
    McDArray<McVec3f> tempD;
    McDArray<double> xz, yz;
    xc.append(McVec3f(0, 0, 0));
    xc.append(McVec3f(0, 1, 0));
    xc.append(McVec3f(1, 0, 0));
    xc.append(McVec3f(10, 10, 0));
    tempD.append(McVec3f(1, 2, 3));
    tempD.append(McVec3f(0, 10, 0));
    tempD.append(McVec3f(4, 0, 0));
    tempD.append(McVec3f(4, 20, 10));
    xc.append(McVec3f(-100, -100, 0));
    tempD.append(McVec3f(4, 0, 10));
    McMat3f R;
    R.setValue(cos(0.3), -sin(0.3), 0, sin(0.3), cos(0.3), 0, 0, 0, 1);
    for (int i = 0; i < tempD.size(); i++) {
        tempD[i].normalize();
        McVec3f transCoord, rotCoord, rotDir;
        transCoord = McVec3f(1.0, 1.0, 0.0) + xc[i];
        R.multMatrixVec(transCoord, rotCoord);
        R.multMatrixVec(tempD[i], rotDir);
        xd.append(tempD[i]);
        yd.append(rotDir);
        yc.append(rotCoord);
    }
    yd.append(yd[2]);
    yc.append(yc[2] + McVec3f(0.01, 0.05, 0.0));
    cpd.convertMcDArraysToMatrix(xc, xd, cpd.Xc, cpd.Xd);
    cpd.convertMcDArraysToMatrix(yc, yd, cpd.Yc, cpd.Yd);
    return cpd;
}

TEST(CoherentPointDriftRigidFisherMises, givenMat_TestTace_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    McDMatrix<double> idMat(10, 10);
    idMat.makeIdentity();
    EXPECT_EQ(CoherentPointDriftRigidFisherMises::trace(idMat), 10);
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDScaledAndRotated_CheckConvergesR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    double rho = 1.8;
    double scale = 1.0;
    CoherentPointDriftRigidFisherMises cpdr =
        initCPDdataVMScaledAndRotatedMoreComplex(scale, rho);
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;
    cpdr.params.maxIterations = 100;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);
    double rhoErg = mtalign::rotationAngle2d(Rc);
    EXPECT_NEAR(rhoErg, -rho, 5.e-3);
    EXPECT_NEAR(s, 1.0 / scale, 1.e-3);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.001);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDScaledAndRotated_CheckAlignVMCoupledRsRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    double rho = 2.0;
    double scale = 20.0;
    CoherentPointDriftRigidFisherMises cpdr =
        initCPDdataVMScaledAndRotatedMoreComplex(scale, rho);
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;
    cpdr.params.maxIterations = 5;

    McDArray<McVec2i> dummy;

    cpdr.align(Rc, s, t, Rd, dummy);
    double rhoErg = mtalign::rotationAngle2d(Rc);
    EXPECT_NEAR(rhoErg, -rho, 5.e-3);
    EXPECT_NEAR(s, 1.0 / scale, 1.e-3);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.001);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyScaled_CheckUncoupledConvergesRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMScaledAndRotated(10.0, 0.5);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 100;
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;

    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyScaled_CheckConverges1R_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMScaledAndRotated(1.0, 0.0);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 10;
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;

    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDRotAndScaled_CheckConverges2R_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMScaledAndRotated(1.0, 0.6);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 100;
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.002);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.002);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyRotated_CheckConverges3R_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMOnlyRotated(0.5);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 10;
    std::cout << "\n maxIterations should be 1!";
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyScaled_CheckConverges8R_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMOnlyScaled(10.0);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 100;
    std::cout << "\n maxIterations should be 1!";
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyScaled_CheckConverges4R_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMOnlyScaled(1.0);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    printMatrix(cpdr.Xd, "Xd");
    printMatrix(cpdr.Yd, "Yd");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 5;
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPDOnlyRotated_CheckAlignVMCoupledRsRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMOnlyRotated(0.3);
    printMatrix(cpdr.Xc, "Xc");
    printMatrix(cpdr.Yc, "Yc");
    cpdr.params.w = 0.1;
    cpdr.params.maxIterations = 10;
    cpdr.params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD1_CheckAlignVMRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdataVMRotated();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(Rc, s, t, Rd, dummy);

    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.001);
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPD2NoDirRotated_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata2NoDirRotated();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(R, s, t, Rd, dummy);
    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
        std::cout << "\n Xc[" << i << "]:" << cpdr.Xc[i][0] << " "
                  << cpdr.Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << cpdr.Yc[i][0] << " "
                  << cpdr.Yc[i][1];
        std::cout << "\n shiftedYc[" << i << "]:" << shiftedYc[i][0] << " "
                  << shiftedYc[i][1];
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPD2Rotated_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata2Rotated();
    for (int i = 0; (i < cpdr.Xc.nRows()); i++) {

        std::cout << "\n Xc[" << i << "]:" << cpdr.Xc[i][0] << " "
                  << cpdr.Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << cpdr.Yc[i][0] << " "
                  << cpdr.Yc[i][1];
    }
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    McDArray<McVec2i> dummy;
    // cpdr.mWithScale=false;
    cpdr.align(R, s, t, Rd, dummy);
    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, R, s, t, shiftedYc);

    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
        std::cout << "\n Xc[" << i << "]:" << cpdr.Xc[i][0] << " "
                  << cpdr.Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << cpdr.Yc[i][0] << " "
                  << cpdr.Yc[i][1];
        std::cout << "\n shiftedYc[" << i << "]:" << shiftedYc[i][0] << " "
                  << shiftedYc[i][1];
    }
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPD2_CheckRotationMatrixIsRotationMatrixR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata2();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, P;

    McDArray<McVec2i> dummy;

    /*******************************************************************************/
    // perform E-step
    // 1. Compute new ys
    // no shift, since we only test one iteration
    // compute the P matrix, conditional P(x|m)
    std::cout << "\n..compute P..";
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 0.783331, 0.1, cpdr.params.w, P,
                  LL);

    /*******************************************************************************/
    // M-Step
    // compute Np
    double Np = cpdr.getNP(P);
    // compue muX and muY
    McDMatrix<double> PT = P;
    PT.transpose();
    McDMatrix<double> muX, muY;
    cpdr.getMu(cpdr.Xc, PT, Np, muX);
    cpdr.getMu(cpdr.Yc, P, Np, muY);
    McDMatrix<double> XcHat, YcHat;
    cpdr.getHat(cpdr.Xc, muX, XcHat);
    cpdr.getHat(cpdr.Yc, muY, YcHat);
    // Initialize all the long expressions
    CoherentPointDriftRigidFisherMises::allTerms terms;
    cpdr.initTerms(terms, XcHat, cpdr.Xd, YcHat, cpdr.Yd, P);

    // solve for R, s and t
    // compute R
    cpdr.computeRc(terms, R);
    std::cout << "\nR is :\n" << R[0][0] << " " << R[0][1] << "\n" << R[1][0]
              << " " << R[1][1] << "\n";
    EXPECT_NEAR(R.det(), 1.0, 0.0001);
    McDMatrix<double> RT = R;
    RT.transpose();
    R *= RT;
    EXPECT_NEAR(R[0][0], 1.0, 0.0001);
    EXPECT_NEAR(R[1][1], 1.0, 0.0001);
    EXPECT_NEAR(R[1][0], 0.0, 0.0001);
    EXPECT_NEAR(R[0][1], 0.0, 0.0001);

    cpdr.computeS(terms, R, s);
    EXPECT_TRUE(s > 0.0);
    std::cout << "...new s is " << s << "...";
    cpdr.computeT(muX, s, R, muY, t);
    std::cout << "...new t is " << t[0] << " " << t[1] << "...";
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPD1NoDir_CheckRotationMatrixIsRotationMatrixR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1NoDir();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    t.fill(0.0);
    double s = 1.0;
    McDMatrix<double> R, P;

    McDArray<McVec2i> dummy;
    /*******************************************************************************/
    // perform E-step
    // 1. Compute new ys
    // no shift, since we only test one iteration
    // compute the P matrix, conditional P(x|m)
    std::cout << "\n..compute P..";
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 0.783331, 0.1, cpdr.params.w, P,
                  LL);

    /*******************************************************************************/
    // M-Step
    // compute Np
    double Np = cpdr.getNP(P);
    std::cout << "...Np is " << Np;
    // compue muX and muY
    McDMatrix<double> PT = P;
    PT.transpose();
    McDMatrix<double> muX, muY;
    cpdr.getMu(cpdr.Xc, PT, Np, muX);
    std::cout << "...muX is " << muX[0][0] << " " << muX[1][0];
    cpdr.getMu(cpdr.Yc, P, Np, muY);
    std::cout << "...muY is " << muY[0][0] << " " << muY[1][0];
    McDMatrix<double> XcHat, YcHat;
    cpdr.getHat(cpdr.Xc, muX, XcHat);
    cpdr.getHat(cpdr.Yc, muY, YcHat);
    // Initialize all the long expressions
    CoherentPointDriftRigidFisherMises::allTerms terms;
    cpdr.initTerms(terms, XcHat, cpdr.Xd, YcHat, cpdr.Yd, P);

    // solve for R, s and t
    // compute R
    cpdr.computeRc(terms, R);
    std::cout << "\nR is :\n" << R[0][0] << " " << R[0][1] << "\n" << R[1][0]
              << " " << R[1][1] << "\n";
    EXPECT_NEAR(R.det(), 1.0, 0.0001);
    McDMatrix<double> RT = R;
    RT.transpose();
    R *= RT;
    EXPECT_NEAR(R[0][0], 1.0, 0.0001);
    EXPECT_NEAR(R[1][1], 1.0, 0.0001);
    EXPECT_NEAR(R[1][0], 0.0, 0.0001);
    EXPECT_NEAR(R[0][1], 0.0, 0.0001);

    cpdr.computeS(terms, R, s);
    EXPECT_TRUE(s > 0.0);
    std::cout << "...new s is " << s << "...";
    cpdr.computeT(muX, s, R, muY, t);
    std::cout << "...new t is " << t[0] << " " << t[1] << "...";
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD1_CheckGetHatRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1();
    cpdr.params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    std::cout << "...next: computeP()...";
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 1, 0.1, 0.1, P, LL);

    P.transpose();
    double Np = cpdr.getNP(P);
    std::cout << "...next: getMu()...";
    cpdr.getMu(cpdr.Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);

    McDMatrix<double> XHat;
    std::cout << "...next: getHat()...";
    cpdr.getHat(cpdr.Xc, muX, XHat);
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD2_CheckGetHatRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata2();
    cpdr.params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 1, 0.1, 0.1, P, LL);

    P.transpose();
    double Np = cpdr.getNP(P);
    cpdr.getMu(cpdr.Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);

    McDMatrix<double> XHat;
    cpdr.getHat(cpdr.Xc, muX, XHat);
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD_CheckGetMuRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1();
    cpdr.params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 1, 0.1, 0.1, P, LL);

    double Np = cpdr.getNP(P);
    P.transpose();
    cpdr.getMu(cpdr.Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);
}

TEST(CoherentPointDriftRigidFisherMises,
     givenTestCPD1NoDir_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1NoDir();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(R, s, t, Rd, dummy);
    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        std::cout << "\n point: " << i << "\n";
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD2_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata2();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(R, s, t, Rd, dummy);
    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD1_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1();
    cpdr.params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    McDArray<McVec2i> dummy;
    cpdr.align(R, s, t, Rd, dummy);
    McDMatrix<double> shiftedYc;
    cpdr.shiftYCoords(cpdr.Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < cpdr.Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], cpdr.Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], cpdr.Xc[i][1], 0.01);
    }
}

TEST(CoherentPointDriftRigidFisherMises, givenTestCPD1_CheckPOKR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1();
    McDMatrix<double> P;
    double LL;
    cpdr.computeP(cpdr.Xc, cpdr.Xd, cpdr.Yc, cpdr.Yd, 1, 0.1, 0.1, P, LL);
    for (int i = 0; i < P.nRows(); i++) {
        for (int j = 0; j < P.nCols(); j++)
            std::cout << P[i][j] << " ";
        std::cout << "\n";
    }
}

TEST(CoherentPointDriftRigidFisherMises, givenNothing_CheckInitCPDOKR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    CoherentPointDriftRigidFisherMises cpdr = initCPDdata1();
    EXPECT_EQ(cpdr.Xc.nRows(), 3);
    EXPECT_EQ(cpdr.Xc.nCols(), 2);
    EXPECT_EQ(cpdr.Xd.nRows(), 3);
    EXPECT_EQ(cpdr.Xd.nCols(), 3);
    EXPECT_EQ(cpdr.Yc.nRows(), 4);
    EXPECT_EQ(cpdr.Yc.nCols(), 2);
    EXPECT_EQ(cpdr.Yd.nRows(), 4);
    EXPECT_EQ(cpdr.Yd.nCols(), 3);
}

TEST(CoherentPointDriftRigidFisherMises,
     givenMtrices_lapackSolvesCorrect_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    double aArray[6] = { 4, 7, 0, 13 };
    McDMatrix<double> A(2, 2, aArray);
    double bArray[9] = { 1, 2, 3, 4 };
    McDMatrix<double> B(2, 2, bArray);

    // pseudoinverse
    McDMatrix<double> AT = A;
    AT.transpose();
    McDMatrix<double> temp = AT * A;
    temp.inverse();
    McDMatrix<double> result = temp * AT * B;
    // same with lapack
    McDMatrix<double> result_blas;
    CoherentPointDriftRigidFisherMises::solve(A, B, result_blas);
    for (int i = 0; i < result.nEntries(); i++) {

        EXPECT_NEAR(result_blas.dataPtr()[i], result.dataPtr()[i], 1.e-7);
    }
    EXPECT_EQ(result_blas.nRows(), result.nRows());
    EXPECT_EQ(result_blas.nCols(), result.nCols());
}

TEST(CoherentPointDriftRigidFisherMises, givenMtrices_blasMultsCorrect_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    double aArray[6] = { 1, 2, 3, 4, 5, 6 };
    McDMatrix<double> A(2, 3, aArray);
    double bArray[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    McDMatrix<double> B(3, 3, bArray);
    double cArray[6] = { 1, 2, 3, 4, 5, 6 };
    McDMatrix<double> C(2, 3, cArray);
    McDMatrix<double> result = A * B + C;
    McDMatrix<double> result_blas;
    CoherentPointDriftRigidFisherMises::multMatrices(A, B, C, 1.0, 1.0,
                                                     result_blas);
    for (int i = 0; i < result_blas.nEntries(); i++) {

        EXPECT_NEAR(result_blas.dataPtr()[i], result.dataPtr()[i], 1.e-7);
    }
    EXPECT_EQ(result_blas.nRows(), result.nRows());
    EXPECT_EQ(result_blas.nCols(), result.nCols());
}

