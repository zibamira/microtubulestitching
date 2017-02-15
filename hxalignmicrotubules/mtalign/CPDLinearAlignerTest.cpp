#include <hxalignmicrotubules/mtalign/CPDLinearAligner.h>

#include <stdio.h>
#include <iostream>

#include <gtest/internal/gtest.h>
#include <hxcore/internal/TestingData.h>
#include <hxcore/internal/TestingObjectPoolCleaner.h>
#include <hxspatialgraph/internal/HxSpatialGraph.h>
#include <mclib/McMat3f.h>
#include <mclib/internal/TestingDevNullRedirect.h>

#include <hxalignmicrotubules/mtalign/rotation.h>
#include <hxalignmicrotubules/mtalign/math.h>
#include <hxalignmicrotubules/mtalign/data.h>

namespace ma = mtalign;

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

// Derive from CPDLinearAligner to access protected members
// in tests.  `this->` is used below to make the access explicit.
class mtalign__CPDLinearAlignerTest : public ::testing::Test,
                                      public ma::CPDLinearAligner {

  protected:
    void initCPDdata1NoDir() {

        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = false;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdata1() {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdata2() {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdata2Rotated() {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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
            std::cout << "\n transCoord: " << transCoord.x << " "
                      << transCoord.y;
            R.multMatrixVec(transCoord, rotCoord);
            std::cout << "\n rotCoord:   " << rotCoord.x << " " << rotCoord.y;
            R.multMatrixVec(tempD[i], rotDir);
            xd.append(tempD[i]);
            yd.append(rotDir);
            yc.append(rotCoord);
        }
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdata2NoDirRotated() {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = false;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdataVMOnlyScaled(const double scale) {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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

        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdataVMOnlyRotated(const double ang) {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdataVMScaledAndRotated(const double scale, const double angle) {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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

        R.setValue(cos(angle), -1 * sin(angle), 0, sin(angle), cos(angle), 0, 0,
                   0, 1);
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdataVMScaledAndRotatedMoreComplex(const double scale,
                                                  const double angle) {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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

        R.setValue(cos(angle), -1 * sin(angle), 0, sin(angle), cos(angle), 0, 0,
                   0, 1);
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }

    void initCPDdataVMRotated() {
        this->params.withScaling = true;
        this->params.usePositions = true;
        this->params.useDirections = true;
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
        ma::FacingPointSets points;
        points.ref.positions = xc;
        points.ref.directions = xd;
        points.trans.positions = yc;
        points.trans.directions = yd;
        this->setPoints(points);
    }
};

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDScaledAndRotated_CheckConvergesR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    double rho = 1.8;
    double scale = 1.0;
    initCPDdataVMScaledAndRotatedMoreComplex(scale, rho);
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;
    this->params.maxIterations = 100;

    this->align(Rc, s, t, Rd);
    double rhoErg = mtalign::rotationAngle2d(Rc);
    EXPECT_NEAR(rhoErg, -rho, 5.e-3);
    EXPECT_NEAR(s, 1.0 / scale, 1.e-3);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.001);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDScaledAndRotated_CheckAlignVMCoupledRsRunsR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    double rho = 2.0;
    double scale = 20.0;
    initCPDdataVMScaledAndRotatedMoreComplex(scale, rho);
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;
    this->params.maxIterations = 5;

    this->align(Rc, s, t, Rd);
    double rhoErg = mtalign::rotationAngle2d(Rc);
    EXPECT_NEAR(rhoErg, -rho, 5.e-3);
    EXPECT_NEAR(s, 1.0 / scale, 1.e-3);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.001);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyScaled_CheckUncoupledConvergesRunsR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMScaledAndRotated(10.0, 0.5);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 100;
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyScaled_CheckConverges1R_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMScaledAndRotated(1.0, 0.0);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 10;
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDRotAndScaled_CheckConverges2R_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMScaledAndRotated(1.0, 0.6);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 100;
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.002);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.002);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyRotated_CheckConverges3R_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMOnlyRotated(0.5);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 10;
    std::cout << "\n maxIterations should be 1!";
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);
    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyScaled_CheckConverges8R_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMOnlyScaled(10.0);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 100;
    std::cout << "\n maxIterations should be 1!";
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyScaled_CheckConverges4R_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMOnlyScaled(1.0);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    printMatrix(this->Xd, "Xd");
    printMatrix(this->Yd, "Yd");
    this->params.w = 0.1;
    this->params.maxIterations = 5;
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPDOnlyRotated_CheckAlignVMCoupledRsRunsR_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMOnlyRotated(0.3);
    printMatrix(this->Xc, "Xc");
    printMatrix(this->Yc, "Yc");
    this->params.w = 0.1;
    this->params.maxIterations = 10;
    this->params.eDiffRelStop = 1.e-6;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    printMatrix(Rc, "newR");
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD1_CheckAlignVMRunsR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdataVMRotated();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> Rc, Rd;

    this->align(Rc, s, t, Rd);

    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, Rc, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.001);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.001);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPD2NoDirRotated_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata2NoDirRotated();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    this->align(R, s, t, Rd);
    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
        std::cout << "\n Xc[" << i << "]:" << this->Xc[i][0] << " "
                  << this->Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << this->Yc[i][0] << " "
                  << this->Yc[i][1];
        std::cout << "\n shiftedYc[" << i << "]:" << shiftedYc[i][0] << " "
                  << shiftedYc[i][1];
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPD2Rotated_CheckAlignRunsR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata2Rotated();
    for (int i = 0; (i < this->Xc.nRows()); i++) {

        std::cout << "\n Xc[" << i << "]:" << this->Xc[i][0] << " "
                  << this->Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << this->Yc[i][0] << " "
                  << this->Yc[i][1];
    }
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    // this->mWithScale=false;
    this->align(R, s, t, Rd);
    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, R, s, t, shiftedYc);

    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
        std::cout << "\n Xc[" << i << "]:" << this->Xc[i][0] << " "
                  << this->Xc[i][1];
        std::cout << "\n Yc[" << i << "]:" << this->Yc[i][0] << " "
                  << this->Yc[i][1];
        std::cout << "\n shiftedYc[" << i << "]:" << shiftedYc[i][0] << " "
                  << shiftedYc[i][1];
    }
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPD2_CheckRotationMatrixIsRotationMatrixR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata2();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, P;

    /*******************************************************************************/
    // perform E-step
    // 1. Compute new ys
    // no shift, since we only test one iteration
    // compute the P matrix, conditional P(x|m)
    std::cout << "\n..compute P..";
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 0.783331, 0.1,
                   this->params.w, P, LL);

    /*******************************************************************************/
    // M-Step
    // compute Np
    double Np = this->getNP(P);
    // compue muX and muY
    McDMatrix<double> PT = P;
    PT.transpose();
    McDMatrix<double> muX, muY;
    this->getMu(this->Xc, PT, Np, muX);
    this->getMu(this->Yc, P, Np, muY);
    McDMatrix<double> XcHat, YcHat;
    this->getHat(this->Xc, muX, XcHat);
    this->getHat(this->Yc, muY, YcHat);
    // Initialize all the long expressions
    ma::CPDLinearAlignerTerms terms(XcHat, this->Xd, YcHat, this->Yd, P);

    // solve for R, s and t
    // compute R
    terms.computeRc(R);
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

    terms.computeS(R, s);
    EXPECT_TRUE(s > 0.0);
    std::cout << "...new s is " << s << "...";
    this->computeT(muX, s, R, muY, t);
    std::cout << "...new t is " << t[0] << " " << t[1] << "...";
}

TEST_F(mtalign__CPDLinearAlignerTest,
       givenTestCPD1NoDir_CheckRotationMatrixIsRotationMatrixR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1NoDir();
    this->params.w = 0.1;
    McDVector<double> t;
    t.fill(0.0);
    double s = 1.0;
    McDMatrix<double> R, P;

    /*******************************************************************************/
    // perform E-step
    // 1. Compute new ys
    // no shift, since we only test one iteration
    // compute the P matrix, conditional P(x|m)
    std::cout << "\n..compute P..";
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 0.783331, 0.1,
                   this->params.w, P, LL);

    /*******************************************************************************/
    // M-Step
    // compute Np
    double Np = this->getNP(P);
    std::cout << "...Np is " << Np;
    // compue muX and muY
    McDMatrix<double> PT = P;
    PT.transpose();
    McDMatrix<double> muX, muY;
    this->getMu(this->Xc, PT, Np, muX);
    std::cout << "...muX is " << muX[0][0] << " " << muX[1][0];
    this->getMu(this->Yc, P, Np, muY);
    std::cout << "...muY is " << muY[0][0] << " " << muY[1][0];
    McDMatrix<double> XcHat, YcHat;
    this->getHat(this->Xc, muX, XcHat);
    this->getHat(this->Yc, muY, YcHat);
    // Initialize all the long expressions
    ma::CPDLinearAlignerTerms terms(XcHat, this->Xd, YcHat, this->Yd, P);

    // solve for R, s and t
    // compute R
    terms.computeRc(R);
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

    terms.computeS(R, s);
    EXPECT_TRUE(s > 0.0);
    std::cout << "...new s is " << s << "...";
    this->computeT(muX, s, R, muY, t);
    std::cout << "...new t is " << t[0] << " " << t[1] << "...";
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD1_CheckGetHatRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1();
    this->params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    std::cout << "...next: computeP()...";
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 1, 0.1, 0.1, P, LL);

    P.transpose();
    double Np = this->getNP(P);
    std::cout << "...next: getMu()...";
    this->getMu(this->Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);

    McDMatrix<double> XHat;
    std::cout << "...next: getHat()...";
    this->getHat(this->Xc, muX, XHat);
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD2_CheckGetHatRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata2();
    this->params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 1, 0.1, 0.1, P, LL);

    P.transpose();
    double Np = this->getNP(P);
    this->getMu(this->Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);

    McDMatrix<double> XHat;
    this->getHat(this->Xc, muX, XHat);
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD_CheckGetMuRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1();
    this->params.w = 0.1;
    McDMatrix<double> muX;
    McDMatrix<double> P;
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 1, 0.1, 0.1, P, LL);

    double Np = this->getNP(P);
    P.transpose();
    this->getMu(this->Xc, P, Np, muX);
    EXPECT_EQ(muX.nRows(), 2);
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD1NoDir_CheckAlignRunsR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1NoDir();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    this->align(R, s, t, Rd);
    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        std::cout << "\n point: " << i << "\n";
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD2_CheckAlignRunsR_E4MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata2();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    this->align(R, s, t, Rd);
    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD1_CheckAlignRunsR_E5MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1();
    this->params.w = 0.1;
    McDVector<double> t;
    double s;
    McDMatrix<double> R, Rd;

    this->align(R, s, t, Rd);
    McDMatrix<double> shiftedYc;
    this->shiftYCoords(this->Yc, R, s, t, shiftedYc);
    std::cout << "\n shifted: \n";
    for (int i = 0; (i < shiftedYc.nRows()) && (i < this->Xc.nRows()); i++) {
        EXPECT_NEAR(shiftedYc[i][0], this->Xc[i][0], 0.01);
        EXPECT_NEAR(shiftedYc[i][1], this->Xc[i][1], 0.01);
    }
}

TEST_F(mtalign__CPDLinearAlignerTest, givenTestCPD1_CheckPOKR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1();
    McDMatrix<double> P;
    double LL;
    this->computeP(this->Xc, this->Xd, this->Yc, this->Yd, 1, 0.1, 0.1, P, LL);
    for (int i = 0; i < P.nRows(); i++) {
        for (int j = 0; j < P.nCols(); j++)
            std::cout << P[i][j] << " ";
        std::cout << "\n";
    }
}

TEST_F(mtalign__CPDLinearAlignerTest, givenNothing_CheckInitCPDOKR_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    initCPDdata1();
    EXPECT_EQ(this->Xc.nRows(), 3);
    EXPECT_EQ(this->Xc.nCols(), 2);
    EXPECT_EQ(this->Xd.nRows(), 3);
    EXPECT_EQ(this->Xd.nCols(), 3);
    EXPECT_EQ(this->Yc.nRows(), 4);
    EXPECT_EQ(this->Yc.nCols(), 2);
    EXPECT_EQ(this->Yd.nRows(), 4);
    EXPECT_EQ(this->Yd.nCols(), 3);
}

TEST_F(mtalign__CPDLinearAlignerTest, givenMtrices_lapackSolvesCorrect_E3MS) {
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
    ma::solve(A, B, result_blas);
    for (int i = 0; i < result.nEntries(); i++) {

        EXPECT_NEAR(result_blas.dataPtr()[i], result.dataPtr()[i], 1.e-7);
    }
    EXPECT_EQ(result_blas.nRows(), result.nRows());
    EXPECT_EQ(result_blas.nCols(), result.nCols());
}

TEST_F(mtalign__CPDLinearAlignerTest, givenMtrices_blasMultsCorrect_E3MS) {
    TestingDevNullRedirect silentout(stdout);

    double aArray[6] = { 1, 2, 3, 4, 5, 6 };
    McDMatrix<double> A(2, 3, aArray);
    double bArray[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    McDMatrix<double> B(3, 3, bArray);
    double cArray[6] = { 1, 2, 3, 4, 5, 6 };
    McDMatrix<double> C(2, 3, cArray);
    McDMatrix<double> result = A * B + C;
    McDMatrix<double> result_blas;
    ma::multMatrices(A, B, C, 1.0, 1.0, result_blas);
    for (int i = 0; i < result_blas.nEntries(); i++) {

        EXPECT_NEAR(result_blas.dataPtr()[i], result.dataPtr()[i], 1.e-7);
    }
    EXPECT_EQ(result_blas.nRows(), result.nRows());
    EXPECT_EQ(result_blas.nCols(), result.nCols());
}
