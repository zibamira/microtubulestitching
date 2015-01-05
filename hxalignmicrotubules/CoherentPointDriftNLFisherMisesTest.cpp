#include <hxalignmicrotubules/CoherentPointDriftNLFisherMises.h>

#include <stdio.h>
#include <iostream>

#include <gtest/gtest.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <mclib/McMat3f.h>

static void checkSameCoord(const McDVector<double> coord1,
                           const McDVector<double> coord2) {
    float eps = 2.e-2;
    for (int i = 0; i < coord1.size(); i++) {
        ASSERT_NEAR(coord1[i], coord2[i], eps);
    }
}

static void dontPrint(const char* msg) {
    // Use the following line to temporarily enable printing.
    // printf("%s", msg);
}

static CoherentPointDriftNLFisherMises makeCPDdataNLVMRotated() {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
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
    cpd.convertCoordsToMatrix(xc, cpd.xs);
    cpd.convertCoordsToMatrix(yc, cpd.ys);
    cpd.convertDirectionsToMatrix(xd, cpd.xDirs);
    cpd.convertDirectionsToMatrix(yd, cpd.yDirs);
    return cpd;
}

TEST(CoherentPointDriftNLFisherMises,
     givenPoints_runCDPNLFMAndCheckResult_E3MS) {
    CoherentPointDriftNLFisherMises cpd = makeCPDdataNLVMRotated();
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;

    cpd.align(G, W, pairs);
    McDMatrix<double> shifted;
    cpd.shiftYs(cpd.ys, G, W, shifted);
    checkSameCoord(cpd.xs.getRowVector(0), shifted.getRowVector(0));
    checkSameCoord(cpd.xs.getRowVector(1), shifted.getRowVector(1));
    checkSameCoord(cpd.xs.getRowVector(2), shifted.getRowVector(2));
}

TEST(CoherentPointDriftNLFisherMises,
     givenAssymPoints43_checkRuns_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));
    xsarray.append(McVec3f(2, 1.1, 0));
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(3, 3, 0));
    ysarray.append(McVec3f(3, 2, 0));
    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);

    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;

    cpd.align(G, W, pairs);

    double res[9] = { 1,                   0.32047626849284516,
                      0.56610623428191031, 0.32047626849284516,
                      1,                   0.56610623428191031,
                      0.56610623428191031, 0.56610623428191031,
                      1 };
    McDMatrix<double> erg(3, 3, res);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(erg[i][j], G[i][j], 1.e-5);

    McDMatrix<double> shifted;
    cpd.shiftYs(cpd.ys, G, W, shifted);
    checkSameCoord(cpd.xs.getRowVector(0), shifted.getRowVector(0));
    checkSameCoord(cpd.xs.getRowVector(1), shifted.getRowVector(1));
    shifted[2][1] -= 0.05;
    checkSameCoord(cpd.xs.getRowVector(2), shifted.getRowVector(2));
}

TEST(CoherentPointDriftNLFisherMises, givenAssymPoints34_checkRuns_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    cpd.xs.resize(3, 3);
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));

    cpd.xs.resize(4, 3);
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(3, 3, 0));
    ysarray.append(McVec3f(3, 2, 0));
    ysarray.append(McVec3f(2, 1.1, 0));
    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;

    cpd.align(G, W, pairs);

    double res[9] = { 1,                   0.32047626849284516,
                      0.56610623428191031, 0.32047626849284516,
                      1,                   0.56610623428191031,
                      0.56610623428191031, 0.56610623428191031,
                      1 };
    McDMatrix<double> erg(3, 3, res);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(erg[i][j], G[i][j], 1.e-5);

    McDMatrix<double> shifted;
    cpd.shiftYs(cpd.ys, G, W, shifted);
    checkSameCoord(cpd.xs.getRowVector(0), shifted.getRowVector(0));
    checkSameCoord(cpd.xs.getRowVector(1), shifted.getRowVector(1));
    checkSameCoord(cpd.xs.getRowVector(2), shifted.getRowVector(2));
}

TEST(CoherentPointDriftNLFisherMises, givenPoints_testG_GUI_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    cpd.xs.resize(3, 2);
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));

    cpd.xs.resize(3, 2);
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(3, 3, 0));
    ysarray.append(McVec3f(3, 2, 0));

    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;
    cpd.align(G, W, pairs);

    double res[9] = { 1,                   0.34686364525689239,
                      0.58895130975055354, 0.34686364525689239,
                      1,                   0.58895130975055354,
                      0.58895130975055354, 0.58895130975055354,
                      1 };
    McDMatrix<double> erg(3, 3, res);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(erg[i][j], G[i][j], 1.e-5);
}

TEST(CoherentPointDriftNLFisherMises,
     givenPoints_testComputeP_GUI_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    cpd.xs.resize(3, 3);
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));

    cpd.xs.resize(3, 3);
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(3, 3, 0));
    ysarray.append(McVec3f(3, 2, 0));
    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;
    double initSigmaSquare =
        (1.0 / (2 * 3 * 3)) * cpd.sumSquaredDistances(cpd.xs, cpd.ys);
    double LL;
    cpd.computeP(cpd.xs, cpd.ys, initSigmaSquare, 0.1, cpd.params.w, P, LL);
    double res[9] = { 0.286168,  0.310922, 0.295566, 0.0358612, 0.155592,
                      0.0740157, 0.101303, 0.219948, 0.209085 };
    McDMatrix<double> erg(3, 3, res);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(erg[i][j], P[i][j], 1.e-5);
}

TEST(CoherentPointDriftNLFisherMises,
     givenPoints_runCDPAndCheckResult_GUI_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    cpd.xs.resize(3, 3);
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));

    cpd.xs.resize(3, 3);
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(3, 3, 0));
    ysarray.append(McVec3f(3, 2, 0));
    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    McDMatrix<double> P;
    cpd.align(G, W, pairs);
    McDMatrix<double> shifted;
    cpd.shiftYs(cpd.ys, G, W, shifted);
    checkSameCoord(cpd.xs.getRowVector(0), shifted.getRowVector(0));
    checkSameCoord(cpd.xs.getRowVector(1), shifted.getRowVector(1));
    checkSameCoord(cpd.xs.getRowVector(2), shifted.getRowVector(2));
}

TEST(CoherentPointDriftNLFisherMises, givenPoints_runCDP_GUI_E3MS) {
    CoherentPointDriftNLFisherMises cpd;
    cpd.setPrint(dontPrint);
    McDArray<McVec3f> xsarray, ysarray;
    cpd.xs.resize(3, 3);
    xsarray.append(McVec3f(1, 1, 0));
    xsarray.append(McVec3f(2, 2, 0));
    xsarray.append(McVec3f(2, 1, 0));

    cpd.xs.resize(3, 3);
    ysarray.append(McVec3f(1, 1, 0));
    ysarray.append(McVec3f(2, 2, 0));
    ysarray.append(McVec3f(2, 1, 0));
    cpd.convertCoordsToMatrix(xsarray, cpd.xs);
    cpd.convertCoordsToMatrix(ysarray, cpd.ys);
    cpd.params.beta = 1;
    cpd.params.lambda = 1;
    cpd.params.w = 0.1;
    McDMatrix<double> W;
    McDMatrix<double> G;
    McDArray<McVec2i> pairs;
    cpd.align(G, W, pairs);
    McDMatrix<double> shifted;
    cpd.shiftYs(cpd.ys, G, W, shifted);
    checkSameCoord(cpd.xs.getRowVector(0), shifted.getRowVector(0));
    checkSameCoord(cpd.xs.getRowVector(1), shifted.getRowVector(1));
    checkSameCoord(cpd.xs.getRowVector(2), shifted.getRowVector(2));
}
