#pragma once

#include <mclib/internal/McDMatrix.h>

namespace mtalign {

/// `CPDLinearAlignerDerivatives` implements the Q derivatives (see supporting
/// information [Weber 2014]) and is used internally by `CPDLinearAligner`.
class CPDLinearAlignerDerivatives {
  public:
    McDMatrix<double> A, BT, C, DT, E;
    double s;
    double sigmaSquare;
    double rho;
    double kappa;
    McDMatrix<double> R, RdRho, RdRhodRho;
    double Np;

    double trace(const McDMatrix<double>& mat) {
        double trace = 0.0;
        for (int i = 0; i < mat.nCols(); i++)
            trace += mat[i][i];
        return trace;
    }

    void initCurParams(const McDVector<double>& at) {
        s = at[0];
        rho = at[1];
        kappa = at[2];
        sigmaSquare = at[3];

        R.resize(2, 2);
        R[0][0] = cos(rho);
        R[0][1] = -sin(rho);
        R[1][0] = sin(rho);
        R[1][1] = cos(rho);

        RdRho.resize(2, 2);
        RdRho[0][0] = -sin(rho);
        RdRho[0][1] = -cos(rho);
        RdRho[1][0] = cos(rho);
        RdRho[1][1] = -sin(rho);

        RdRhodRho.resize(2, 2);
        RdRhodRho[0][0] = -cos(rho);
        RdRhodRho[0][1] = sin(rho);
        RdRhodRho[1][0] = -sin(rho);
        RdRhodRho[1][1] = -cos(rho);
    }

    double getInitialSigmaSquare() {
        McDMatrix<double> secondTerm = BT;
        secondTerm *= R;
        double valueTerm1 =
            (trace(A) - 2.0 * s * trace(secondTerm) + s * s * trace(C)) / 2.0;
        double valueTerm2 = Np;
        return valueTerm1 / valueTerm2;
    }

    void initStartValues(McDVector<double>& at) {
        at[3] = getInitialSigmaSquare();
    }

    double getFuncValue(McDVector<double>& at) {
        initCurParams(at);
        McDMatrix<double> BTR = BT;
        BTR *= R;
        mcassert(kappa > 0);
        McDMatrix<double> DTR = DT;
        DTR *= R;
        return 1.0 / (2.0 * sigmaSquare) *
                   (trace(A) - 2.0 * s * trace(BTR) + s * s * trace(C)) +
               Np * log(sigmaSquare) - kappa * (trace(DTR) + trace(E)) -
               Np * log(kappa / (2.0 * M_PI * (exp(kappa) - exp(-kappa))));
    }

    double dS() {
        mcassert(s > 0.0);
        McDMatrix<double> firstTerm = BT;
        firstTerm *= R;
        return 1.0 / sigmaSquare * (-trace(firstTerm) + s * trace(C));
    }

    double dRho() {
        McDMatrix<double> firstTerm = BT;
        firstTerm *= RdRho;
        McDMatrix<double> secondTerm = DT;
        secondTerm *= RdRho;
        return -s / sigmaSquare * trace(firstTerm) - kappa * trace(secondTerm);
    }

    double dKappa() {
        mcassert(kappa > 0);
        McDMatrix<double> firstTerm = DT;
        firstTerm *= R;
        return -trace(firstTerm) - trace(E) -
               Np * (1.0 / kappa - 1.0 / tanh(kappa));
    }

    double dSigmaSquare() {
        McDMatrix<double> secondTerm = BT;
        secondTerm *= R;
        return -1.0 / (2.0 * sigmaSquare * sigmaSquare) *
                   (trace(A) - 2.0 * s * trace(secondTerm) + s * s * trace(C)) +
               Np / sigmaSquare;
    }

    double dSdS() { return 1.0 / sigmaSquare * trace(C); }

    double dSdRho() {
        McDMatrix<double> firstTerm = BT;
        firstTerm *= RdRho;
        return -1.0 / sigmaSquare * (trace(firstTerm));
    }

    double dSdKappa() { return 0; }

    double dSdSigmaSquare() {
        McDMatrix<double> firstTerm = BT;
        firstTerm *= R;
        return 1.0 / (sigmaSquare * sigmaSquare) *
               (trace(firstTerm) - s * trace(C));
    }

    double dRhodRho() {
        McDMatrix<double> firstTerm = BT;
        firstTerm *= RdRhodRho;
        McDMatrix<double> secondTerm = DT;
        secondTerm *= RdRhodRho;
        return -s / sigmaSquare * trace(firstTerm) - kappa * trace(secondTerm);
    }

    double dRhodKappa() {
        McDMatrix<double> firstTerm = DT;
        firstTerm *= RdRho;
        return -trace(firstTerm);
    }

    double dRhodSigmaSquare() {
        McDMatrix<double> firstTerm = BT;
        firstTerm *= RdRho;
        return s / (sigmaSquare * sigmaSquare) * trace(firstTerm);
    }

    double dKappadKappa() {
        return -Np *
               (-1.0 / (kappa * kappa) + 1.0 / (sinh(kappa) * sinh(kappa)));
    }

    double dKappadSigmaSquare() { return 0.0; }

    double dSigmaSquaredSigmaSquare() {
        McDMatrix<double> secondTerm = BT;
        secondTerm *= R;
        return 1.0 / (sigmaSquare * sigmaSquare * sigmaSquare) *
                   (trace(A) - 2.0 * s * trace(secondTerm) + s * s * trace(C)) -
               Np / (sigmaSquare * sigmaSquare);
    }

    void gradient(const McDVector<double>& at, McDVector<double>& gradient) {
        gradient.resize(at.size());
        initCurParams(at);
        McDVector<double> atB, atF;
        atB = at;
        atF = at;
        double eps = 1.e-10;
        atB[1] -= eps;
        atF[1] += eps;

        gradient[0] = dS();
        gradient[1] = dRho();
        atB = at;
        atF = at;
        atB[2] -= eps;
        atF[2] += eps;
        gradient[2] = dKappa();
        gradient[3] = dSigmaSquare();
    }

    void hessian(const McDVector<double>& at, McDMatrix<double>& hessian) {
        hessian.resize(at.size(), at.size());
        initCurParams(at);
        hessian[0][0] = dSdS();
        hessian[0][1] = dSdRho();
        hessian[1][0] = hessian[0][1];
        hessian[0][2] = dSdKappa();
        hessian[2][0] = hessian[0][2];
        hessian[0][3] = dSdSigmaSquare();
        hessian[3][0] = hessian[0][3];
        hessian[1][1] = dRhodRho();
        hessian[1][2] = dRhodKappa();
        hessian[2][1] = hessian[1][2];
        hessian[1][3] = dRhodSigmaSquare();
        hessian[3][1] = hessian[1][3];
        hessian[2][2] = dKappadKappa();
        hessian[2][3] = dKappadSigmaSquare();
        hessian[3][2] = hessian[2][3];
        hessian[3][3] = dSigmaSquaredSigmaSquare();
    }
};

}  // namespace mtalign
