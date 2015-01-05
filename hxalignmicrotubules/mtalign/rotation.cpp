#include <hxalignmicrotubules/mtalign/rotation.h>

namespace ma = mtalign;

double ma::rotationAngle2d(const McDMatrix<double>& R) {
    double rho = 0.0;
    const double rhoFromSin = asin(R[1][0]);
    const double rhoFromCos = acos(R[0][0]);
    const double sinRho = R[1][0];
    const double cosRho = R[0][0];
    if (sinRho > 0.0 && cosRho > 0.0) {
        rho = rhoFromSin;
    }
    if (sinRho > 0.0 && cosRho <= 0.0) {
        rho = rhoFromCos;
    }
    if (sinRho <= 0.0 && cosRho <= 0.0) {
        rho = -rhoFromSin + M_PI;
    }
    if (sinRho <= 0.0 && cosRho > 0.0) {
        rho = -rhoFromCos + 2.0 * M_PI;
    }
    if (rho > M_PI)
        rho = -(2.0 * M_PI - rho);
    if (rho < -M_PI)
        rho = -(-2.0 * M_PI - rho);
    return rho;
}
