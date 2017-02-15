#pragma once

#include <mclib/internal/McDMatrix.h>

namespace mtalign {

/// `rotationAngle2d()` computes the angle (in the range `[-pi, pi]`) for a
/// two-dimensional rotation from the top-left 2x2 sub-matrix of `R`.
double rotationAngle2d(const McDMatrix<double>& R);

}  // namespace mtalign
