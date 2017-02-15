#pragma once

#include <mclib/internal/McDMatrix.h>
#include <mclib/internal/McDVector.h>

namespace mtalign {

void multMatrices(const McDMatrix<double>& A, const McDMatrix<double>& B,
                  const McDMatrix<double>& C, const double alpha,
                  const double beta, McDMatrix<double>& result);

void multThreeMatrices(const McDMatrix<double>& A, const McDMatrix<double>& B,
                       const McDMatrix<double>& C, McDMatrix<double>& result);

double trace(const McDMatrix<double>& mat);

void solve(const McDMatrix<double>& A, const McDMatrix<double>& B,
               McDMatrix<double>& X);

double gauss(const McDVector<double>& mean, const double sigmaSquare,
             const McDVector<double> x);

double fisherMises(const McDVector<double>& mean, const double kappa,
                   const McDVector<double> x);

}  // namespace mtalign
