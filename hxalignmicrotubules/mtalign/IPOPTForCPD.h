#pragma once

#include <coin/IpTNLP.hpp>
#include <hxalignmicrotubules/mtalign/CPDLinearAlignerDerivatives.h>
#include <mclib/internal/McDVector.h>

namespace mtalign {

/// `IPOPTForCPD` is the base class to implement derivatives objects for Ipopt.
/// It is used internally by `CPDLinearAligner`.
class IPOPTForCPD : public Ipopt::TNLP {
  public:
    double s;
    double rho;
    double kappa;
    double sigmaSquare;
    CPDLinearAlignerDerivatives gradAndHessAndFunc;
    McDVector<double> resultValues;
};

}  // namespace mtalign
