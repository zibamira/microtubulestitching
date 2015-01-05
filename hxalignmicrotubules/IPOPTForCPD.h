#pragma once

#include <coin/IpTNLP.hpp>
#include <hxalignmicrotubules/DerivativesForRigidRegistration.h>
#include <mclib/McDVector.h>

class IPOPTForCPD : public Ipopt::TNLP {
  public:
    double s;
    double rho;
    double kappa;
    double sigmaSquare;
    DerivativesForRigidRegistration gradAndHessAndFunc;
    McDVector<double> resultValues;
};
