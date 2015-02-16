#pragma once

namespace mtalign {

class IPOPTForCPD;

/// `createIPOPTForCPDLinearAlignerWithoutScale()` returns an instance of
/// `IPOPTForCPD` that implements derivatives without scaling for Ipopt.  It is
/// used internally by `CPDLinearAligner`.
IPOPTForCPD* createIPOPTForCPDLinearAlignerWithoutScale();

}  // namespace mtalign
