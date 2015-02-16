#pragma once

namespace mtalign {

class IPOPTForCPD;

/// `createIPOPTForCPDLinearAligner()` returns an instance of `IPOPTForCPD`
/// that implements derivatives that include scaling for Ipopt.  It is used
/// internally by `CPDLinearAligner`.
IPOPTForCPD* createIPOPTForCPDLinearAligner();

}  // namespace mtalign
