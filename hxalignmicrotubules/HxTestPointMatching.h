#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFloatTextN.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortToggleList.h>

#include <hxalignmicrotubules/api.h>

class HXALIGNMICROTUBULES_API HxTestPointMatching : public HxCompModule {
    HX_HEADER(HxTestPointMatching);

  public:
    HxTestPointMatching();

    ~HxTestPointMatching();

    void compute();

    void update();

    int parse(Tcl_Interp* t, int argc, char** argv);

    HxPortMultiMenu portMatchingAlgorithmType;
    HxPortIntTextN portProjectionType;
    HxPortFloatTextN portThresholds;
    HxPortFloatTextN portParams;
    HxPortToggleList portUseParams;
    HxPortFloatTextN portPGMPairParam;
    HxPortFloatTextN portPGMDummieSignificance;
    HxPortFloatTextN portPairDummy;
    HxPortDoIt mDoIt;
};
