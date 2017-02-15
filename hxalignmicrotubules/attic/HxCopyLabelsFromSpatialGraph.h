#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxConnection.h>
#include <hxcore/HxPortMultiMenu.h>

#include <hxalignmicrotubules/api.h>

class HXALIGNMICROTUBULES_API HxCopyLabelsFromSpatialGraph
    : public HxCompModule {

    HX_HEADER(HxCopyLabelsFromSpatialGraph);

  public:

    void compute();
    /// for tcl commands
    int parse(Tcl_Interp* t, int argc, char** argv);

    void update();

    HxPortMultiMenu portPairLabel;

    HxConnection portOtherSG;

    HxPortDoIt mDoIt;

    void updatePairLabelPort();
};
