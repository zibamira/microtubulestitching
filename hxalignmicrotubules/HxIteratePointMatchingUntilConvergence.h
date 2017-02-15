#pragma once

#include <hxcore/HxCompModule.h>
#include <hxcore/HxConnection.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntTextN.h>
#include <hxcore/HxPortMultiMenu.h>

#include <hxalignmicrotubules/api.h>

class HXALIGNMICROTUBULES_API HxIteratePointMatchingUntilConvergence
    : public HxCompModule {
    HX_HEADER(HxIteratePointMatchingUntilConvergence);

  public:

    void compute();

    int parse(Tcl_Interp* t, int argc, char** argv);

    void update();

    HxPortMultiMenu portEvidenceHeuristic;
    HxPortIntTextN portNumIterations;
    HxConnection connectionToFEModule;
    HxConnection connectionToPMEvalModule;
    HxPortDoIt mDoIt;

  private:
    void updateEvidenceToAssignPort();
    int getNumberOfPointsToAssign();
    void addEvidenceForNodes(const McDArray<int>& nodesToAssign);
    void getListOfNeededEvidenceNodes(McDArray<int>& nodesToAssign);
    void assignNeededEvidence();
    void addEvidence(const int node1, const int node2);
};
