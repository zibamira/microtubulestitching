# Microtubule filament editor MicrotubuleAlignTool plugin calls expected algorithms

 - Purpose: Verify that the GUI of the microtubule aligner plugin dispatches to
   the expected low-level functions.
 - Preconditions: A spatial graph stack for testing;
   `data/spatialgraph/fullp0p1.am` should work.
 - Platforms: LinuxAMD64, Win64VC9.

## Steps

 1. Execute an alignment with transform type `Rigid fisher mises` as explained
    in the step 'Registration (rigid)' of the protocol 'Non-rigid microtubule
    alignment'[^protocol].

    Expected:

    - A transformation should be computed.
    - The result spreadsheet should contain reasonable values.
    - For `fullp0p1.am`, the computation takes minutes (optimize build).

 2. Execute an alignment with transform type `non-linear` as explained in the
    step 'Registration (non-linear)' of the protocol 'Non-rigid microtubule
    alignment'[^protocol].

    Expected:

    - A transformation should be computed.
    - The result spreadsheet should contain reasonable values.
    - For `fullp0p1.am`, the computation takes less than a minute (optimize
      build).


[^protocol]: Repo `2013-04_microtubuli_non-rigid-alignment-protocol`, also
    available at <http://www.zib.de/kratz/hidden/2014/XfnsdjE/latest/non-rigid-alignment-protocol/latest/Non-rigid-alignment-protocol.html>
