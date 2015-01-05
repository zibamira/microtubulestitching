# Package Structure

The supplementary software uses several packages in `zib-amira`.

The first group of packages, which are contained in this folder, only depends on
each other and can be easily made available together:

 - `hxalignmicrotubules`: Algorithms described in the paper.
 - `dai`: Inference in graphical models.  It is a slightly modified version of
   the upstream at `git://git.tuebingen.mpg.de/libdai.git`.  We should probably
   change how we track the upstream: we should import entires tree snapshots in
   order to avoid tainting our history with authors that we don't know; or we
   could track libdai as a submodule that contains a fork of libdai with our
   modifications and fixes.

The first group of packages depends on a second group of packages that are used
by other unrelated packages in zib-amira.git:

 - `pointmatching`: Clique-based alignment.  Also used by `hxalignspatialgraph`.
 - `hxgraphalgorithms`: Also used by various other packages.
 - `ipopt`: Interior point optimization library from
   <https://projects.coin-or.org/Ipopt>.

# Classes, functions, call paths, ..

Some of the code has been cleaned up, restructured and documented.  Such code is
usually found in the namespace `mtalign`.  Use the following command to create
a local Doxygen documentation in `product/share/devreflocal/index.html`; then
start from the namespace documentation:

    tclsh src/amira/devtools/genDevRef.scro --preset ZIBAmira \
        --dirs src/zib-amira/microtubulestitching/hxalignmicrotubules/

`class QxMicrotubuleAlignSpatialGraphTool` contains mainly boilerplate code to
map GUI events to the algorithm classes.  It uses `HxManualMTAlign` to implement
interactive manipulation of slice transformations.

The enum `MicrotubuleSpatialGraphAligner::CPDType` describes coherent point
drift transformation types.

`class MicrotubuleSpatialGraphAligner` contains a mixture of things.  Primarily,
it receives calls from the GUI class and calls to the low-level algorithms.  But
it also includes some algorithmic details.

Algorithm 1 'linear alignment from orientation' and Algorithm 2 'linear
alignment from position and orientation' are implemented in `class
CoherentPointDriftRigidFisherMises`.  Algorithm 3 'elastic alignment' is
implemented in `class CoherentPointDriftNLFisherMises`.  The algorithms are now
available through `mtalign::cpd()`.

Call path from the GUI to `CoherentPointDriftNLFisherMises` (Algorithm 3):

 - Qt `ui.applyCPD SIGNAL(clicked())` is processed in
   `QxMicrotubuleAlignSpatialGraphTool::applyCPD()`, which calls
   `MicrotubuleSpatialGraphAligner::warpAll()`.

 - `MicrotubuleSpatialGraphAligner::warpAll()`: Iterates over slices.  Calls
   `warpSlices()` for each pair.  Collects results.  Stores sigma, kappa, ... in
   a spread sheet.  Applies transformations.  Non-linear transform will be
   applied only if there is no previous non-linear transform.  There is no
   option to restrict CPD processing to a section pair.

 - `warpSlices()`: Calls `mPointRepresentationCreatorForMicrotubules` to produce
   point representation from spatial graph.  Dispatches to specific algorithm:
   `warpSlicesRigid()` for Algorithm 1 and 2; `cpdElastic()` for Algorithm 3.

 - `mtalign::cpd()`: Instantiates `CoherentPointDriftNLFisherMises`.  Configures
   it.  Calls `align()` on it to execute Algorithm 3.  Retrieves result and
   stores them as `WarpResult`.

The class `mtalign::SliceSelector` partitions a `HxSpatialGraph` into slices
based on the value of an attribute.

The functions in `mtalign/project.h` compute feature points from an
`HxSpatialGraph`.  `mtalign::projectEndPoints()` implements the different ways
how to project points to the plane (see `projectionType` and
`MicrotubuleSpatialGraphAligner::ProjectionTypes`).

Functions in `mtalign/matching.h` depend on the package `pointmatching`.  They
provide the building blocks for clique-based alignment.  The building blocks are
used by `MicrotubuleSpatialGraphAligner`.  The algorithms from `pointmatching`
expect arguments of the interface type `PointRepresentation`, which is
implemented as small helper classes in the implementation files
`mtalign/matching*.cpp`.  The following algorithms from `pointmatching` are
used:

    StartTransformationGenerator3d
    PointMatchingScoringFunction
    GreedyPointMatchingAlgorithm
    PointMatchingDataStruct
    ExactPointMatchingAlgorithm

`Auto-Align` always computes a matching and then uses the matched points to
compute a transform unless `None` is specified in `Transform`.  The matching
type is controlled as follows:

 - `Compute transform = Initial`: Clique-based initialization,
   `matchingCliqueTransforms()`, followed by `matchingExact()` or
   `matchingGreedy()` as specified in `PM algorithm`.   `PM algorithm = PGM` is
   invalid.

 - `Compute transform = Optimum`: `matchingExact()`, `matchingGreedy()`, or
   `MicrotubulePGMPointMatcher` as specified in `PM algorithm`.

The call path is as follows:

 - `alignAllOrPair()`: call `alignPairAndTestScaling()`, collect matrices, apply
   transforms.

 - `alignPairAndTestScaling()`: for each scale, for each gap, `align()`; find
   best transform (based on `score`); fill spreadsheet; set attributes for
   `matchedRefPointIds` and `matchedTransPointIds`; return matrix.

 - `align()`: either call `alignNonPGM()` or `alignPGM()`.

 - `alignNonPGM()`: setup parameters and call functions `mtalign::matching*()`.

 - `MicrotubulePGMPointMatcher`: compute PGM matching; set PGM-specific
   attributes.

`DerivativesForRigidRegistration` contains the Q derivatives from the supporting
information.

Be careful to maintain compatibility with the supplementary data deposited at
Dryad <http://dx.doi.org/10.5061/dryad.v8j20>.  It contains Amira scripts that
use the script object `TryStitchingParametersAndIterate.scro` and the following
classes:

    HxCPDSpatialGraphWarp
    HxComparePointMatchings
    HxIteratePointMatchingUntilConvergence
    HxRotateSpatialGraphStackSliceAndCDP
    HxTestPointMatching
    SpreadSheetWrapper

## Transformations

`Mat4f` transforms are applied to the points and multiplied to matrices that are
already stored in the parameters of the spatial graph.
`MicrotubuleSpatialGraphAligner::applyTransform()` creates a selection of the
relevant slices (either all above the current pair or a specified number of
slices).  It then constructs a `MicrotubuleTransformOperation`, which transforms
the selected points and multiplies the matrix to the matrices stored in the
parameters (see `MicrotubuleTransformOperation::appendTransform()`).

Points are transformed forward and backwards when searching parameters in
`MicrotubuleSpatialGraphAligner::alignPairAndTestScaling()`.  Finally, the
points should be back to their original position.  This should perhaps be
changed to keep the original spatial graph stack unmodified during a search and
only apply the final transform.

`MovingLeastSquares` can be applied only once.  See
`MicrotubuleSpatialGraphAligner::warpAll()`.

# Testing data

`hxalignmicrotubules` uses testing data from `hxalignspatialgraph`.  The
following files might be relevant when testing the alignment:

    fullp0p1.am
    sgs_align.am
    sgs_aligned.am
    sgs_many-pairs.am
    spatialgraphstack_2slices_1kvertices_upl.am
    spatialgraphstack_2slices_1kvertices_upl-sn.am
    spatialgraphstack_3slices_2kvertices.am
    spatialgraphstack_3slices_2kvertices_matchedGreedy_withGap.am
    spatialgraphstack_3slices_2kverticesWithLabels.am

# Changelog

## June 2014

The normalization factor has been changed [^c4c577] to match the definition of
the Fisher-Mises distribution.  The change caused minor differences in the test
results (see below).  The differences should not matter in practice.

    CoherentPointDriftNLFisherMisesAccTest.cpp:161:
    Value of: cpd.mKappaOut
         New: is < 90.903
         Old: 94.4124 (of type double)
    CoherentPointDriftNLFisherMisesAccTest.cpp:163:
    Value of: cpd.mEOut
         New: is < 1492.2
         Old: 1966.41 (of type double)
    CoherentPointDriftNLFisherMisesAccTest.cpp:164:
    Value of: cpd.mNumIterationsOut
         New: 45
         Old: 42

[^c4c577]: zib-amira:c4c57770265640acdcd47e7e025a8f8b9ecdaaf2
    CoherentPointDriftNLFisherMises: Fix fisherMises() normalization factor
