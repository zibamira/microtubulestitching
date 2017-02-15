# Introduction

This document contains various technical details about the microtubule stitching
implementation that is described in [Weber 2014].

# Package structure

The supplementary software uses several packages in `zib-amira`.

The first group of packages, which are contained in this folder, only depends on
each other.  They contain the main algorithms from [Weber 2014] an can be easily
made available outside ZIB:

 - `hxalignmicrotubules`: Algorithms described in the paper.
 - `dai`: Inference in graphical models.  It is a slightly modified version of
   the upstream from `git://git.tuebingen.mpg.de/libdai.git`.  See below for
   details about the git history.

`hxalignmicrotubules` depends on a second group of packages that are used by
other unrelated packages in zib-amira.git:

 - `pointmatching`: clique-based alignment.  Also used by `hxalignspatialgraph`.
 - `hxgraphalgorithms`: used by various other packages.
 - `ipopt`: Interior point optimization library from
   <https://projects.coin-or.org/Ipopt>.

# Classes, functions, call paths, ..

## mtalign

Code that has been cleaned up, restructured and documented is usually found in
the namespace `mtalign`.  Use the following command to create a local Doxygen
documentation in `product/share/devreflocal/index.html`; then start from the
namespace documentation:

    tclsh src/amira/devtools/genDevRef.scro --preset ZIBAmira \
        --dirs src/zib-amira/microtubulestitching/hxalignmicrotubules/

    open product/share/devreflocal/hxalignmicrotubules/namespacemtalign.html

The folder `hxalignmicrotubules/examples/` and the Doxygen documentation contain
an example how to use some of the low-level functions that implement the main
algorithms of [Weber 2014] and are mostly independent of Amira.

## Amira integration

This section contains some unsorted notes about the integration of the low-level
functions with Amira.

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

`class QxMicrotubuleAlignSpatialGraphTool` contains mainly boilerplate code to
map GUI events to the algorithm classes.  It uses `HxManualMTAlign` to implement
interactive manipulation of slice transformations.

The enum `MicrotubuleSpatialGraphAligner::CPDType` describes coherent point
drift transformation types.

`class MicrotubuleSpatialGraphAligner` is a mixed bag.  Primarily, it receives
calls from the GUI class and calls down to the low-level algorithms.  But it
also includes some algorithmic details.

Algorithm 1 'linear alignment from orientation' and Algorithm 2 'linear
alignment from position and orientation' are implemented in `class
CPDLinearAligner`.  Algorithm 3 'elastic alignment' is implemented in `class
CPDElasticAligner`.  The algorithms are available through `mtalign::cpd()`.

Call path from the GUI to `CPDElasticAligner` (Algorithm 3):

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

 - `mtalign::cpd()`: Instantiates `CPDElasticAligner`.  Configures it.  Calls
   `align()` on it to execute Algorithm 3.  Retrieves result and stores them as
   `WarpResult`.

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

# Git history

## libdai

In the past, we imported the full git history from the upstream into the git
repository `zib-amira.git` at ZIB.  We should consider changing this in the
future and import only entire tree snapshots in order to avoid tainting our
history with authors that we don't know; or we could track libdai as a submodule
that contains a fork of libdai with our modifications and fixes.

## How to export the microtubulestitching source to a separate git repo?

The sub-directory `zib-amira/microtubulestitching` is exported to a separate git
repo.  The full history is maintained in `zib-amira`.  The export is done as
squashed commits from release tags `zibamira-*.*` using the low-level `git
commit-tree`.

Prepare a squashed commit with:

```bash
./export-to-github.sh
```

Then review and push to GitHub.

# Useful commands

You can call the GitHub API to test how GitHub will translate a markdown file:

    file=README-technical.md
    curl --data-binary @$file -H "Content-Type:text/plain" -s \
        https://api.github.com/markdown/raw >${file}.html

# Changelog

2015-02: The sign of the threshold for filtering lines that are nearly parallel
to the section boundary (parameter `Remove angel lower` in GUI) has been fixed
[8f93cb].  The change caused minor differences in the test results that should
not matter in practice.

 - [8f93cb]: zib-amira@8f93cbf9ffa9e9808562677e74eef2a1f3c56191 'mtalign: fix
   sign in angleToPlaneFilter; polish EndPointParams 2/n'


2015-02: The Fisher-Mises normalization factor for the linear CPD has been
changed to match the definition in the paper [d04065].  The change did not
affect any tests.

 - [d04065]: zib-amira@d040657d21079cf70bcedc9cb21a9c8947d01ac0 'mtalign: change
   Fisher-Mises normalization factor in linear CPD to match paper'


2014-06: The normalization factor for the elastic CPD has been changed to match
the definition of the Fisher-Mises distribution [c4c577].  The change caused
minor differences in the test results (see below).  The differences should not
matter in practice.

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

 - [c4c577]: zib-amira:c4c57770265640acdcd47e7e025a8f8b9ecdaaf2
   CoherentPointDriftNLFisherMises: Fix fisherMises() normalization factor
