#pragma once

#include <hxalignmicrotubules/mtalign/Context.h>
#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>
#include <hxalignmicrotubules/mtalign/SliceSelector.h>
#include <hxalignmicrotubules/mtalign/cpd.h>
#include <hxalignmicrotubules/mtalign/data.h>
#include <hxalignmicrotubules/mtalign/fitTransform.h>
#include <hxalignmicrotubules/mtalign/matching.h>
#include <hxalignmicrotubules/mtalign/matchingClique.h>
#include <hxalignmicrotubules/mtalign/matchingExact.h>
#include <hxalignmicrotubules/mtalign/matchingGreedy.h>
#include <hxalignmicrotubules/mtalign/matchingPGM.h>
#include <hxalignmicrotubules/mtalign/matchingPGM_sg.h>
#include <hxalignmicrotubules/mtalign/project.h>
#include <hxalignmicrotubules/mtalign/rotation.h>

// Do not include the following files, because they are used only internally.
#if 0
#include <hxalignmicrotubules/mtalign/PGMMatcher.h>
#endif

/// Namespace `mtalign` contains functions for alignment of microtubules across
/// section boundaries as published in [Weber 2014].
///
/// Some functions are illustrated in the example below.
///
/// The endpoints on the two sides of a section boundary are stored as
/// `FacingPointSets`.
///
/// `projectEndPoints()` and variants compute `FacingPointSets` for two
/// consecutive sections from an `HxSpatialGraph` (as defined by the attribute
/// `TransformInfo`).
///
/// `matchingDirect()` and `matchingTwoStepBest()` compute matchings using the
/// lower-level functions `matchingCliqueTransforms()`, `matchingGreedy()`, and
/// `matchingExact()`.
///
/// `matchingPGM()` computes a special type of matching using a probabilistic
/// graphical model.
///
/// `fitTransform()` fits a transform to a matching.
///
/// `cpd()` and related functions implement variants of coherent point drift to
/// align consecutive sections.
///
/// **References**
///
///  - [Weber 2014] Britta Weber, Erin M. Tranfield, Johanna L. Hoeoeg, Daniel
///    Baum, Claude Antony, Tony Hyman, Jean-Marc Verbavatz, Steffen Prohaska
///    (2014).  Automated stitching of microtubule centerlines across serial
///    electron tomograms.  PLoS ONE, e113222.
///    <http://dx.doi.org/10.1371/journal.pone.0113222>.
///
/// **Example**
///
/// `hxalignmicrotubules/examples/HxCPDAlignerExample.cpp`:
///
/// \include HxCPDAlignerExample.cpp
///
namespace mtalign { }

/// \example hxalignmicrotubules/examples/HxCPDAlignerExample.cpp
///
/// `HxCPDAlignerExample` illustrates how to use the CPD functions in `mtalign`.
