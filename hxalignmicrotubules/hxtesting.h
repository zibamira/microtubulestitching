#pragma once

template <typename T> class McHandle;
class HxSpatialGraph;

namespace hxtesting {

/// `makeSpatialGraphThreeStackedSections()` creates an `HxSpatialGraph` that
/// has three stacked sections as indicated by the int attribute `section` with
/// values 1, 2, 3.  The sections are ordered, meaning that all nodes of each
/// section have consecutive node indices.  The spatial graph has a second int
/// attribute `unordered` that cannot be used to partition the graph into
/// sections, because it is not ordered.
McHandle<HxSpatialGraph> makeSpatialGraphThreeStackedSections();

}  // hxtesting
