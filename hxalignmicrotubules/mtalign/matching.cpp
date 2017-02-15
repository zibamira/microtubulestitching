#include <hxalignmicrotubules/mtalign/matching.h>

#include <mclib/McDArray.h>
#include <mclib/McMath.h>
#include <mclib/McMat4.h>
#include <mclib/McVec3.h>
#include <mclib/McVec3.h>
#include <mclib/McRot.h>
#include <mclib/internal/McAssert.h>

#include <hxalignmicrotubules/mtalign/matchingClique.h>
#include <hxalignmicrotubules/mtalign/matchingExact.h>
#include <hxalignmicrotubules/mtalign/matchingGreedy.h>
#include <hxalignmicrotubules/mtalign/fitTransform.h>

namespace ma = mtalign;

#ifndef NDEBUG
static bool hasZShift(const McMat4f& mat) {
    McVec3f scale, translation;
    McRotation rotation, scaleOrientation;
    mat.getTransform(translation, rotation, scale, scaleOrientation);
    return fabs(translation.z) > 1.e-2;
}
#endif

ma::Matching ma::matchingDirect(const ma::FacingPointSets& pts,
                                McMat4f startTransformation,
                                const ma::PGMPairWeightsParams& weightConfig,
                                const ma::MatchingAlgorithm algo) {
    mcassert(!hasZShift(startTransformation));
    ma::Matching matching;
    switch (algo) {
    case ma::MA_EXACT:
        matching = ma::matchingExact(pts, startTransformation, weightConfig);
        break;
    case ma::MA_GREEDY:
        matching = ma::matchingGreedy(pts, startTransformation, weightConfig);
        break;
    }
    return matching;
}

/// Checks whether the transform is acceptable (e.g. no flipping allowed).
static bool isTransformationOK(const McMat4f& matrix,
                               const McDArray<McVec3f>& transVertices) {
    // At least 3 points used to find unambiguous matrix and no flipping in z.
    return (transVertices.size() >= 3) && (matrix[2][2] > 0.0f);
}

ma::Matching
ma::matchingTwoStepBest(const ma::FacingPointSets& pts,
                        const ma::MatchingParams& params,
                        const ma::PGMPairWeightsParams& weightConfig,
                        const ma::MatchingAlgorithm algo) {

    const McDArray<McMat4f> startTransformations =
        ma::matchingCliqueTransforms(pts, params);

    ma::Matching greatestMatching;
    for (int i = 0; i < startTransformations.size(); ++i) {
        ma::Matching matching =
            matchingDirect(pts, startTransformations[i], weightConfig, algo);

        const McMat4f tf =
            fitTransformRigid(pts, matching, params.transformType);
        mcassert(!hasZShift(tf));
        if (!isTransformationOK(tf, pts.trans.positions)) {
            continue;
        }
        if (matching.matchedRefPointIds.size() <=
            greatestMatching.matchedRefPointIds.size()) {
            continue;
        }
        greatestMatching = matching;
    }

    return greatestMatching;
}
