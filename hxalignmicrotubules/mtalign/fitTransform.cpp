#include <hxalignmicrotubules/mtalign/fitTransform.h>

#include <hxcore/HxMessage.h>
#include <mclib/internal/McAlignPointSets.h>
#include <mclib/McDArray.h>
#include <mclib/McMat4.h>
#include <mclib/McVec3.h>
#include <mclib/internal/McDMatrix.h>

namespace ma = mtalign;

McMat4f ma::fitTransformRigid(const FacingPointSets& pts,
                              const Matching& matching, TransformType tfType) {
    mcassert(tfType == TF_RIGID || tfType == TF_RIGID_ISO_SCALE);

    McDArray<McVec3f> coords1Z0 = pts.ref.positions;
    McDArray<McVec3f> coords2Z0 = pts.trans.positions;
    for (int i = 0; i < coords1Z0.size(); i++)
        coords1Z0[i].z = 0;
    for (int i = 0; i < coords2Z0.size(); i++)
        coords2Z0[i].z = 0;

    // mcAlignPointSets expects transform type as integers:
    const int rigid = 0;
    const int rigidIsoScale = 1;

    McMat4f tf;
    mcAlignPointSets(tf, coords1Z0.dataPtr(), coords2Z0.dataPtr(),
                     matching.matchedRefPointIds, matching.matchedTransPointIds,
                     matching.matchedRefPointIds.size(),
                     tfType == TF_RIGID ? rigid : rigidIsoScale);
    return tf;
}

McMat4f ma::fitTransformAffine(const FacingPointSets& pts,
                               const Matching& matching) {
    if (!matching.matchedRefPointIds.size() ||
        !matching.matchedTransPointIds.size()) {
        return McMat4f::IDENTITY;
    }

    // Do a final affine transform for the best point matching.
    theMsg->printf("----------------");
    int numRows = 2 * matching.matchedRefPointIds.size();
    int numCols = 6;
    McDArray<float> a(numCols);
    McDArray<float> b(numRows);
    McDMatrix<float> M(numRows, numCols);
    M.fill(0.0f);
    for (int i = 0; i < (numRows / 2); ++i) {
        M[2 * i][0] = pts.trans.positions[matching.matchedTransPointIds[i]][0];
        M[2 * i][1] = pts.trans.positions[matching.matchedTransPointIds[i]][1];
        M[2 * i][2] = 1.0f;
        M[2 * i + 1][3] =
            pts.trans.positions[matching.matchedTransPointIds[i]][0];
        M[2 * i + 1][4] =
            pts.trans.positions[matching.matchedTransPointIds[i]][1];
        M[2 * i + 1][5] = 1.0f;
        b[2 * i] = pts.ref.positions[matching.matchedRefPointIds[i]][0];
        b[2 * i + 1] = pts.ref.positions[matching.matchedRefPointIds[i]][1];
    }
    McDMatrix<float> Q1, Q2;
    McDArray<float> S(numRows);
    M.SVD(Q1, &S[0], Q2);

    McDMatrix<float> Sinv(numCols, numCols);
    Sinv.makeIdentity();
    for (int i = 0; i < numCols; ++i) {
        Sinv[i][i] = 1.0f / S[i];
    }
    McDMatrix<float> Minv(numCols, numRows);
    Minv = Q2;
    Minv *= Sinv;
    Minv *= (Q1.transpose());

    Minv.multVec(&b[0], &a[0]);

    McMat4f m = McMat4f::IDENTITY;

    m[0][0] = a[0];
    m[1][0] = a[1];
    m[3][0] = a[2];
    m[0][1] = a[3];
    m[1][1] = a[4];
    m[3][1] = a[5];

    for (int i = 0; i < 4; ++i) {
        theMsg->printf("%f\t%f\t%f\t%f", m[i][0], m[i][1], m[i][2], m[i][3]);
    }
    return m;
}

McMat4f ma::fitTransform(const FacingPointSets& pts, const Matching& matching,
                         TransformType tfType) {
    McMat4f tf;
    switch (tfType) {

    case TF_NONE:
        tf = McMat4f::IDENTITY;
        break;

    case TF_RIGID:
    case TF_RIGID_ISO_SCALE:
        tf = fitTransformRigid(pts, matching, tfType);
        break;

    case TF_AFFINE:
        tf = fitTransformAffine(pts, matching);
        break;
    }

    tf[2][2] = 1.0f;  // no z-scaling
    tf[3][2] = 0.0f;  // no z-translation
    return tf;
}
