#pragma once

#include <mclib/McDArray.h>
#include <mclib/McVec2d.h>

#include <hxalignmicrotubules/api.h>

/// Moving least squares warping with rigid transform as described in the paper
/// Schaefer et al. (2006) "Image deformation using moving least squares".
class HXALIGNMICROTUBULES_API MovingLeastSquares {
  public:
    MovingLeastSquares();

    /// Source landmarks.
    const McDArray<McVec2d>& ps() const { return mPs; }

    /// Destination landmarks.
    const McDArray<McVec2d>& qs() const { return mQs; }

    /// See paper.
    double alpha() const { return mAlpha; }

    ///
    void setAlpha(const double a);

    /// \pre `src.size() == dst.size()`.
    void setLandmarks(const McDArray<McVec2d>& src,
                      const McDArray<McVec2d>& dst);

    /// \pre `src.size() == dst.size()`.
    void setLandmarks(const McDArray<McVec3f>& src,
                      const McDArray<McVec3f>& dst);

    /// Interpolate destination position for `point`.
    McVec2d interpolate(const McVec2d& point) const;

  private:
    double mAlpha;
    McDArray<McVec2d> mPs;
    McDArray<McVec2d> mQs;
};
