#pragma once

#include <mclib/McBox2f.h>
#include <mclib/McDArray.h>
#include <mclib/McVec2.h>
#include <mclib/McVec2i.h>

#include <hxalignmicrotubules/api.h>

/// Moving least squares warping with rigid transform as described in the paper
/// Schaefer et al. (2006) "Image deformation using moving least squares".
class HXALIGNMICROTUBULES_API MovingLeastSquares {

    private:

        class Grid2D {

            private:

                class Cell2D {

                    public:

                        Cell2D();

                    public:

                        McDArray< int > mElements;
                };

            public:

                Grid2D();

            public:

                const McVec2d& getCellSize() const;

                int getNeighbors(McDArray< int >& neighbors, const McVec2d& position, double radius, int minNeighbors) const;

                void initialize(const McDArray< McVec2d >& points);

            private:

                void clear();

                int coordToCell(const McVec2i& cellCoord) const;

                int coordToCell(const McVec2d& position) const;

                void initializeEmptyGrid();

            private:


                McDArray< Cell2D > mCells;
                McVec2d            mCellSize;
                McVec2i            mDimensions;
                McVec2d            mGridStart;
                McVec2d            mGridEnd;
                double             mMaxRadius;
    };


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
    /// Sets an additional box. All points are scaled
    /// and translated such that the box lies in the
    /// origin and the longest side is 1. This reduces
    /// numerical problems with large data sets.
    void setLandmarks(const McDArray<McVec2d>& src,
                      const McDArray<McVec2d>& dst,
                      const McBox2f& box);

    /// \pre `src.size() == dst.size()`.
    void setLandmarks(const McDArray<McVec3f>& src,
                      const McDArray<McVec3f>& dst);

    /// Swaps source and destination landmarks to get the
    /// inverse transformation.
    void swapLandmarks();

    /// Initializes a grid for interpolateFast.
    /// This call shall be done only once after setting
    /// landmarks or swapping landmarks.
    void initializeGrid();

    /// Interpolate destination position for `point`.
    McVec2d interpolate(const McVec2d& point) const;

    /// Interpolates destination position for `point`.
    /// A grid is used to select only Landmarks in the
    /// neighborhood.
    McVec2d interpolateFast(const McVec2d& point, int minPoints = 3) const;

  private:

    /// Computes a transformation such that all points lie in the
    /// unit square (increasing numerical stability)
    void updateTransformation();

  private:

    double mAlpha;
    McDArray<McVec2d> mPs;
    McDArray<McVec2d> mQs;

    McBox2f           mBox;
    bool              mBoxTransform;
    Grid2D            mGrid;
    double            mScale;
    McVec2d           mTranslate;
};
