#pragma once

#include <hxspatialgraph/internal/SpatialGraphSelection.h>
#include <mclib/McHandable.h>

#include <hxspreadsheet/internal/HxSpreadSheet.h>

#include <hxalignmicrotubules/mtalign.h>

class SpreadSheetWrapper;
class MovingLeastSquares;
class HxSpatialGraph;

namespace mtalign {
struct CPDParams;
struct WarpResult;
struct MatchingParams;
struct Matching;
}

class HxNeuronEditorSubApp;

class MicrotubuleSpatialGraphAligner : public McHandable {
    struct AlignmentResults {
        float duration_s;
        float score;
        McDArray<int> matchedRefPointIds;
        McDArray<int> matchedTransPointIds;
        SpatialGraphSelection refVertexSelection;
        SpatialGraphSelection transVertexSelection;
    };

  public:
    MicrotubuleSpatialGraphAligner(HxSpatialGraph* graph);
    ~MicrotubuleSpatialGraphAligner();

    /// Return a list of attribute names which can be used for alignment.
    McDArray<McString> getPossibleAlignmentAttributeNames() const;

    /// Return a list of projection types which can be used for alignment.
    static McDArray<McString> getPossiblePointMatchingAlgorithms();

    /// Return a list of projection types which can be used for alignment.
    static McDArray<McString> getPossibleProjectionTypes();

    /// Return a list of transform types which can be used for alignment.
    static McDArray<McString> getPossibleTransformTypes();

    McString getNameOfIthSlice(const int sliceNum) const;

    /// Set/get whether to automatically estimate gap size.
    void setAutoEstimateGapSize(const int autoEstimate);
    int getAutoEstimateGapSize() const;

    void setEndPointRegion(const float slicePercentage);
    void setMaxDistanceDifference(const float maxDistDiff);
    void setMinCliqueSizeFraction(const float minCliqueSizeFraction);
    void setMaxNumCliques(const int maxCliques);
    void setAlpha(const float alpha);
    void setPointMatchingAlgorithm(const McString pmAlg);
    void setGapSize(const float gapSize);
    void setDMin(const float dMin);
    void setDMax(const float dMax);
    void setDeltaH(const float deltaH);
    void setProjectionType(const McString pType);
    void setTransformType(const McString tType);
    void setMinAngle(const float angle);
    void setCreateMatchingLabels(const bool enabled);
    void setScaleTestMin(const float scaleTestMin);
    void setScaleTestMax(const float scaleTestMax);
    void setScaleTestIncrement(const float scaleTestIncrement);

    void
    setMaxAngleDiffForInitMatching(const float maxAngleDiffForInitMatching);

    void setMinAngleForOptMatching(const float minAngleForOptMatching);
    void setMaxDistForAngle(const float maxDistForAngle);
    void setUseAbsouteValueForEndPointRegion(
        const bool useAbsouteValueForEndPointRegionLabel);

    void setNumMaxPointsForInitialTransform(
        const int numMaxPointsForInitialTransform) {
        mNumMaxPointsForInitialTransform = numMaxPointsForInitialTransform;
    }

    void setAlignType(const bool performOptAlign) {
        mAlignType = performOptAlign;
    }

    void setWeightConfig(const mtalign::PGMPairWeightsParams config) {
        mWeightConfig = config;
    }

    void setPairFactorParam(double pairFactorParam) {
        mPairFactorParam = pairFactorParam;
    }

    void setAngleToPlaneFilter(double angleToPlaneFilter) {
        mAngleToPlaneFilter = angleToPlaneFilter;
    }

    void alignAllOrPair(HxNeuronEditorSubApp* editor = 0,
                        const int refSliceNum = -1,
                        const int transSliceNum = -1);

    void getSliceSelection(const int sliceNum,
                           SpatialGraphSelection& slice) const;

    /// Return the current transform for slice \a sliceNum.
    McMat4f getTransform(const int sliceNum) const;

    /// Return the z-component of the current transformation of slice \a
    /// sliceNum
    float getSliceZPosition(const int sliceNum) const;

    /// Move all slices starting from \a sliceNum such that the z-component of
    /// \a sliceNum becomes \a newZ.
    void setSliceZPosition(const int sliceNum, const float newZ,
                           HxNeuronEditorSubApp* editor = 0);

    /// Apply transform \a transform to \a num slices starting with \a start
    /// (from start to last slice, when num<0);
    void applyTransform(const McMat4f& transform, const int start = 0,
                        const int num = -1, HxNeuronEditorSubApp* editor = 0);

    void removeIntermediateAndAddLabelForUnmatched(const char* matchingLabel);

    void joinMatchedEnds(const char* matchingLabel);

    void addLabelForUnmatched(const char* matchingLabel);

    /// `getNumSlices()` returns the number of sections detected in the spatial
    /// graph based on the attribute `TransformInfo`.
    int getNumSlices() const;

  private:
    const mtalign::SliceSelector mSelectionHelper;

  public:
    static const int NumProjectionTypes;
    static const McString ProjectionTypes[];
    static const int NumTransformTypes;
    static const McString TransformTypes[];

    void warpAll(HxNeuronEditorSubApp* editor,
                 const mtalign::CPDParams& params);

    static void applyMLSToSelection(MovingLeastSquares& mls, HxSpatialGraph* sg,
                                    const SpatialGraphSelection& sel);

  protected:
    HxSpatialGraph* mGraph;  // The spatialgraph to be aligned.

    // Parameters for initial matching, suffix
    float mMaxDistanceDifference;
    float mMinCliqueSizeFraction;
    int mMaxNumCliques;
    float mAlpha;
    float mMaxAngleDiffForInitMatching;
    int mNumMaxPointsForInitialTransform;

    // Global parameters, suffix.
    float mEndPointRegion;
    bool mUseAbsoluteValueForEndPointRegion;
    McString mProjectionType;
    McString mTransformType;
    bool mCreateMatchingLabels;
    bool mAlignType;  // initial or opt using all points.

    float mMaxDistForAngle;
    float mAngleToPlaneFilter;

    mtalign::PGMPairWeightsParams mWeightConfig;

    // Parameters for opt matching
    static const int NumPointMatchingAlgorithmTypes;
    static const McString PointMatchingAlgorithmTypes[];
    McString mPMAlgorithm;
    double mPairFactorParam;  // PGM only

    // Parameters for scaling and gap testing, suffix
    float mGapSize;
    float mDMin;
    float mDMax;
    float mDeltaH;
    float mScaleTestMin;
    float mScaleTestMax;
    float mScaleTestIncrement;

    // Automatically find gap size If 0, gap size is taken as fixed.
    //
    //  - 1: gap size is estimated by varying trans plane from min to max.
    //  - 2: Gapsize and mid plane is varied.
    //
    int mAutoEstimateGapSize;

    // Label for vertices not taken into account for matching.
    int mNotUsedLabelValue;

    // Label value for vertices taken into account, but not matched.
    int mNotMatchedLabelValue;

    void alignPairAndTestScaling(const int refSliceNum, const int transSliceNum,
                                 McMat4f& transMat,
                                 HxNeuronEditorSubApp* editor,
                                 McHandle<SpreadSheetWrapper> spreadSheet);

    /// Creates a new attribute defined on the graph vertices for holding the
    /// found point correspondences.
    void createPointMatchingAttribute(HxNeuronEditorSubApp* editor);

    /// Assign colors to all vertices used for matching. Corresponding points
    /// get the same color.
    void setPointMatchingAttribute(const SpatialGraphSelection& refVertexSel,
                                   const SpatialGraphSelection& transVertexSel,
                                   const McDArray<int>& matchedRefPoints,
                                   const McDArray<int>& matchedTransPoints,
                                   const int refSliceNum,
                                   const int transSliceNum,
                                   HxNeuronEditorSubApp* editor);

    /// Computes and prints the average point distances (to help estimating the
    /// threshold for the d-bounded matching).  Create extra labels for the
    /// pointmatching attribute.
    void rewriteLabels(HxNeuronEditorSubApp* editor);

    void rewriteBinaryLabel(const McDArray<int>& points, McString labelName);

    void setUsedForAlignmentLabel(const SpatialGraphSelection& refSelection,
                                  const SpatialGraphSelection& transSelection);
    void rewriteLabel(const SpatialGraphSelection& selection,
                      const McString& labelName);

    /// Returns the parameter bundle containing the transformation information
    /// for `slice`.
    HxParameter* getTransformParameter(const int slice) const;

    /// Returns whether the spatialgraph has a TransformInfo parameter bundle.
    bool hasTransformInfo() const;

    /// Create or reuse the evaluation spreadsheet.
    McHandle<SpreadSheetWrapper> getEvaluationSpreadsheet();

    /// Returns the attribute to sort the slices by from the TransformInfo
    /// parameter bundle.
    McString getTransformInfoAttribute() const;

    /// Sets all vertices in slice refSliceNum and transSliceNum to
    /// mNotUsedLabelValue if they are 50 percent from the slices that should
    /// be matched.
    void clearPointMatchingAttributeForSlices(const int refSliceNum,
                                              const int transSliceNum,
                                              HxNeuronEditorSubApp* editor);

    mtalign::EndPointParams
    makePointRepresentationParams(int refSliceNum, int transSliceNum,
                                  float projectionPlane);

    mtalign::MatchingParams makeMatchingParams() const;

    struct AlignmentInputParams {
        AlignmentInputParams(float projPlane, int refSlNum, int transSlNum) {
            projectionPlane = projPlane;
            refSliceNum = refSlNum;
            transSliceNum = transSlNum;
        };
        float projectionPlane;
        int refSliceNum;
        int transSliceNum;
    };

    void writeResultsToSpreadSheet(const AlignmentInputParams& inputParams,
                                   const AlignmentResults& alignmentResults,
                                   const float gap, const float scale,
                                   McHandle<SpreadSheetWrapper> spreadSheet);

    McMat4f alignNonPGM(const AlignmentInputParams& inputParams,
                        AlignmentResults& alignmentResults);

    McMat4f alignPGM(const AlignmentInputParams& inputParams,
                     AlignmentResults& alignmentResults);

    McMat4f align(const AlignmentInputParams& inputParams,
                  AlignmentResults& alignmentResults,
                  HxNeuronEditorSubApp* editor);

    void getPointMatchingString(McString& curName);

    void warpSlices(const int refSliceNum, const int transSliceNum,
                    const mtalign::CPDParams& params,
                    mtalign::WarpResult& deformation);

    void applyMLS(MovingLeastSquares& mls, const int start,
                  HxNeuronEditorSubApp* editor);

    void addWarpPointsToParams(const mtalign::MLSParams& mlsParams,
                               const int sliceNum);
};
