#pragma once

#include <hxneuroneditor/internal/HxNeuronEditorSubApp.h>
#include <hxneuroneditor/internal/QxNeuronEditorToolBox.h>
#include <mclib/McBitfield.h>

#include <hxalignmicrotubules/mtalign/PGMPairWeights.h>

#include <hxalignmicrotubules/ui_QxMicrotubuleAlignSpatialGraphTool.h>

#include <hxalignmicrotubules/api.h>

class MicrotubuleSpatialGraphAligner;
class HxManualMTAlign;
class HxSpatialGraph;
class QWidget;
class QObject;
class QButtonGroup;
class QShortcut;

class HXALIGNMICROTUBULES_API QxMicrotubuleAlignSpatialGraphTool
    : public QObject,
      public QxNeuronEditorToolBox {

    Q_OBJECT

  public:
    QxMicrotubuleAlignSpatialGraphTool(HxNeuronEditorSubApp* editor);
    ~QxMicrotubuleAlignSpatialGraphTool();

    virtual QWidget* toolcard();

    /// Set the row'th slice to be the one to be transformed and the (row-1)'th
    /// slice as the reference slice.
    void selectSlice(const int row);

    /// Show/hide row-th slice.
    void showSlice(const int row, const bool show);

    /// Show/hide all slices.
    void showAllSlices(const bool show);

    int getRefSliceNum() const;
    int getTransSliceNum() const;

  public slots:
    /// Handle click on one of the cells in the table.
    void cellClicked(const int row, const int col);

    /// Show all slices in 3D viewer.
    void showAllClicked();

    /// Hide all slices in 3D viewer.
    void hideAllClicked();

    /// Toggle manual transform.
    void transformClicked();

    /// Notification if item in table changed.
    void tableItemChanged(QTableWidgetItem* item);

    // This is called whenever the alignment toolbox icon is clicked.
    virtual void updateToolcard();

  protected slots:
    /// Do the alignment.
    void align();

    /// Update toolbox after graph has changed outside toolbox
    void updateAfterSpatialGraphChange(
        HxNeuronEditorSubApp::SpatialGraphChanges changes);

    /// Set text of end point region label when state changed.
    void endPointRegionOptionChanged();

    void addEvidenceLabel();

    void addLandmark();

    void joinEnds();

    void removeAllIntermediate();

    void applyCPD();

  protected:
    // Toolbox GUI.
    QWidget* uiParent;
    Ui::QxMicrotubuleAlignSpatialGraphTool ui;
    McHandle<MicrotubuleSpatialGraphAligner> mAligner;
    HxSpatialGraph* mGraph;

    McBitfield mVisibleSlices;

    // Manual alignment.
    bool mManualTransforming;
    McHandle<HxSpatialGraph> mManualSlice;
    McHandle<HxManualMTAlign> mManualAligner;

    SpatialGraphSelection mOriginalVisibleSelection;
    SpatialGraphSelection mManualSliceSelection;

    void startManualTransform();
    void stopManualTransform();

    void updateSliceTableZCoords();

    int mRefSlice;
    int mTransSlice;

    /// Button group for the pair-wise or entire stack alignment radio buttons.
    QButtonGroup* mButtonGroup;

    void updateGraph();
    void initAligner();

    /// Fills the sorting attribute combobox.
    void fillComboBoxes();

    /// Fills the slice table of the GUI.
    void fillTable();

    QShortcut* mAddEvidenceForPGMMatching;
    QShortcut* mAddLandmarksForWarping;

    void createParamsForPGMMatching(mtalign::PGMPairWeightsParams& config);

    bool shouldAddVertexAttributeToVertexColoringBox(const char* item);

    void copyViewerSettingsForManualTransform();
};

#ifdef __cplusplus
extern "C" {
#endif

/// Hook to register the `AlignSpatialGraphTool` at the `HxNeuronEditorSubApp`.
/// `extern "C"` to avoid name mangling, so that `HxDSO` can find it.
void HXALIGNMICROTUBULES_API
QxMicrotubuleAlignSpatialGraphTool_init_plugin(HxNeuronEditorSubApp* editor);

#ifdef __cplusplus
}
#endif
