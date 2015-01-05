#include <hxalignmicrotubules/QxMicrotubuleAlignSpatialGraphTool.h>

#include <QShortcut>

#include <hxcore/HxResource.h>
#include <hxneuroneditor/HxNeuronEditorSubApp.h>
#include <hxspatialgraph/HierarchicalLabels.h>
#include <hxspatialgraph/HxSpatialGraph.h>
#include <hxspatialgraph/HxSpatialGraphView.h>
#include <mclib/McException.h>

#include <hxalignmicrotubules/MicrotubuleSpatialGraphAligner.h>
#include <hxalignmicrotubules/HxManualMTAlign.h>
#include <hxalignmicrotubules/mtalign.h>

namespace ma = mtalign;

QxMicrotubuleAlignSpatialGraphTool::QxMicrotubuleAlignSpatialGraphTool(
    HxNeuronEditorSubApp* editor)
    : QxNeuronEditorToolBox(editor),
      uiParent(0),
      mAligner(0),
      mGraph(0),
      mManualTransforming(false),
      mButtonGroup(0) {
    HxViewer* viewer = editor->get3DViewer();
    QxViewerPanel* panel = viewer ? viewer->getPanel() : NULL;
    if (panel) {
        mAddEvidenceForPGMMatching =
            new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_E), (QWidget*)panel);
        mAddLandmarksForWarping =
            new QShortcut(QKeySequence(Qt::CTRL + Qt::Key_W), (QWidget*)panel);
    }
}

QxMicrotubuleAlignSpatialGraphTool::~QxMicrotubuleAlignSpatialGraphTool() {
    delete mButtonGroup;
    delete uiParent;
}

QWidget* QxMicrotubuleAlignSpatialGraphTool::toolcard() {
    if (!uiParent) {
        uiParent = new QWidget;
        uiParent->setObjectName("MicrotubuleAlignSpatialGraphTool");
        ui.setupUi(uiParent);

        // Automatic alignment: pair/all.
        mButtonGroup = new QButtonGroup(uiParent);
        mButtonGroup->addButton(ui.pairRadioButton, 0);
        mButtonGroup->addButton(ui.allRadioButton, 1);
        ui.pairRadioButton->setChecked(true);

        // Automatic alignment: general options.
        mRefSlice = 0;
        mTransSlice = 1;
        ui.includeLineEdit->setText(QString("40"));
        ui.maxCliquesLineEdit->setText("1000");

        // Automatic alignment: starting transformation options.
        ui.maxPointDistLineEdit->setText("600");
        ui.cliqueSizeLineEdit->setText("0.3");
        ui.numPointsForStartTransform->setText("50");

        ui.dMinLineEdit->setText(QString("0"));
        ui.dMaxLineEdit->setText(QString("0"));
        ui.deltaHLineEdit->setText(QString("100"));

        connect(ui.alignButton, SIGNAL(clicked()), this, SLOT(align()));

        // Manual transform GUI.
        ui.sliceTable->setColumnCount(3);
        QStringList labels;
        labels.append(QString("Slice value"));
        labels.append(QString("Show"));
        labels.append(QString("Z-Coord"));
        ui.sliceTable->setHorizontalHeaderLabels(labels);
        ui.sliceTable->verticalHeader()->hide();
        ui.sliceTable->setSelectionMode(QAbstractItemView::SingleSelection);
        ui.sliceTable->setSelectionBehavior(QAbstractItemView::SelectRows);
        ui.scaleTestMin->setText(QString("1.0"));
        ui.scaleTestMax->setText(QString("1.0"));
        ui.scaleTestIncrement->setText(QString("0.01"));
        ui.dMaxAngleDiffForInitialMatching->setText(QString("20.0"));
        ui.angularMatchingSmoothing->setText(QString("2000"));
        ui.useAbsouteValueForEndPointRegionCheckBox->setChecked(false);
        ui.radioButtonInitialAlign->setChecked(true);
        connect(ui.sliceTable, SIGNAL(cellClicked(int, int)), this,
                SLOT(cellClicked(int, int)));
        connect(ui.sliceTable, SIGNAL(itemChanged(QTableWidgetItem*)), this,
                SLOT(tableItemChanged(QTableWidgetItem*)));

        connect(ui.showAllButton, SIGNAL(clicked()), this,
                SLOT(showAllClicked()));
        connect(ui.hideAllButton, SIGNAL(clicked()), this,
                SLOT(hideAllClicked()));
        connect(ui.transformButton, SIGNAL(clicked()), this,
                SLOT(transformClicked()));

        connect(mEditor, SIGNAL(spatialGraphChanged(
                             HxNeuronEditorSubApp::SpatialGraphChanges)),
                this, SLOT(updateAfterSpatialGraphChange(
                          HxNeuronEditorSubApp::SpatialGraphChanges)));

        connect(ui.useAbsouteValueForEndPointRegionCheckBox,
                SIGNAL(toggled(bool)), this,
                SLOT(endPointRegionOptionChanged()));

        connect(mAddEvidenceForPGMMatching, SIGNAL(activated()), this,
                SLOT(addEvidenceLabel()));
        connect(mAddLandmarksForWarping, SIGNAL(activated()), this,
                SLOT(addLandmark()));

        // Init opt parameters.
        ui.use3dDistanceThres->setChecked(true);
        ui.threeDDistanceThres->setText("1117.6");
        ui.useAngleThres->setChecked(true);
        ui.angleThres->setText("26.8");
        ui.useProjectedDistanceThres->setChecked(true);
        ui.projectedDistanceThres->setText("784");
        ui.use3dDistanceWeight->setChecked(true);
        ui.threeDDistanceWeight->setText("242.7");
        ui.useAngleWeight->setChecked(true);
        ui.angleWeight->setText("5.8");
        ui.useProjectedDistanceWeight->setChecked(true);
        ui.projectedDistanceWeight->setText("170.2");
        ui.exponentialWeightButton->setChecked(true);
        ui.pairFactorParam->setText("150");
        ui.createMatchingLabelsCheckBox->setChecked(false);

        connect(ui.joinButton, SIGNAL(clicked()), this, SLOT(joinEnds()));
        connect(ui.finalButton, SIGNAL(clicked()), this,
                SLOT(removeAllIntermediate()));

        connect(mEditor, SIGNAL(spatialGraphChanged(
                             HxNeuronEditorSubApp::SpatialGraphChanges)),
                this, SLOT(updateToolcard()));

        ui.dummySignificance->setText("0.01");

        ui.projectionComboBox->clear();
        ui.transformComboBox->clear();
        ui.pmAlgorithmComboBox->clear();
        McDArray<McString> projTypes =
            MicrotubuleSpatialGraphAligner::getPossibleProjectionTypes();
        for (int i = 0; i < projTypes.size(); ++i) {
            ui.projectionComboBox->addItem(QString(projTypes[i].getString()));
        }
        McDArray<McString> transTypes =
            MicrotubuleSpatialGraphAligner::getPossibleTransformTypes();
        for (int i = 0; i < transTypes.size(); ++i) {
            ui.transformComboBox->addItem(QString(transTypes[i].getString()));
        }
        McDArray<McString> algoTypes = MicrotubuleSpatialGraphAligner::
            getPossiblePointMatchingAlgorithms();
        for (int i = 0; i < algoTypes.size(); ++i) {
            ui.pmAlgorithmComboBox->addItem(QString(algoTypes[i].getString()));
        }

        ui.cpdBetaLineEdit->setText("10");
        ui.cpdLambdaLineEdit->setText("1");
        ui.cpdWLineEdit->setText("0.1");
        connect(ui.applyCPD, SIGNAL(clicked()), this, SLOT(applyCPD()));
        ui.CPDTypeComboBox->addItem("NA");
        ui.CPDTypeComboBox->addItem("RigidFisherMises");
        ui.CPDTypeComboBox->addItem("non-linear");
        ui.CPDTypeComboBox->setCurrentIndex(1);
        ui.cpdWithCoordButton->setChecked(true);
        ui.cpdUseScalingButton->setChecked(true);
        ui.cpdUseDirectionButton->setChecked(true);
        ui.angleToPlaneFilter->setText("0.01");
        ui.cpdLandmarkSampleDist->setText("0");
        ui.mlsAlpha->setText("2");
        ui.cpdMaxSigmaSquare->setText("0.001");
    }

    return uiParent;
}

void QxMicrotubuleAlignSpatialGraphTool::fillComboBoxes() {
    ui.sortComboBox->clear();
    if (mAligner) {
        McDArray<McString> attNames =
            mAligner->getPossibleAlignmentAttributeNames();
        for (int i = 0; i < attNames.size(); ++i) {
            ui.sortComboBox->addItem(QString(attNames[i].getString()));
        }
    }
}

void QxMicrotubuleAlignSpatialGraphTool::initAligner() {
    mcassert(mAligner);
    QString attName = ui.sortComboBox->currentText();
    McString attName2 = McString(attName.toAscii().constData());
    if (attName2 == "") {
        theMsg->printf("Error: no valid attribute selected");
        return;
    }
    Qt::CheckState gapTestState = ui.GapSizeAutoCheckBox->checkState();
    int gapTest = (gapTestState == Qt::Unchecked) ? 0 : 1;
    mAligner->setAutoEstimateGapSize(gapTest);

    mAligner->setDMin(ui.dMinLineEdit->text().toFloat());
    mAligner->setDMax(ui.dMaxLineEdit->text().toFloat());
    mAligner->setDeltaH(ui.deltaHLineEdit->text().toFloat());

    mAligner->setEndPointRegion(ui.includeLineEdit->text().toFloat());
    mAligner->setMaxDistanceDifference(
        ui.maxPointDistLineEdit->text().toFloat());
    mAligner->setMaxNumCliques(ui.maxCliquesLineEdit->text().toInt());
    mAligner->setMinCliqueSizeFraction(ui.cliqueSizeLineEdit->text().toFloat());
    mAligner->setCreateMatchingLabels(
        ui.createMatchingLabelsCheckBox->isChecked());

    McString projType =
        McString(ui.projectionComboBox->currentText().toAscii().constData());
    mAligner->setProjectionType(projType);
    McString transType =
        McString(ui.transformComboBox->currentText().toAscii().constData());
    mAligner->setTransformType(transType);
    McString pmAlg =
        McString(ui.pmAlgorithmComboBox->currentText().toAscii().constData());
    mAligner->setPointMatchingAlgorithm(pmAlg);
    mAligner->setScaleTestMin(ui.scaleTestMin->text().toFloat());
    mAligner->setScaleTestMax(ui.scaleTestMax->text().toFloat());
    mAligner->setScaleTestIncrement(ui.scaleTestIncrement->text().toFloat());
    mAligner->setMaxAngleDiffForInitMatching(
        ui.dMaxAngleDiffForInitialMatching->text().toFloat());
    mAligner->setUseAbsouteValueForEndPointRegion(
        ui.useAbsouteValueForEndPointRegionCheckBox->isChecked());
    mAligner->setMaxDistForAngle(ui.angularMatchingSmoothing->text().toFloat());
    mAligner->setNumMaxPointsForInitialTransform(
        ui.numPointsForStartTransform->text().toInt());
    mAligner->setAlignType(ui.radioButtonOptAlign->isChecked());
    mtalign::PGMPairWeightsParams config;
    createParamsForPGMMatching(config);
    mAligner->setWeightConfig(config);
    mAligner->setPairFactorParam(ui.pairFactorParam->text().toDouble());
    mAligner->setAngleToPlaneFilter(ui.angleToPlaneFilter->text().toDouble());
}

void QxMicrotubuleAlignSpatialGraphTool::createParamsForPGMMatching(
    mtalign::PGMPairWeightsParams& config) {
    config.distanceThreshold3d = ui.threeDDistanceThres->text().toDouble();
    config.useDistanceThreshold3d = ui.use3dDistanceThres->isChecked();
    config.distanceThresholdProjected =
        ui.projectedDistanceThres->text().toDouble();
    config.useDistanceThresholdProjected =
        ui.useProjectedDistanceThres->isChecked();
    config.angleThreshold = ui.angleThres->text().toDouble();
    config.useAngleThreshold = ui.useAngleThres->isChecked();
    // Weight parameters.
    config.angleWeightParam = ui.angleWeight->text().toDouble();
    config.useAngleWeight = ui.useAngleWeight->isChecked();
    config.dist3dParam = ui.threeDDistanceWeight->text().toDouble();
    config.useDist3dWeight = ui.use3dDistanceWeight->isChecked();
    config.distProjectedParam = ui.projectedDistanceWeight->text().toDouble();
    config.useProjectedDistWeight = ui.useProjectedDistanceWeight->isChecked();
    config.weightType = ui.linareWeightButton->isChecked()
                            ? mtalign::PGMPairWeightsParams::LINEAR
                            : mtalign::PGMPairWeightsParams::EXPONENTIAL;
    config.dummySignificance = ui.dummySignificance->text().toDouble();
}

void QxMicrotubuleAlignSpatialGraphTool::updateGraph() {
    // Check whether spatialgraph changed.
    if (mEditor->getSpatialGraph() != mGraph) {
        mGraph = mEditor->getSpatialGraph();
        if (mGraph) {
            mAligner = new MicrotubuleSpatialGraphAligner(mGraph);
            fillComboBoxes();
            fillTable();
        } else {
            mAligner = 0;
            ui.sliceTable->clear();
        }
    }
    if (mGraph) {
        updateSliceTableZCoords();
    }
}

void QxMicrotubuleAlignSpatialGraphTool::updateToolcard() {
    updateGraph();
    // Node coloring (currently only labels are supported).
    ui.wNodeColoringComboBox->blockSignals(true);
    ui.wNodeColoringComboBox->clear();
    HxSpatialGraphView* view = mEditor->getView();
    int numNodeColorItems = view->portVertexColoring.getNum();
    for (int i = 0; i < numNodeColorItems; ++i) {
        QString item = view->portVertexColoring.getLabel(i);
        if (shouldAddVertexAttributeToVertexColoringBox(qPrintable(item)) ||
            (i == 0))  // constant color
        {
            ui.wNodeColoringComboBox->addItem(item);
        }
    }

    const int currentIdx = view->portVertexColoring.getIndex(0);
    const QString currentItem = view->portVertexColoring.getLabel(currentIdx);

    int comboBoxIndex = ui.wNodeColoringComboBox->findText(currentItem);
    if (comboBoxIndex >= 0) {
        ui.wNodeColoringComboBox->setCurrentIndex(comboBoxIndex);
    } else {
        ui.wNodeColoringComboBox->setCurrentIndex(0);
    }

    ui.wNodeColoringComboBox->blockSignals(false);

    copyViewerSettingsForManualTransform();
}

void QxMicrotubuleAlignSpatialGraphTool::align() {
    initAligner();
    if (mAligner) {
        const bool manualTrans = mManualTransforming;
        if (manualTrans) {
            stopManualTransform();
        }
        if (ui.pairRadioButton->isChecked()) {
            // Align pair.
            int refSliceNum = getRefSliceNum();
            int transSliceNum = getTransSliceNum();
            mAligner->alignAllOrPair(mEditor, refSliceNum, transSliceNum);
        } else if (ui.allRadioButton->isChecked()) {
            // Align all.
            mAligner->alignAllOrPair(mEditor);
        }
        if (manualTrans) {
            startManualTransform();
        }
    }
}
void QxMicrotubuleAlignSpatialGraphTool::applyCPD() {
    initAligner();
    ma::CPDParams cpdParams;
    QString cpdType = ui.CPDTypeComboBox->currentText();

    if (QString::compare(cpdType, QString("Rigid")) == 0) {
        mcthrow("'Rigid' is not implemented!");
    } else if (QString::compare(cpdType, QString("RigidFisherMises")) == 0) {
        cpdParams.type = ma::CPD_RIGID_FISHER_MISES;
        cpdParams.rigid.withScaling = ui.cpdUseScalingButton->isChecked();
        cpdParams.rigid.usePositions = ui.cpdWithCoordButton->isChecked();
        cpdParams.rigid.useDirections = ui.cpdUseDirectionButton->isChecked();
        cpdParams.rigid.w = ui.cpdWLineEdit->text().toDouble();
        cpdParams.rigid.maxIterations = 200;
        cpdParams.rigid.eDiffRelStop = 1.e-5;
        cpdParams.rigid.sigmaSquareStop = 1.e-7;
    } else if (QString::compare(cpdType, QString("non-linear")) == 0) {
        cpdParams.type = ma::CPD_ELASTIC;
        cpdParams.elastic.beta = ui.cpdBetaLineEdit->text().toDouble();
        cpdParams.elastic.lambda = ui.cpdLambdaLineEdit->text().toDouble();
        cpdParams.elastic.useDirections = ui.cpdUseDirectionButton->isChecked();
        cpdParams.elastic.w = ui.cpdWLineEdit->text().toDouble();
        cpdParams.elastic.maxIterations = 200;
        cpdParams.elastic.eDiffRelStop = 1.e-5;
        cpdParams.elastic.sigmaSquareStop = 1.e-7;
        cpdParams.elastic.sampleDistForWarpingLandmarks =
            ui.cpdLandmarkSampleDist->text().toFloat();
    } else {
        mcerror("Must not reach this point.");
        mcthrow("Internal error: received invalid cpdType from GUI.");
    }
    cpdParams.alphaForMLS = ui.mlsAlpha->text().toDouble();
    cpdParams.maxAcceptableSigmaSquare =
        ui.cpdMaxSigmaSquare->text().toDouble();
    mAligner->warpAll(mEditor, cpdParams);
}

void QxMicrotubuleAlignSpatialGraphTool::fillTable() {
    ui.sliceTable->blockSignals(true);
    ui.sliceTable->clear();

    if (mAligner) {
        const int numValues = mAligner->getNumSlices();
        mVisibleSlices.resize(numValues, McBitfield::INIT_SET);
        ui.sliceTable->setRowCount(numValues);
        QStringList labels;
        labels.append(QString("Slice number (value)"));
        labels.append(QString("Show"));
        labels.append(QString("Z-Coord"));

        for (int i = 0; i < numValues; ++i) {
            // Display value in first column.
            McString filename = mAligner->getNameOfIthSlice(i);
            QString str(filename.dataPtr());

            QTableWidgetItem* item = new QTableWidgetItem(str);
            QSize bestSize(200, 20);
            item->setSizeHint(bestSize);
            ui.sliceTable->setItem(i, 0, item);
            item->setTextAlignment(Qt::AlignLeft);
            item->setFlags(Qt::ItemIsEnabled);

            // Display visibility icon in second column.
            QTableWidgetItem* visItem = new QTableWidgetItem();
            visItem->setIcon(
                QIcon(HxResource::readIcon("alignSpatialGraphEyeIcon.png")));
            ui.sliceTable->setItem(i, 1, visItem);

            // Display z-coordinate of slice in third column.
            QTableWidgetItem* zCoordItem = new QTableWidgetItem();
            ui.sliceTable->setItem(i, 2, zCoordItem);
            zCoordItem->setTextAlignment(Qt::AlignHCenter);
        }
        updateSliceTableZCoords();
        ui.sliceTable->resizeRowsToContents();
        ui.sliceTable->resizeColumnsToContents();
        ui.sliceTable->setMinimumHeight((40) * mAligner->getNumSlices());
    }

    ui.sliceTable->blockSignals(false);
}

void QxMicrotubuleAlignSpatialGraphTool::updateSliceTableZCoords() {
    ui.sliceTable->blockSignals(true);
    const int numValues = mAligner->getNumSlices();
    for (int i = 0; i < numValues; ++i) {
        QString zStr;
        zStr.sprintf("%f", mAligner->getSliceZPosition(i));
        QTableWidgetItem* zCoordItem = ui.sliceTable->item(i, 2);
        zCoordItem->setText(zStr);
    }
    ui.sliceTable->blockSignals(false);
}

void
QxMicrotubuleAlignSpatialGraphTool::tableItemChanged(QTableWidgetItem* item) {
    if (item->column() == 2) {
        int sliceNum = item->row();
        float newZ = item->text().toFloat();
        mAligner->setSliceZPosition(sliceNum, newZ, mEditor);
        updateSliceTableZCoords();
    }
}

void QxMicrotubuleAlignSpatialGraphTool::cellClicked(const int row,
                                                     const int col) {
    switch (col) {
    case 0:
        selectSlice(row);
        break;
    case 1:
        showSlice(row, !(mVisibleSlices[row]));
        break;
    case 2:
        selectSlice(row);
        break;
    }
}

void QxMicrotubuleAlignSpatialGraphTool::selectSlice(const int row) {
    if (row == 0)
        return;

    const bool manualTransforming = mManualTransforming;
    if (manualTransforming) {
        stopManualTransform();
    }

    mRefSlice = row - 1;
    mTransSlice = row;

    // Show only slices to be aligned.
    showAllSlices(false);
    showSlice(mRefSlice, true);
    showSlice(mTransSlice, true);

    if (manualTransforming) {
        startManualTransform();
    }
}

void QxMicrotubuleAlignSpatialGraphTool::showSlice(const int row,
                                                   const bool show) {
    if (show && !mVisibleSlices[row]) {
        SpatialGraphSelection sel;
        mAligner->getSliceSelection(row, sel);
        SpatialGraphSelection visible = mEditor->getVisibleElements();
        visible.addSelection(sel);
        mEditor->setVisibleElements(visible);
        mVisibleSlices.set(row, true);
        QTableWidgetItem* visItem = ui.sliceTable->item(row, 1);
        visItem->setIcon(
            QIcon(HxResource::readIcon("alignSpatialGraphEyeIcon.png")));
    } else if (!show && mVisibleSlices[row]) {
        SpatialGraphSelection sel;
        mAligner->getSliceSelection(row, sel);
        SpatialGraphSelection visible = mEditor->getVisibleElements();
        visible.subtractSelection(sel);
        mEditor->setVisibleElements(visible);
        mVisibleSlices.set(row, false);
        QTableWidgetItem* visItem = ui.sliceTable->item(row, 1);
        visItem->setIcon(
            QIcon(HxResource::readIcon("alignSpatialGraphClosedEyeIcon.png")));
    }
}

void QxMicrotubuleAlignSpatialGraphTool::showAllSlices(const bool show) {
    for (int i = 0; i < mAligner->getNumSlices(); ++i) {
        showSlice(i, show);
    }
}

void QxMicrotubuleAlignSpatialGraphTool::showAllClicked() {
    showAllSlices(true);
}

void QxMicrotubuleAlignSpatialGraphTool::hideAllClicked() {
    showAllSlices(false);
}

void QxMicrotubuleAlignSpatialGraphTool::transformClicked() {
    mManualTransforming = !mManualTransforming;
    if (mManualTransforming) {
        startManualTransform();
    } else {
        stopManualTransform();
    }
}

void QxMicrotubuleAlignSpatialGraphTool::startManualTransform() {
    mcassert(mGraph);

    ui.transformButton->setText("Stop trans");

    mOriginalVisibleSelection = mEditor->getVisibleElements();
    mAligner->getSliceSelection(mTransSlice, mManualSliceSelection);
    mManualSlice = mGraph->getSubgraph(mManualSliceSelection);

    mManualAligner = new HxManualMTAlign();
    mManualAligner->setViewerMask(1 << 14);
    mManualAligner->portData.connect(mManualSlice);
    mManualAligner->fire();

    copyViewerSettingsForManualTransform();

    // Hide slice being transformed in 3D viewer of editor.
    SpatialGraphSelection sel = mEditor->getVisibleElements();
    sel.subtractSelection(mManualSliceSelection);
    mEditor->setVisibleElements(sel);

    mManualTransforming = true;
    mManualAligner->startTransform();
}

void QxMicrotubuleAlignSpatialGraphTool::stopManualTransform() {
    mcassert(mGraph);

    ui.transformButton->setText("Transform");

    SbMatrix trans = mManualAligner->stopTransform();
    McMat4f mat = McMat4f(*reinterpret_cast<McMat4f*>(&trans));

    mAligner->applyTransform(mat, mTransSlice, -1, mEditor);

    mEditor->setVisibleElements(mOriginalVisibleSelection);
    SpatialGraphSelection highSel = mEditor->getHighlightedElements();
    highSel.clear();
    mEditor->setHighlightedElements(highSel);

    mManualAligner = 0;
    mManualTransforming = false;
}

int QxMicrotubuleAlignSpatialGraphTool::getRefSliceNum() const {
    return mRefSlice;
}

int QxMicrotubuleAlignSpatialGraphTool::getTransSliceNum() const {
    return mTransSlice;
}

void QxMicrotubuleAlignSpatialGraphTool::updateAfterSpatialGraphChange(
    HxNeuronEditorSubApp::SpatialGraphChanges changes) {
    if (HxNeuronEditorSubApp::SpatialGraphDataSetChange & changes) {
        updateGraph();
    }
    return;
}

void QxMicrotubuleAlignSpatialGraphTool::endPointRegionOptionChanged() {
    const bool state = ui.useAbsouteValueForEndPointRegionCheckBox->isChecked();
    if (state) {
        ui.includeLabel->setText(QString("Boundary    :"));
    } else {
        ui.includeLabel->setText(QString("Boundary (%):"));
    }
}

void
QxMicrotubuleAlignSpatialGraphTool_init_plugin(HxNeuronEditorSubApp* editor) {
    QxMicrotubuleAlignSpatialGraphTool* alignTool =
        new QxMicrotubuleAlignSpatialGraphTool(editor);
    alignTool->setText("MicrotubuleAlignTool");
    editor->registerToolBox(alignTool);
}

void QxMicrotubuleAlignSpatialGraphTool::addLandmark() {
    HxSpatialGraph* graph = mEditor->getSpatialGraph();
    SpatialGraphSelection highSel = mEditor->getHighlightedElements();
    const int numSelected = highSel.getNumSelectedVertices();
    if (numSelected != 2) {
        theMsg->printf("Can only add landmark if two nodes are selected.");
        return;
    }
    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(graph->addAttribute(
            "WarpPairs", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    EdgeVertexAttribute* trafoAtt =
        dynamic_cast<EdgeVertexAttribute*>(graph->addAttribute(
            "TransformInfo", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    if (!trafoAtt) {
        theMsg->printf("Landmarks are only useful if you have a TransformInfo "
                       "attribute. I'm not doing anything.");
        return;
    }

    // Make sure both nodes are in different and adjacent slices.
    McDArray<int> trafoVals;
    for (int i = 0; i < numSelected; ++i) {
        trafoVals.append(
            trafoAtt->getIntDataAtIdx(highSel.getSelectedVertex(i)));
    }
    if (trafoVals.size() > 1) {
        if (abs(trafoVals[0] - trafoVals[1]) != 1) {
            theMsg->printf("Landmarks assignment must be of nodes that are in "
                           "two different slices. I'm not doing anything.");
            return;
        }
    }

    // Check if the nodes are currently assigned to anything.
    McDArray<int> curAttVals;
    for (int i = 0; i < numSelected; ++i) {
        curAttVals.append(att->getIntDataAtIdx(highSel.getSelectedVertex(i)));
    }

    // Create label group or get it if it exists already.
    mGraph->addNewLabelGroup("WarpPairs", false, true);
    HierarchicalLabels* labelGroup = mGraph->getLabelGroup("WarpPairs");
    // Set everything that has these values to 0.
    for (int i = 0; i < numSelected; ++i) {
        if (curAttVals[i] > 0) {
            for (int j = 0; j < graph->getNumVertices(); j++) {
                int curVal = att->getIntDataAtIdx(j);
                if (curVal == curAttVals[i])
                    att->setIntDataAtIdx(j, 0);
            }
            // Remove the label, too.
            McDArray<int> dummy;
            labelGroup->removeLabel(curAttVals[i], dummy);
        }
    }

    const int labelId = labelGroup->getMaxLabelId();
    McString s;
    s.printf("Assignment%d", labelId + 1);
    SbColor color;
    color[0] = float(rand()) / float(RAND_MAX);
    color[1] = float(rand()) / float(RAND_MAX);
    color[2] = float(rand()) / float(RAND_MAX);
    int val = mGraph->addLabel("WarpPairs", 0, s.getString(), color);

    for (int i = 0; i < numSelected; ++i) {
        att->setIntDataAtIdx(highSel.getSelectedVertex(i), val);
    }
    if (mEditor) {
        mGraph->touch(HxData::NEW_PARAMETERS);
        const HxNeuronEditorSubApp::SpatialGraphChanges changes =
            HxNeuronEditorSubApp::SpatialGraphLabelTreeChange |
            HxNeuronEditorSubApp::SpatialGraphAttributeValueChange;
        mEditor->updateAfterSpatialGraphChange(changes);
    }
}

void QxMicrotubuleAlignSpatialGraphTool::addEvidenceLabel() {
    HxSpatialGraph* graph = mEditor->getSpatialGraph();
    SpatialGraphSelection highSel = mEditor->getHighlightedElements();
    const int numSelected = highSel.getNumSelectedVertices();
    if (numSelected != 1 && numSelected != 2) {
        theMsg->printf("Can only add evidence if at least one node and at most "
                       "two are selected. Select either one or two nodes.");
        return;
    }
    EdgeVertexAttribute* att =
        dynamic_cast<EdgeVertexAttribute*>(graph->addAttribute(
            "Evidence", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    EdgeVertexAttribute* trafoAtt =
        dynamic_cast<EdgeVertexAttribute*>(graph->addAttribute(
            "TransformInfo", HxSpatialGraph::VERTEX, McPrimType::mc_int32, 1));
    if (!trafoAtt) {
        theMsg->printf("Evidence is only useful if you have a TransformInfo "
                       "attribute. I'm not doing anything.");
        return;
    }

    // Make sure both nodes are in different and adjacent slices.
    McDArray<int> trafoVals;
    for (int i = 0; i < numSelected; ++i) {
        trafoVals.append(
            trafoAtt->getIntDataAtIdx(highSel.getSelectedVertex(i)));
    }
    if (trafoVals.size() > 1) {
        if (abs(trafoVals[0] - trafoVals[1]) != 1) {
            theMsg->printf("Evidence assignment must be of nodes that are in "
                           "two different slices. I'm not doing anything.");
            return;
        }
    }

    // Check if the nodes are currently assigned to anything.
    McDArray<int> curAttVals;
    for (int i = 0; i < numSelected; ++i) {
        curAttVals.append(att->getIntDataAtIdx(highSel.getSelectedVertex(i)));
    }

    // Create label group or get it if it exists already.
    mGraph->addNewLabelGroup("Evidence", false, true);
    HierarchicalLabels* labelGroup = mGraph->getLabelGroup("Evidence");
    // Set everything that has these values to 0.
    for (int i = 0; i < numSelected; ++i) {
        if (curAttVals[i] > 0) {
            for (int j = 0; j < graph->getNumVertices(); j++) {
                int curVal = att->getIntDataAtIdx(j);
                if (curVal == curAttVals[i])
                    att->setIntDataAtIdx(j, 0);
            }
            // Remove the label, too.
            McDArray<int> dummy;
            labelGroup->removeLabel(curAttVals[i], dummy);
        }
    }

    const int labelId = labelGroup->getMaxLabelId();
    McString s;
    s.printf("Assignment%d", labelId + 1);
    SbColor color;
    color[0] = float(rand()) / float(RAND_MAX);
    color[1] = float(rand()) / float(RAND_MAX);
    color[2] = float(rand()) / float(RAND_MAX);
    const int val = mGraph->addLabel("Evidence", 0, s.getString(), color);

    for (int i = 0; i < numSelected; ++i) {
        att->setIntDataAtIdx(highSel.getSelectedVertex(i), val);
    }
    if (mEditor) {
        mGraph->touch(HxData::NEW_PARAMETERS);
        HxNeuronEditorSubApp::SpatialGraphChanges changes =
            HxNeuronEditorSubApp::SpatialGraphLabelTreeChange |
            HxNeuronEditorSubApp::SpatialGraphAttributeValueChange;
        mEditor->updateAfterSpatialGraphChange(changes);
    }
}

bool
QxMicrotubuleAlignSpatialGraphTool::shouldAddVertexAttributeToVertexColoringBox(
    const char* item) {
    HxSpatialGraphView* view = mEditor->getView();
    // Check if the current item is a label.
    if (view->isLabelGroupAttribute(item, HxSpatialGraph::VERTEX))
        return true;

    // Check if the current item is a vertex attribute with type int_32.
    // If so, we also need the labelcolor field to display it.
    HxSpatialGraph* graph = mEditor->getSpatialGraph();
    if (!graph)
        return false;
    EdgeVertexAttribute* attribute = dynamic_cast<EdgeVertexAttribute*>(
        mEditor->getSpatialGraph()
            ->findAttribute(HxSpatialGraph::VERTEX, item));
    if (!attribute)
        return false;
    McPrimType type = attribute->primType();
    if (McPrimType::mc_int32 == type)
        return true;
    return false;
}

void QxMicrotubuleAlignSpatialGraphTool::joinEnds() {
    McString matchingLabel(qPrintable(ui.wNodeColoringComboBox->currentText()));
    mAligner->joinMatchedEnds(matchingLabel.dataPtr());
    SpatialGraphSelection selectAll(mGraph);
    selectAll.selectAllVerticesAndEdges();
    selectAll.selectAllPoints();
    mEditor->setHighlightedElements(selectAll);
    mEditor->setVisibleElements(selectAll);
    if (mEditor) {
        const HxNeuronEditorSubApp::SpatialGraphChanges changes =
            HxNeuronEditorSubApp::SpatialGraphLabelTreeChange |
            HxNeuronEditorSubApp::SpatialGraphAttributeValueChange |
            HxNeuronEditorSubApp::SpatialGraphGeometryChange;
        mEditor->updateAfterSpatialGraphChange(changes);
    }
}

void QxMicrotubuleAlignSpatialGraphTool::removeAllIntermediate() {
    const McString matchingLabel(
        qPrintable(ui.wNodeColoringComboBox->currentText()));
    mAligner->removeIntermediateAndAddLabelForUnmatched(
        matchingLabel.dataPtr());
    SpatialGraphSelection selectAll(mGraph);
    selectAll.selectAllVerticesAndEdges();
    selectAll.selectAllPointsOnEdgesFromSelection(selectAll);
    mEditor->setHighlightedElements(selectAll);
    mEditor->setVisibleElements(selectAll);
    if (mEditor) {
        mGraph->touch(HxData::NEW_PARAMETERS);
        const HxNeuronEditorSubApp::SpatialGraphChanges changes =
            HxNeuronEditorSubApp::SpatialGraphLabelTreeChange |
            HxNeuronEditorSubApp::SpatialGraphAttributeValueChange |
            HxNeuronEditorSubApp::SpatialGraphGeometryChange;
        mEditor->updateAfterSpatialGraphChange(changes);
    }
}

void
QxMicrotubuleAlignSpatialGraphTool::copyViewerSettingsForManualTransform() {
    if (!mManualSlice) {
        return;
    }

    const HxSpatialGraphView* viewer = mEditor->getView();
    mcassert(viewer);

    mManualAligner->showNodes(viewer->portItemsToShow.getValue(0));
    mManualAligner->showSegments(viewer->portItemsToShow.getValue(1));

    const int nodeColorIndex = viewer->portVertexColoring.getIndex();
    const QString nodeColorAttributeName =
        viewer->portVertexColoring.getLabel(nodeColorIndex);
    const EdgeVertexAttribute* nodeColorAttribute =
        mManualSlice->getSpatialGraphInterface()->findVertexAttribute(
            qPrintable(nodeColorAttributeName));
    mManualAligner->setNodeColorAttribute(nodeColorAttribute);

    const int nodeScaleIndex = viewer->portVertexScaling.getIndex();
    const QString nodeScaleAttributeName =
        viewer->portVertexScaling.getLabel(nodeScaleIndex);
    const EdgeVertexAttribute* scaleAttrib =
        mManualSlice->getSpatialGraphInterface()->findVertexAttribute(
            qPrintable(nodeScaleAttributeName));
    mManualAligner->setNodeScaleAttribute(scaleAttrib);

    const float scaleFactor = viewer->portVertexScaleValue.getValue();
    mManualAligner->setNodeScaleFactor(scaleFactor);

    const int edgeColorIndex = viewer->portEdgeColoring.getIndex();
    const QString edgeColorAttributeName =
        viewer->portEdgeColoring.getLabel(edgeColorIndex);
    const EdgeVertexAttribute* edgeColorAttribute =
        mManualSlice->getSpatialGraphInterface()->findEdgeAttribute(
            qPrintable(edgeColorAttributeName));
    mManualAligner->setEdgeColorAttribute(edgeColorAttribute);

    const float segmentWidth = viewer->portSegmentWidth.getValue();
    mManualAligner->setSegmentWidth(segmentWidth);
}
