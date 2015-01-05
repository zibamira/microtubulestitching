#include <hxalignmicrotubules/HxTestPointMatching.h>

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QLineEdit>
#include <QString>
#include <QToolBar>

#include <hxcore/HxMain.h>
#include <hxspatialgraph/HxSpatialGraph.h>

namespace SteerFE {

void qxtSelectFilamentEditor() {
    QList<QToolBar*> subAppToolBar =
        theMain->findChildren<QToolBar*>("subApplicationToolBar");
    mcassert(1 == subAppToolBar.size());
    QList<QAction*> allActions = subAppToolBar[0]->findChildren<QAction*>();
    QList<QAction*> filamentEditorAction;
    for (int i = 0; i < allActions.size(); i++) {
        if (allActions[i]->text() == QString("Filament Editor"))
            filamentEditorAction.append(allActions[i]);
    }
    mcassert(1 == filamentEditorAction.size());
    filamentEditorAction[0]->trigger();
}

QAction* findNamedAction(QToolBar* toolbar) {
    foreach (QAction* a, toolbar->actions()) {
        if (a->text() == "MicrotubuleAlignTool") {
            return a;
        }
    }
    return 0;
}

void qxtSelectAlignTool() {
    QList<QToolBar*> toolbar =
        theMain->findChildren<QToolBar*>("FilamentEditorToolboxToolbar");
    mcassert(1 == toolbar.size());
    QAction* alignToolAction = findNamedAction(toolbar[0]);
    mcassert(NULL != alignToolAction);
    alignToolAction->trigger();
}

void qxtActivateAlignTool() {
    qxtSelectFilamentEditor();
    qxtSelectAlignTool();
}

void qxtDeactivateFilamentEditor() {
    QList<QToolBar*> subAppToolBar =
        theMain->findChildren<QToolBar*>("subApplicationToolBar");
    mcassert(1 == subAppToolBar.size());
    QList<QAction*> allActions = subAppToolBar[0]->findChildren<QAction*>();
    QList<QAction*> opEditorAction;
    for (int i = 0; i < allActions.size(); i++) {
        if (allActions[i]->text() == QString("Object Pool"))
            opEditorAction.append(allActions[i]);
    }
    mcassert(1 == opEditorAction.size());
    opEditorAction[0]->trigger();
}

class QxtButton {
  public:
    QAbstractButton* const q;

    explicit QxtButton(QAbstractButton* w) : q(w) {}

    void check() {
        mcassert(q->isCheckable());
        q->setChecked(true);
    }

    void uncheck() {
        mcassert(q->isCheckable());
        q->setChecked(false);
    }

    void click() { q->click(); }
};

class QxtCheckBox {
  public:
    QCheckBox* const q;

    explicit QxtCheckBox(QCheckBox* w) : q(w) {}

    void check() {
        mcassert(q->isCheckable());
        q->setChecked(true);
    }

    void uncheck() {
        mcassert(q->isCheckable());
        q->setChecked(false);
    }

    void checkPartially() {
        mcassert(q->isTristate());
        q->setCheckState(Qt::PartiallyChecked);
    }

    void click() { q->click(); }
};

class QxtLineEdit {
  public:
    QLineEdit* const q;

    explicit QxtLineEdit(QLineEdit* w) : q(w) {}

    void setText(const char* text) { q->setText(text); }
};

class QxtComboBox {
  public:
    QComboBox* const q;

    explicit QxtComboBox(QComboBox* w) : q(w) {}

    void setCurrentIndex(const int theIndex) { q->setCurrentIndex(theIndex); }
};

class QxtWidget {
  public:
    QWidget* const q;

    explicit QxtWidget(QWidget* w) : q(w) {}

    QxtWidget find(const char* name) {
        QList<QWidget*> widgets = q->findChildren<QWidget*>(name);
        mcassert(1 == widgets.size());
        return QxtWidget(widgets[0]);
    }

    QxtButton findButton(const char* name) {
        QList<QAbstractButton*> button =
            q->findChildren<QAbstractButton*>(name);
        mcassert(1 == button.size());
        return QxtButton(button[0]);
    }

    QxtLineEdit findLineEdit(const char* name) {
        QList<QLineEdit*> te = q->findChildren<QLineEdit*>(name);
        mcassert(1 == te.size());
        return QxtLineEdit(te[0]);
    }
    QxtComboBox findComboBox(const char* name) {
        QList<QComboBox*> te = q->findChildren<QComboBox*>(name);
        mcassert(1 == te.size());
        return QxtComboBox(te[0]);
    }
    QxtCheckBox findCheckBox(const char* name) {
        QList<QCheckBox*> te = q->findChildren<QCheckBox*>(name);
        mcassert(1 == te.size());
        return QxtCheckBox(te[0]);
    }
};

QxtWidget qxtRoot() { return QxtWidget(theMain); }

class QxtAlignSpatialGraphTool {
  public:
    QxtWidget tool;

    QxtAlignSpatialGraphTool(QxtWidget w) : tool(w) {}

    QxtAlignSpatialGraphTool defaults() {
        tool.findButton("allRadioButton").check();
        tool.findButton("startTMinAngleCheckBox").check();
        tool.findButton("showAdvancedOptionsCheckBox").check();
        tool.findLineEdit("maxPointDistLineEdit").setText("1000");
        tool.findLineEdit("alphaLineEdit").setText("0.0001");
        tool.findLineEdit("minAngleLineEdit").setText("60");
        tool.findButton("angleForMatchingAutoCheckBox").check();
        tool.findCheckBox("GapSizeAutoCheckBox").uncheck();
        tool.findLineEdit("scaleTestMin").setText("1");
        tool.findLineEdit("scaleTestMax").setText("1");
        tool.findLineEdit("angularMatchingSmoothing").setText("100000");
        tool.findLineEdit("dMinAngleForOptMatching").setText("170");
        tool.findCheckBox("usePGMCheckBox").uncheck();
        return *this;
    }
    QxtAlignSpatialGraphTool PGMdefaults() {
        tool.findButton("allRadioButton").check();
        tool.findButton("showAdvancedOptionsCheckBox").check();
        tool.findLineEdit("maxPointDistLineEdit").setText("600");
        tool.findLineEdit("alphaLineEdit").setText("0.0001");
        tool.findLineEdit("minAngleLineEdit").setText("60");
        tool.findButton("angleForMatchingAutoCheckBox").check();
        tool.findCheckBox("GapSizeAutoCheckBox").uncheck();
        tool.findLineEdit("scaleTestMin").setText("1");
        tool.findLineEdit("scaleTestMax").setText("1");
        tool.findLineEdit("angularMatchingSmoothing").setText("2000");
        tool.findLineEdit("dMinAngleForOptMatching").setText("160");
        tool.findCheckBox("usePGMCheckBox").check();
        tool.findComboBox("projectionComboBox").setCurrentIndex(7);
        tool.findLineEdit("angularMatchingSmoothing").setText("2000");
        return *this;
    }

    QxtAlignSpatialGraphTool allAngles() {
        tool.findButton("minAngleCheckBox").uncheck();
        return *this;
    }

    QxtAlignSpatialGraphTool limitAngles() {
        tool.findButton("minAngleCheckBox").check();
        return *this;
    }
    QxtAlignSpatialGraphTool pairWise() {
        tool.findButton("pairRadioButton").check();
        return *this;
    }

    QxtAlignSpatialGraphTool withScaling() {
        tool.findLineEdit("scaleTestMin").setText("0.95");
        tool.findLineEdit("scaleTestMax").setText("1.05");
        tool.findLineEdit("scaleTestIncrement").setText("0.01");
        return *this;
    }

    QxtAlignSpatialGraphTool noScaling() {
        tool.findLineEdit("scaleTestMin").setText("1");
        tool.findLineEdit("scaleTestMax").setText("1");
        return *this;
    }

    QxtAlignSpatialGraphTool withLabels() {
        tool.findButton("createMatchingLabelsCheckBox").check();
        return *this;
    }

    QxtAlignSpatialGraphTool noLabels() {
        tool.findButton("createMatchingLabelsCheckBox").uncheck();
        return *this;
    }

    QxtAlignSpatialGraphTool tangentApprox(const char* approxDist) {
        tool.findComboBox("projectionComboBox").setCurrentIndex(7);
        tool.findLineEdit("angularMatchingSmoothing").setText(approxDist);
        return *this;
    }

    QxtAlignSpatialGraphTool varyGap(const char* minGap, const char* maxGap,
                                     const char* stepSize) {
        tool.findCheckBox("GapSizeAutoCheckBox").checkPartially();
        tool.findLineEdit("dMinLineEdit").setText(minGap);
        tool.findLineEdit("dMaxLineEdit").setText(maxGap);
        tool.findLineEdit("deltaHLineEdit").setText(maxGap);
        return *this;
    }

    QxtAlignSpatialGraphTool align() {
        tool.findButton("alignButton").click();
        return *this;
    }
    QxtAlignSpatialGraphTool usePGM() {
        tool.findCheckBox("usePGMCheckBox").check();
        return *this;
    }
};

QxtAlignSpatialGraphTool qxtAlignSpatialGraphTool() {
    return QxtAlignSpatialGraphTool(
        qxtRoot().find("MicrotubuleAlignSpatialGraphTool"));
}

}  // namespace

HX_INIT_CLASS(HxTestPointMatching, HxCompModule);

HxTestPointMatching::HxTestPointMatching()
    : HxCompModule(HxSpatialGraph::getClassTypeId()),
      portMatchingAlgorithmType(this, "algorithmType", 3),
      portProjectionType(this, "projectionType", 1),
      portThresholds(this, "thresholds", 3),
      portParams(this, "parameters", 3),
      portUseParams(this, "useParamsForWeights", 3),
      portPGMPairParam(this, "pgmPairParam", 1),
      portPGMDummieSignificance(this, "pgmDummieSignificance", 1),
      portPairDummy(this, "pairDummy", 1),
      mDoIt(this, "apply") {

    portMatchingAlgorithmType.setLabel(0, "Greedy");
    portMatchingAlgorithmType.setLabel(1, "Exact");
    portMatchingAlgorithmType.setLabel(2, "PGM");

    portThresholds.setLabel(0, "3dDistance");
    portThresholds.setLabel(1, "projectedDistance");
    portThresholds.setLabel(2, "angle");
    portThresholds.setValue(0, 1500);
    portThresholds.setValue(1, 600);
    portThresholds.setValue(2, 20);

    portUseParams.setLabel(0, "3dDistance");
    portUseParams.setLabel(1, "projectedDistance");
    portUseParams.setLabel(2, "angle");

    portParams.setLabel(0, "3dDistance");
    portParams.setLabel(1, "projectedDistance");
    portParams.setLabel(2, "angle");
    portParams.setValue(0, 1000);
    portParams.setValue(1, 300);
    portParams.setValue(2, 10);
    portPGMPairParam.setValue(0, 100);

    portUseParams.setValue(0, 0);
    portUseParams.setValue(1, 1);
    portUseParams.setValue(2, 1);

    portPGMDummieSignificance.setValue(0, 0.01);
}

HxTestPointMatching::~HxTestPointMatching() {}

void HxTestPointMatching::compute() {

    if (!mDoIt.wasHit()) {
        return;
    }

    if (portData.source() == NULL) {
        return;
    }

    // Recompute all
    int alogithmIdx = portMatchingAlgorithmType.getValue(0);
    McString algorithmName =
        qPrintable(QString(portMatchingAlgorithmType.getLabel(alogithmIdx)));
    HxSpatialGraph* origGraph =
        dynamic_cast<HxSpatialGraph*>(portData.source());
    HxSpatialGraph* sg = origGraph->duplicate();

    if (sg == NULL) {
        return;
    }

    SteerFE::qxtActivateAlignTool();
    SteerFE::QxtWidget mainWidget(SteerFE::qxtRoot().find("QxNeuronEditor"));
    SteerFE::QxtComboBox box = mainWidget.findComboBox("wGraphDataComboBox");
    int ind = box.q->findText(QString(origGraph->getLabel().getString()));
    box.q->setCurrentIndex(ind);
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findButton("allRadioButton")
        .check();
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findButton("radioButtonOptAlign")
        .check();
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findComboBox("projectionComboBox")
        .setCurrentIndex(portProjectionType.getValue());

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findComboBox("transformComboBox")
        .setCurrentIndex(3);
    if (portUseParams.getValue(0))
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("use3dDistanceWeight")
            .check();
    else
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("use3dDistanceWeight")
            .uncheck();

    if (portUseParams.getValue(1))
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("useProjectedDistanceWeight")
            .check();
    else
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("useProjectedDistanceWeight")
            .uncheck();

    if (portUseParams.getValue(2))
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("useAngleWeight")
            .check();
    else
        SteerFE::qxtAlignSpatialGraphTool()
            .tool.findButton("useAngleWeight")
            .uncheck();

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("threeDDistanceWeight")
        .setText(QString::number(portParams.getValue(0)).toLocal8Bit().data());
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("projectedDistanceWeight")
        .setText(QString::number(portParams.getValue(1)).toLocal8Bit().data());
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("angleWeight")
        .setText(QString::number(portParams.getValue(2)).toLocal8Bit().data());

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("threeDDistanceThres")
        .setText(
            QString::number(portThresholds.getValue(0)).toLocal8Bit().data());
    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("projectedDistanceThres")
        .setText(
            QString::number(portThresholds.getValue(1)).toLocal8Bit().data());
    SteerFE::qxtAlignSpatialGraphTool().tool.findLineEdit("angleThres").setText(
        QString::number(portThresholds.getValue(2)).toLocal8Bit().data());

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("pairFactorParam")
        .setText(
            QString::number(portPGMPairParam.getValue(0)).toLocal8Bit().data());

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findComboBox("pmAlgorithmComboBox")
        .setCurrentIndex(portMatchingAlgorithmType.getValue());

    SteerFE::qxtAlignSpatialGraphTool()
        .tool.findLineEdit("dummySignificance")
        .setText(QString::number(portPGMDummieSignificance.getValue(0))
                     .toLocal8Bit()
                     .data());

    SteerFE::qxtAlignSpatialGraphTool().align();
    SteerFE::qxtDeactivateFilamentEditor();
}

int HxTestPointMatching::parse(Tcl_Interp* t, int argc, char** argv) {
    return HxCompModule::parse(t, argc, argv);
}

void HxTestPointMatching::update() {
    if (portParams.isNew() || portPGMPairParam.isNew() ||
        portPGMDummieSignificance.isNew()) {
        const double dist3dParam = portParams.getValue(0);
        const double distProjectedParam = portParams.getValue(1);
        const double angleWeightParam = portParams.getValue(2);
        const double pairParam = portPGMPairParam.getValue(0);
        const double dummySignificance = portPGMDummieSignificance.getValue(0);
        const double dist3dCutoff = -1.0 * dist3dParam * log(dummySignificance);
        const double distProjectedCutOff =
            -1.0 * distProjectedParam * log(dummySignificance);
        const double angleCutoff =
            -1.0 * angleWeightParam * log(dummySignificance);
        const double pairCutoff = -1.0 * pairParam * log(dummySignificance);
        // compute thresholds from parameters
        portThresholds.setValue(0, dist3dCutoff);
        portThresholds.setValue(1, distProjectedCutOff);
        portThresholds.setValue(2, angleCutoff);
        portPairDummy.setValue(0, pairCutoff);
        portThresholds.setSensitivity(0, 0);
        portThresholds.setSensitivity(1, 0);
        portThresholds.setSensitivity(2, 0);
        portPairDummy.setSensitivity(0, 0);
    }
}
