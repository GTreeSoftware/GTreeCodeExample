/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/


#ifndef NEURONTRACER_H
#define NEURONTRACER_H

#include <QMainWindow>
#include <QtGui>
#include <QtWidgets>
#include <QTimer>
#include <ctime>
#include <memory>

#include "Dialog/creathdf5dialog.h"///////////////////////////////////////////////////

#include "ngtypes/basetypes.h"
#include "ngtypes/ineurondataobject.h"
#include "ngtypes/ParamPack.h"
#include "Thread/binarythread.h"
#include "Function/IO/Prestore.h"
class BinaryFilter;
class BinarySettingsDockWidget;
class NeuroGLWidget;
class SparseTraceFilter;
class LargeSparseTraceFilter;
class INeuronBigReader;
class TraceStackReader;
class TraceStackWriter;
class HDF5Writer;
NG_SMART_POINTER_TYPEDEF(BinaryFilter, NGBinaryFilter);
NG_SMART_POINTER_TYPEDEF(SparseTraceFilter, NGSparseTraceFilter);
NG_SMART_POINTER_TYPEDEF(LargeSparseTraceFilter, NGLargeSparseTraceFilter);
NG_SMART_POINTER_TYPEDEF(INeuronBigReader, NGNeuronBigReader);
NG_SMART_POINTER_TYPEDEF(std::vector<VectorVec5d>, SMARTTREE);
NG_SMART_POINTER_TYPEDEF(TraceStackReader, NGTraceStackReader);
NG_SMART_POINTER_TYPEDEF(TraceStackWriter, NGTraceStackWriter);
NG_SMART_POINTER_TYPEDEF(HDF5Writer, NGHDF5Writer);

namespace Ui {
class SparseTracer;
}

class SparseTracer : public QMainWindow
{
    Q_OBJECT

public:
    explicit SparseTracer(QWidget *parent = 0);
    ~SparseTracer();

    bool SaveSeperateTreeAsOne(QString savePath);

signals:
//	void readyToHis();
	void renderCaliber();

private slots:
    void on_actionOpenImage_triggered();
    void on_actionSaveTree_triggered();
    void on_actionClear_triggered();
    void on_actionRun_triggered();
    void EndRunThread();
    //void on_actionStep_triggered();
    void on_actionStop_triggered();
    void UpdateStatus_Slot();
    void ApplyReadNewImage_Slot();
    void on_actionSaveImage_triggered();
    void on_actionSaveBack_triggered();
    //for big data
    void SetROIByBoxWidget_Slot();
    //
    void on_actionChoose_Line_toggled(bool);
    void on_actionChoose_Vertex_toggled(bool);
    void on_actionDelete_Line_triggered();
    void on_actionCut_Vertex_triggered();
    void on_actionDraw_Line_toggled(bool);
    void on_action2DView_toggled(bool);
    void on_actionZoom_in_triggered();
    void on_actionZoom_out_triggered();
    void on_actionVisible_toggled(bool);
    void on_actionSaveSVM_triggered();
    void ChangeMoveCursor_Slot(bool);//
    void TimingSave_Slot();
    void on_actionPick_Soma_toggled(bool arg);
    void on_actionSelect_Tree_toggled(bool arg);
    void on_actionNGTree_toggled(bool arg);
    void on_actionSBWT_toggled(bool arg);
    void on_actionNGTree_Trace_triggered();
    void on_actionDelete_Soma_triggered();
    void on_actionDelete_Tree_triggered();
    void on_actionTest_triggered();
    void EditActionGroup_triggered(QAction *);
    void SetProjectMaxThickness_Slot(int);
    void MaxBoundNumChanged_Slot(int);
    void Set2DViewStatus_Slot(bool);
    void SetDrawStatus_Slot(bool);
    void on_actionTrain_triggered();
    void on_actionTreeChecker_triggered();
    void UpdateGLBoxWidget_Slot();
    void on_actionLocal_Run_triggered();

    void on_actionSaveLayerImage_triggered();
    void on_actionSaveLayerSwc_triggered();
    void ClearMOSTDCache_Slot();//action_ClearCache_Button
    //tree widget
    void ActivateTree_Slot(int arg);
    void CompareTree_Slot(int arg);
    void ToggleTreeVisible_Slot(int);
    void GotoDiff_Slot(int);

    void ToggleTraverse_Slot(bool);
    void BackTraverse_Slot();
    void NextTraverse_Slot();
    void ToggleShowDiffArrow_Slot(bool);

    void StartTraverseFromHere_Slot(int, int);
    void ResetTraverse_Slot();
    //timer
    void on_actionStart_triggered();
    void on_actionHalt_triggered();
    void on_actionReset_triggered();
    void AddRunningTime_Slot();
    //
    void CacheComplete_Slot();
    void on_actionNeuroGPS_triggered();
    void on_actionSaveSoma_triggered();
    void on_actionOpenSoma_triggered();

	void on_actionSaveCaliber_triggered();//2018-3-22
    //
    void OpacAdjustApply_Slot();
    void PreviewOpacAdjust_Slot(bool);
    void PreviewBinary_Slot(bool);
    void PreviewAxonBinary_Slot(bool);

    void on_actionSnapshot_triggered();
    void on_actionCreate_HDF5_triggered();
    void EndCreateHDF5_Slot();
    void UpdateHDF5ProgressBar_Slot(int);
    void StopCreateHDF5_Slot();

protected:
    bool GetDrawLineIsChecked();
    bool CheckReadImageRangeValid();
    virtual void closeEvent(QCloseEvent * event);
    QString TranslateLayerListIntoText(const std::vector<int>& layerList);
    bool TranslateTextIntoLayerList(const QString&, std::vector<int>& layerList);
    bool ReadImage();
    int SetScale(double res);
    void ClearAll();

private:
    Ui::SparseTracer *ui;
    QGridLayout *grid;
    QFrame *frame;//style, for beautiful
    QHBoxLayout *frameLayout;
    BinarySettingsDockWidget *settingDock;
    QSplitter* splitter;
    NeuroGLWidget *glbox;
    QActionGroup *toolBarButtonGroup;
    QProgressDialog *progressDialog;
    
    //thread
    BinaryThread *worker;
    //status
    enum EDITTOOLBARTAG{ NONE, LINE, VERTEX, BOX, SOMA, TREE, CONNECT, DRAW } oldEditToolbarTag, newEditToolbarTag;
    bool hasSoma_ = true;
    int layerRadiusDisplay_ = 7;
    int sumSteps_;
    int statusTip_;
    int timerForSaveTime_ = 300000;
    QTimer *timer;
    QTimer *runningTimer_;
    clock_t beg, end;
    QString statusString_;
    int xScale_, yScale_, zScale_;
    //Data
    NGParamPack paramPack;
    //int xRangeMax_, yRangeMax_, zRangeMax_;
    double boxWidgetBound[6];
    QString imageFileName_;
    int traceStepLength_;
    IDataPointer soma;
    IDataPointer resTree_;
    IDataPointer treeConInfo;
    
    SMARTTREE manualLabelCurve_;
    VectorVec3i binPtSet;
    IDataPointer m_Source;//test

    NGBinaryFilter filter;
    NGSparseTraceFilter sparseFilter;
    NGLargeSparseTraceFilter largeFilter;
    NGNeuronBigReader mostdReader;
    NGNeuronBigReader mostdTraverseReader;
    //
    NGHDF5Writer writer;
    //date time
    QStringList previousSaveTime_;
    bool traverseMode_ = false;
    NGPrestore preStore;
    //
    QDialog *snapDialog;
};

#endif // NEURONTRACER_H
