
/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef NEUROGLWIDGET_H
#define NEUROGLWIDGET_H

//----- VTK

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingContextOpenGL2)
VTK_MODULE_INIT(vtkRenderingFreeType)
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2)
#include <QVTKWidget.h>
#include <QVTKOpenGLWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkImageData.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkSmartPointer.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkLODProp3D.h>
#include <vtkProperty.h>
#include <vtkRegularPolygonSource.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkFollower.h>
#include <vtkProp3DFollower.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeRayCastCompositeFunction.h>
#include <vtkVolumeRayCastMapper.h>
#include <vtkVolumeRayCastIsosurfaceFunction.h>
#include <vtkVolumeRayCastMIPFunction.h>
#include <vtkImageImport.h>
#include <vtkCamera.h>
#include <vtkTexture.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkWorldPointPicker.h>
#include <vtkCellPicker.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkBoxWidget.h>
#include <vtkAreaPicker.h>
#include <vtkImageTracerWidget.h>
#include <vtkImageActor.h>
#include <vtkPointPicker.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkPointPicker.h>
#include <vtkEvent.h>
#include <vtkIdFilter.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkPlanes.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyLine.h>
#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <QMouseEvent>
#include <QAction>
#include <QCursor>
#include <QColor>
#include <vector>
#include <set>
#include <memory>
#include "ngtypes/volume.h"
#include "ngtypes/basetypes.h"
#include "ngtypes/ineurondataobject.h"
#include "nginteractorstyletrackballcamera.h"
#include "ngtypes/ParamPack.h"
#include "ngtypes/tree.h"
#include "Function/Trace/SemiAutoTracer.h"
#include <vtkIndent.h>
#define LINEWIDTH_NOR 2
#define POINTSIZE_NOR 3
#define LINEWIDTH_CHOOSEN 5
#define VTK_NEW(type, var) var = vtkSmartPointer<type>::New()
#define VTK_CREATE(type, var) vtkSmartPointer<type> var = vtkSmartPointer<type>::New()
NG_SMART_POINTER_TYPEDEF(std::vector<VectorVec5d>, SMARTTREE);

class Soma;
//struct CurveInBox;

class NeuroGLWidget :public QVTKOpenGLWidget//QVTKWidget//public QGLWidget
{
    Q_OBJECT
public:
    explicit NeuroGLWidget(QWidget *parent = 0);
    ~NeuroGLWidget();
    enum EDITMODE{ FORBIDDEN, NORMALMODE, BRANCHMODE, POINTMODE, SOMAMODE, TREEMODE, DRAWLINEMODE, TRAVERSEMODE };
    virtual void SetInputNeuronPopulation(IDataPointer);
    virtual void SetInputVolume(IDataPointer, DATATYPE = DATATYPE::IMAGE16);
    void UpdateSomaList();
    void BuildCurrentBoundaryCollection();//only display collected active tree boundary info
    void BuildHistoryBoundaryActorCollection();//only display active tree history boundary info
    void RebuildNeuronTreeActor(std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>&, bool = false);
    void RebuildActiveTreeActor();
    void RebuildLayerTreeActor(std::vector<int>& activeLayer);
    void SetParam(NGParamPack arg){ paramPack = arg; semiTracer_->SetParam(paramPack); }
    //big data viewer
    bool Is2DView()const{ return is2DView_; }
    bool IsIn2DDrawMode()const{ return mode_ == DRAWLINEMODE; }
    double GetBoundaryDistanceThreshold() const { return paramPack->boundaryDistanceThreshold_; }
    //int GetLevel()const{ return paramPack->mostdLevel_; }
    virtual void RemoveAllActor();
    void BuildCurInitTracePt();
    void BuildSphereActor(vtkSmartPointer<vtkActor>, const Vec3d&);
    void ClearAllActorRenderActor();
    void GetBoxWidgetBound(double bound[])const;
    void On2DMouseWheel(bool);
    void RemoveInitPtActor(){ if (initPtActor_) actorRenderer_->RemoveActor(initPtActor_); }
    void SentBoxWidgetChangedSignal();
    void SentChangeMoveCursorSignal_Slot(bool);
    void SetNormalMode(){ SetEditMode(NORMALMODE); }//choose edit mode
    void SetScale(int, int, int);
    void TwoDViewZoom(bool zoomIn);
    //void UpdateOffset();
    std::vector<int>& GetLayerDisplay()  { return layerDisplay_; }
    std::vector<int> GetLayerDisplay() const  { return layerDisplay_; }
    std::vector<int>& GetBranchDisplay()  { return branchDisplay_; }
    std::vector<int> GetBranchDisplay() const  { return branchDisplay_; }
    //void SetLayerDisplay(const std::vector<int>& val) { layerDisplay_ = val; }
    vtkSmartPointer<vtkRenderer> GetActiveRenderer();
    NeuroGLWidget::EDITMODE GetMode(){ return mode_; }
    void DeactiveAllTrees();
    bool IsMovePointMode(){ return isMovePointMode_; }
    void SetMovePointPos(double x, double y, double z);
    //sulei
    bool AddCellId2PickedArray(vtkSmartPointer<vtkActor> pickActor, vtkIdType idArg);
    void ChangePickedCellColor( vtkIdType , bool isChanged);
    void ClearPickedCellsArray();
    void EnterBranchMode();
    void EnterPointMode();
    void ExitBranchMode();
    void ExitPointMode();
    void InitImageData(int w, int h, VTK_DEFINE(vtkImageData, imArg), int wM, int hM);
    void MouseClick();
    void ShowVertex( vtkIdType arg, bool flag);
    //modified 20161205
    bool CollectionBoundaryPt(BOUNDINFOLIST&);
    const BOUNDINFOLIST& CurBoundInfo() const { return curBoundInfo; }
    void ResetCurForbiddenBoundInfo(){ curForbiddenFoundBoundInfo.clear(); }
    std::vector<std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer> >& GetTreeActorCollection(){ return neuronActorList_; }
    void ChangePickedVertColor(std::pair<vtkSmartPointer<vtkActor>, vtkIdType>& arg, bool isChanged);
    vtkSmartPointer<vtkAreaPicker> GetAreaPicker(){ return areaPicker_; }
    vtkSmartPointer<vtkRenderer> GetVolumeRenderer(){ return volumeRenderer_; }
    vtkSmartPointer<vtkRenderer> GetActorRenderer(){ return actorRenderer_; }
    void ToggledForbidden(bool arg){ if (arg){ bakmode_ = mode_; mode_ = FORBIDDEN; } else mode_ = bakmode_; }
    bool IsForbidden()const{ return mode_ == FORBIDDEN; }
    //tree widget
    void RebuildActiveTreeChangedPopulation(int arg); 
    void ToggleTargetTreeVisible(int);
    void FocusOnPoint(Vec5d &pointArg);
    void ToggleTraverseMode(bool);
    void ToggleDiffArrows(bool);
    void LightDiffArrow(int);
    void UpdateTraverseFlagActor();

    //qt5 action
    QAction *pInteract, *pXYReset, *pXZReset, *pYZReset, *pChooseBranchMode, *pDeleteBranch,
        *pChooseNodeMode, *pCutLine, *pDrawLineMode, *pChooseSomaMode, *pDeleteSoma,
        *pChooseTreeMode, *pDeleteTrees, *pUndoDraw, *pMergeTree,
        *pSetHeadAsTreeRoot, *pSetTailAsTreeRoot, *pSaveTree, *pSetActiveTree, *pUnDel, *pShowRoot,
        *pShowDirection, *p2DView, *pShowPosition, *pSeperateBranch, *pConjChild, *pToggleEnableTrace,
        *pToggleCurInitPtShow, *pUndelTrees, *pBoundCheck, *pForceBoundCheckAll, *pSmartSemiTrace, //*pActiveTreeLayerDisplay,
        *pAddFlag, *pDelFlag, *pShowLevelBranchId, *pStartTraverseHere_, *pSmoothCurve, *pSaveImageAndSwcForTest,
        *pInsertPt_, *pDelPt_, *pStartMovePt_, *pReconnectParent_, *pRefreshTree_, *pCaliberBuild_,*pAllBuilder_, *pSetDendrite_, *pSetAxon_,
        *pHideDendrite_, *pHideAxon_, *pHideCaliber_;
    //
    EDITMODE mode_, bakmode_ = NORMALMODE;

    void Snapshot(QString& path, int vx);//int vy, 

protected:
    virtual void keyPressEvent(QKeyEvent* event);
    void ReadMemoryImage(std::shared_ptr<CVolume>);
    void ReadMemoryImage(std::shared_ptr<SVolume>, DATATYPE);
    void ReadMemoryImage(std::shared_ptr<IVolume>);
    int GetSampleRate();
    void contextMenuEvent(QContextMenuEvent* event);
    void CreateActions();
    void PickPoint(double *xyz);
    void SetVolumeRenderData(vtkSmartPointer<vtkImageImport>);
    vtkSmartPointer<vtkActor> CreateArrowActor(const Vec3d&, const Vec3d&, bool = false);
    vtkSmartPointer<vtkDoubleArray> GetXYZScale();
    //interactive
    void ClearActorCollection(std::vector<vtkSmartPointer<vtkActor> > &arg);
    void ClearActorCollection(vtkSmartPointer<vtkActorCollection> arg);
    void CreatePointChooseActor(double *p = 0);
    void EnterDrawLineMode();
    void EnterSomaMode();
    void EnterTreeMode();
    void ExitDrawLineMode();
    void ExitSomaMode();
    void ExitTreeMode();
    void SetEditMode(NeuroGLWidget::EDITMODE);
    // new version!!!!! sulei
    double PickXYSliceMaxPoint(double pickPoint[3]);
    double PickXZSliceMaxPoint(double pickPoint[3]);
    double PickYZSliceMaxPoint(double pickPoint[3]);
    void AppendNeuronTreeActor(NeuronTreePointer, bool = false);
    void BuildLineActor(vtkSmartPointer<vtkActor> actor, const Line5d& line, bool flag);
    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir);
    void ChangePickedVertColor(vtkIdType idArg, bool isChanged);
    void NGResetCamera(bool);
    void ProjectPick(double pickPoint[3]);
    void ResetVTKFollowerCamera();
    void Set2dCamera(int);
    void SetActorCollectionVisible(std::vector<vtkSmartPointer<vtkActor> >, bool);
    vtkSmartPointer<vtkActor> BuildTreeActor(NeuronTree &nTree, bool active = false, bool layer = false);
    int BuildSubTree(NeuronTreePointer newNeuronTree);
    bool MergeTwoTrees(NeuronTreePointer arg1, NeuronTreePointer arg2);
    void errorcheck();
    void ForceBoundCollectManualCurves();
    bool CheckValidBeforeBoundaryCollect(const Vec3d&);
    bool IsInCurrentImageRange(const Vec3d&);
    void  BuildVertexActor(const VectorVec5d& tmpPtSet, vtkSmartPointer<vtkActor>);
    bool Draw3DPoint(double pickPoint[3]);
    bool GetCurveLevelAndBranch(int &, int&, NEURITETYPE &);
    void  BuildSomaActor(const Cell& tmpSoma, vtkSmartPointer<vtkFollower>);
    void BuildLocalCaliberActor(const CellVec3d &localCaliber, vtkSmartPointer<vtkActor> actor);
    void HideDenAxonBase(bool, NEURITETYPE);

signals:
    void BoxWidgetChanged_Signal();
    void ChangeMoveCursor_Signal(bool);
    void clippingDistanceChanged_Signal(int);
    void InteractModeChanged_Signal(bool);
    void Set2DProject_Signal(int);
    void Set2DView_Signal(bool);
    void SetDraw_Signal(bool);
    void SetTraversePos_Signal(int, int);
    void NeedUpdateTreeWidget();
    void NeedUpdateCheckWidget();
    void NextTraverse_Signal();

public slots:
    void DeleteBranch();
    void SeparateBranch();
    void DeleteSoma();
    void DeleteTrees();
    void DivideLine();
    void MergeTwoPickedTree();
    void SaveSelectedTrees();
    void SetActiveTree_Slot();
    void SetClippingDistance_Slot(int arg);
    void SetHeadAsTreeRoot();
    void SetInteractEnable_Slot(bool);
    void SetTailAsTreeRoot();
    void Toggle2D3DView(bool arg);
    void Toggle2DDraw(bool arg){ SetEditMode(arg ? DRAWLINEMODE : NORMALMODE); }
    void ToggleLineMode(bool arg){ SetEditMode(arg ? BRANCHMODE : NORMALMODE); }
    void TogglePointMode(bool arg){ SetEditMode(arg ? POINTMODE : NORMALMODE); }
    void ToggleShowCurrentLinesDirection(bool);
    void ToggleShowTreeRoot(bool arg);
    void ToggleSomaMode(bool arg){ SetEditMode(arg ? SOMAMODE : NORMALMODE); }
    void ToggleTreeMode(bool arg){ SetEditMode(arg ? TREEMODE : NORMALMODE); }
    void ToggleTreesVisible(bool arg);
    void UnDeleteOneActiveTreeBranch();
    void UndoDraw();
    void UpdateBoxWidget_Slot();
    void UpdateBoxWidgetWithoutReset_Slot(double *);
    void UpdateVolumeOpac_Slot();
    void XYReset_Slot();
    void XZReset_Slot();
    void YZReset_Slot();
    /*merge or connect the first line to the second. move the children of the first to the second. donot care about the parent of the first line.*/
    void ShowPosition_Slot();
    void ConjChildCurve();
    void ToggleEnableCurveTrace();
    void ToggleShowInit(bool arg);
    void UndeleteTrees();
    void ToggleBoundCheck(bool arg){ isBoundCheck_ = arg; BuildCurrentBoundaryCollection(); }
    void ForceBoundCheckAll();
    void ToggleSmartSemiTracer(bool arg){ isSmartSemiTracer_ = arg; }
    void ToggleActiveTreeLayerDisplay(bool arg);
    void Annotation_Slot();
    void Unannotation_Slot();
    void ShowLevelBranchID();
    void StartTraverseFromHere_Slot();
    void InsertPointBefore_Slot();
    void RemovePoint_Slot();
    bool ToggleMovePoint_Slot(bool);
    void ToggleReconParent_Slot(bool);
    void RefreshTree_Slot();
    void BuildCaliber_Slot();
    void SetBranchAsDenrite_Slot();
    void SetBranchAsAxon_Slot();
    void HideDendrite_Slot(bool);
    void HideAxon_Slot(bool);
    void HideCaliber_Slot(bool);
	//khtao
	void SmoothCurrentCurve_Slot();
	void SaveImageAndSwcForTest_Slot();
    void ToggleShowTraverseFlag_Slot();

	//huale
	void renderCaliber_Slot();
	void AllCaliberBulider();
    
private:
    //flag
    bool bHasStartPoint_;
    bool is2DView_;
    bool isDrawHasStartPoint_;
    bool isTreeVisible_;
    double *drawOldPoint_;
    bool isBoundCheck_ = true;
    bool isForceBoundCheckAll_ = false;
    bool isSmartSemiTracer_ = true;
    bool isActiveTreeLayerDisplay_ = false;
    bool isMovePointMode_ = false;
    bool isMoveTreeRootMode_ = false;
    bool isReconnectParent_ = false;
    int eventPos_[2];
    //draw line id
    size_t drawLinePID_, drawVertexPID_;
    int pickLineNum = 0;
    VectorVec3d pick3DLineList_;
    vtkSmartPointer<vtkActor> pickLineActor_;
    vtkSmartPointer<vtkActor> pickVertexActor_;
    //volume info
    int volumeX_;
    int volumeY_;
    int volumeZ_;
    int xScale_, yScale_, zScale_;
    int ratio_;
    //image_type
    NGParamPack paramPack;
    //volume data
    IDataPointer volumeData;
    //Sample Rate
    int sampleRate_;//0 : no sample ; 1 : 1/4 rate ; 2: 1/16 rate
    int sampleRateRatio_[3];
    //m_Trees
    IDataPointer m_Source;
    NG_SMART_POINTER_DEFINE(Soma, manualSoma_);
    Vec3d manualStartPoint_;
    //layer display
    std::vector<int> layerDisplay_;
    std::vector<int> branchDisplay_;
    std::vector<size_t> levelCurveID_;
    //std::vector<size_t > allVisibleLayerCurveID_;
    //old color
    std::deque<uchar *> oldColor_;
    //offset
    //double xOffset_, yOffset_, zOffset_;
    struct Camera2DParam{
        int xyz_;
        double dFromCamera2BlockEnd_;
        double clippingDistance_;
        double parallelScale_;
        double zoomRatio_;
        double moveStep_;
        double rangeStep_;
    } camera2dParam_;
    //boundary
    BOUNDINFOLIST curBoundInfo;
    VectorVec3d curForbiddenFoundBoundInfo;
    //annotation object
    std::vector < std::pair<Vec3d, vtkSmartPointer<vtkFollower> > > annotationFlagList_;
    typedef  std::pair<Vec3d, vtkSmartPointer<vtkFollower> > annotationFlagClass;
    //for parent layer display and deactivate tree 
    vtkSmartPointer<vtkActor> otherLayerActor_;
    //vtk display object
    vtkSmartPointer<vtkTransform> sliceTransform_;
    vtkSmartPointer<vtkCamera> camera2d_;
    vtkSmartPointer<vtkCamera> camera3d_;
    vtkSmartPointer<vtkCamera> cameraBox_;
    vtkSmartPointer<vtkRenderer> actorRenderer_;
    vtkSmartPointer<vtkRenderer> volumeRenderer_;
    vtkSmartPointer<vtkOrientationMarkerWidget> marker_;
    vtkSmartPointer<NGInteractorStyleTrackballCamera> interactorStyle_;
    vtkSmartPointer<NGInteractorObjectStyle> interactorObjectStyle_;
    vtkSmartPointer<vtkBoxWidget> boxWidget_;
    vtkSmartPointer<vtkActor> imageBoxActor_;
    vtkSmartPointer<vtkActor> initPtActor_;
    vtkSmartPointer<vtkActor> treeRootActor_;
    std::vector<vtkSmartPointer<vtkActor> > lineDirectionActor_;
    vtkSmartPointer<vtkActor> localCurveCaliberActor_;
    //soma object
    std::vector<vtkSmartPointer<vtkActor> > somaActorCollection_;
    //pick data
    vtkIdType pickedCellId_;
    std::vector<vtkIdType> pickedLineList_;//active tree cell id
    std::pair<vtkIdType, vtkIdType> pickedVertexID_;//active tree cell id point id
    std::vector<std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer> > pickedTreeList_;//actor and NeuronTree
    //pick object
    vtkSmartPointer<vtkActor> pickVertActor_;
    vtkSmartPointer<vtkFollower> choosedPointActor_;
    vtkSmartPointer<vtkFollower> choosedSomaActor_;
    std::vector<vtkSmartPointer<vtkActor> > choosedSomaActorCollection_;
    std::vector<vtkSmartPointer<vtkActor> > pickRayActorCollection_;
    vtkSmartPointer<vtkActor> drawLineActor_;
    vtkSmartPointer<vtkActor> drawVertActor_;
    vtkActor* activeTreeActor_ = nullptr;
    std::vector<std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer> > neuronActorList_;
    //boundary actor
    std::vector<vtkSmartPointer<vtkActor> > currentBoundaryActorCollection_;
    std::vector<vtkSmartPointer<vtkActor> > historyBoundaryActorCollection_;
    //picker
    vtkSmartPointer<vtkCellPicker> cellPicker_;
    vtkSmartPointer<vtkWorldPointPicker> worldPicker_;
    vtkSmartPointer<vtkAreaPicker> areaPicker_;
    vtkSmartPointer<vtkImageTracerWidget> imageTracer_;
    vtkSmartPointer<vtkImageActor> traceImageActor_;
    //volume mapper;
    vtkSmartPointer<vtkVolume> oldVolume_;
    vtkSmartPointer<vtkActor> somaActor_;
    vtkSmartPointer<vtkActor> outlineActor_;
    vtkSmartPointer<vtkActor> sliceOutlineActor_;
    //transform
    vtkSmartPointer<vtkTransform> transformForVolume;
    vtkSmartPointer<vtkTransformFilter> transformModelForVolume;
    //for new operator version sulei
    VTK_DEFINE(vtkActor, vertsActor_);
    VTK_DEFINE(vtkPointPicker, pointPicker_);
    //bool canPickPoint_;// isDrawing_;
    NGSemiAutoTracer semiTracer_;
    //to check golden standard
    std::vector<vtkSmartPointer<vtkActor> > diffArrowList_;
    vtkSmartPointer<vtkActor> traverseFlagActor_;
    //move point
    vtkSmartPointer<vtkActor> movePointActor_;
};

class vtkWidgetsCallback : public vtkCommand
{
public:
    vtkWidgetsCallback(){ pVTKWidget = NULL; t = vtkTransform::New(); }
    ~vtkWidgetsCallback(){ t->Delete(); }
    static vtkWidgetsCallback *New()
    {
        return new vtkWidgetsCallback;
    }
    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
        widget->GetTransform(t);
        widget->GetProp3D()->SetUserTransform(t);
        if (pVTKWidget) pVTKWidget->SentBoxWidgetChangedSignal();
    }
    void SetQVTKWidget(NeuroGLWidget* arg){ pVTKWidget = arg; }
private:
    NeuroGLWidget* pVTKWidget;
    vtkTransform *t;
};

#endif // NEUROGLWIDGET_H
