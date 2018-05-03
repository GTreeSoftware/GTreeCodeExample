/* SparseTracer
*	function: trace single neuron
*	Author : zhou hang
*	2015-10-28
*/
#include <vtkAutoInit.h>  
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageShrink3D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCubeSource.h>
#include <vtkPropPicker.h>
#include <vtkLineSource.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkOutlineFilter.h>
#include <vtkMath.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkArrowSource.h>
#include <vtkAppendPolyData.h>
#include <vtkOutlineSource.h>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QMenu>
#include <QKeySequence>
#include <vtkPngReader.h>
#include <vtkProperty.h>
#include <vtkPlatonicSolidSource.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyVertex.h>
#include <vtkImageMapper.h>
#include <vtkPlaneSource.h>
#include <vtkTextureMapToPlane.h>
#include <vtkActor2D.h>
#include <vtkRenderLargeImage.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <vector>
#include <QDateTime>
#include <QFileDialog>
#include <QString>
#include <QInputDialog>
#include <QApplication>
#include "neuroglwidget.h"
#include "Function/volumealgo.h"
#include "Function/Trace/traceutil.h"
#include "../ngtypes/soma.h"
#include "../ngtypes/volume.h"
#include "../ngtypes/tree.h"
#include "../Function/IO/treewriter.h"
#include "../Function/Trace/SparseTrace/SVMTraceFilter.h"
#include "../Function/NGUtility.h"
#include "../Function/Trace/CorrectTrace.h"
#include "../Function/Trace/CaliberBuilder.h"

static uchar traceColor[45] = {0,255,255,
                            255,255,0,
                            255,0,255,
                            255, 125, 0,
                            255, 0, 125,
                            125, 255, 0,
                            0, 255, 125,
                            125, 0, 255,
                            0, 125, 255,
                            255, 125, 125,
                            125, 125, 255,
                            125, 255, 125,
                            0, 0, 255,
                            0, 255, 0,
                            255, 0, 0};

static uchar NORMALCELLCOLOR[3] = { 0, 255, 0 };//green
static uchar PICKEDCELLCOLOR[3] = { 255, 0, 0 };//red
static uchar NORMALVERTEXCOLOR[3] = { 0, 0, 255 };//blue
static uchar PICKEDVERTEXCOLOR[3] = { 255, 255, 0 };//yellow
static uchar DEACTIVECOLOR[3] = { 240, 150, 200 };//pink


NeuroGLWidget::NeuroGLWidget(QWidget *parent) :QVTKOpenGLWidget(parent)//QVTKWidget(parent)
{
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
    VTK_CREATE(vtkGenericOpenGLRenderWindow, window);
    this->SetRenderWindow(&*window);

    is2DView_ = false;
    isTreeVisible_ = true;
    sampleRate_ = 0;
    sampleRateRatio_[0]=1;
    sampleRateRatio_[1]=2;
    sampleRateRatio_[2]=4;
    xScale_ = yScale_ = zScale_ = 1;
    //xOffset_ = yOffset_ = zOffset_ = 0;
    //resample_ = 1.0;
    isDrawHasStartPoint_ = false;
    //2d view
    camera2dParam_.xyz_ = -1;
    camera2dParam_.clippingDistance_ = 10;
    camera2dParam_.dFromCamera2BlockEnd_ = 1000;
    camera2dParam_.parallelScale_ = 0;
    camera2dParam_.zoomRatio_ = 0.9;
    camera2dParam_.moveStep_ = 5.0;
    camera2dParam_.rangeStep_ = 2.0;
    drawOldPoint_ = 0;
    //flag
    mode_ = NORMALMODE;
    bHasStartPoint_ = false;
    //
    CreateActions();
    VTK_NEW(vtkCamera, camera3d_);
    camera3d_->SetViewUp(0, 1, 0);
    VTK_NEW(vtkCamera, camera2d_);camera2d_->ParallelProjectionOn();
    actorRenderer_ = vtkSmartPointer<vtkRenderer>::New();
    actorRenderer_->SetLayer(1);
    volumeRenderer_ = vtkSmartPointer<vtkRenderer>::New();
    volumeRenderer_->SetLayer(0);
    GetRenderWindow()->SetNumberOfLayers(2);
    GetRenderWindow()->AddRenderer(actorRenderer_);
    GetRenderWindow()->AddRenderer(volumeRenderer_);
    volumeRenderer_->SetActiveCamera(camera3d_);
    //
    VTK_NEW(vtkCamera, cameraBox_);
    cameraBox_->ParallelProjectionOn();
    actorRenderer_->SetActiveCamera(volumeRenderer_->GetActiveCamera());/*must set camera*/
    actorRenderer_->GetActiveCamera()->ParallelProjectionOn();
    volumeRenderer_->GetActiveCamera()->ParallelProjectionOn();
    interactorStyle_ = vtkSmartPointer<NGInteractorStyleTrackballCamera>::New();
    interactorStyle_->SetWidget(this);
    interactorStyle_->AutoAdjustCameraClippingRangeOn();
    GetInteractor()->SetInteractorStyle(interactorStyle_);
    //Add Observer
    VTK_NEW(vtkCellPicker, cellPicker_);
    VTK_NEW(vtkWorldPointPicker, worldPicker_);
    VTK_NEW(vtkAreaPicker, areaPicker_);
    cellPicker_->SetTolerance(0.005);
    cellPicker_->SetPickFromList(1);
    areaPicker_->PickFromListOn();
    VTK_NEW(vtkFollower, choosedPointActor_);
    VTK_NEW(vtkFollower, choosedSomaActor_);
    VTK_NEW(vtkActor, initPtActor_);
    //-- Set OrientationMarker like Amira
    VTK_CREATE(vtkAxesActor, axes);
    VTK_NEW(vtkOrientationMarkerWidget, marker_);
    marker_->SetOutlineColor(1, 1, 1);
    marker_->SetOrientationMarker(axes);
    marker_->SetInteractor(GetRenderWindow()->GetInteractor());
    marker_->SetViewport(0.0, 0.0, 0.1, 0.1);
    marker_->SetEnabled(1);
    marker_->InteractiveOn();
    //for big data
    VTK_NEW(vtkActor, otherLayerActor_);
    VTK_NEW(vtkActor,sliceOutlineActor_);
    VTK_CREATE(vtkCubeSource, cubeData);
    cubeData->SetBounds(0, 1, 0, 1, 0, 1);
    VTK_CREATE(vtkPolyDataMapper, cubeMapper);
    cubeMapper->SetInputConnection(cubeData->GetOutputPort());
    VTK_NEW(vtkActor, imageBoxActor_);
    imageBoxActor_->SetMapper(cubeMapper);
    imageBoxActor_->GetProperty()->SetRepresentationToWireframe();
    imageBoxActor_->GetProperty()->SetColor(1.0,1.0,0.0);//yellow
    imageBoxActor_->GetProperty()->SetLineWidth(1.0);
    actorRenderer_->AddActor(imageBoxActor_);
    VTK_NEW(vtkBoxWidget, boxWidget_);
    boxWidget_->SetInteractor(GetInteractor());
    boxWidget_->SetPlaceFactor(1.0);
    actorRenderer_->AddActor(imageBoxActor_);
    boxWidget_->SetProp3D(imageBoxActor_);
    boxWidget_->PlaceWidget();
    boxWidget_->SetCurrentRenderer(actorRenderer_);
    VTK_CREATE(vtkWidgetsCallback, callback);
    callback->SetQVTKWidget(this);
    boxWidget_->AddObserver(vtkCommand::InteractionEvent, callback);
    boxWidget_->SetRotationEnabled(false);
    boxWidget_->Off();
    volumeRenderer_->ResetCamera();
   
    pickedVertexID_.first = -1;
    pickedVertexID_.second = -1;

    VTK_NEW(vtkActor, vertsActor_);
    VTK_NEW(vtkPointPicker, pointPicker_);
    pointPicker_->PickFromListOn();
    semiTracer_ = SemiAutoTracer::New();
}

NeuroGLWidget::~NeuroGLWidget()
{
    if (drawOldPoint_) {
        delete[] drawOldPoint_; drawOldPoint_ = 0;
    }
}

void NeuroGLWidget::CreateActions()
{
    pInteract = new QAction("Interact", this);
    pInteract->setCheckable(true);
    pInteract->setChecked(false);
    connect(pInteract, SIGNAL(toggled(bool)), this, SLOT(SetInteractEnable_Slot(bool)));

    pXYReset = new QAction("XY", this);
    pXZReset = new QAction("XZ", this);
    pYZReset = new QAction("YZ", this);
    connect(pXYReset, SIGNAL(triggered()), this, SLOT(XYReset_Slot()));
    connect(pXZReset, SIGNAL(triggered()), this, SLOT(XZReset_Slot()));
    connect(pYZReset, SIGNAL(triggered()), this, SLOT(YZReset_Slot()));

    pChooseBranchMode = new QAction(tr("Choose_Branch   C"), this); pChooseBranchMode->setCheckable(true);
    connect(pChooseBranchMode, SIGNAL(toggled(bool)), this, SLOT(ToggleLineMode(bool)));
    pDeleteBranch = new QAction(tr("Delete_Branch   D"), this);
    connect(pDeleteBranch, SIGNAL(triggered()), this, SLOT(DeleteBranch()));
    pChooseNodeMode = new QAction(tr("Choose_Node   Z"), this); pChooseNodeMode->setCheckable(true);
    connect(pChooseNodeMode, SIGNAL(toggled(bool)), this, SLOT(TogglePointMode(bool)));
    pCutLine = new QAction("Cut_Line   X", this);
    connect(pCutLine, SIGNAL(triggered()), this, SLOT(DivideLine()));
    pDrawLineMode = new QAction(tr("Draw_Line   V"), this); pDrawLineMode->setCheckable(true);
    connect(pDrawLineMode, SIGNAL(toggled(bool)), this, SLOT(Toggle2DDraw(bool)));
    pChooseSomaMode = new QAction("Choose_Soma", this); pChooseSomaMode->setCheckable(true);
    connect(pChooseSomaMode, SIGNAL(toggled(bool)), this, SLOT(ToggleSomaMode(bool)));
    pDeleteSoma = new QAction("Delete_Soma", this);
    connect(pDeleteSoma, SIGNAL(triggered()), this, SLOT(DeleteSoma()));
    pChooseTreeMode = new QAction("Choose_Tree", this); pChooseTreeMode->setCheckable(true);
    connect(pChooseTreeMode, SIGNAL(toggled(bool)), this, SLOT(ToggleTreeMode(bool)));
    pDeleteTrees = new QAction("Delete_Trees", this);
    connect(pDeleteTrees, SIGNAL(triggered()), this, SLOT(DeleteTrees()));
    pUndoDraw = new QAction(tr("Undo_Draw   U"), this);
    connect(pUndoDraw, SIGNAL(triggered()), this, SLOT(UndoDraw()));
    pMergeTree = new QAction("Merge_2_Trees", this);
    connect(pMergeTree, SIGNAL(triggered()), this, SLOT(MergeTwoPickedTree()));
    pSetHeadAsTreeRoot = new QAction("Set_Head_As_Root", this);
    connect(pSetHeadAsTreeRoot, SIGNAL(triggered()), this, SLOT(SetHeadAsTreeRoot()));
    pSetTailAsTreeRoot = new QAction("Set_Tail_As_Root", this);
    connect(pSetTailAsTreeRoot, SIGNAL(triggered()), this, SLOT(SetTailAsTreeRoot()));
    pSaveTree = new QAction("Save_Tree", this);
    connect(pSaveTree, SIGNAL(triggered()), this, SLOT(SaveSelectedTrees()));
    pSetActiveTree = new QAction("Toggle_Active", this);
    connect(pSetActiveTree, SIGNAL(triggered()), this, SLOT(SetActiveTree_Slot()));
    pUnDel = new QAction("Undo", this);
    connect(pUnDel, SIGNAL(triggered()), this, SLOT(UnDeleteOneActiveTreeBranch()));
    pShowRoot = new QAction("Show_Root", this); pShowRoot->setCheckable(true);
    connect(pShowRoot, SIGNAL(toggled(bool)), this, SLOT(ToggleShowTreeRoot(bool)));
    pShowDirection = new QAction("Show_Direction", this); pShowDirection->setCheckable(true);
    connect(pShowDirection, SIGNAL(toggled(bool)), this, SLOT(ToggleShowCurrentLinesDirection(bool)));
    p2DView = new QAction("2DView", this); p2DView->setCheckable(true);
    connect(p2DView, SIGNAL(toggled(bool)), this, SLOT(Toggle2D3DView(bool)));
    pShowPosition = new QAction("Show_Position", this); 
    connect(pShowPosition, SIGNAL(triggered()), this, SLOT(ShowPosition_Slot()));
    pSeperateBranch = new QAction("Sepetate_Branch", this);
    connect(pSeperateBranch, SIGNAL(triggered()), this, SLOT(SeparateBranch()));
    pConjChild = new QAction("Conj_Child_Curve", this);
    connect(pConjChild, SIGNAL(triggered()), this, SLOT(ConjChildCurve()));
    pToggleEnableTrace = new QAction("Toggle_Trace   E", this);
    connect(pToggleEnableTrace, SIGNAL(triggered()), this, SLOT(ToggleEnableCurveTrace()));
    pToggleCurInitPtShow = new QAction("Toggle_Show_Init", this);
    pToggleCurInitPtShow->setCheckable(true);
    pToggleCurInitPtShow->setChecked(true);
    connect(pToggleCurInitPtShow, SIGNAL(toggled(bool)), this, SLOT(ToggleShowInit(bool)));
    pUndelTrees = new QAction("Undel_Trees", this);
    connect(pUndelTrees, SIGNAL(triggered()), this, SLOT(UndeleteTrees()));
    pBoundCheck = new QAction("BoundCheck  B", this); pBoundCheck->setCheckable(true); pBoundCheck->setChecked(true);
    connect(pBoundCheck, SIGNAL(toggled(bool)), this, SLOT(ToggleBoundCheck(bool)));
    pForceBoundCheckAll = new QAction("Force_Bound_Check", this);
    connect(pForceBoundCheckAll, SIGNAL(triggered()), this, SLOT(ForceBoundCheckAll()));
    pSmartSemiTrace = new QAction("Smart_Semi_Trace", this); pSmartSemiTrace->setCheckable(true); pSmartSemiTrace->setChecked(true);
    connect(pSmartSemiTrace, SIGNAL(toggled(bool)), this, SLOT(ToggleSmartSemiTracer(bool)));

    pAddFlag = new QAction("Annotate", this); 
    connect(pAddFlag, SIGNAL(triggered()), this, SLOT(Annotation_Slot()));
    pDelFlag = new QAction("Undo_Annotate", this);
    connect(pDelFlag, SIGNAL(triggered()), this, SLOT(Unannotation_Slot()));
    pShowLevelBranchId = new QAction("Level_Branch_ID", this);
    connect(pShowLevelBranchId, SIGNAL(triggered()), this, SLOT(ShowLevelBranchID()));

    pStartTraverseHere_ = new QAction("Start_Here", this);
    connect(pStartTraverseHere_, SIGNAL(triggered()), this, SLOT(StartTraverseFromHere_Slot()));
    /*khtao*/
    pSmoothCurve = new QAction("Smooth_Curve", this);
    connect(pSmoothCurve, SIGNAL(triggered()), this, SLOT(SmoothCurrentCurve_Slot()));
    pSaveImageAndSwcForTest = new QAction("SaveImageAndSwcForTest", this);
    connect(pSaveImageAndSwcForTest, SIGNAL(triggered()), this, SLOT(SaveImageAndSwcForTest_Slot()));

    pInsertPt_ = new QAction("Insert_Point", this);
    connect(pInsertPt_, SIGNAL(triggered()), this, SLOT(InsertPointBefore_Slot()));
    pDelPt_ = new QAction("Remove_Point", this);
    connect(pDelPt_, SIGNAL(triggered()), this, SLOT(RemovePoint_Slot()));
    pStartMovePt_ = new QAction("Move_Point M", this); pStartMovePt_->setCheckable(true); pStartMovePt_->setChecked(false);
    connect(pStartMovePt_, SIGNAL(toggled(bool)), this, SLOT(ToggleMovePoint_Slot(bool)));

    pReconnectParent_ = new QAction("Recon_Parent", this); pReconnectParent_->setCheckable(true); pReconnectParent_->setChecked(false);
    connect(pReconnectParent_, SIGNAL(toggled(bool)), this, SLOT(ToggleReconParent_Slot(bool)));

    pRefreshTree_ = new QAction("Refresh_Tree", this);
    connect(pRefreshTree_, SIGNAL(triggered()), this, SLOT(RefreshTree_Slot()));

    pCaliberBuild_ = new QAction("Build_Caliber", this);
    connect(pCaliberBuild_, SIGNAL(triggered()), this, SLOT(BuildCaliber_Slot()));

	pAllBuilder_ = new QAction("Build_AllCaliber", this);
	connect(pAllBuilder_, SIGNAL(triggered()), this, SLOT(AllCaliberBulider()));//20180328

    pSetDendrite_ = new QAction("Set_Dendrite", this);
    connect(pSetDendrite_, SIGNAL(triggered()), this, SLOT(SetBranchAsDenrite_Slot()));
    pSetAxon_ = new QAction("Set_Axon", this);
    connect(pSetAxon_, SIGNAL(triggered()), this, SLOT(SetBranchAsAxon_Slot()));

    pHideDendrite_ = new QAction("Hide_Dendrite", this); pHideDendrite_->setCheckable(true); pHideDendrite_->setChecked(false);
    connect(pHideDendrite_, SIGNAL(toggled(bool)), this, SLOT(HideDendrite_Slot(bool)));
    pHideAxon_ = new QAction("Hide_Axon", this); pHideAxon_->setCheckable(true); pHideAxon_->setChecked(false);
    connect(pHideAxon_, SIGNAL(toggled(bool)), this, SLOT(HideAxon_Slot(bool)));
    //pHideCaliber_
    pHideCaliber_ = new QAction("Hide_Caliber", this); pHideCaliber_->setCheckable(true); pHideCaliber_->setChecked(false);
    connect(pHideCaliber_, SIGNAL(toggled(bool)), this, SLOT(HideCaliber_Slot(bool)));
}

void NeuroGLWidget::contextMenuEvent( QContextMenuEvent* )
{
    if (mode_ == FORBIDDEN) {
        printf("neuroGPS-Tree2 is running...\n");
        return;
    }
	QApplication::restoreOverrideCursor();
    QMenu *menu = new QMenu(this);
    pInteract->setChecked(boxWidget_->GetEnabled());
    menu->addAction(pInteract);
    menu->addAction(pXYReset);
    menu->addAction(pXZReset);
    menu->addAction(pYZReset);
    menu->addAction(p2DView);
    menu->addSeparator();
    switch (mode_)
    {
    case NeuroGLWidget::NORMALMODE:
        menu->addAction(pChooseBranchMode);
        menu->addAction(pChooseNodeMode);
        menu->addAction(pDrawLineMode);
        menu->addAction(pChooseTreeMode);
        menu->addAction(pChooseSomaMode);
        break;
    case NeuroGLWidget::BRANCHMODE:
        menu->addAction(pChooseBranchMode);
        menu->addAction(pToggleEnableTrace);
        menu->addAction(pDeleteBranch);
        menu->addAction(pSeperateBranch);
        menu->addAction(pUnDel);
        menu->addAction(pShowDirection);
        menu->addAction(pConjChild);
        menu->addAction(pSetHeadAsTreeRoot);
        menu->addAction(pSetTailAsTreeRoot);
        menu->addAction(pToggleCurInitPtShow);
        menu->addAction(pForceBoundCheckAll);
        menu->addAction(pShowLevelBranchId);
        menu->addAction(pSetDendrite_);
        menu->addAction(pSetAxon_);
        menu->addAction(pSmoothCurve);
        menu->addAction(pCaliberBuild_);
		menu->addAction(pAllBuilder_);
        menu->addAction(pHideCaliber_);
        //menu->addAction(pSaveImageAndSwcForTest);
        break;
    case NeuroGLWidget::POINTMODE:
        menu->addAction(pChooseNodeMode);
        menu->addAction(pShowPosition);
        menu->addAction(pCutLine);
        menu->addAction(pAddFlag);
        menu->addAction(pDelFlag);
        menu->addAction(pInsertPt_);
        menu->addAction(pDelPt_);
        menu->addAction(pStartMovePt_);
        menu->addAction(pReconnectParent_);
        break;
    case NeuroGLWidget::DRAWLINEMODE:
        menu->addAction(pDrawLineMode);
        menu->addAction(pUndoDraw);
        menu->addAction(pSmartSemiTrace);
        break;
    case NeuroGLWidget::SOMAMODE:
        menu->addAction(pChooseSomaMode);
        menu->addAction(pDeleteSoma);
        break;
    case NeuroGLWidget::TREEMODE:
        menu->addAction(pChooseTreeMode);
        menu->addAction(pShowRoot);
        menu->addAction(pHideDendrite_);
        menu->addAction(pHideAxon_);
        //menu->addAction(pActiveTreeLayerDisplay);
        menu->addAction(pSetActiveTree);
        menu->addAction(pDeleteTrees);
        menu->addAction(pMergeTree);
        menu->addAction(pRefreshTree_);
        menu->addAction(pUndelTrees);
        menu->addAction(pSaveTree);
        break;
    case NeuroGLWidget::TRAVERSEMODE:
        menu->addAction(pStartTraverseHere_);
        menu->addAction(pShowPosition);
        break;
    default:
        break;
    }
    menu->addSeparator();
    menu->addAction(pBoundCheck);
	menu->exec(QCursor::pos());
    delete menu;
}

void NeuroGLWidget::ClearActorCollection(vtkSmartPointer<vtkActorCollection> arg)
{
    arg->InitTraversal();
    vtkSmartPointer<vtkActor> iterActor = arg->GetNextActor();
    while(iterActor) {
        volumeRenderer_->RemoveActor(iterActor);
        actorRenderer_->RemoveActor(iterActor);
        iterActor = arg->GetNextActor();
    }
    arg->RemoveAllItems();
    arg->Modified();
}

void NeuroGLWidget::ClearActorCollection(std::vector<vtkSmartPointer<vtkActor> > &arg)
{
    for (auto &it : arg) {
        volumeRenderer_->RemoveActor(it);
        actorRenderer_->RemoveActor(it);
    }
    arg.clear();
}

void NeuroGLWidget::ClearAllActorRenderActor()
{
    for (auto &it : neuronActorList_) {
        volumeRenderer_->RemoveActor(it.first);
        actorRenderer_->RemoveActor(it.first);
    }
    pickedLineList_.clear();
    neuronActorList_.clear();
    ClearActorCollection(currentBoundaryActorCollection_);
    ClearActorCollection(historyBoundaryActorCollection_);
    ClearActorCollection(somaActorCollection_);
    actorRenderer_->RemoveActor(otherLayerActor_);
    if (initPtActor_) actorRenderer_->RemoveActor(initPtActor_);
    if (somaActor_) actorRenderer_->RemoveActor(somaActor_);
}

void NeuroGLWidget::CreatePointChooseActor(double *p)
{
    VTK_CREATE(vtkRegularPolygonSource, polygonSource);
    polygonSource->SetNumberOfSides(4);
    polygonSource->SetRadius(0.2);
    polygonSource->SetCenter(0, 0, 0);
    VTK_CREATE(vtkPolyDataMapper, mapper);
    mapper->SetInputConnection(polygonSource->GetOutputPort());
    VTK_NEW(vtkFollower, choosedPointActor_);
    choosedPointActor_->GetProperty()->SetLineWidth(3.0);
    choosedPointActor_->GetProperty()->SetRepresentationToWireframe();
    choosedPointActor_->SetMapper(mapper);
    choosedPointActor_->SetPosition(p[0] * paramPack->xRes_, p[1] * paramPack->yRes_, p[2] * paramPack->zRes_);
    choosedPointActor_->GetProperty()->SetColor(0, 127, 255);
    choosedPointActor_->SetCamera(volumeRenderer_->GetActiveCamera());
    choosedPointActor_->Modified();
}

int NeuroGLWidget::GetSampleRate()
{
    sampleRate_ = 0;
    int ratio = sampleRateRatio_[sampleRate_];//sample ratio
    while( (volumeX_ * volumeY_ * volumeZ_ / (ratio * ratio)) > 1073741824 //|| volumeX_ > 1024
        ){//|| volumeY_ > 1024 || volumeZ_ > 1024
            if(sampleRate_ < 2) {//sampleRate increase and sample ratio decrease
                ratio = sampleRateRatio_[++sampleRate_];
            }
            else return sampleRate_;//image texture is too big
    }
    return sampleRate_;
}

vtkSmartPointer<vtkDoubleArray> NeuroGLWidget::GetXYZScale()
{
    vtkSmartPointer<vtkDoubleArray> scale_array = vtkSmartPointer<vtkDoubleArray>::New();
    scale_array->SetNumberOfComponents(3);
    scale_array->InsertNextTuple3(paramPack->xRes_, paramPack->yRes_, paramPack->zRes_ );
    return scale_array;
}

void NeuroGLWidget::SetInputVolume(IDataPointer pImg, DATATYPE type)
{
    if(!pImg){
        QMessageBox::warning(this, "NULL Image","this image is null");
        return;
    }
    volumeData = pImg;
    camera2dParam_.xyz_ = -1;
    switch(type){
    case DATATYPE::IMAGE1:
            {
                NG_CREATE_DYNAMIC_CAST(CVolume, tmpImg, pImg);
                ReadMemoryImage(tmpImg);
            }
            break;
    case DATATYPE::IMAGE8:
    case DATATYPE::IMAGE16:
            {
                NG_CREATE_DYNAMIC_CAST(SVolume, tmpImg, pImg);
                ReadMemoryImage(tmpImg, type);
            }
            break;

    case DATATYPE::IMAGE32:
            {
                NG_CREATE_DYNAMIC_CAST(IVolume, tmpImg, pImg);
                ReadMemoryImage(tmpImg);
            }
            break;
        default:
            return;
    }
    actorRenderer_->SetBackground(0,0,0);
    volumeRenderer_->SetBackground(0,0,0);
    NGResetCamera(false);
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::SetInputNeuronPopulation(IDataPointer arg)
{
    m_Source = arg;
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
    // check whether arg is soma file
    if (!tmpTree){
        printf( "empty dendrites.\n");
        return;
    }
    // clear actor
    ClearAllActorRenderActor();
    ResetCurForbiddenBoundInfo();
    // Create points and cells for the spiral
    std::vector<NeuronTreePointer> &m_Tree = tmpTree->m_pop;
    if (m_Tree.empty()) {
        return;
    }
    //create dendrites
    //
    //auto &curTopo = m_Tree[0]->m_Topo;
    for (size_t i = 0; i < m_Tree.size(); ++i){//directly create dendrite. Separate Tree by topo.
        vtkSmartPointer<vtkActor> treeActor;
        if (i == 0lu) treeActor = BuildTreeActor(*(m_Tree[i]), true);
        else treeActor = BuildTreeActor(*(m_Tree[i]));
        treeActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
        actorRenderer_->AddActor(treeActor);
        neuronActorList_.push_back(std::make_pair(treeActor, m_Tree[i]));
    }
    //set active tree
    paramPack->activeTree = m_Tree[0];//force to as active
    activeTreeActor_ = neuronActorList_[0].first.GetPointer();
    BuildHistoryBoundaryActorCollection();
    BuildCurrentBoundaryCollection();
    BuildCurInitTracePt();
    UpdateSomaList();
    NGResetCamera(false);
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::RemoveAllActor()
{
    ClearAllActorRenderActor();
    if(oldVolume_) volumeRenderer_->RemoveVolume(oldVolume_);
    if(outlineActor_)  volumeRenderer_->RemoveActor(outlineActor_);
}

void NeuroGLWidget::ReadMemoryImage(std::shared_ptr<CVolume> volume)
{
    // turn to 1d data array;
    volumeX_ = volume->x();
    volumeY_ = volume->y();
    volumeZ_ = volume->z();
    int newX = volumeX_;
    int newY = volumeY_;
    int newZ = volumeZ_;

    //turn to vtk data
    double den = std::pow(2.0, paramPack->mostdLevel_ - 1);
    VTK_CREATE(vtkImageImport, importer);
    importer->SetWholeExtent(0,newX - 1,0,newY-1,0,newZ-1);
    importer->SetDataExtentToWholeExtent();
    importer->SetDataScalarTypeToUnsignedChar();
    importer->SetNumberOfScalarComponents(1);
    importer->CopyImportVoidPointer(volume->GetPointer(), newX * newY * newZ);
    importer->SetDataSpacing(volume->XResolution() * den, volume->YResolution() * den, volume->ZResolution() * den);
    //importer->SetDataSpacing(volume->XResolution() * ratio_,volume->YResolution() * ratio_, volume->ZResolution());
    SetVolumeRenderData(importer);
}

void NeuroGLWidget::ReadMemoryImage(std::shared_ptr<SVolume> volume, DATATYPE arg)
{
    // turn to 1d data array;
    volumeX_ = volume->x();
    volumeY_ = volume->y();
    volumeZ_ = volume->z();
    vtkIdType newX = vtkIdType(volumeX_);
    vtkIdType newY = vtkIdType(volumeY_);
    vtkIdType newZ = vtkIdType(volumeZ_);

    //turn to vtk data
    //double den = std::pow(2.0, paramPack->mostdLevel_ - 1);
    vtkSmartPointer<vtkImageImport> importer = vtkSmartPointer<vtkImageImport>::New();
    //importer->SetWholeExtent(0,newX - 1,0,newY-1,0,newZ-1);
    //importer->SetDataExtentToWholeExtent();
    importer->SetNumberOfScalarComponents(1);
    int scale = 1;
    switch(arg){
    case DATATYPE::IMAGE8:
            {
                if (newX > 2048) scale = int(std::ceil(double(newX) / 2048.0));
                if (newY > 2048) scale = int(std::ceil(double(newX) / 2048.0)) > scale ? int(std::ceil(double(newX) / 2048.0)) : scale;
                if (scale >1) 
                    printf("image size is larger than 2048, will be scale to 1/%d times\n", scale);
                uchar *tmp = new uchar[newX/scale * newY/scale * newZ];
                int nx = newX / scale; int ny = newY / scale;
                vtkIdType xy = newX * newY/(scale *scale);
                for (int i = 0; i < nx; ++i) {
                    for (int j = 0; j < ny; ++j) {
                        for (int ij = 0; ij < newZ; ++ij) {
                            tmp[i+j * nx + ij * xy] = volume->operator()(i*scale, j*scale, ij) / 4;
                        }
                    }
                }
                importer->SetWholeExtent(0, nx - 1, 0, ny - 1, 0, newZ - 1);
                importer->SetDataExtentToWholeExtent();
                importer->SetDataScalarTypeToUnsignedChar();
                importer->CopyImportVoidPointer(tmp, nx * ny * newZ);
                //importer->CopyImportVoidPointer(volume->GetPointer(), newX * newY * newZ*2);
                importer->SetDataSpacing(volume->XResolution()*scale, volume->YResolution()*scale, volume->ZResolution());
                importer->Update();
                SetVolumeRenderData(importer);
                delete[] tmp;
            }
            break;
    case DATATYPE::IMAGE16:
            {
                if (newX > 2048) scale = int(std::ceil(double(newX) / 2048.0));
                if (newY > 2048) scale = int(std::ceil(double(newX) / 2048.0)) > scale ? int(std::ceil(double(newX) / 2048.0)) : scale;
                if (scale > 1){
                    printf("image size is larger than 2048, will be scale to 1/%d times\n", scale);
                    unsigned short *tmp = new unsigned short[newX / scale * newY / scale * newZ];
                    int nx = newX / scale; int ny = newY / scale;
                    vtkIdType xy = newX * newY / (scale *scale);
                    for (int i = 0; i < nx; ++i) {
                        for (int j = 0; j < ny; ++j) {
                            for (int ij = 0; ij < newZ; ++ij) {
                                tmp[i + j * nx + ij * xy] = volume->operator()(i*scale, j*scale, ij);
                            }
                        }
                    }
                    importer->SetWholeExtent(0, nx - 1, 0, ny - 1, 0, newZ - 1);
                    importer->SetDataExtentToWholeExtent();
                    importer->SetDataScalarTypeToUnsignedChar();
                    importer->CopyImportVoidPointer(tmp, nx * ny * newZ * 2);
                    importer->SetDataSpacing(volume->XResolution()*scale, volume->YResolution()*scale, volume->ZResolution());
                    importer->Update();
                    SetVolumeRenderData(importer);
                    delete[] tmp;
                }
                else{
                    importer->SetWholeExtent(0, newX - 1, 0, newY - 1, 0, newZ - 1);
                    importer->SetDataExtentToWholeExtent();
                    importer->SetDataScalarTypeToUnsignedShort();
                    importer->SetImportVoidPointer(volume->GetPointer(), newX * newY * newZ * 2);
                    importer->SetDataSpacing(volume->XResolution(), volume->YResolution(), volume->ZResolution());//volume resolution has mostlevel
                    importer->Update();
                    SetVolumeRenderData(importer);
                }
            }
            break;
        default:
            break;
    }
}

void NeuroGLWidget::ReadMemoryImage(std::shared_ptr<IVolume> volume)
{
    //TODO:
}

void NeuroGLWidget::UpdateVolumeOpac_Slot()
{
    //2015-8-13
    if(oldVolume_){
        vtkVolumeProperty *oldProperty;
        oldProperty = oldVolume_->GetProperty();
        VTK_CREATE(vtkPiecewiseFunction, newOpacFunction);
        VTK_CREATE(vtkColorTransferFunction, colorTransfer);

        newOpacFunction->AddPoint(paramPack->lowOpac_,0.99);
        newOpacFunction->AddPoint(paramPack->highOpac_,1.0);
        newOpacFunction->ClampingOn();
        colorTransfer->AddRGBPoint(paramPack->lowOpac_,0.0,0.0,0.0);
        colorTransfer->AddRGBPoint(paramPack->highOpac_,1.0,1.0,1.0);

        oldProperty->SetScalarOpacity(newOpacFunction);
        oldProperty->SetColor(colorTransfer);
    }
    //this->update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::SetVolumeRenderData(vtkSmartPointer<vtkImageImport> importer)
{
    //LOD data for low resolution view
    /* Part3 high, low resolution mapper*/
    VTK_CREATE(vtkSmartVolumeMapper, hrVolmapper);
    //hrVolmapper->AutoAdjustSampleDistancesOn();
    //VTK_CREATE(vtkFixedPointVolumeRayCastMapper, hrVolmapper);
    //VTK_CREATE(vtkGPUVolumeRayCastMapper, hrVolmapper);
    //hrVolmapper->SetAutoAdjustSampleDistances(0);
    //hrVolmapper->SetInteractiveUpdateRate(interactiveUpdateRate);
    //hrVolmapper->SetSampleDistance(0.1);
    hrVolmapper->SetBlendModeToMaximumIntensity();
    //hrVolmapper->SetRequestedRenderModeToGPU();
    /* Part4 color trans function, only gray*/
    VTK_CREATE(vtkColorTransferFunction, colorTransfer);
    colorTransfer->AddRGBPoint(paramPack->lowOpac_, 0.0, 0.0, 0.0);
    colorTransfer->AddRGBPoint(paramPack->highOpac_, 1.0, 1.0, 1.0);
    /* Part5 opaque Function*/
    VTK_CREATE(vtkPiecewiseFunction, opacFunction);
    opacFunction->AddPoint(paramPack->lowOpac_, 0.99);
    opacFunction->AddPoint(paramPack->highOpac_, 1.0);
    opacFunction->ClampingOn();//
    /* Part6 set properties*/
    VTK_CREATE(vtkVolumeProperty, volproperty);
    volproperty->SetColor(colorTransfer);
    volproperty->SetScalarOpacity(opacFunction);
    volproperty->SetScalarOpacityUnitDistance(0.2);
    //volproperty->ShadeOn();
    //volproperty->SetAmbient(0.1);
    //volproperty->SetDiffuse(0.9);
    //volproperty->SetSpecular(0.2);
    //volproperty->SetSpecularPower(10);
    volproperty->SetInterpolationType(2);
    /* Part7 set high resolution mapper and low resolution mapper*/
    hrVolmapper->SetInputData(importer->GetOutput());
    //hrVolmapper->SetInterpolationMode(2);
    //lrVolmapper->SetInput(shrinkFilter->GetOutput());
    if(!transformForVolume) transformForVolume = vtkSmartPointer<vtkTransform>::New();
    transformForVolume->Identity();
    //transformForVolume->Translate(xOffset_, yOffset_ , zOffset_ );
    transformForVolume->Translate(paramPack->xMin_ * paramPack->xRes_, paramPack->yMin_* paramPack->yRes_, paramPack->zMin_ *paramPack->zRes_);

    /* Part 8*/
    VTK_CREATE(vtkVolume, volume);
    volume->SetProperty(volproperty);
    volume->SetMapper(hrVolmapper);
    volume->SetUserTransform(transformForVolume);
    /* Part9 if at low frame per second, use high resolution.*/
    //-- set renderer
    if(oldVolume_) volumeRenderer_->RemoveVolume(oldVolume_);
    if(outlineActor_)  actorRenderer_->RemoveActor(outlineActor_);
    volumeRenderer_->AddActor(volume);
    oldVolume_ = volume;
    VTK_CREATE(vtkOutlineFilter, outlineData);
    VTK_CREATE(vtkPolyDataMapper, mapOutLine);
    VTK_NEW(vtkActor, outlineActor_);
    outlineData->SetInputConnection(importer->GetOutputPort());
    mapOutLine->SetInputConnection(outlineData->GetOutputPort());
    outlineActor_->SetMapper(mapOutLine);
    outlineActor_->GetProperty()->SetColor(1,1,0);//yellow
    outlineActor_->SetUserTransform(transformForVolume);
    actorRenderer_->AddActor(outlineActor_);
}

void NeuroGLWidget::SetScale( int x, int y, int z )
{
    xScale_ = x;
    yScale_ = y;
    zScale_ = z;
}

void NeuroGLWidget::XYReset_Slot()
{
    if (is2DView_) {
        //camera2dParam_.xyz_ = 0;
        Set2dCamera(0);
        volumeRenderer_->SetActiveCamera(camera2d_);
        actorRenderer_->SetActiveCamera(camera2d_);
        std::for_each(annotationFlagList_.begin(), annotationFlagList_.end(), [&](annotationFlagClass &arg){arg.second->SetCamera(camera2d_); });
    }
    else{
        camera3d_ = vtkSmartPointer<vtkCamera>::New();
        camera3d_->ParallelProjectionOn();
        actorRenderer_->SetActiveCamera(camera3d_);
        volumeRenderer_->SetActiveCamera(camera3d_);
        if (isTreeVisible_) {
            actorRenderer_->ResetCamera();
        }
        else {
            volumeRenderer_->ResetCamera();
        }
        //NGResetCamera(isTreeVisible_);
    }
    ResetVTKFollowerCamera();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::YZReset_Slot()
{
    if (is2DView_) {
        Set2dCamera(1);
        volumeRenderer_->SetActiveCamera(camera2d_);
        actorRenderer_->SetActiveCamera(camera2d_);
        std::for_each(annotationFlagList_.begin(), annotationFlagList_.end(), [&](annotationFlagClass &arg){arg.second->SetCamera(camera2d_); });
    }
    else{
        camera3d_ = vtkSmartPointer<vtkCamera>::New();
        camera3d_->ParallelProjectionOn();
        camera3d_->Azimuth(90);
        actorRenderer_->SetActiveCamera(camera3d_);
        volumeRenderer_->SetActiveCamera(camera3d_);
        if (isTreeVisible_) {
            actorRenderer_->ResetCamera();
        }
        else {
            volumeRenderer_->ResetCamera();
        }
        //NGResetCamera(isTreeVisible_);
    }
    ResetVTKFollowerCamera();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::XZReset_Slot()
{
    if (is2DView_) {
        Set2dCamera(2);
        volumeRenderer_->SetActiveCamera(camera2d_);
        actorRenderer_->SetActiveCamera(camera2d_);
        std::for_each(annotationFlagList_.begin(), annotationFlagList_.end(), [&](annotationFlagClass &arg){arg.second->SetCamera(camera2d_); });
    }
    else{
        camera3d_ = vtkSmartPointer<vtkCamera>::New();
        camera3d_->ParallelProjectionOn();
        camera3d_->Elevation(45);
        camera3d_->OrthogonalizeViewUp();
        camera3d_->Elevation(45);
        camera3d_->OrthogonalizeViewUp();
        actorRenderer_->SetActiveCamera(camera3d_);
        volumeRenderer_->SetActiveCamera(camera3d_);
        if (isTreeVisible_) {
            actorRenderer_->ResetCamera();
        }
        else {
            volumeRenderer_->ResetCamera();
        }
        //NGResetCamera(isTreeVisible_);
    }
    ResetVTKFollowerCamera();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::DeleteSoma()
{
    if (SOMAMODE == mode_){      
        if (!choosedSomaActorCollection_.empty()){
            std::vector<size_t> delList;
            for (size_t i = 0; i < choosedSomaActorCollection_.size(); ++i) {
                vtkSmartPointer<vtkActor> iterActor = choosedSomaActorCollection_[i];
                actorRenderer_->RemoveActor(iterActor);
                auto actor = std::find(somaActorCollection_.begin(), somaActorCollection_.end(), iterActor);
                size_t id = std::distance(somaActorCollection_.begin(), actor);
                somaActorCollection_.erase(actor);
                delList.push_back(id);
                cellPicker_->GetPickList()->RemoveItem(iterActor);
            }
            for (auto &it : delList)   manualSoma_->GetAllCell()[it].ID = -1;
            auto &somaList = manualSoma_->GetAllCell();
            somaList.erase(std::remove_if(somaList.begin(), somaList.end(), [&](const Cell& arg){return arg.ID == -1; }), somaList.end());
        }
        choosedSomaActorCollection_.clear();
        //update();
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::DeleteTrees()
{
    if (mode_ == TREEMODE) {
        if (!pickedTreeList_.empty()) {
            //r u sure to delete tree?
            if (QMessageBox::No == QMessageBox::question(this, "Delete Tree","Are you sure to delete these trees?",
                QMessageBox::Yes | QMessageBox::No, QMessageBox::No)) {
                return;
            }
            //
            NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
            std::vector<NeuronTreePointer> &treeList = tmpTree->m_pop;
            auto &delTreeList = tmpTree->deletedPop; delTreeList.clear();//clear last deleted trees
            for (auto &it : pickedTreeList_) {
                //remove actor from renderer
                actorRenderer_->RemoveActor(it.first);
                //find actor and data from neuronActorList_
                auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == it.first.GetPointer(); });
                //remove data
                auto iter = std::find_if(treeList.begin(), treeList.end(), [&](const NeuronTreePointer &arg){return &(*arg) == &(*corNeuronTree->second); });
                //(*corNeuronTree).second->Clear();
                std::move(iter, iter + 1, std::back_inserter(delTreeList));//move to deleted tree list, but the iterator still exist
                treeList.erase(iter);
                //remove active
                if (corNeuronTree->second.get() == paramPack->activeTree.get()) {
                    if (tmpTree->m_pop.size()>0lu) {
                        paramPack->activeTree = tmpTree->m_pop[0];
                        auto newActive = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.second.get() == paramPack->activeTree.get(); });
                        RebuildNeuronTreeActor(*newActive, true);
                        RemoveInitPtActor();
                        activeTreeActor_ = newActive->first.GetPointer();
                    }
                    else{
                        if (!paramPack->UserRecord.empty()){
                            paramPack->activeTree->m_ChangedBranch.clear();
                            paramPack->UserRecord.clear();
                            paramPack->activeTree->m_DelBranch.clear();
                        }
                        paramPack->activeTree.reset();//set as null
                        activeTreeActor_ = nullptr;
                    }
                }
                //remove actor from neuronActorList_
                neuronActorList_.erase(corNeuronTree);
            }
            pickedTreeList_.clear();
            printf("delete tree data.\n");
            BuildHistoryBoundaryActorCollection();
            BuildCurrentBoundaryCollection();
            emit NeedUpdateTreeWidget();
            //update();
            GetRenderWindow()->Render();
        }
    }
}

void NeuroGLWidget::TwoDViewZoom(bool zoomIn)
{
    double curParallelScale = camera2d_->GetParallelScale();
    if(zoomIn){
        curParallelScale *= camera2dParam_.zoomRatio_;
    }else{
        curParallelScale /= camera2dParam_.zoomRatio_;
    }
    if(curParallelScale < 1){
        cout<<"Can not be larger!"<<endl;
    }else if(curParallelScale > 2000){
        cout<<"Can not be smaller!"<<endl;
    }else{
        camera2d_->SetParallelScale(curParallelScale);
        camera2d_->Modified();
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::On2DMouseWheel( bool arg)
{
    double levelRes = std::pow(2.0, paramPack->mostdLevel_ - 1);
    int Beg3, End3;
    double incre = 0.0;
    double xOffset_ = paramPack->xMin_  * paramPack->xRes_;
    double yOffset_ = paramPack->yMin_  * paramPack->yRes_;
    double zOffset_ = paramPack->zMin_  * paramPack->zRes_;

    switch (camera2dParam_.xyz_)
    {
    case 0:
        Beg3 = zOffset_;
        End3 = zOffset_ + volumeZ_ * paramPack->zRes_ * levelRes - 1;
        incre = camera2dParam_.clippingDistance_ * paramPack->zRes_;
        break;
    case 1:
        Beg3 = xOffset_;
        End3 = xOffset_ + volumeX_ * paramPack->xRes_ * levelRes - 1;
        incre = camera2dParam_.clippingDistance_ * paramPack->xRes_;
        break;
    case 2:
        Beg3 = yOffset_;
        End3 = yOffset_ + volumeY_ * paramPack->yRes_ * levelRes - 1;
        incre = camera2dParam_.clippingDistance_ * paramPack->yRes_;
        break;
    default:
        printf("error in On2DMouseWheel\n");
        return;
        //break;
    }
    double cameraPostion = End3 + camera2dParam_.dFromCamera2BlockEnd_;
    double position[3];
    double focalPoint[3];
    double clippingRange[2];
    camera2d_->GetPosition(position);
    camera2d_->GetFocalPoint(focalPoint);
    camera2d_->GetClippingRange(clippingRange);
    double clipStart = clippingRange[0];
    double newClipStart(0), newClipEnd(0);
    newClipStart = clipStart + (arg == true ? 0.5 : -0.5) * incre;//camera2dParam_.clippingDistance_
    newClipEnd = newClipStart + incre;//camera2dParam_.clippingDistance_
    
    if (newClipEnd > camera2dParam_.dFromCamera2BlockEnd_ + End3 - Beg3
        || newClipStart > camera2dParam_.dFromCamera2BlockEnd_ + End3 - Beg3){
        newClipEnd = camera2dParam_.dFromCamera2BlockEnd_ + End3 - Beg3;
        newClipStart = newClipEnd - incre;//camera2dParam_.clippingDistance_
    }
    else if (newClipStart < camera2dParam_.dFromCamera2BlockEnd_ || newClipEnd < camera2dParam_.dFromCamera2BlockEnd_){
        newClipStart = camera2dParam_.dFromCamera2BlockEnd_ - 0.5;
        newClipEnd = newClipStart + incre;//camera2dParam_.clippingDistance_
    }
    if(newClipEnd - newClipStart < 1.0) return;
    double focalPoint3 = cameraPostion - (newClipStart + newClipEnd) / 2;
    camera2d_->SetClippingRange(newClipStart, newClipEnd);
    //if (NGUtility::Round(newClipEnd) - NGUtility::Round(newClipStart) < 0.5) {
        //camera2d_->SetClippingRange(newClipStart, newClipStart+0.5);
    //}else camera2d_->SetClippingRange(newClipStart, newClipEnd);
    switch (camera2dParam_.xyz_)
    {
    case 0:
        camera2d_->SetFocalPoint(focalPoint[0], focalPoint[1], focalPoint3);
        break;
    case 1:
        camera2d_->SetFocalPoint(focalPoint3, focalPoint[1], focalPoint[2]);
        break;
    case 2:
        camera2d_->SetFocalPoint(focalPoint[0], focalPoint3, focalPoint[2]);
        break;
    default:
        break;
    }
    camera2d_->Modified();
    GetRenderWindow()->Render();
    //update();
}

void NeuroGLWidget::SetEditMode( NeuroGLWidget::EDITMODE arg)
{
    if (mode_ == arg) {
        return;
    }
    if (isMovePointMode_) {
        QMessageBox::information(this, "Move Point Mode", "Please exit move point mode first.");
        return;
    }
    switch(mode_){
        case NeuroGLWidget::NORMALMODE:
            break;
        case NeuroGLWidget::BRANCHMODE:
            ExitBranchMode();
            break;
        case NeuroGLWidget::POINTMODE:
            ExitPointMode();
            break;
        case NeuroGLWidget::DRAWLINEMODE:
            ExitDrawLineMode();
            break;
        case NeuroGLWidget::SOMAMODE:
            ExitSomaMode();
            break;
        case NeuroGLWidget::TREEMODE:
            ExitTreeMode();
            break;
        default:
            break;
    }

    
    switch (arg){
        case NeuroGLWidget::NORMALMODE:
            break;
        case NeuroGLWidget::BRANCHMODE:
            EnterBranchMode();
            break;
        case NeuroGLWidget::POINTMODE:
            EnterPointMode();
            break;
        case NeuroGLWidget::DRAWLINEMODE:
            EnterDrawLineMode();
            break;
        case NeuroGLWidget::SOMAMODE:
            EnterSomaMode();
            break;
        case NeuroGLWidget::TREEMODE:
            EnterTreeMode();
            break;
        default:
            break;
    }
    mode_ = arg;
    GetRenderWindow()->Render();
    //update();
}

vtkSmartPointer<vtkActor> NeuroGLWidget::CreateArrowActor( const Vec3d& p1, const Vec3d& p2, bool isCur)
{
    //pickRayActor_
    VTK_CREATE(vtkArrowSource, arrowSource);
    VTK_CREATE(vtkTransform, tmpTransform);
    VTK_CREATE(vtkActor, pickRayActor);
    arrowSource->SetShaftRadius(0.08);
    arrowSource->SetTipRadius(0.2);
    double startPoint[3], endPoint[3];
    startPoint[0] = p1(0) * paramPack->xRes_;
    startPoint[1] = p1(1) * paramPack->yRes_;
    startPoint[2] = p1(2) * paramPack->zRes_;
    endPoint[0] = p2(0) * paramPack->xRes_;
    endPoint[1] = p2(1) * paramPack->yRes_;
    endPoint[2] = p2(2) * paramPack->zRes_;
    // Compute a basis
    double normalizedX[3];
    double normalizedY[3];
    double normalizedZ[3];
    // The X axis is a vector from start to end
    vtkMath::Subtract(endPoint, startPoint, normalizedX);
    double length = vtkMath::Norm(normalizedX);
    vtkMath::Normalize(normalizedX);
    // The Z axis is an arbitrary vector cross X
    double arbitrary[3];
    arbitrary[0] = 1.0;
    arbitrary[1] = 0.0;
    arbitrary[2] = 0.0;
    vtkMath::Cross(normalizedX, arbitrary, normalizedZ);
    vtkMath::Normalize(normalizedZ);
    // The Y axis is Z cross X
    vtkMath::Cross(normalizedZ, normalizedX, normalizedY);
    vtkSmartPointer<vtkMatrix4x4> matrix =
        vtkSmartPointer<vtkMatrix4x4>::New();
    // Create the direction cosine matrix
    matrix->Identity();
    for (unsigned int i = 0; i < 3; i++) {
        matrix->SetElement(i, 0, normalizedX[i]);
        matrix->SetElement(i, 1, normalizedY[i]);
        matrix->SetElement(i, 2, normalizedZ[i]);
    }   
    tmpTransform->Translate(startPoint);
    tmpTransform->Concatenate(matrix);
    tmpTransform->Scale(length, length, length);
    // Transform the polydata
    vtkSmartPointer<vtkTransformPolyDataFilter> transformPD = 
        vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    transformPD->SetTransform(tmpTransform);
    transformPD->SetInputConnection(arrowSource->GetOutputPort());
    //Create a mapper and actor for the arrow
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(transformPD->GetOutputPort());
    pickRayActor->SetMapper(mapper);
    if (!isCur) pickRayActor->GetProperty()->SetColor(0, 0, 255);
    else pickRayActor->GetProperty()->SetColor(255, 200, 0);
    pickRayActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    return pickRayActor;
}

void NeuroGLWidget::BuildCurrentBoundaryCollection()
{
    //destroy old actor
    ClearActorCollection(currentBoundaryActorCollection_);
    if (isBoundCheck_) {
        if (paramPack->activeTree) {
            curBoundInfo.clear();
            CollectionBoundaryPt(curBoundInfo);
            //create boundary actor
            for (size_t i = 0; i < curBoundInfo.size(); ++i) {
                vtkSmartPointer<vtkActor> tmpActor = CreateArrowActor(curBoundInfo[i].first, curBoundInfo[i].first + 8.0 *curBoundInfo[i].second, true);
                currentBoundaryActorCollection_.push_back(tmpActor);
                actorRenderer_->AddActor(tmpActor);
            }
        }
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir)
{
    int nx = int(ptCurve.size());
    Vec3d tmpFir;
    tmpFir.setZero();
    double tmp(0.0);
    Vec3d tmpVec;
    for (int i = 0; i < nx - 1; ++i){
        tmpVec = ptCurve[i + 1] - ptCurve[i];
        tmp = tmpVec.norm();
        tmpFir += tmp * tmpVec;
    }
    fir = tmpFir.normalized();
}

void NeuroGLWidget::SetInteractEnable_Slot( bool arg )
{
    boxWidget_->SetEnabled(arg);
    boxWidget_->SetCurrentRenderer(actorRenderer_);
}

void NeuroGLWidget::UpdateBoxWidget_Slot( )
{
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();
    transform->Translate(paramPack->xMin_ * paramPack->xRes_ ,
        paramPack->yMin_ * paramPack->yRes_ , paramPack->zMin_ * paramPack->zRes_ );

    transform->Scale((paramPack->xMax_ - paramPack->xMin_)*paramPack->xRes_ ,
        (paramPack->yMax_ - paramPack->yMin_)*paramPack->yRes_,
        (paramPack->zMax_ - paramPack->zMin_)*paramPack->zRes_);
    boxWidget_->SetTransform(transform);
    boxWidget_->GetProp3D()->SetUserTransform(transform);
    NGResetCamera(false);
}

void NeuroGLWidget::UpdateBoxWidgetWithoutReset_Slot(double *arg)
{
    if (!arg) return;
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform->Identity();
    transform->Translate(arg[0] * paramPack->xRes_,
        arg[2] * paramPack->yRes_, arg[4] * paramPack->zRes_);

    transform->Scale((arg[1] - arg[0])*paramPack->xRes_,
        (arg[3] - arg[2])*paramPack->yRes_,
        (arg[5] - arg[4])*paramPack->zRes_);
    boxWidget_->SetTransform(transform);
    boxWidget_->GetProp3D()->SetUserTransform(transform);
}

void NeuroGLWidget::SentBoxWidgetChangedSignal()
{
    emit BoxWidgetChanged_Signal();
}

void NeuroGLWidget::GetBoxWidgetBound( double bound[] ) const
{
    boxWidget_->GetProp3D()->GetBounds(bound);
}

void NeuroGLWidget::BuildCurInitTracePt( )
{
    const Vec3d &initPt = paramPack->activeTree->m_curInitPt;
    if (initPt.sum() < 1.0) {
        printf("there is no curInitTrace Pt.\n");
        return;
    }
    RemoveInitPtActor();
    BuildSphereActor(initPtActor_, initPt);
    initPtActor_->GetProperty()->SetColor(0,125,255);
    actorRenderer_->AddActor(initPtActor_);
}

void NeuroGLWidget::Toggle2D3DView(bool arg)
{
    is2DView_ = arg;
    if(is2DView_){
        actorRenderer_->RemoveActor(outlineActor_);
        Set2dCamera(camera2dParam_.xyz_);
        volumeRenderer_->SetActiveCamera(camera2d_);
        actorRenderer_->SetActiveCamera(camera2d_);
        oldVolume_->GetProperty()->SetInterpolationTypeToNearest();
        for (auto &it : annotationFlagList_) 
            it.second->SetCamera(camera2d_);
        emit Set2DView_Signal(arg);
    }else{
        actorRenderer_->AddActor(outlineActor_);
        volumeRenderer_->SetActiveCamera(camera3d_);
        actorRenderer_->SetActiveCamera(camera3d_);
        oldVolume_->GetProperty()->SetInterpolationTypeToLinear();
        for (auto &it : annotationFlagList_)
            it.second->SetCamera(camera3d_);
        emit Set2DView_Signal(arg);
    }
    ResetVTKFollowerCamera();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::Set2dCamera(int arg)
{
    double clippingStart, clippingEnd;
    double cameraPostion, focalPoint3;
    double Beg1, End1, Beg2, End2, Beg3, End3;
    double levelRes = std::pow(2.0, paramPack->mostdLevel_ - 1);
    if (arg != camera2dParam_.xyz_ || arg < 0) {
        camera2dParam_.xyz_ = arg < 0? 0:arg;
        VTK_NEW(vtkCamera, camera2d_);
        camera2d_->ParallelProjectionOn();
        double xOffset = paramPack->xMin_ * paramPack->xRes_;
        double yOffset = paramPack->yMin_ * paramPack->yRes_;
        double zOffset = paramPack->zMin_ * paramPack->zRes_;
        switch (camera2dParam_.xyz_){
        case 0:
            Beg1 = xOffset; End1 = xOffset + volumeX_ * paramPack->xRes_*levelRes - 1;
            Beg2 = yOffset; End2 = yOffset + volumeY_ * paramPack->yRes_*levelRes - 1;
            Beg3 = zOffset; End3 = zOffset + volumeZ_ * paramPack->zRes_*levelRes - 1;
            if (camera2d_->GetFocalPoint()[2] < Beg3 || camera2d_->GetFocalPoint()[2] > End3) {
                clippingStart = camera2dParam_.dFromCamera2BlockEnd_ - camera2dParam_.clippingDistance_ * paramPack->zRes_*levelRes / 2;
                clippingEnd = clippingStart + camera2dParam_.clippingDistance_ * paramPack->zRes_*levelRes;
            }
            else{
                camera2d_->GetClippingRange(clippingStart, clippingEnd);
            }
            cameraPostion = End3 + camera2dParam_.dFromCamera2BlockEnd_;
            focalPoint3 = cameraPostion - (clippingStart + clippingEnd) / 2;
            camera2dParam_.parallelScale_ = (End2 - Beg2) / 2 + 5;
            camera2d_->SetParallelScale(camera2dParam_.parallelScale_);
            camera2d_->SetPosition((Beg1 + End1) / 2, (Beg2 + End2) / 2, cameraPostion);
            camera2d_->SetClippingRange(clippingStart, clippingEnd);
            camera2d_->SetFocalPoint((Beg1 + End1) / 2, (Beg2 + End2) / 2, focalPoint3);
            break;
        case 1:
            Beg3 = xOffset; End3 = xOffset + volumeX_ * paramPack->xRes_ - 1;
            Beg1 = yOffset; End1 = yOffset + volumeY_ * paramPack->yRes_ - 1;
            Beg2 = zOffset; End2 = zOffset + volumeZ_ * paramPack->zRes_ - 1;
            if (camera2d_->GetFocalPoint()[0] < Beg3 || camera2d_->GetFocalPoint()[0] > End3) {
                clippingStart = camera2dParam_.dFromCamera2BlockEnd_ - camera2dParam_.clippingDistance_ * paramPack->xRes_*levelRes / 2;
                clippingEnd = clippingStart + camera2dParam_.clippingDistance_ * paramPack->xRes_*levelRes;
            }
            else{
                camera2d_->GetClippingRange(clippingStart, clippingEnd);
            }
            cameraPostion = End3 + camera2dParam_.dFromCamera2BlockEnd_;
            focalPoint3 = cameraPostion - (clippingStart + clippingEnd) / 2;
            camera2dParam_.parallelScale_ = (End2 - Beg2) / 2 + 5;
            camera2d_->Azimuth(90);//XZ
            camera2d_->OrthogonalizeViewUp();
            camera2d_->SetParallelScale(camera2dParam_.parallelScale_);
            camera2d_->SetPosition(cameraPostion, (Beg1 + End1) / 2, (Beg2 + End2) / 2);
            camera2d_->SetClippingRange(clippingStart, clippingEnd);
            camera2d_->SetFocalPoint(focalPoint3, (Beg1 + End1) / 2, (Beg2 + End2) / 2);
            break;
        case 2:
            Beg1 = xOffset; End1 = xOffset + volumeX_ * paramPack->xRes_ - 1;
            Beg3 = yOffset; End3 = yOffset + volumeY_ * paramPack->yRes_ - 1;
            Beg2 = zOffset; End2 = zOffset + volumeZ_ * paramPack->zRes_ - 1;
            if (camera2d_->GetFocalPoint()[1] < Beg3 || camera2d_->GetFocalPoint()[1] > End3) {
                clippingStart = camera2dParam_.dFromCamera2BlockEnd_ - camera2dParam_.clippingDistance_ * paramPack->yRes_*levelRes / 2;
                clippingEnd = clippingStart + camera2dParam_.clippingDistance_ * paramPack->yRes_*levelRes;
            }
            else{
                camera2d_->GetClippingRange(clippingStart, clippingEnd);
            }
            cameraPostion = End3 + camera2dParam_.dFromCamera2BlockEnd_;
            focalPoint3 = cameraPostion - (clippingStart + clippingEnd) / 2;
            camera2dParam_.parallelScale_ = (End2 - Beg2) / 2 + 5;
            camera2d_->Elevation(45);
            camera2d_->OrthogonalizeViewUp();
            camera2d_->Elevation(45);
            camera2d_->OrthogonalizeViewUp();
            camera2d_->SetParallelScale(camera2dParam_.parallelScale_);
            camera2d_->SetPosition((Beg1 + End1) / 2, cameraPostion, (Beg2 + End2) / 2);
            camera2d_->SetClippingRange(clippingStart, clippingEnd);
            camera2d_->SetFocalPoint((Beg1 + End1) / 2, focalPoint3, (Beg2 + End2) / 2);//
            break;
        default:
            break;
        }
    }
    emit Set2DProject_Signal(arg);
}

void NeuroGLWidget::EnterDrawLineMode()
{
    printf("enter draw mode.\n");
    VTK_CREATE(vtkPolyDataMapper, mapper1);
    VTK_CREATE(vtkPolyDataMapper, mapper2);
    VTK_NEW(vtkActor, drawLineActor_);
    VTK_NEW(vtkActor, drawVertActor_);
    drawLineActor_->SetMapper(mapper1);
    drawVertActor_->SetMapper(mapper2);
    drawVertActor_->GetProperty()->SetPointSize(4.0);
    drawVertActor_->GetProperty()->SetColor(0.0, 1.0, 0.0);
    drawLineActor_->GetProperty()->SetLineWidth(3.0);
    drawLineActor_->GetProperty()->SetColor(1.0, 0, 1.0);
    actorRenderer_->AddActor(drawLineActor_);
    actorRenderer_->AddActor(drawVertActor_);
    pickLineNum = 0;
}

void NeuroGLWidget::ExitDrawLineMode()
{
    printf("exit Draw mode.\n");
    //drawn lines is tree. 
    isDrawHasStartPoint_ = false;
    VTK_DEFINE(vtkMapper, mapper) = drawLineActor_->GetMapper();
    actorRenderer_->RemoveActor(drawLineActor_);
    actorRenderer_->RemoveActor(drawVertActor_);
    actorRenderer_->RemoveActor(pickLineActor_);
    if (mapper == NULL) return;
    if (paramPack->manualLabelCurve_.size() < 2lu){
        QMessageBox::information(this, "invalid drawn line", "You may draw a line of which the size is below 2, please redraw.");
        return;
    }
    //save manual tree
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, paramPack->separateTree);
    std::vector<NeuronTreePointer> &treeList = tmpTree->m_pop;
    if (paramPack->activeTree) {//connect to active tree
        auto &lineList = paramPack->activeTree->m_curveList;
        auto &topo = paramPack->activeTree->m_Topo;
        auto & it = paramPack->manualLabelCurve_;
        if (it.size() >= 2lu) {
            if ((it.front().block(0, 0, 3, 1) - lineList[drawLinePID_].back().block(0,0,3,1)).norm() < 1.0) {
                std::copy(it.begin() + 1, it.end(), std::back_inserter(lineList[drawLinePID_]));
            }
            else if ((it.front().block(0, 0, 3, 1) - lineList[drawLinePID_].front().block(0, 0, 3, 1)).norm() < 1.0) {
                auto treeIt = std::find(topo.begin() + 1, topo.end(), LineID(drawLinePID_));
                auto newIt = topo.append_child(treeIt, LineID(lineList.size()));
                topo.move_after(treeIt, newIt);
                lineList.push_back(it);
            }
            else{
                auto treeIt = std::find(topo.begin() + 1, topo.end(), LineID(drawLinePID_));
                topo.append_child(treeIt, LineID(lineList.size()));
                lineList.push_back(it);
            }
        }
        //delete connect curves initial trace arrow
        for (size_t i = 0; i < paramPack->activeTree->m_traceInitInfo.size(); ++i) {
            const auto &pos = paramPack->activeTree->m_traceInitInfo[i].first;
            if ( (it.front().block(0,0,3,1) - pos).norm() < 2.0 ) {
                paramPack->activeTree->m_traceInitInfo.erase(paramPack->activeTree->m_traceInitInfo.begin() + i);
                break;
            }
        }
        //post procedure
        RebuildActiveTreeActor();
        if (isBoundCheck_) ForceBoundCollectManualCurves();
    }
    else{//create individual tree
        auto & it = paramPack->manualLabelCurve_;
        if (it.size() >= 2lu){
            //const Line5d& curTree = it;
            treeList.push_back(std::make_shared<NeuronTree>());
            treeList.back()->m_curveList.push_back(it);
            treeList.back()->m_type = 2;
            treeList.back()->m_Topo.set_head(std::numeric_limits<size_t>::max());
            treeList.back()->m_Topo.append_child(treeList.back()->m_Topo.begin(), LineID(0lu));
            VTK_CREATE(vtkActor, treeActor); 
            treeActor = BuildTreeActor(*(treeList.back()), true);
            actorRenderer_->AddActor(treeActor);
            neuronActorList_.push_back(std::make_pair(treeActor, treeList.back()));
            delete[] drawOldPoint_; drawOldPoint_ = 0;
            paramPack->activeTree = treeList.back();
            activeTreeActor_ = treeActor; treeActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
            BuildCurrentBoundaryCollection();
            emit NeedUpdateTreeWidget();
        }
    }
    //force to bound check manual curves when enable bound check
    paramPack->manualLabelCurve_.clear();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ForceBoundCollectManualCurves()
{
    ////collect
    auto &curCurve = paramPack->manualLabelCurve_;
    if (curCurve.size() < 2lu){
        QMessageBox::warning(this, "please debug", "there is curves of which size is 1, please debug.\n");
        return;
    }
    Vec3d pt = curCurve.back().block(0, 0, 3, 1);
    if (!CheckValidBeforeBoundaryCollect(pt)) return;
    //delete initial point of  parent curve
    Vec3d pPt = paramPack->activeTree->m_curveList[drawLinePID_].back().block(0, 0, 3, 1);
    auto &globalTraceInfo = paramPack->activeTree->m_traceInitInfo;
    if ((curCurve.front().block(0, 0, 3, 1) - pPt).norm() < 1.0) {
        auto foundIt2 = std::find_if(globalTraceInfo.begin(), globalTraceInfo.end(),
            [&](const std::pair<Vec3d, Vec3d>& arg){return (arg.first - pPt).norm() < 2.0; });
        if (foundIt2 != globalTraceInfo.end()) globalTraceInfo.erase(foundIt2);
    }
    // get initial point
    int sz = int(curCurve.size());
    Vec3d tmpDir;
    VectorVec3d tmpTail;
    NGUtility::GetPartVectorVec3d(curCurve, std::max(0, sz - 5), sz - 1, tmpTail);// <= maxid
    CalcPtCurveDirection(tmpTail, tmpDir);
    globalTraceInfo.push_back(std::make_pair(pt, tmpDir));
    //create boundary actor
    vtkSmartPointer<vtkActor> tmpActor = CreateArrowActor(pt, pt + 8.0 *tmpDir, false);//global
    currentBoundaryActorCollection_.push_back(tmpActor);
    actorRenderer_->AddActor(tmpActor);
}

double NeuroGLWidget::PickXYSliceMaxPoint(double pickPoint[3] )
{
    double xOffset_ = paramPack->xMin_  * paramPack->xRes_;
    double yOffset_ = paramPack->yMin_  * paramPack->yRes_;
    double zOffset_ = paramPack->zMin_  * paramPack->zRes_;
    double levelRes = std::pow(2.0, paramPack->mostdLevel_-1);
    int x = vtkMath::Round( (pickPoint[0] - xOffset_) / paramPack->xRes_ / levelRes);
    int y = vtkMath::Round( (pickPoint[1] - yOffset_) / paramPack->yRes_ / levelRes);
    int maxVal = 0, ind = 0;
    double pos[3];
    camera2d_->GetPosition(pos);
    double tMin, tMax;   
    camera2d_->GetClippingRange(tMin, tMax);
    std::vector<int> pixList;
    std::vector<std::vector<size_t > > maxIndList;
    if (volumeData->DataType() == DATATYPE::IMAGE8){
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[2] - tMax) / paramPack->zRes_ - paramPack->zMin_) / levelRes));
        int rMax = std::min(volume->z() - 1, NGUtility::Round(((pos[2] - tMin) / paramPack->zRes_ - paramPack->zMin_ - 1) / levelRes));
        for(int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(x, y, k));
        }
        if (pixList.empty()) return double(zOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    } else{
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[2] - tMax) / paramPack->zRes_ - paramPack->zMin_) / levelRes));
        int rMax = std::min(volume->z() - 1, NGUtility::Round(((pos[2] - tMin) / paramPack->zRes_ - paramPack->zMin_ - 1) / levelRes));
        for (int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(x, y, k));
        }
        if (pixList.empty()) return double( zOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    }
    return double(ind * levelRes * paramPack->zRes_ + zOffset_);
}

double NeuroGLWidget::PickYZSliceMaxPoint(double pickPoint[3])
{
    double xOffset_ = paramPack->xMin_  * paramPack->xRes_;
    double yOffset_ = paramPack->yMin_  * paramPack->yRes_;
    double zOffset_ = paramPack->zMin_  * paramPack->zRes_;
    double levelRes = std::pow(2.0, paramPack->mostdLevel_ - 1);
    int z = vtkMath::Round((pickPoint[2] - zOffset_) / paramPack->zRes_ / levelRes);
    int y = vtkMath::Round((pickPoint[1] - yOffset_) / paramPack->yRes_ / levelRes);
    int maxVal = 0, ind = 0;
    double pos[3];
    camera2d_->GetPosition(pos);
    double tMin, tMax;
    camera2d_->GetClippingRange(tMin, tMax);
    std::vector<int> pixList;
    std::vector<std::vector<size_t > > maxIndList;
    if (volumeData->DataType() == DATATYPE::IMAGE8){
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[0] - tMax) / paramPack->xRes_ - paramPack->xMin_) / levelRes));
        int rMax = std::min(volume->x() - 1, NGUtility::Round(((pos[0] - tMin) / paramPack->xRes_ - paramPack->xMin_ - 1) / levelRes));
        for (int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(k, y, z));
        }
        if (pixList.empty()) return double(xOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    }
    else{
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[0] - tMax) / paramPack->xRes_ - paramPack->xMin_) / levelRes));
        int rMax = std::min(volume->x() - 1, NGUtility::Round(((pos[0] - tMin) / paramPack->xRes_ - paramPack->xMin_ - 1) / levelRes));
        for (int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(k, y, z));
        }
        if (pixList.empty()) return double(xOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    }
    return double(ind * levelRes * paramPack->xRes_ + xOffset_);
}

double NeuroGLWidget::PickXZSliceMaxPoint(double pickPoint[3])
{
    double xOffset_ = paramPack->xMin_  * paramPack->xRes_;
    double yOffset_ = paramPack->yMin_  * paramPack->yRes_;
    double zOffset_ = paramPack->zMin_  * paramPack->zRes_;
    double levelRes = std::pow(2.0, paramPack->mostdLevel_ - 1);
    int x = vtkMath::Round((pickPoint[0] - xOffset_) / paramPack->xRes_ / levelRes);
    int z = vtkMath::Round((pickPoint[2] - zOffset_) / paramPack->zRes_ / levelRes);
    int maxVal = 0, ind = 0;
    double pos[3];
    camera2d_->GetPosition(pos);
    double tMin, tMax;
    camera2d_->GetClippingRange(tMin, tMax);
    std::vector<int> pixList;
    std::vector<std::vector<size_t > > maxIndList;
    if (volumeData->DataType() == DATATYPE::IMAGE8){
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[1] - tMax) / paramPack->yRes_ - paramPack->yMin_) / levelRes));
        int rMax = std::min(volume->y() - 1, NGUtility::Round(((pos[1] - tMin) / paramPack->yRes_ - paramPack->yMin_ - 1) / levelRes));
        for (int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(x, k, z));
        }
        if (pixList.empty()) return double(yOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    }
    else{
        NG_CREATE_DYNAMIC_CAST(SVolume, volume, volumeData);
        int rMin = std::max(0, NGUtility::Round(((pos[1] - tMax) / paramPack->yRes_ - paramPack->yMin_) / levelRes));
        int rMax = std::min(volume->y() - 1, NGUtility::Round(((pos[1] - tMin) / paramPack->yRes_ - paramPack->yMin_ - 1) / levelRes));
        for (int k = rMin; k <= rMax; ++k){
            pixList.push_back(volume->operator ()(x, k, z));
        }
        if (pixList.empty()) return double(yOffset_);
        maxVal = *(std::max_element(pixList.begin(), pixList.end()));
        int lastValue = -1;
        for (size_t k = 0; k < pixList.size(); ++k) {
            if (pixList[k] == maxVal) {
                if (lastValue < maxVal) maxIndList.push_back(std::vector<size_t >());
                maxIndList.back().push_back(k);
            }
            lastValue = pixList[k];
        }
        auto maxLenIt = std::max_element(maxIndList.begin(), maxIndList.end(), [](std::vector<size_t >& lhs, std::vector<size_t > & rhs)->bool{ return lhs.size() < rhs.size(); });
        ind = int((*maxLenIt)[(*maxLenIt).size() / 2lu]) + rMin;
    }
    return double(ind * levelRes * paramPack->yRes_ + yOffset_);
}

void NeuroGLWidget::SetClippingDistance_Slot( int arg )
{
    camera2dParam_.clippingDistance_ = double(arg);//camera2dParam use pixel instead of micrometer
    //out of range check
    double incre;
    if (camera2dParam_.clippingDistance_ < 1.0)
        camera2dParam_.clippingDistance_ = 1.0;
    if (camera2dParam_.xyz_ == 0){
        if (camera2dParam_.clippingDistance_ > volumeZ_ ) camera2dParam_.clippingDistance_ = volumeZ_ ;
        incre = camera2dParam_.clippingDistance_ * paramPack->zRes_ / 2.0;
    }else if (camera2dParam_.xyz_ == 1){
        if (camera2dParam_.clippingDistance_ > volumeX_ ) camera2dParam_.clippingDistance_ = volumeX_ ;
        incre = camera2dParam_.clippingDistance_ * paramPack->xRes_ / 2.0;
    }
    else if (camera2dParam_.xyz_ == 2){
        if (camera2dParam_.clippingDistance_ > volumeY_ ) camera2dParam_.clippingDistance_ = volumeY_ ;
        incre = camera2dParam_.clippingDistance_ * paramPack->yRes_ / 2.0;
    }

    if(is2DView_){
        double clippingRange[2];
        camera2d_->GetClippingRange(clippingRange);
        double center = (clippingRange[0] + clippingRange[1]) / 2.0;
        clippingRange[0] = center - incre;
        clippingRange[1] = center + incre;
        camera2d_->SetClippingRange(clippingRange);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ToggleTreesVisible(bool arg)
{
    isTreeVisible_ = arg;
    initPtActor_->SetVisibility(isTreeVisible_);
    interactorStyle_->SetCurrentRenderer(isTreeVisible_ ? actorRenderer_ : volumeRenderer_);
    if (layerDisplay_.empty() || layerDisplay_[0]==0) {
        for (auto &it : neuronActorList_) {
            //it.first->SetVisibility(isTreeVisible_);
            if(!isTreeVisible_) actorRenderer_->RemoveActor(it.first);
            else actorRenderer_->AddActor(it.first);
        }
    }
    else{
        if (!isTreeVisible_) {
            actorRenderer_->RemoveActor(activeTreeActor_);
            actorRenderer_->RemoveActor(otherLayerActor_);
        }
        else {
            actorRenderer_->AddActor(activeTreeActor_);
            actorRenderer_->AddActor(otherLayerActor_);
        }
    }
    if (traverseFlagActor_) {
        traverseFlagActor_->SetVisibility(isTreeVisible_);
    }
    
    SetActorCollectionVisible(currentBoundaryActorCollection_, isTreeVisible_);
    SetActorCollectionVisible(historyBoundaryActorCollection_, isTreeVisible_);
    pToggleCurInitPtShow->setChecked(isTreeVisible_);//signal
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::SetActorCollectionVisible(std::vector<vtkSmartPointer<vtkActor> > collector, bool arg)
{
    for (auto & iter : collector)  {
        iter->SetVisibility(arg);
        if (!arg) actorRenderer_->RemoveActor(iter);
        else actorRenderer_->AddActor(iter);
    }
}

void NeuroGLWidget::SentChangeMoveCursorSignal_Slot( bool arg)
{
    emit ChangeMoveCursor_Signal(arg);
}

void NeuroGLWidget::NGResetCamera(bool arg)
{
    if (arg) {
        actorRenderer_->ResetCamera();
    } else {
        //volumeRenderer_->ResetCamera();
        if(!camera3d_) camera3d_ = vtkSmartPointer<vtkCamera>::New();
        camera3d_->ParallelProjectionOn();
        actorRenderer_->SetActiveCamera(camera3d_);
        volumeRenderer_->SetActiveCamera(camera3d_);
        actorRenderer_->ResetCamera();
        if (paramPack->OrigImage){
            double oldfocus[3];  camera3d_->GetFocalPoint(oldfocus);
            double oldpos[3];  camera3d_->GetPosition(oldpos);
            std::vector<int> curBlockRange{ paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_ };
            double focalPoint[3] = { double(curBlockRange[0] + curBlockRange[1])*paramPack->xRes_ / 2.0, 
                double(curBlockRange[2] + curBlockRange[3])*paramPack->yRes_ / 2.0, double(curBlockRange[4] + curBlockRange[5])*paramPack->zRes_ / 2.0 };
            double posOff[3]; posOff[0] = focalPoint[0] - oldfocus[0]; posOff[1] = focalPoint[1] - oldfocus[1]; posOff[2] = focalPoint[2] - oldfocus[2];
            double cameraPos[3] = { oldpos[0] + posOff[0], oldpos[1] + posOff[1], oldpos[2] + posOff[2] };
            double paraScale = double(curBlockRange[3] - curBlockRange[2]) * paramPack->xRes_ / 2.0 + 30;
            camera3d_->SetFocalPoint(focalPoint);
            camera3d_->SetPosition(cameraPos);
            camera3d_->SetParallelScale(paraScale);
        }
        //update();
        GetRenderWindow()->Render();
    }
}

vtkSmartPointer<vtkRenderer> NeuroGLWidget::GetActiveRenderer()
{
    if (isTreeVisible_) return actorRenderer_;
    else return volumeRenderer_;
}

void NeuroGLWidget::EnterSomaMode()
{
    cout << "Enter Soma Mode!" << endl;
    cellPicker_->InitializePickList();
    for (auto & iterActor : somaActorCollection_) {
        cellPicker_->AddPickList(iterActor);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::EnterTreeMode()
{
    printf("Enter Tree Mode!\n");
    //lineActorCollection save , do not need manual tree.
    if (!neuronActorList_.empty()) {
        cellPicker_->InitializePickList();
        for (auto &iterActor : neuronActorList_) {
            if (iterActor.first->GetVisibility() == 0) continue;
            cellPicker_->AddPickList(iterActor.first);
        }
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ExitSomaMode()
{
    printf("Leave Soma Mode!\n");
    if (!choosedSomaActorCollection_.empty()) {
        for (auto & iter:choosedSomaActorCollection_) {
            iter->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
            iter->GetProperty()->SetColor(255, 255, 0);
        }
        choosedSomaActorCollection_.clear();
    }
    
    cellPicker_->InitializePickList();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ExitTreeMode()
{
    printf("Exit Tree Mode!\n");
    if (!pickedTreeList_.empty()){
        for (auto &it : pickedTreeList_) {
            vtkSmartPointer<vtkActor> actor = it.first;
            actor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);//just change thick
        }
    }
    pickedTreeList_.clear();
    cellPicker_->InitializePickList();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ProjectPick(double pickPoint[3])
{
    double focalPoint[3];
    camera2d_->GetFocalPoint(focalPoint);
    int eventPosition[2];
    this->GetInteractor()->GetEventPosition(eventPosition);
    volumeRenderer_->SetDisplayPoint(eventPosition[0], eventPosition[1], 0);
    volumeRenderer_->DisplayToWorld();
    double worldPoint[4];
    volumeRenderer_->GetWorldPoint(worldPoint);
    pickPoint[0] = worldPoint[0];
    pickPoint[1] = worldPoint[1];
    pickPoint[2] = worldPoint[2];
    switch (camera2dParam_.xyz_){
    case 0:
        pickPoint[2] = PickXYSliceMaxPoint(pickPoint);
        break;
    case 1://YZ
        pickPoint[0] = PickYZSliceMaxPoint(pickPoint);
        break;
    case 2://XZ
        pickPoint[1] = PickXZSliceMaxPoint(pickPoint);
        break;
    default:
        break;
    }
}

void NeuroGLWidget::ResetVTKFollowerCamera()
{
    choosedPointActor_->SetCamera(volumeRenderer_->GetActiveCamera());
    for (auto & iterActor : somaActorCollection_) {
        vtkFollower* iterF = dynamic_cast<vtkFollower*>(iterActor.GetPointer());
        iterF->SetCamera(volumeRenderer_->GetActiveCamera());
    }
    
}

void NeuroGLWidget::EnterBranchMode()
{
    if (!neuronActorList_.empty() && activeTreeActor_) {
        cellPicker_->InitializePickList();
        cellPicker_->AddPickList(activeTreeActor_);//just need active tree
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ExitBranchMode()
{
    if (!neuronActorList_.empty()){
        for (auto &it : neuronActorList_) cellPicker_->DeletePickList(it.first);
        cellPicker_->InitializePickList();
        ClearPickedCellsArray();
        pickedLineList_.clear();
    }
}

// modified 20161205
void NeuroGLWidget::ChangePickedCellColor( vtkIdType arg, bool isSelect)
{
    auto findIt = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& lhs){ return lhs.second.get() == paramPack->activeTree.get(); });
    VTK_DEFINE(vtkPolyData, treePoly) = vtkPolyData::SafeDownCast(findIt->first->GetMapper()->GetInput());
    VTK_CREATE(vtkUnsignedCharArray, colorScalars) = vtkUnsignedCharArray::SafeDownCast(treePoly->GetCellData()->GetScalars());
    if (isSelect){
        unsigned char* old = new unsigned char[3];
        //uchar *old = new uchar[3];
        //colorScalars->GetTupleValue(arg, old);
        colorScalars->GetTypedTuple(arg, old);
        oldColor_.push_back(old);
        //colorScalars->SetTupleValue(arg, PICKEDCELLCOLOR);
        colorScalars->SetTypedTuple(arg, PICKEDCELLCOLOR);
    }
    else{
        //colorScalars->SetTupleValue(arg, NORMALCELLCOLOR);
        if (oldColor_.empty()) {
            NG_ERROR_MESSAGE("debug old color");
            colorScalars->SetTypedTuple(arg, NORMALCELLCOLOR);
            //system("pause");
        }
        uchar *restoreColor = oldColor_.front();
        colorScalars->SetTypedTuple(arg, restoreColor);
        delete[] restoreColor; restoreColor = 0;
        oldColor_.pop_front();
    }
    colorScalars->Modified();
    treePoly->GetCellData()->Update();
    treePoly->Modified();
    findIt->first->Modified();
}

void NeuroGLWidget::MouseClick()
{
    if (mode_ == NORMALMODE) return;
    memset(eventPos_, 0, 2 * sizeof(int));
    GetInteractor()->GetEventPosition(eventPos_);
    if (BRANCHMODE == mode_){//choose branch
        cellPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
        vtkSmartPointer<vtkActor> pickActor = cellPicker_->GetActor();
        pickedCellId_ = cellPicker_->GetCellId();

        if (pickActor.GetPointer() == activeTreeActor_){
            if (neuronActorList_.empty() || pickedCellId_ < 0)  return;//exception
            auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == pickActor.GetPointer(); });
            if (corNeuronTree == neuronActorList_.end()) {
                printf("no active tree.\n");
                return;
            }
            if (&(*(corNeuronTree->second)) != &(*(paramPack->activeTree))) return;//must be active tree
            AddCellId2PickedArray(pickActor, pickedCellId_);
        }
    }else if (POINTMODE == mode_ || TRAVERSEMODE == mode_){// && !canPickPoint_
        cellPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
        vtkSmartPointer<vtkActor> pickActor = cellPicker_->GetActor();
        vtkIdType pickedCellId = cellPicker_->GetCellId();//line id
        if (pickActor.GetPointer() == activeTreeActor_ && pickedCellId > -1){
            if (neuronActorList_.empty() || pickedCellId < 0){
                printf("MouseClick2 point mode error.\n"); return;
            }
            auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == pickActor.GetPointer(); });
            if (&(*(corNeuronTree->second)) != &(*(paramPack->activeTree))) return;
            if (isReconnectParent_) {
                if (pickedVertexID_.first != -1 && pickedVertexID_.first != pickedCellId){
                    //check whether there are child curve connect to front part the current curve
                    auto &curTopo = paramPack->activeTree->m_Topo;
                    auto findIt = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedVertexID_.first));
                    std::vector<tree<LineID>::iterator> childrenList;
                    if (findIt.number_of_children() != 0) {
                        //find children
                        auto childSiblingIt = curTopo.child(findIt, 0);//if null number_of_children must be 0
                        auto endIt = curTopo.end();
                        for (; childSiblingIt.node != NULL; childSiblingIt = curTopo.next_sibling(childSiblingIt))
                            childrenList.push_back(childSiblingIt);
                        //find corresponding parent line.
                        for (auto &it : childrenList) {
                            auto &allLine = paramPack->activeTree->m_curveList;
                            auto &curLine = paramPack->activeTree->m_curveList[pickedVertexID_.first];
                            size_t conId = KPTREE::FindNearestID(allLine[it->id].front().block(0, 0, 3, 1), curLine);
                            if (conId <= size_t(pickedVertexID_.second)) {
                                QMessageBox::warning(this, "recon parent warning!", "There are some curves connect to previous part of this curve,\n please process these curves first.");
                                return;
                            }
                        }
                    }

                    auto pickedVertexIDCp = pickedVertexID_;
                    ShowVertex(pickedCellId, true);
                    pointPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
                    vtkSmartPointer<vtkActor> pickVertsActor = pointPicker_->GetActor();
                    vtkIdType pickedPointId = pointPicker_->GetPointId();
                    ShowVertex(pickedCellId, false);
                    pickedVertexID_ = pickedVertexIDCp;
                    if (pickedPointId != -1) {
                        auto &allLine = paramPack->activeTree->m_curveList;
                        auto &curLine = paramPack->activeTree->m_curveList[pickedVertexID_.first];
                        auto &addPos = allLine[pickedCellId][pickedPointId];
                        auto iter = curLine.begin() + pickedVertexID_.second;
                        curLine.erase(curLine.begin(), iter);
                        curLine.insert(curLine.begin(), addPos);
                        //topology
                       /* auto &curTopo = paramPack->activeTree->m_Topo;
                        auto findIt = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedVertexID_.first));*/
                        auto newParentIt = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedCellId));
                        if (pickedPointId == 0) {
                            curTopo.move_before(newParentIt, findIt);
                        } else {
                            if (curTopo.number_of_children(newParentIt) > 0) {
                                curTopo.move_before(curTopo.child(newParentIt, 0), findIt); 
                            }
                            else{
                                curTopo.append_child(newParentIt, LineID(100000));
                                auto eraseIt = std::find(curTopo.begin() + 1, curTopo.end(), LineID(100000));
                                curTopo.move_ontop(eraseIt, findIt);
                            }
                        }
                        RebuildActiveTreeActor();
                        ShowVertex(pickedVertexIDCp.first, true);
                        //update();
						pickedVertexID_.second = 1;
                        GetRenderWindow()->Render();
                    }
                }
                return;
            }
            
            //if (isActiveTreeLayerDisplay_) {//level display cannot pick other branch
            //    auto iter = std::find(levelCurveID_.begin(), levelCurveID_.end(), pickedCellId);
            //    if (iter == levelCurveID_.end()) return;
            //}
            
            if (pickedVertexID_.first == pickedCellId  ){//cancel select cur line
            }
            else{
                if (pickedVertexID_.first != -1){// cancel old line
                    ChangePickedCellColor(pickedVertexID_.first, false);
                    ShowVertex(pickedVertexID_.first, false);
                    pointPicker_->InitializePickList();
                    pickedVertexID_.second = -1;
                }
                pickedVertexID_.first = pickedCellId;
                ChangePickedCellColor(pickedVertexID_.first, true);
                ShowVertex(pickedVertexID_.first, true);
            }
        } else return;
        pointPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
        vtkSmartPointer<vtkActor> pickVertsActor = pointPicker_->GetActor();
        vtkIdType pickedPointId = pointPicker_->GetPointId();
        if (pickVertsActor.GetPointer() != vertsActor_.GetPointer()) return;
        if (pickedPointId >= 0){
            if (pickedVertexID_.second == pickedPointId){ //&& pickedVertexID_.first == pickedCellId
                ChangePickedVertColor(pickedVertexID_.second, false);
                pickedVertexID_.second = -1;
                if (pickedVertexID_.first == pickedCellId) {
                    ChangePickedCellColor(pickedVertexID_.first, false);
                    ShowVertex(pickedVertexID_.first, false);
                    pointPicker_->InitializePickList();
                    pickedVertexID_.first = -1;
                    pickedVertexID_.second = -1;
                }
                /*ChangePickedCellColor(pickedVertexID_.first, false);
                ShowVertex(pickedVertexID_.first, false);
                pickedVertexID_.first = -1;*/
            }
            else{
                if (pickedVertexID_.second != -1)
                    ChangePickedVertColor(pickedVertexID_.second, false);
                ChangePickedVertColor( pickedPointId, true);
                pickedVertexID_.second = pickedPointId;
            }
        }
    }
   
    if (DRAWLINEMODE == mode_){
        double pickPoint[3];
        if (!is2DView_){
            if (interactorStyle_->numberOfClicks ==2){//(GetRenderWindow()->GetInteractor()->GetRepeatCount() == 1){
                if (!Draw3DPoint(pickPoint)) return;   
            }else return;
        }
        else if (is2DView_) {
            ProjectPick(pickPoint);
        }
        if (isDrawHasStartPoint_){
            if (drawOldPoint_){
                Vec5d tmp = NGUtility::MakeVec5d(pickPoint[0] / paramPack->xRes_,
                    pickPoint[1] / paramPack->yRes_,
                    pickPoint[2] / paramPack->zRes_);
                if ((Vec3d(drawOldPoint_) - tmp.block(0, 0, 3, 1)).norm() < 0.2) {
                    printf("you draw too close to last point.\n");  return;
                }
                if (!IsInCurrentImageRange(tmp.block(0, 0, 3, 1))){
                    printf("you draw near image boundary.\n");  return;
                }
                Line5d &line = paramPack->manualLabelCurve_;
                if (isSmartSemiTracer_) {
                    //semi tracer
					semiTracer_->SetInput(paramPack->OrigImage);
                    semiTracer_->SetBeginEnd(line.back().block(0, 0, 3, 1), tmp.block(0, 0, 3, 1));
                    semiTracer_->SetGlobalOffset(paramPack->xMin_, paramPack->yMin_, paramPack->zMin_);
                    semiTracer_->Update();
                    auto &tmpLine = semiTracer_->GetSemiTraceLine();
					if (!tmpLine.empty())
					{//smooth //20180409
						NGCorrectTrace mycorrector = CorrectTrace::New();
						mycorrector->SetParam(paramPack);
						mycorrector->SetLinePoints(tmpLine, 1);
						if (!mycorrector->Update()){
							QMessageBox::warning(this, "error!", "maybe the connect point is not in the current range!");
						}
						VectorVec5d n_tmpLine = mycorrector->GetSmoSemiLine();
						mycorrector->SemiFlag = 0;
						for (auto &it : n_tmpLine)  line.push_back(it);
					}
					else
					{
						for (auto &it : tmpLine)  line.push_back(it);
					}
                   
					
                } else line.push_back(tmp);
                //line.push_back(tmp);
                BuildLineActor(drawLineActor_, line, true);
                BuildLineActor(drawVertActor_, line, false);
                drawLineActor_->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
                drawOldPoint_[0] = pickPoint[0];
                drawOldPoint_[1] = pickPoint[1];
                drawOldPoint_[2] = pickPoint[2];
            }
        }
        else{
            if (!drawOldPoint_)  drawOldPoint_ = new double[3];
            drawOldPoint_[0] = pickPoint[0] / paramPack->xRes_;
            drawOldPoint_[1] = pickPoint[1] / paramPack->yRes_;
            drawOldPoint_[2] = pickPoint[2] / paramPack->zRes_;
            paramPack->manualLabelCurve_.clear();
            Line5d &line = paramPack->manualLabelCurve_;
            if (paramPack->activeTree) {//connect to active tree line
                Vec3d curPt(drawOldPoint_[0], drawOldPoint_[1], drawOldPoint_[2]);
                std::vector<Line5d> &lineList = paramPack->activeTree->m_curveList;
                drawLinePID_ = KPTREE::FastSearchNearestLine(curPt, lineList, drawVertexPID_, 3lu, 10.0);//threv = 10.0
                if (drawLinePID_ < std::numeric_limits<size_t>::max()) {
                    line.push_back(NGUtility::MakeVec5d(lineList[drawLinePID_][drawVertexPID_](0),
                        lineList[drawLinePID_][drawVertexPID_](1),
                        lineList[drawLinePID_][drawVertexPID_](2)));
                }   else{
                    printf("please draw a closer point near active tree.\n");
                    return;
                }
            }
            else{
                line.push_back(NGUtility::MakeVec5d(drawOldPoint_[0],
                    drawOldPoint_[1],
                    drawOldPoint_[2]));
            }
            BuildLineActor(drawLineActor_, line, true);
            BuildLineActor(drawVertActor_, line, false);
            drawLineActor_->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
            isDrawHasStartPoint_ = true;
        }
        //printf("manual size: %d\n",paramPack->manualLabelCurve_.size());
    }
    if (mode_ == SOMAMODE) {
        cellPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
        vtkSmartPointer<vtkActor> pickActor = cellPicker_->GetActor();
        if (pickActor){
            float lineWidth = pickActor->GetProperty()->GetLineWidth();
            if (static_cast<int>(lineWidth) == LINEWIDTH_NOR){
                pickActor->GetProperty()->SetLineWidth(LINEWIDTH_CHOOSEN);
                pickActor->GetProperty()->SetColor(0, 255, 0);//green
                choosedSomaActorCollection_.push_back(pickActor);
            }
            else if (static_cast<int>(lineWidth) == LINEWIDTH_CHOOSEN){
                pickActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
                auto it = std::find(somaActorCollection_.begin(), somaActorCollection_.end(), pickActor);
                if (it != somaActorCollection_.end()) pickActor->GetProperty()->SetColor(255, 0, 0);
                auto it2 = std::find(choosedSomaActorCollection_.begin(), choosedSomaActorCollection_.end(), pickActor);
                choosedSomaActorCollection_.erase(it2);
            }
        }
        else if (is2DView_) {
            if (!paramPack->SomaList) paramPack->SomaList = std::make_shared<Soma>();
            NG_DYNAMIC_CAST(Soma, manualSoma_, paramPack->SomaList);
            double pickPoint[3];
            ProjectPick(pickPoint);
            auto tmpSoma = Soma::MakeCell(pickPoint[0] / paramPack->xRes_,
                pickPoint[1] / paramPack->yRes_,
                pickPoint[2] / paramPack->zRes_);
            manualSoma_->push_back(tmpSoma);
            auto &cell = manualSoma_->GetAllCell().back();
            if (cell.x < paramPack->xMin_ || cell.x > paramPack->xMax_ ||
                cell.y < paramPack->yMin_ || cell.y > paramPack->yMax_ ||
                cell.z < paramPack->zMin_ || cell.z > paramPack->zMax_) {
                manualSoma_->GetAllCell().pop_back();
                return;
            }
            VTK_NEW(vtkFollower, choosedSomaActor_);
            BuildSomaActor(manualSoma_->GetAllCell().back(), choosedSomaActor_);
            //VTK_CREATE(vtkRegularPolygonSource, polygonSource);
            //polygonSource->SetCenter(0, 0, 0);
            //polygonSource->SetNumberOfSides(3);
            ////polygonSource->SetRadius();
            //VTK_CREATE(vtkPolyDataMapper, mapper);
            //mapper->SetInputConnection(polygonSource->GetOutputPort());
            //VTK_NEW(vtkFollower, choosedSomaActor_);
            //choosedSomaActor_->GetProperty()->SetRepresentationToSurface();
            //choosedSomaActor_->SetMapper(mapper);
            //choosedSomaActor_->SetPosition(pickPoint[0], pickPoint[1], pickPoint[2]);
            //choosedSomaActor_->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
            //choosedSomaActor_->GetProperty()->SetColor(255, 0, 0);
            somaActorCollection_.push_back(choosedSomaActor_);
            cellPicker_->AddPickList(choosedSomaActor_);
            //choosedSomaActor_->Modified();
            actorRenderer_->AddActor(choosedSomaActor_);
        }
    }
    if (mode_ == TREEMODE) {
        cellPicker_->Pick(eventPos_[0], eventPos_[1], 0, actorRenderer_);
        vtkSmartPointer<vtkActor> pickActor = cellPicker_->GetActor();
        if (pickActor){
            auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == pickActor.GetPointer(); });
            if (corNeuronTree == neuronActorList_.end()) {
                printf("please debug tree mode.\n");
                return;
            }
            //find whether it is inside dendrites
            float lineWidth = pickActor->GetProperty()->GetLineWidth();
            if (static_cast<int>(lineWidth) == LINEWIDTH_NOR){
                pickActor->GetProperty()->SetLineWidth(LINEWIDTH_CHOOSEN);
                pickedTreeList_.push_back(*corNeuronTree);
            }
            else if (static_cast<int>(lineWidth) == LINEWIDTH_CHOOSEN){
                pickActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
                pickedTreeList_.erase(std::remove_if(pickedTreeList_.begin(), pickedTreeList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == pickActor.GetPointer(); }), pickedTreeList_.end());
            }
        }
    }
    //update();
    GetRenderWindow()->Render();
}

bool NeuroGLWidget::IsInCurrentImageRange(const Vec3d&pt)
{
    int nxx[2] = { paramPack->xMin_, paramPack->xMax_ - 1 };
    int nyy[2] = { paramPack->yMin_, paramPack->yMax_ - 1 };
    int nzz[2] = { paramPack->zMin_, paramPack->zMax_ - 1 };
    if (NGUtility::Round(pt(0)) < nxx[0] || NGUtility::Round(pt(0)) > nxx[1] ||
        NGUtility::Round(pt(1)) < nyy[0] || NGUtility::Round(pt(1)) > nyy[1] ||
        NGUtility::Round(pt(2)) < nzz[0] || NGUtility::Round(pt(2)) > nzz[1])
        return false;
    return true;
}

bool NeuroGLWidget::AddCellId2PickedArray(vtkSmartPointer<vtkActor> pickActor, vtkIdType idArg)
{
    auto pickedObj = std::make_pair(pickActor, idArg);
    auto &topo = paramPack->activeTree->m_Topo;
    auto topIt = std::find_if(topo.begin(), topo.end(), [&](const LineID &arg){ return arg.id == idArg; });
    if (topIt == topo.end()) {
        QMessageBox::warning(this, "Please debug choose branch.", "AddCellId2PickedArray:you choose invalid branch.");
        return false;
    }
    auto subtree = paramPack->activeTree->m_Topo.subtree(topIt, topo.next_sibling(topIt));

    if (pickedLineList_.empty()) {//add
        for (auto it = subtree.begin(); it != subtree.end(); ++it) {
            pickedLineList_.push_back(it->id);
            ChangePickedCellColor(it->id, true);//set as picked color
        }
    }else{
        for(auto &it : pickedLineList_) ChangePickedCellColor(it, false);//restore normal or original color
        if (idArg == pickedLineList_[0]) {
            pickedLineList_.clear();
        }
        else{
            pickedLineList_.clear();
            for (auto it = subtree.begin(); it != subtree.end(); ++it) {
                pickedLineList_.push_back(it->id);
                ChangePickedCellColor(it->id, true);
            }
        }
    } 
    return true;
}

void NeuroGLWidget::InitImageData(int w, int h, VTK_DEFINE(vtkImageData, imArg), int wM, int hM)
{
    imArg->SetOrigin(-w / 2, -h / 2, 0);
    imArg->SetExtent(0, w - 1, 0, h - 1, 0, 0);
    imArg->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
    imArg->Modified();
    for (int i = 0; i < w; ++i){
        for (int j = 0; j < h; ++j){
            uchar *pixel = reinterpret_cast<uchar *>(imArg->GetScalarPointer(i, j, 0));
            if ((i >= 0 && i <= wM) || (i <= w - 1 && i >= w - wM) || (j >= 0 && j <= hM) || (j <= h - 1 && j >= h - hM)){
                pixel[0] = 100;
                pixel[1] = 0;
                pixel[2] = 100;
                pixel[3] = 255;
            }
            else{
                pixel[0] = 100;
                pixel[1] = 100;
                pixel[2] = 100;
                pixel[3] = 10;
            }
        }
    }
    imArg->Modified();
}

void NeuroGLWidget::ClearPickedCellsArray()
{
    size_t idN = pickedLineList_.size();
    for (size_t i = 0; i < idN; ++i)
        ChangePickedCellColor(pickedLineList_[i], false);
    pickedLineList_.clear();
}

void NeuroGLWidget::DeleteBranch()
{
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            if (pickedIdNum == paramPack->activeTree->m_curveList.size()) {
                QMessageBox::information(this, "Cannot Delete Branch", "You will delete a whole tree. Please turn to tree mode to delete whole tree.");
                return;
            }
            auto newNeuronTree = std::make_shared<NeuronTree>();//backup
            int pid = BuildSubTree(newNeuronTree);//delete branch id list
            if (pid < -1) { GetRenderWindow()->Render();  return; }//update();
            if (!paramPack->UserRecord.empty()){
                paramPack->activeTree->m_ChangedBranch.clear();
                paramPack->UserRecord.clear();
                paramPack->activeTree->m_DelBranch.clear();
            }
            paramPack->activeTree->m_DelBranch.push_back(std::make_pair(pid, newNeuronTree));
            paramPack->UserRecord.push_back(Record::DeleteBranch);
            pUnDel->setEnabled(true);
            //
            if (!layerDisplay_.empty() && layerDisplay_[0] > 0) {//remove delete branch id
                int level, branchID; NEURITETYPE type;
                if (GetCurveLevelAndBranch(level, branchID, type)){
                    auto delIt = std::remove(branchDisplay_.begin(), branchDisplay_.end(), branchID);
                    if (delIt != branchDisplay_.end()) {
                        branchDisplay_.erase(delIt);
                    }
                }
            }
            pickedLineList_.clear();
            while (!oldColor_.empty()) {
                delete[] oldColor_.back();
                oldColor_.back() = NULL;
                oldColor_.pop_back();
            }
            //
            RebuildActiveTreeActor();
            ExitBranchMode();
            EnterBranchMode();
            //update suspicious point
            if ( !paramPack->activeTree->m_suspiciousPositions_.empty()) {
                auto &delTree = paramPack->activeTree->m_DelBranch[0].second;
                delTree->BuildHash();
                Vec3i tmp;
                auto &sp = paramPack->activeTree->m_suspiciousPositions_;
                auto &sc = paramPack->activeTree->m_suspiciousCheckState_;
                for (size_t k = 0; k < sp.size(); ++k){
                    auto &pt = sp[k];
                    tmp << NGUtility::Round(pt(0)), NGUtility::Round(pt(1)), NGUtility::Round(pt(2));
                    if (delTree->m_hashList.isExistInRange(tmp, 10)) {//remove invalid suspicious point
                        pt << 0.0, 0.0, 0.0;
                        sc[k] = -1;
                    }
                }
                sp.erase(std::remove_if(sp.begin(), sp.end(), [](const Vec3d& arg){ return arg.sum() < 1.0; }), sp.end());
                sc.erase(std::remove(sc.begin(), sc.end(), -1), sc.end());
            }
        }
        emit NeedUpdateCheckWidget();
        if (!diffArrowList_.empty()) {//It has been shown.
            ToggleDiffArrows(true);
        }
        //update();
        GetRenderWindow()->Render();
        //
        int kk = activeTreeActor_->GetMapper()->GetInput()->GetNumberOfCells();
        if (kk != paramPack->activeTree->m_curveList.size() ) {
            QMessageBox::warning(this, "Caught bug!", "If you use debug mode, call zhouhang ");
            system("pause");
        }
    }
}

void NeuroGLWidget::SeparateBranch()
{
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            auto newNeuronTree = std::make_shared<NeuronTree>();
            int pid = BuildSubTree(newNeuronTree);
            if (pid < -1) { GetRenderWindow()->Render(); return; }//update(); 
            if (!paramPack->activeTree->m_DelBranch.empty()) {
                paramPack->activeTree->m_DelBranch[0].second.reset();
                paramPack->activeTree->m_DelBranch.clear();
            }
            NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
            tmpTree->m_pop.push_back(newNeuronTree);
            pickedLineList_.clear();
            while (!oldColor_.empty()){
                uchar *restoreColor = oldColor_.front();
                delete[] restoreColor; restoreColor = 0;
                oldColor_.pop_front();
            }
            RebuildActiveTreeActor();
            AppendNeuronTreeActor(newNeuronTree);
            ExitBranchMode();
            EnterBranchMode();
            emit NeedUpdateTreeWidget();
        }
    }
}

int NeuroGLWidget::BuildSubTree(NeuronTreePointer newNeuronTree)
{
    auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
    auto topoIt = std::find(corNeuronTree->second->m_Topo.begin() + 1, corNeuronTree->second->m_Topo.end(), LineID(pickedLineList_[0]));
    std::vector<Line5d> &oldLines = corNeuronTree->second->m_curveList;
        
    int pid = 0;
    if (paramPack->activeTree->m_Topo.parent(topoIt)->id > 1000000llu && paramPack->activeTree->m_Topo.parent(topoIt).number_of_children() < 2) {
        pid = -2;
    }
    else if (paramPack->activeTree->m_Topo.parent(topoIt)->id > 1000000llu){
        pid = -1;
    }
    else{
        pid = int(paramPack->activeTree->m_Topo.parent(topoIt)->id);
    }
    if (pid == -2) {
        QMessageBox::warning(this, "You Will Choose Whole Tree to delete or separate", "Please Delete whole tree using tree mode Separation is not allowed.\n");
        //update();
        GetRenderWindow()->Render();
        return pid;
    }
    newNeuronTree->m_type = 2;
    //move in tree topo, root is child, branch is the children of child
    newNeuronTree->m_Topo.set_head(LineID(std::numeric_limits<size_t>::max()));
    auto subtree = paramPack->activeTree->m_Topo.subtree(topoIt, corNeuronTree->second->m_Topo.next_sibling(topoIt));
    //collection ID
    std::vector<size_t> deletedIDList; 
    Vec3d endPt;
    for (auto it = subtree.begin(); it != subtree.end(); ++it) { // no head
        newNeuronTree->m_curveList.push_back(Line5d());
        endPt = oldLines[(*it).id].back().block(0, 0, 3, 1);
        newNeuronTree->m_curveList.back().swap(oldLines[(*it).id]);//move lines according to topo id.
        deletedIDList.push_back((*it).id);
        //collection history bound info
        auto boundIt = std::find_if(corNeuronTree->second->m_traceInitInfo.begin(), corNeuronTree->second->m_traceInitInfo.end(), [&](const std::pair<Vec3d, Vec3d>& arg){ return (arg.first - endPt).cwiseAbs().sum() < 1.0; });
        while (boundIt != corNeuronTree->second->m_traceInitInfo.end()) {
            newNeuronTree->m_traceInitInfo.push_back(std::make_pair(Vec3d(0, 0, 0), Vec3d(0, 0, 0)));
            newNeuronTree->m_traceInitInfo.back().swap(*boundIt);
            boundIt = std::find_if(corNeuronTree->second->m_traceInitInfo.begin(), corNeuronTree->second->m_traceInitInfo.end(), [&](const std::pair<Vec3d, Vec3d>& arg){ return (arg.first - endPt).cwiseAbs().sum() < 1.0; });
        }
    }
    //collect deleted line as new tree and update id
    size_t tmpID = 0lu;
    newNeuronTree->m_Topo.append_child(newNeuronTree->m_Topo.begin(), subtree.begin());
    for (auto it = newNeuronTree->m_Topo.begin() + 1; it != newNeuronTree->m_Topo.end(); ++it) {
        it->id = tmpID;//curves save according to topo deep-first search
        ++tmpID;
    }
    //erase corresponding history bound info from old tree. they have bee move to new tree
    BOUNDINFOLIST &oldInfo = corNeuronTree->second->m_traceInitInfo;
    Vec3d zeroTmp; zeroTmp.setZero();
    oldInfo.erase(std::remove_if(oldInfo.begin(), oldInfo.end(), [&](const std::pair<Vec3d, Vec3d>& arg){return (arg.first - zeroTmp).cwiseAbs().sum() < 1.0; }), oldInfo.end());
    //delete line
    auto removeIt = std::remove_if(oldLines.begin(), oldLines.end(), [](const Line5d& arg){return arg.empty(); });//no change the sort
    oldLines.erase(removeIt, oldLines.end());
    if (oldLines.empty()) {//delete empty tree
        NG_ERROR_MESSAGE("It should be tree mode, how could the whole tree deleted? Please Debug");
    }
    //delete topo it
    corNeuronTree->second->m_Topo.erase(topoIt);
    //update ID
    KPTREE::UpdateDeleteTreeIDList(*(corNeuronTree->second), deletedIDList);
    return pid;
}

void NeuroGLWidget::ConjChildCurve()
{
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            if (!KPTREE::ConjunctChildCurve(*(paramPack->activeTree), pickedLineList_[0])){
                QMessageBox::warning(this, "operation failed", "cannot conjunct child curve.");
                return;
            }
            //auto topoIt = std::find(paramPack->activeTree->m_Topo.begin() + 1, paramPack->activeTree->m_Topo.end(), LineID(pickedLineList_[0]));
            //if (paramPack->activeTree->m_Topo.number_of_children(topoIt) != 0lu) {
            //    auto childIt = paramPack->activeTree->m_Topo.child(topoIt, 0);
            //    if (childIt.node == NULL) return;
            //    std::vector<Line5d> &oldLines = paramPack->activeTree->m_curveList;
            //    Line5d& curLine = oldLines[topoIt->id];
            //    Vec3d curLineEnd = curLine.back().block(0, 0, 3, 1);
            //    for (; childIt != paramPack->activeTree->m_Topo.end() && childIt.node!=NULL; childIt = paramPack->activeTree->m_Topo.next_sibling(childIt)) {
            //        Line5d& childLine = oldLines[childIt->id];
            //        if (( childLine[0].block(0,0,3,1) - curLineEnd).norm() < 0.01 ) {
            //            size_t id = childIt->id;
            //            std::vector<size_t> deleteID; deleteID.push_back(id);
            //            //modify lines
            //            curLine.reserve(curLine.size() + childLine.size());
            //            std::copy(childLine.begin(), childLine.end(), std::back_inserter(curLine));
            //            oldLines.erase(oldLines.begin() + id);
            //            //modify topo
            //            auto subtree = paramPack->activeTree->m_Topo.subtree(childIt, paramPack->activeTree->m_Topo.next_sibling(childIt));
            //            paramPack->activeTree->m_Topo.reparent(topoIt, subtree.begin());
            //            paramPack->activeTree->m_Topo.erase(childIt);
            //            KPTREE::UpdateDeleteTreeIDList(*(paramPack->activeTree), deleteID);
            //            break;
            //        }
            //    }
                RebuildActiveTreeActor();
                ExitBranchMode();
                EnterBranchMode();
            //}
        }
    }
}

//only display 
void NeuroGLWidget::BuildHistoryBoundaryActorCollection()
{
    ClearActorCollection(historyBoundaryActorCollection_);
    if (paramPack->activeTree) {
        BOUNDINFOLIST &bound = paramPack->activeTree->m_traceInitInfo;
        for (size_t i = 0; i < bound.size(); ++i) {
            vtkSmartPointer<vtkActor> tmpActor = CreateArrowActor(std::get<0>(bound[i]), std::get<0>(bound[i]) + 8.0 * std::get<1>(bound[i]));
            historyBoundaryActorCollection_.push_back(tmpActor);
            actorRenderer_->AddActor(tmpActor);
        }
    }
}

void NeuroGLWidget::EnterPointMode()
{
    if (!neuronActorList_.empty()) {
        cellPicker_->InitializePickList();
        for (auto &it : neuronActorList_) {
            if (it.first.GetPointer() == activeTreeActor_) {
                cellPicker_->AddPickList(it.first);
                break;
            }
        }
    }
}

void NeuroGLWidget::ExitPointMode()
{//just select one line and one node
    if (neuronActorList_.empty()) return;
    
    for (auto &it : neuronActorList_) {
        if (it.first.GetPointer() == activeTreeActor_) {
            cellPicker_->DeletePickList(it.first);
            break;
        }
    }
    cellPicker_->InitializePickList();
    if ( pickedVertexID_.first != -1){
        if (pickedVertexID_.second != -1) {
            ChangePickedVertColor(pickedVertexID_.second, false);
        }
        auto tmp = pickedVertexID_.first;
        ShowVertex(tmp, false);
        pointPicker_->InitializePickList();
        ChangePickedCellColor(tmp, false);
        pickedVertexID_.first = -1;
        pickedVertexID_.second = -1;
    }
    
}

void NeuroGLWidget::ShowVertex(vtkIdType arg, bool flag)
{
    if (POINTMODE == mode_ || TRAVERSEMODE == mode_){
        if (arg < 0lu){
            return;
        }
        pointPicker_->InitializePickList();
        if (flag){
            VTK_DEFINE(vtkPolyData, treePoly) = vtkPolyData::SafeDownCast(activeTreeActor_->GetMapper()->GetInput());
            VTK_CREATE(vtkIdList, pointList);
            treePoly->GetCellPoints(arg, pointList);
            vtkIdType pointNum = pointList->GetNumberOfIds();
            if (pointNum > 0){
                VTK_CREATE(vtkPolyData, vertsPoly);
                VTK_CREATE(vtkPolyDataMapper, vertsMapper);
                VTK_CREATE(vtkCellArray, vertsCells);
                VTK_CREATE(vtkPoints, vertsPoints);
                VTK_CREATE(vtkUnsignedCharArray, vertScalars);
                vertScalars->SetNumberOfComponents(3);
                vertsCells->InsertNextCell(pointNum);
                for (vtkIdType i = 0; i < pointNum; ++i){
                    double *p = treePoly->GetPoint(pointList->GetId(i));
                    vertsPoints->InsertNextPoint(p);
                    vertsCells->InsertCellPoint(i);
                    //vertScalars->InsertNextTupleValue(NORMALVERTEXCOLOR);
                    vertScalars->InsertNextTypedTuple(NORMALVERTEXCOLOR);
                }
                vertsPoly->SetPoints(vertsPoints);
                vertsPoly->SetVerts(vertsCells);
                vertsPoly->GetPointData()->SetScalars(vertScalars);
                vertsMapper->SetInputData(vertsPoly);
                vertsMapper->SetScalarModeToUsePointData();
                vertsActor_->SetMapper(vertsMapper);
                vertsActor_->GetProperty()->SetPointSize(10.0);
                vertsActor_->Modified();
                actorRenderer_->AddActor(vertsActor_);
                pointPicker_->AddPickList(vertsActor_);
                pointPicker_->Modified();
            }
        }  else{
            actorRenderer_->RemoveViewProp(vertsActor_);
            /*pointPicker_->InitializePickList();
            pickedVertexID_.second = -1;*/
        }
    }
    //update();
    GetRenderWindow()->Render();
}

// modified 20161205
void NeuroGLWidget::DivideLine()
{
    if (POINTMODE == mode_ &&pickedVertexID_.first != -1 && pickedVertexID_.second != -1){// canPickPoint_ && 
        /*cut line into 2 lines, the second is the child of the first. search the children connect to original line, and identify which connect to which.*/
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        size_t idsNum = corNeuronTree->second->m_curveList[pickedVertexID_.first].size();
        if (idsNum <= 2llu || pickedVertexID_.second < 1 || pickedVertexID_.second >= int(idsNum) - 1){
            printf("cut line is too short, or the cut dot is the end node.\n");
            return;
        }
        tree<LineID>&topo = corNeuronTree->second->m_Topo;
        std::vector<Line5d>& lines = corNeuronTree->second->m_curveList;
        tree<LineID>::iterator curIt = std::find(topo.begin()+1, topo.end(), LineID(pickedVertexID_.first));
        std::vector<size_t> connectNodeId;
        std::vector<tree<LineID>::iterator> childrenList;
        if (curIt.number_of_children() != 0) {
            //find children
            auto childSiblingIt = topo.child(curIt, 0);//if null number_of_children must be 0
            auto endIt = topo.end();
            for (; childSiblingIt.node != NULL; childSiblingIt = topo.next_sibling(childSiblingIt))
                childrenList.push_back(childSiblingIt);
            //find corresponding parent line.
            for (auto &it : childrenList) {
                connectNodeId.push_back(KPTREE::FindNearestID(lines[it->id].front().block(0, 0, 3, 1), lines[curIt->id]));
            }
        }
        //cut line. second half is save in the end.
        lines.push_back(Line5d());
        std::copy(lines[pickedVertexID_.first].begin() + pickedVertexID_.second, lines[pickedVertexID_.first].end(), std::back_inserter(lines.back()));
        lines[pickedVertexID_.first].erase(lines[pickedVertexID_.first].begin() + pickedVertexID_.second +1, lines[pickedVertexID_.first].end());
        auto secondHalf = topo.append_child(curIt, LineID(lines.size() - 1lu));
        //reparent corresponding line
        if (!connectNodeId.empty()) {//curIt.number_of_children() != 0
            for (size_t i = 0; i < connectNodeId.size(); ++i) {
                size_t id = connectNodeId[i];
                if (vtkIdType(id) > pickedVertexID_.second) {
                    auto it = secondHalf;
                    if (topo.number_of_children(secondHalf) == 0lu) {
                         topo.move_ontop(topo.append_child(secondHalf, LineID(100000000)) , childrenList[i]);
                    }
                    else{
                        topo.move_before(topo.child(secondHalf, 0), childrenList[i]);
                    }
                }
            }
        }
        if (!paramPack->UserRecord.empty()){
            paramPack->activeTree->m_ChangedBranch.clear();
            paramPack->UserRecord.clear();
            paramPack->activeTree->m_DelBranch.clear();
        }
        //old color should be empty
        while (!oldColor_.empty()) {
            delete[] oldColor_.back();
            oldColor_.back() = NULL;
            oldColor_.pop_back();
        }
        //display
        ShowVertex(pickedVertexID_.first, false);
        pointPicker_->InitializePickList();
        pickedVertexID_.first = -1;
        pickedVertexID_.second = -1;
        pickedLineList_.clear();
        RebuildActiveTreeActor();
        ExitPointMode();
        EnterPointMode();
        //update();
        GetRenderWindow()->Render();
        //
        int kk = activeTreeActor_->GetMapper()->GetInput()->GetNumberOfCells();
        if (kk != paramPack->activeTree->m_curveList.size() ) {
            QMessageBox::warning(this, "Caught bug!", "If you use debug mode, call zhouhang ");
            system("pause");
        }
    }
}

void NeuroGLWidget::ChangePickedVertColor( vtkIdType arg, bool isChanged)
{
    if (arg < 0){
        return;
    }

    vtkPolyData *vertsPoly = vtkPolyData::SafeDownCast(vertsActor_->GetMapper()->GetInput());
    vtkUnsignedCharArray *vertsColors = vtkUnsignedCharArray::SafeDownCast(vertsPoly->GetPointData()->GetScalars());
    if (!isChanged){
        vertsColors->SetTypedTuple(arg, NORMALVERTEXCOLOR);
    }
    else{
        vertsColors->SetTypedTuple(arg, PICKEDVERTEXCOLOR);
    }
    vertsColors->Modified();
    vertsPoly->GetPointData()->Update();
    vertsPoly->Modified();
    vertsActor_->Modified();
    //update();
    GetRenderWindow()->Render();
}

bool NeuroGLWidget::CollectionBoundaryPt(BOUNDINFOLIST& arg)
{
    if (paramPack->activeTree){
        if (!(paramPack->activeTree->m_curveList.empty())) {
            size_t oldLen = arg.size();
            const std::vector<Line5d> &curCurve = paramPack->activeTree->m_curveList;
            int nxx[2] = { paramPack->xMin_, paramPack->xMax_ - 1 };
            int nyy[2] = { paramPack->yMin_, paramPack->yMax_ - 1 };
            int nzz[2] = { paramPack->zMin_, paramPack->zMax_ - 1 };
            BOUNDINFOLIST tmp;
            auto &tmpTree = paramPack->activeTree->m_Topo;
            bool flag = false; 
            int rpx, rpy, rpz;
            int innerBoundX, outerBoundX, innerBoundY, outerBoundY, innerBoundZ, outerBoundZ;
            innerBoundX = nxx[0] - paramPack->boundaryDistanceThreshold_;
            outerBoundX = nxx[1] + paramPack->boundaryDistanceThreshold_;
            innerBoundY = nyy[0] - paramPack->boundaryDistanceThreshold_;
            outerBoundY = nyy[1] + paramPack->boundaryDistanceThreshold_;
            innerBoundZ = nzz[0] - paramPack->boundaryDistanceThreshold_;
            outerBoundZ = nzz[1] + paramPack->boundaryDistanceThreshold_;
            for (size_t j = 0; j < curCurve.size(); ++j){
                if (curCurve[j].size() < 2lu){
                    QMessageBox::warning(this, "please debug", "there is curves of which size is 1, please debug.\n");
                    continue;
                }
                Vec3d pt = curCurve[j].back().block(0, 0, 3, 1);
                rpx = NGUtility::Round(pt(0));
                rpy = NGUtility::Round(pt(1));
                rpz = NGUtility::Round(pt(2));
                if (rpx < innerBoundX || rpx > outerBoundX ||
                    rpy < innerBoundY || rpy > outerBoundY ||
                    rpz < innerBoundZ || rpz > outerBoundZ)
                    continue;//must be in bounding box.

                if (AbsDiff(rpx, nxx[0]) <= paramPack->boundaryDistanceThreshold_ ||
                    AbsDiff(rpx, nxx[1]) <= paramPack->boundaryDistanceThreshold_ ||
                    AbsDiff(rpy, nyy[0]) <= paramPack->boundaryDistanceThreshold_ ||
                    AbsDiff(rpy, nyy[1]) <= paramPack->boundaryDistanceThreshold_ ||
                    AbsDiff(rpz, nzz[0]) <= paramPack->boundaryDistanceThreshold_ ||
                    AbsDiff(rpz, nzz[1]) <= paramPack->boundaryDistanceThreshold_){
                     //whether have been put into bound info
                     if (!CheckValidBeforeBoundaryCollect(pt)) continue;
                     if (!isForceBoundCheckAll_) {//check whether is it the traced curve last.
                         auto checkIt = std::find_if(paramPack->allowableBoundCheckPtList_.begin(), paramPack->allowableBoundCheckPtList_.end(),
                             [&](const Vec3d& arg){ return (arg - pt).norm() < 2.0; });
                         if (checkIt == paramPack->allowableBoundCheckPtList_.end()) continue;
                     }
                     //check if its end connect curve's front.
                     flag = false;
                     auto iter = std::find_if(tmpTree.begin() + 1, tmpTree.end(), [&](const LineID& arg){ return arg.id == j; });
                     if (tmpTree.number_of_children(iter) != 0lu) {
                         for (auto childIt = tmpTree.child(iter, 0); childIt.node != NULL; childIt = tmpTree.next_sibling(childIt)) {
                             if (childIt->id > 10000000lu) {
                                 QMessageBox::warning(this, "error in CollectionBoundaryPt", "the id error!");
                             }
                             if ((curCurve[iter->id].back().block(0, 0, 3, 1) - curCurve[childIt->id].front().block(0, 0, 3, 1)).norm() < 1.0) {
                                 flag = true;
                                 break;
                             }
                         }
                     }
                     if (flag) continue;
                     //collect
                     int sz = int(curCurve[j].size());
                     Vec3d tmpDir;
                     VectorVec3d tmpTail;
                     NGUtility::GetPartVectorVec3d(*(curCurve.begin() + j), std::max(0, sz - 5), sz - 1, tmpTail);// <= maxid
                     CalcPtCurveDirection(tmpTail, tmpDir);
                     tmp.push_back(std::make_pair(pt, tmpDir));
                 }
                 else  continue;//must be in bounding box.
            }
            if (tmp.empty()) {
                //printf("Debug\n");
            }
            if (tmp.size() <= paramPack->maxBoundNum_) {
                arg.reserve(arg.size() + tmp.size());
                std::copy(tmp.begin(), tmp.end(), std::back_inserter(arg));
            }
            else printf("the boundary point number is large.\n");
            return ((arg.size() - oldLen) != 0u);
        }
        else{
            printf("error in CollectionBoundaryPt.\n");
            return false;
        }
    }
    return false;
}

bool NeuroGLWidget::CheckValidBeforeBoundaryCollect(const Vec3d& pt)
{
    if ((pt - paramPack->activeTree->m_curInitPt).norm() < 10.0){
        return false;
    }
    //do not find the init trace point from forbidden bound
    auto foundIt1 = std::find_if(curForbiddenFoundBoundInfo.begin(), curForbiddenFoundBoundInfo.end(), [&](const Vec3d &arg){return (arg - pt).norm() < 2.0; });
    if (foundIt1 != curForbiddenFoundBoundInfo.end()) {
        return false;
    }
    //do not find the init trace point from history trace init point.
    auto foundIt2 = std::find_if(paramPack->activeTree->m_traceInitInfo.begin(), paramPack->activeTree->m_traceInitInfo.end(),
        [&](const std::pair<Vec3d, Vec3d>& arg){return (arg.first - pt).norm() < 2.0; });
    if (foundIt2 != paramPack->activeTree->m_traceInitInfo.end()) return false;
    return true;
}

void NeuroGLWidget::MergeTwoPickedTree()
{
    if (mode_ != TREEMODE ) return;
    if (pickedTreeList_.size()!=2) {
        printf("please choose two tree to merge.\n");
        return;
    }
    NeuronTreePointer arg1 = pickedTreeList_[0].second;
    NeuronTreePointer arg2 = pickedTreeList_[1].second;
    if (!MergeTwoTrees(arg1, arg2)) return;
    //rebuild 2nd tree. rebuild bound and set active if possible
    auto corNeuronTree1 = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(arg1)); });
    auto corNeuronTree2 = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(arg2)); });
    if (&(*paramPack->activeTree) == &(*(arg1)) || &(*paramPack->activeTree) == &(*(arg2))) {//update active tree
        //set active
        if (&(*paramPack->activeTree) == &(*(arg1))) {
            paramPack->activeTree = arg2;
            activeTreeActor_ = corNeuronTree2->first;
        }
        BuildHistoryBoundaryActorCollection();
        BuildCurrentBoundaryCollection();
        RebuildActiveTreeActor();
    } else{
        RebuildNeuronTreeActor(*corNeuronTree2, false);//set as active
    }
    corNeuronTree2->first->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    //clear pick list
    pickedTreeList_.clear();
    //delete tree1
    actorRenderer_->RemoveActor(corNeuronTree1->first);
    neuronActorList_.erase(corNeuronTree1);
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
    tmpTree->m_pop.erase(std::remove_if(tmpTree->m_pop.begin(), tmpTree->m_pop.end(), [](const NeuronTreePointer& arg){return arg->m_curveList.empty(); }), tmpTree->m_pop.end());//arg1->m_curveList is empty
    emit NeedUpdateTreeWidget();
    GetRenderWindow()->Render();
}

bool NeuroGLWidget::MergeTwoTrees(NeuronTreePointer arg1, NeuronTreePointer arg2)
{
    if (mode_ != TREEMODE&& mode_ != BRANCHMODE) return  false;
    /*merge the first tree as subtree of the second tree.*/
    Vec3d root1 = arg1->m_curveList[(arg1->m_Topo.begin() + 1)->id].front().block(0, 0, 3, 1);
    //find connect line and vertex
    const std::vector<Line5d> &lines2 = arg2->m_curveList;
    size_t lineID, vertID; lineID = vertID = std::numeric_limits<size_t>::max();
    double minDist = 1000000.0, tmpDist;
    for (size_t i = 0; i < lines2.size(); ++i) {
        for (size_t j = 0; j < lines2[i].size(); ++j){
            if (std::abs(root1(0) - lines2[i][j](0)) > 50.0 || std::abs(root1(1) - lines2[i][j](1)) > 50.0 || std::abs(root1(2) - lines2[i][j](2)) > 50.0) continue;
            tmpDist = (root1 - lines2[i][j].block(0, 0, 3, 1)).norm();
            if (tmpDist < minDist) {
                minDist = tmpDist;
                lineID = i;
                vertID = j;
            }
        }
    }
    if (minDist >20.0) {//can not merge two trees far away
        QMessageBox::warning(this, "Cannot merge trees!", "The distance of two trees is so far, please draw lines to make them closer.");
        return false;
    }
    //if far, add node to connect lines
    if (minDist > 0.01){
        auto it = arg1->m_Topo.begin()+1;
        for (; it != arg1->m_Topo.end() && it.node != NULL ; it = arg1->m_Topo.next_sibling(it))
            arg1->m_curveList[it->id].insert(arg1->m_curveList[it->id].begin(), arg2->m_curveList[lineID][vertID]);
    }
    //find the iterator of 2nd tree
    tree<LineID>::iterator pIter2 = std::find(arg2->m_Topo.begin() + 1, arg2->m_Topo.end(), LineID(lineID));
    //update the first tree id
    size_t addID = arg2->m_curveList.size();
    for (auto it = arg1->m_Topo.begin() + 1; it != arg1->m_Topo.end(); ++it)
        it->id += addID;
    //move bound info
    arg2->m_traceInitInfo.reserve(arg2->m_traceInitInfo.size() + arg1->m_traceInitInfo.size());
    std::copy(arg1->m_traceInitInfo.begin(), arg1->m_traceInitInfo.end(), std::back_inserter(arg2->m_traceInitInfo));
    //move curves to arg2
    size_t oldSz = arg2->m_curveList.size();
    arg2->m_curveList.resize(oldSz + arg1->m_curveList.size());
    for (size_t i = oldSz; i < arg2->m_curveList.size(); ++i) {
        arg2->m_curveList[i].swap(arg1->m_curveList[i - oldSz]);
    }
    if (!paramPack->UserRecord.empty()){
        paramPack->activeTree->m_ChangedBranch.clear();
        paramPack->UserRecord.clear();
        paramPack->activeTree->m_DelBranch.clear();
    }
    arg1->m_curveList.clear();
    //If connect to root, move as sibling of arg2, else move 1st tree as the subtree of the 2nd tree
    if (vertID < 1 && arg2->m_Topo.depth(pIter2) == 1lu) arg2->m_Topo.reparent(arg2->m_Topo.parent(pIter2), arg1->m_Topo.begin());
    else  arg2->m_Topo.reparent(pIter2, arg1->m_Topo.begin());
    return true;
}

void NeuroGLWidget::UndoDraw()
{
    if (mode_ == DRAWLINEMODE) {
        if (paramPack->manualLabelCurve_.empty()){
            if (drawOldPoint_){
                delete[] drawOldPoint_; drawOldPoint_ = 0; isDrawHasStartPoint_ = false;
            }
            return;
        }
        paramPack->manualLabelCurve_.pop_back();
        BuildLineActor(drawLineActor_, paramPack->manualLabelCurve_, true);
        BuildLineActor(drawVertActor_, paramPack->manualLabelCurve_, false);
        if (paramPack->manualLabelCurve_.empty()){
            if (drawOldPoint_){
                delete[] drawOldPoint_; drawOldPoint_ = 0; isDrawHasStartPoint_ = false;
            }
            return;
        }
        //update();
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::BuildLineActor(vtkSmartPointer<vtkActor> actor, const Line5d& line, bool flag)
{
    if (flag) {//build line actor
        VTK_CREATE(vtkPoints, linePoints);
        VTK_CREATE(vtkCellArray, lineCells);
        VTK_CREATE(vtkPolyData, linePoly);
        VTK_CREATE(vtkPolyDataMapper, lineMapper);

        lineCells->InsertNextCell(int(line.size()));
        for (size_t k = 0; k < line.size(); ++k) {
            linePoints->InsertNextPoint(line[k](0) * paramPack->xRes_, line[k](1) * paramPack->yRes_, line[k](2) * paramPack->zRes_);
            lineCells->InsertCellPoint(k);
        }
        linePoly->SetPoints(linePoints);
        linePoly->SetLines(lineCells);
        linePoly->Modified();
        VTK_DEFINE(vtkPolyDataMapper, mapper) = vtkPolyDataMapper::SafeDownCast(actor->GetMapper());
        mapper->SetInputData(linePoly);
        mapper->Modified();
        mapper->Update();
        actor->Modified();
    }
    else{//build vertex actor
        VTK_CREATE(vtkPoints, vertPoints);
        VTK_CREATE(vtkCellArray, vertCells);
        VTK_CREATE(vtkPolyData, vertPoly);
        VTK_CREATE(vtkPolyDataMapper, vertMapper);
        vertCells->InsertNextCell(int(line.size()));
        for (size_t k = 0; k < line.size(); ++k) {
            vertPoints->InsertNextPoint(line[k](0) * paramPack->xRes_, line[k](1) * paramPack->yRes_, line[k](2) * paramPack->zRes_);
            vertCells->InsertCellPoint(k);
        }
        vertPoly->SetPoints(vertPoints);
        vertPoly->SetVerts(vertCells);
        vertPoly->Modified();
        VTK_DEFINE(vtkPolyDataMapper, mapper) = vtkPolyDataMapper::SafeDownCast(actor->GetMapper());
        mapper->SetInputData(vertPoly);
        mapper->Modified();
        mapper->Update();
        actor->Modified();
    }
}

vtkSmartPointer<vtkActor> NeuroGLWidget::BuildTreeActor(NeuronTree &nTree, bool active, bool layer)
{
    const std::vector<Line5d>& curTree = nTree.m_curveList;
    VTK_CREATE(vtkPoints, treePoints);
    VTK_CREATE(vtkCellArray, treeCells);
    VTK_CREATE(vtkPolyDataMapper, treeMapper);
    VTK_CREATE(vtkUnsignedCharArray, treeColorArray);
    treeColorArray->SetNumberOfComponents(3);
    for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
        VTK_CREATE(vtkPolyLine, jPolyLine);
        jPolyLine->GetPointIds()->SetNumberOfIds(int(curTree[j].size()));
        for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
            const Vec5d &kPoint = curTree[j].at(k);
            vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
            jPolyLine->GetPointIds()->SetId(k, xId);
        }
        treeCells->InsertNextCell(jPolyLine);
        if (active) {
            if (layer){
                auto &curTopo = nTree.m_Topo;
                auto it = std::find(curTopo.begin() + 1, curTopo.end(), LineID(j));
                int level = curTopo.depth(it);
                //treeColorArray->InsertNextTupleValue(traceColor + level * 3 % 15);
                treeColorArray->InsertNextTypedTuple(traceColor + level * 3 % 15);
            }
            else treeColorArray->InsertNextTypedTuple(NORMALCELLCOLOR);
        }
        else treeColorArray->InsertNextTypedTuple(active ? NORMALCELLCOLOR: DEACTIVECOLOR);
    }
    treeColorArray->Modified();
    VTK_CREATE(vtkPolyData, treePoly);
    treePoly->SetPoints(treePoints);
    treePoly->SetLines(treeCells);
    treePoly->GetCellData()->SetScalars(treeColorArray);
    treePoly->BuildCells();
    treePoly->Modified();
    treeMapper->SetInputData(treePoly);
    treeMapper->Update();
    VTK_CREATE(vtkActor, treeActor);
    treeActor->SetMapper(treeMapper);
    return treeActor;
}

void NeuroGLWidget::SetHeadAsTreeRoot()
{
    if (mode_ == BRANCHMODE) {
        if (pickedLineList_.empty()) return;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        KPTREE::SetTreeRoot(*(corNeuronTree->second), corNeuronTree->second->m_curveList[pickedLineList_[0]].front().block(0,0,3,1));
        RebuildActiveTreeActor();
        if (!corNeuronTree->second->m_suspiciousPositions_.empty()) {
            printf("traverse flag will be invalid.\n");
            corNeuronTree->second->m_suspiciousPositions_.clear();
            corNeuronTree->second->m_suspiciousCheckState_.clear();
            emit NeedUpdateCheckWidget();
        }
        if (treeRootActor_ && treeRootActor_->GetVisibility() != 0) {
            ToggleShowTreeRoot(false);
            ToggleShowTreeRoot(true);
        }
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::SetTailAsTreeRoot()
{
    if (mode_ == BRANCHMODE) {
        if (pickedLineList_.empty()) return;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        KPTREE::SetTreeRoot(*(corNeuronTree->second), corNeuronTree->second->m_curveList[pickedLineList_[0]].back().block(0, 0, 3, 1));
        RebuildActiveTreeActor();
        if (!corNeuronTree->second->m_suspiciousPositions_.empty()) {
            printf("traverse flag will be invalid.\n");
            corNeuronTree->second->m_suspiciousPositions_.clear();
            corNeuronTree->second->m_suspiciousCheckState_.clear();
            emit NeedUpdateCheckWidget();
        }
        if (treeRootActor_ && treeRootActor_->GetVisibility() != 0) {
            ToggleShowTreeRoot(false);
            ToggleShowTreeRoot(true);
        }
        GetRenderWindow()->Render();
    }
}

//just save one tree as swc file, no sta file.
void NeuroGLWidget::SaveSelectedTrees()
{
    if (mode_ == TREEMODE) {
        if (pickedTreeList_.size() != 1lu){
            printf("please select one tree for saving.\n");
            return;
        }
        QString savePath = QFileDialog::getSaveFileName(this, "save Tree",
            paramPack->defaultDir, tr("swc (*.swc)"));
        std::vector<std::string> fileList = { savePath.toStdString() };
        auto pop = std::make_shared<NeuronPopulation>();
        pop->m_pop.push_back(pickedTreeList_[0].second);
        auto writer = TreeWriter::New();
        writer->SetInput(pop);
        writer->SetParam(paramPack);
        writer->SetOutputFileName(fileList);
        if (!writer->Update()) {
            printf("nima\n");
        }
    }
}

void NeuroGLWidget::SetActiveTree_Slot()
{
    if (mode_ == TREEMODE) {
        if (pickedTreeList_.size() != 1lu){
            printf("please select one tree to active.\n");
            return;
        }
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(pickedTreeList_[0].second)); });
        int arg = std::distance(neuronActorList_.begin(), corNeuronTree);
        RebuildActiveTreeChangedPopulation(arg);
        emit NeedUpdateTreeWidget();
    }
}

void NeuroGLWidget::RebuildActiveTreeChangedPopulation(int arg)
{
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
    if (!tmpPop) {
        printf("there is no tree.\n");
        return;
    }
    bool isoldactive = false;
    for (size_t k = 0; k < tmpPop->m_pop.size(); ++k) {
        if (tmpPop->m_pop[k].get() == paramPack->activeTree.get()) {
            if (arg == int(k)) {
                isoldactive = true;
                break;
            }
        }
    }
    if (isoldactive) { // do nothing
        auto oldActive = paramPack->activeTree;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(oldActive)); });
        RebuildNeuronTreeActor(*corNeuronTree, false);
        pickedTreeList_.clear();
        paramPack->activeTree.reset();
        activeTreeActor_ = nullptr;
        corNeuronTree->first->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    }
    else{
        if (paramPack->activeTree) {
            pickedLineList_.clear();//remove pick information
            auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(paramPack->activeTree)); });
            RebuildNeuronTreeActor(*corNeuronTree, false);
            corNeuronTree->first->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
        }
        paramPack->activeTree = tmpPop->m_pop[arg];//pickedTreeList_[0].second;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(paramPack->activeTree)); });
        activeTreeActor_ = corNeuronTree->first;//pickedTreeList_[0].first.GetPointer();
        //pickedTreeList_.clear();
        //RebuildNeuronTreeActor(*corNeuronTree, true);
        RebuildActiveTreeActor();
        corNeuronTree->first->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::RebuildActiveTreeActor()
{
    if (paramPack->activeTree) {
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.second.get() == paramPack->activeTree.get(); });
        if (layerDisplay_.empty() || layerDisplay_[0] == 0) {//orig tree
            actorRenderer_->RemoveActor(otherLayerActor_);
            RebuildNeuronTreeActor(*corNeuronTree, true);
            BuildCurInitTracePt();
            BuildCurrentBoundaryCollection();
            BuildHistoryBoundaryActorCollection();
            //restore
            std::for_each(annotationFlagList_.begin(), annotationFlagList_.end(), [](annotationFlagClass &arg){arg.second->SetVisibility(true); });
            std::for_each(currentBoundaryActorCollection_.begin(), currentBoundaryActorCollection_.end(), [](vtkSmartPointer<vtkActor> arg){arg->SetVisibility(true); });
            //for (auto &iter : neuronActorList_) {
            //    if (iter.first.GetPointer() == activeTreeActor_) continue;
            //    iter.first->SetVisibility(true);
            //}
        }
        else{//layer display
            RebuildLayerTreeActor(layerDisplay_);
            BuildCurInitTracePt();
            BuildCurrentBoundaryCollection();
            BuildHistoryBoundaryActorCollection();
        }
        //update();
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::UnDeleteOneActiveTreeBranch()
{
    
    if (paramPack->activeTree) {
        bool flag = false;
        if (!paramPack->UserRecord.empty()){
            switch (paramPack->UserRecord.back())
            {
            case Record::DeleteBranch:
                if (paramPack->activeTree->m_DelBranch.empty()) {
                    printf("no deleted branch record.\n");
                }
                else
                {
                    flag = true;
                    if (!MergeTwoTrees(paramPack->activeTree->m_DelBranch.back().second, paramPack->activeTree)) return;
                    if (!paramPack->activeTree->m_DelBranch.empty()) paramPack->activeTree->m_DelBranch.clear();
                }
                break;
            case Record::SmoothCurve:
                if (paramPack->activeTree->m_ChangedBranch.empty()){
                    printf("no changed branch record.\n");
                }
                else
                {
                    flag = true;
                    std::vector<std::pair<size_t, Line5d> > ChangedBranch = paramPack->activeTree->m_ChangedBranch.back();
                    for each(auto arg in ChangedBranch){
                        paramPack->activeTree->m_curveList[arg.first].swap(arg.second);
                    }
                    if (!paramPack->activeTree->m_ChangedBranch.empty()) 
                        paramPack->activeTree->m_ChangedBranch.clear(); 
                }
                break;
            default:
                break;
            }
            paramPack->UserRecord.clear();
        }
        if (paramPack->UserRecord.empty()){
            pUnDel->setDisabled(true);
        }
        if (flag){
            RebuildActiveTreeActor();
            BuildHistoryBoundaryActorCollection();
            BuildCurrentBoundaryCollection();
        }
        else
        {
            QMessageBox::information(this, "Notice", "there is no work to do!");
        }

    }
}

void NeuroGLWidget::AppendNeuronTreeActor(NeuronTreePointer arg, bool flag)
{
    auto treeActor = BuildTreeActor(*arg, flag);
    treeActor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    neuronActorList_.push_back(std::make_pair(treeActor, arg));
    actorRenderer_->AddActor(treeActor);
}

void NeuroGLWidget::ToggleShowTreeRoot(bool arg)
{
    if (arg) {
        if (treeRootActor_) {
            actorRenderer_->RemoveActor(treeRootActor_);
        }
        if (!treeRootActor_) treeRootActor_ = vtkActor::New();
        if (paramPack->activeTree) {
            BuildSphereActor(treeRootActor_, paramPack->activeTree->m_curveList[(*(paramPack->activeTree->m_Topo.begin() + 1)).id].front().block(0, 0, 3, 1));
            treeRootActor_->GetProperty()->SetColor(255, 0, 0);
            actorRenderer_->AddActor(treeRootActor_);
        }
        treeRootActor_->SetVisibility(1);
    }
    else
    {
        if (treeRootActor_) {
            actorRenderer_->RemoveActor(treeRootActor_);
            treeRootActor_->SetVisibility(0);
        }
    }
}

void NeuroGLWidget::BuildSphereActor(vtkSmartPointer<vtkActor> actor, const Vec3d& pos)
{
    VTK_CREATE(vtkSphereSource, tmpSphere);
    tmpSphere->SetCenter(pos(0) * paramPack->xRes_, pos(1) * paramPack->yRes_, pos(2) * paramPack->zRes_);
    double minRes = paramPack->xRes_ < paramPack->yRes_ ?
        (paramPack->xRes_ < paramPack->zRes_ ? paramPack->xRes_ : paramPack->zRes_) :
        (paramPack->yRes_ < paramPack->zRes_ ? paramPack->yRes_ : paramPack->zRes_);
    tmpSphere->SetRadius(minRes * 5.0);
    VTK_CREATE(vtkPolyDataMapper, tmpSphereMapper);
    tmpSphereMapper->SetInputConnection(tmpSphere->GetOutputPort());
    actor->SetMapper(tmpSphereMapper);
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::RebuildNeuronTreeActor(std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer> &arg, bool active)
{
    //VTK_DEFINE(vtkPolyData, treePoly) = vtkPolyData::SafeDownCast(arg.first->GetMapper()->GetInput());
    VTK_CREATE(vtkPolyData, treePoly);
    VTK_CREATE(vtkPoints, treePoints);
    VTK_CREATE(vtkCellArray, treeCells);
    VTK_CREATE(vtkPolyDataMapper, treeMapper);
    //VTK_DEFINE(vtkPolyDataMapper, treeMapper) = vtkPolyDataMapper::SafeDownCast(arg.first->GetMapper());
    VTK_CREATE(vtkUnsignedCharArray, treeColorArray);
    treeColorArray->SetNumberOfComponents(3);
    //save the pick line info
    if (!pickedLineList_.empty()) {
        for (auto &it : pickedLineList_)
            ChangePickedCellColor(it, false);
    }
    const std::vector<Line5d>& curTree = (*arg.second).m_curveList;
    for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
        VTK_CREATE(vtkPolyLine, jPolyLine);
        jPolyLine->GetPointIds()->SetNumberOfIds(int(curTree[j].size()));
        for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
            const Vec5d &kPoint = curTree[j].at(k);
            vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
            jPolyLine->GetPointIds()->SetId(k, xId);
        }
        treeCells->InsertNextCell(jPolyLine);
        //treeColorArray->InsertNextTupleValue(active ? NORMALCELLCOLOR : DEACTIVECOLOR);//InsertNextTypedTuple
        treeColorArray->InsertNextTypedTuple(active ? NORMALCELLCOLOR : DEACTIVECOLOR);//InsertNextTypedTuple
    }
    treeColorArray->Modified();
    treePoly->SetPoints(treePoints);
    treePoly->SetLines(treeCells);
    treePoly->GetCellData()->SetScalars(treeColorArray);
    treePoly->BuildCells();
    treePoly->Modified();
    treeMapper->SetInputData(treePoly);
    treeMapper->Update();
    treeMapper->Modified();
    arg.first->SetMapper(treeMapper);
    //arg.first->Modified();
    //restore the pick line info
    if (!pickedLineList_.empty()) {
        for (auto &it : pickedLineList_)
            ChangePickedCellColor(it, true);
    }
    int kk = arg.first->GetMapper()->GetInput()->GetNumberOfCells();
    if (kk != paramPack->activeTree->m_curveList.size() ) {
        QMessageBox::warning(this, "Caught bug!", "If you use debug mode, call zhouhang ");
        system("pause");
    }
}

void NeuroGLWidget::ToggleShowCurrentLinesDirection(bool arg)
{
    if (mode_ == BRANCHMODE && paramPack->activeTree) {
        if (arg) {
            size_t len = 0lu;
            size_t pickedIdNum = pickedLineList_.size();
            ClearActorCollection(lineDirectionActor_);
            lineDirectionActor_.clear();
            if (pickedIdNum != 0){
                for (size_t i = 0; i < pickedIdNum; ++i){
                    len = paramPack->activeTree->m_curveList[pickedLineList_[i]].size();
                    if (len < 2lu){
                        printf("why there is a dot ? debug.\n");
                        return;
                    }
                    else if (len == 2lu || len == 3lu){
                        auto actor = CreateArrowActor(paramPack->activeTree->m_curveList[pickedLineList_[i]].front().block(0, 0, 3, 1), paramPack->activeTree->m_curveList[pickedLineList_[i]].back().block(0, 0, 3, 1));
                        actorRenderer_->AddActor(actor);
                        lineDirectionActor_.push_back(actor);
                    }
                    else if (len > 3lu) {
                        auto actor = CreateArrowActor(paramPack->activeTree->m_curveList[pickedLineList_[i]][len - 3].block(0,0,3,1), paramPack->activeTree->m_curveList[pickedLineList_[i]][len - 1].block(0,0,3,1));
                        actorRenderer_->AddActor(actor);
                        lineDirectionActor_.push_back(actor);
                    }
                }
            }
            else{
                pShowDirection->setChecked(false);
            }
        }
        else
        {
            ClearActorCollection(lineDirectionActor_);
            lineDirectionActor_.clear();
        }
    }
}

void NeuroGLWidget::ShowPosition_Slot()
{
    if (POINTMODE == mode_ && pickedVertexID_.first != -1 && pickedVertexID_.second != -1){//canPickPoint_ && 
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        Vec3d curPos = corNeuronTree->second->m_curveList[pickedVertexID_.first][pickedVertexID_.second].block(0, 0, 3, 1);
        NG_CREATE_DYNAMIC_CAST(SVolume, tmpImg, paramPack->OrigImage);
        Vec3d tmp; tmp << curPos(0) - paramPack->xMin_, curPos(1) - paramPack->yMin_, curPos(2) - paramPack->zMin_;
        int depth, branchNum; NEURITETYPE type;
        bool restoreFlag = false;
        if (pickedLineList_.empty()){
            pickedLineList_.push_back(pickedVertexID_.first);
            restoreFlag = true;
        }
        if (!GetCurveLevelAndBranch(depth, branchNum, type))
            return;
        if (restoreFlag) pickedLineList_.clear();
        QMessageBox::information(this, "picked point position",
            QString(tr("global:%1 %2 %3\nlocal:%4 %5 %6\nbranchdepth:%7, id: %8 fore:%9\nType: %10"
            ).arg(curPos(0)).arg(curPos(1)).arg(curPos(2)).arg(curPos(0) - paramPack->xMin_).arg(curPos(1) - paramPack->yMin_).arg(curPos(2) - paramPack->zMin_
            ).arg(depth).arg(branchNum
            ).arg(NGUtility::WeighRayValue(tmp, *tmpImg)).arg(
            type == NEURITETYPE::UNKNOWN ? "Unknown" : (
            type == NEURITETYPE::SOMA ? "Soma" : (
            type == NEURITETYPE::AXON ? "Axon" : (
            type == NEURITETYPE::DENDRITE ? "basal Dendrite" : (
            type == NEURITETYPE::APICAL ? "apical Dendrite" : "Special")))))));
    }
}

void NeuroGLWidget::ToggleEnableCurveTrace()
{
    if (mode_ == BRANCHMODE) {
        if (!pickedLineList_.empty()) {
            std::vector<Line5d> &oldLines = paramPack->activeTree->m_curveList;
            Vec3d endPt = oldLines[pickedLineList_[0]].back().block(0, 0, 3, 1);
            //search if have trace info
            bool flag = false;
            BOUNDINFOLIST& boundInfo = paramPack->activeTree->m_traceInitInfo;
            for (auto it = boundInfo.begin(); it != boundInfo.end(); ++it) {
                if ( ((*it).first - endPt).norm() < 0.1 ) {
                    flag = true;
                    curForbiddenFoundBoundInfo.push_back((*it).first);;
                    boundInfo.erase(it);
                    break;
                }
            }
            //search curBoundInfo
            if (!flag) {
                for (auto it = curBoundInfo.begin(); it != curBoundInfo.end(); ++it) {
                    if (((*it).first - endPt).norm() < 0.1) {
                        flag = true;
                        curForbiddenFoundBoundInfo.push_back((*it).first);;
                        curBoundInfo.erase(it);
                        break;
                    }
                }
            }
            //there is no trace info then add it
            if (!flag) {
                Vec3d tmpDir;
                VectorVec3d tmpTail;
                int sz = int(oldLines[pickedLineList_[0]].size());
                NGUtility::GetPartVectorVec3d(oldLines[pickedLineList_[0]], std::max(0, sz - 5), sz - 1, tmpTail);// <= maxid
                CalcPtCurveDirection(tmpTail, tmpDir);
                boundInfo.push_back(std::make_pair(endPt, tmpDir));
                //
                auto it = std::find_if(curForbiddenFoundBoundInfo.begin(), curForbiddenFoundBoundInfo.end(), [&](const Vec3d& arg){return (endPt - arg).norm() < 2.0; });
                if (it != curForbiddenFoundBoundInfo.end()) {
                    curForbiddenFoundBoundInfo.erase(it);
                }
            }
            BuildHistoryBoundaryActorCollection();
            BuildCurrentBoundaryCollection();
        }
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::DeactiveAllTrees()
{
    if (paramPack->activeTree) {
        auto oldActive = paramPack->activeTree;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return &(*(arg.second)) == &(*(oldActive)); });
        RebuildNeuronTreeActor(*corNeuronTree, false);
        pickedTreeList_.clear();
        paramPack->activeTree.reset();
        activeTreeActor_ = nullptr;
        corNeuronTree->first->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    }
}

void NeuroGLWidget::ToggleShowInit(bool arg)
{
    if (initPtActor_)
        initPtActor_->SetVisibility(arg);
    SetActorCollectionVisible(currentBoundaryActorCollection_, arg);
    SetActorCollectionVisible(historyBoundaryActorCollection_,arg);
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::UndeleteTrees()
{
    if (mode_ == TREEMODE) {
        NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
        auto &delTreeList = tmpTree->deletedPop;
        if (delTreeList.empty()) {
            QMessageBox::information(this, "no more deleted trees", "there is no more deleted trees. no operation work.");
            return;
        }
        for (auto & it : delTreeList) AppendNeuronTreeActor(it);
        std::move(delTreeList.begin(), delTreeList.end(), std::back_inserter(tmpTree->m_pop));
        delTreeList.clear();
        ExitTreeMode();
        EnterTreeMode();
    }
}

void NeuroGLWidget::errorcheck()
{
    for (auto it = paramPack->activeTree->m_Topo.begin() + 2; it != paramPack->activeTree->m_Topo.end(); ++it) {
        auto curLineID = it->id;
        const Line5d &line = paramPack->activeTree->m_curveList[curLineID];
        if (paramPack->activeTree->m_Topo.parent(it) == paramPack->activeTree->m_Topo.begin()){
        }
        else{
            size_t parentLineID = paramPack->activeTree->m_Topo.parent(it)->id;
            const Vec3d &headNode = line.front().block(0, 0, 3, 1);
            size_t nearestID = KPTREE::FindNearestID(headNode, paramPack->activeTree->m_curveList[parentLineID]);
            if (nearestID > 1000000){
                printf("error\n");
            }
        }
    }
}

void NeuroGLWidget::keyPressEvent(QKeyEvent* e)
{
    if (IsForbidden()) return;
    switch (e->key()){
    case Qt::Key_C:
        if (mode_ == TRAVERSEMODE || isMovePointMode_ || isReconnectParent_) return;
        if (!pChooseBranchMode->isChecked()) {
            pChooseNodeMode->setChecked(false);
            pChooseSomaMode->setChecked(false);
            pChooseTreeMode->setChecked(false);
            pDrawLineMode->setChecked(false);
            pChooseBranchMode->setChecked(true);
        }
        else pChooseBranchMode->setChecked(false);
        break;
    case Qt::Key_D:
        pDeleteBranch->trigger();
        break;
    case Qt::Key_E:
        ToggleEnableCurveTrace();
        break;
    case Qt::Key_Z:
        if (mode_ == TRAVERSEMODE || isMovePointMode_ || isReconnectParent_) return;
        if (!pChooseNodeMode->isChecked()) {
            pChooseSomaMode->setChecked(false);
            pChooseTreeMode->setChecked(false);
            pDrawLineMode->setChecked(false);
            pChooseBranchMode->setChecked(false);
            pChooseNodeMode->setChecked(true);
        }
        else pChooseNodeMode->setChecked(false);
        break;
    case Qt::Key_X:
        pCutLine->trigger();
        break;
    case Qt::Key_V:
        if (mode_ == TRAVERSEMODE || isMovePointMode_ || isReconnectParent_) return;
        if (!pDrawLineMode->isChecked()) {
            pChooseSomaMode->setChecked(false);
            pChooseTreeMode->setChecked(false);
            pChooseBranchMode->setChecked(false);
            pChooseNodeMode->setChecked(false);
            pDrawLineMode->setChecked(true);
        }
        else pDrawLineMode->setChecked(false);
        break;
    case Qt::Key_Q:
        pInteract->toggle();
        break;
    case Qt::Key_U:
        pUndoDraw->trigger();
        break;
    case Qt::Key_B:
        pBoundCheck->toggle();
        break;
    case Qt::Key_M:
        if (mode_ == POINTMODE) 
            pStartMovePt_->toggle();
        break;
    case Qt::Key_N:
        emit NextTraverse_Signal();
        break;
    case Qt::Key_I:
        if (mode_ == POINTMODE) {
            pInsertPt_->trigger();
        }
        break;
	case Qt::Key_S:
		pSmoothCurve->trigger();
		break;
    default:
        break;
    }
    QVTKOpenGLWidget::keyPressEvent(e);
}

void NeuroGLWidget::ForceBoundCheckAll()
{
    isForceBoundCheckAll_ = true;
    BuildCurrentBoundaryCollection();
    isForceBoundCheckAll_ = false;
}

void NeuroGLWidget::ToggleActiveTreeLayerDisplay(bool arg)
{
    if (arg) {
    }
    else{
        for (auto &iter : neuronActorList_) {
            if (iter.first.GetPointer() == activeTreeActor_) continue;
            RebuildNeuronTreeActor(iter, false);
        }
    }
    isActiveTreeLayerDisplay_ = arg;
    RebuildActiveTreeActor();
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::RebuildLayerTreeActor(std::vector<int>& activeLayer)
{
    if (activeLayer.empty()) return;
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpImg, paramPack->LayerImage);
    if (!tmpImg) {
        NG_ERROR_MESSAGE("layerImage is null.");
        return;
    }
    
    actorRenderer_->RemoveActor(otherLayerActor_);
    VTK_CREATE(vtkPolyData, treePoly);
    VTK_CREATE(vtkPoints, treePoints);
    VTK_CREATE(vtkCellArray, treeCells);
    VTK_DEFINE(vtkPolyDataMapper, treeMapper) = vtkPolyDataMapper::SafeDownCast(activeTreeActor_->GetMapper());
    VTK_CREATE(vtkUnsignedCharArray, treeColorArray);
    treeColorArray->SetNumberOfComponents(3);
    const std::vector<Line5d>& curTree = paramPack->activeTree->m_curveList;
    //save the pick line info
    if (!pickedLineList_.empty()) {
        for (auto &it : pickedLineList_) 
            ChangePickedCellColor(it, false);
    }
    //auto &topo = paramPack->activeTree->m_Topo;
    size_t xMin = size_t(paramPack->xMin_); size_t yMin = size_t(paramPack->yMin_); size_t zMin = size_t(paramPack->zMin_);
    size_t xLen = size_t(paramPack->xMax_) - xMin, yLen = size_t(paramPack->yMax_) - yMin, zLen = size_t(paramPack->zMax_) - zMin;
    size_t x, y, z;
    //find the special curve id
    auto &curTopo = paramPack->activeTree->m_Topo;
    std::vector<size_t> parentID;//search parent id
    
    std::vector<decltype(curTopo.begin())> firstLevel;
    for (auto it = curTopo.begin() + 1; it != curTopo.end(); ++it) {//
        if (layerDisplay_[0] == curTopo.depth(it)){
            firstLevel.push_back(it);
            auto pit = curTopo.parent(it);
            if (pit.node != NULL) {
                parentID.push_back(pit->id);
            }
        }
    }
    //search level subtree 
    levelCurveID_.clear();
    for (size_t i = 0; i < branchDisplay_.size(); ++i) {
        if (branchDisplay_[i] >= firstLevel.size()) break;
        auto subTree = curTopo.subtree(firstLevel[branchDisplay_[i]], curTopo.next_sibling(firstLevel[branchDisplay_[i]]));//level minus layerList [0]
        for (auto iter = subTree.begin(); iter != subTree.end(); ++iter)
            if (NGUtility::IsInVector(layerDisplay_, curTopo.depth(iter) + layerDisplay_[0]))
                levelCurveID_.push_back(iter->id);
    }
    //
    VectorVec5d otherPtSet;
    //allVisibleLayerCurveID_.clear();
    //construct line and parent curve dot
    for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
        VTK_CREATE(vtkPolyLine, jPolyLine);
        if(NGUtility::IsInVector(levelCurveID_, j)){//visible curve
            jPolyLine->GetPointIds()->SetNumberOfIds(int(curTree[j].size()));
            for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
                const Vec5d &kPoint = curTree[j][k];
                vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
                jPolyLine->GetPointIds()->SetId(k, xId);
            }
            //allVisibleLayerCurveID_.push_back(j);
        }
        else if (NGUtility::IsInVector(parentID, j)) {//half visible curve
            //collect valid point
            //VectorVec5d tmpPtSet;
            for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
                const Vec5d &kPoint = curTree[j][k];
                x = size_t(kPoint[0]) - xMin; y = size_t(kPoint[1]) - yMin; z = size_t(kPoint[2]) - zMin;
                if (x < 0 || x > xLen || y < 0 || y > yLen || z < 0 || z > zLen || tmpImg->operator()(x, y, z) < 1) continue;
                otherPtSet.push_back(kPoint);
            }
            //allVisibleLayerCurveID_.push_back(j);
            continue;
        }
        else {//at least two point, avoid vtk pointer error when choose lines. put into original point
            VectorVec5d tmpLine;
            for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
                const Vec5d &kPoint = curTree[j][k];
                x = size_t(kPoint[0]) - xMin; y = size_t(kPoint[1]) - yMin; z = size_t(kPoint[2]) - zMin;
                if (x < 0 || x > xLen || y < 0 || y > yLen || z < 0 || z > zLen || tmpImg->operator()(x, y, z) < 1) break;
                tmpLine.push_back(kPoint);
            }
            if (tmpLine.empty()){ 
                tmpLine.push_back(curTree[j][0]);//soma
                tmpLine.push_back(curTree[j][0]);
            }
            else if (tmpLine.size() == 1){
                tmpLine.push_back(tmpLine[0]);//first node
            }
            jPolyLine->GetPointIds()->SetNumberOfIds(int(tmpLine.size()));
            for (VectorVec5d::size_type k = 0; k < tmpLine.size(); ++k){
                const Vec5d &kPoint = tmpLine[k];
                vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
                jPolyLine->GetPointIds()->SetId(k, xId);
            }
        }
        treeCells->InsertNextCell(jPolyLine);
        auto &curTopo = paramPack->activeTree->m_Topo;
        auto it = std::find(curTopo.begin() + 1, curTopo.end(), LineID(j));
        int level = curTopo.depth(it);
        //treeColorArray->InsertNextTupleValue(traceColor + (level-1) * 3 % 11);//no RGB
        treeColorArray->InsertNextTypedTuple(traceColor + (level - 1) * 3 % 11);//no RGB
    }
    treeColorArray->Modified();
    treePoly->SetPoints(treePoints);
    treePoly->SetLines(treeCells);
    treePoly->GetCellData()->SetScalars(treeColorArray);
    treePoly->BuildCells();
    treePoly->Modified();
    treeMapper->SetInputData(treePoly);
    treeMapper->Update();
    treeMapper->Modified();
    activeTreeActor_->Modified();
    BuildVertexActor(otherPtSet, otherLayerActor_);
    actorRenderer_->AddActor(otherLayerActor_);
    //restore the pick line info
    if (!pickedLineList_.empty()) {
        for (auto &it : pickedLineList_)
            ChangePickedCellColor(it, true);
    }
    //set other trees
    for (auto &iter : neuronActorList_) {//invisible other trees
        if (iter.first.GetPointer() == activeTreeActor_) continue;
        VectorVec5d otherTreePtSet;
        const std::vector<Line5d>& curTree = iter.second->m_curveList;
        for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
            for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
                const Vec5d &kPoint = curTree[j].at(k);
                x = size_t(kPoint[0]) - xMin; y = size_t(kPoint[1]) - yMin; z = size_t(kPoint[2]) - zMin;
                if (x < 0 || x > xLen || y < 0 || y > yLen || z < 0 || z > zLen || tmpImg->operator()(x, y, z) < 1) continue;
                otherTreePtSet.push_back(kPoint);
            }
        }
        //actorRenderer_->RemoveActor(iter.first);
        BuildVertexActor(otherTreePtSet, iter.first);
        actorRenderer_->AddActor(iter.first);
    }
    //set annotation flag visible
    Vec3i flagPos;
    for (auto &it:annotationFlagList_) {
        flagPos << NGUtility::Round(it.first(0)) - paramPack->xMin_, NGUtility::Round(it.first(1)) - paramPack->yMin_, NGUtility::Round(it.first(2)) - paramPack->zMin_;
        if (flagPos(0) < 0 || flagPos(0) > xLen || flagPos(1) < 0 || flagPos(1) > yLen || flagPos(2) < 0 || flagPos(2) > zLen ||
            tmpImg->operator()(flagPos(0), flagPos(1), flagPos(2)) < 1)
            it.second->SetVisibility(false);
    }
    //set arrow flag visible
    for (auto &it : currentBoundaryActorCollection_) {
        double* pos = it->GetPosition();
        flagPos << NGUtility::Round(pos[0]/paramPack->xRes_) - paramPack->xMin_, 
            NGUtility::Round(pos[1]/paramPack->yRes_) - paramPack->yMin_, NGUtility::Round(pos[2]/paramPack->zRes_) - paramPack->zMin_;
        if (flagPos(0) < 0 || flagPos(0) > xLen || flagPos(1) < 0 || flagPos(1) > yLen || flagPos(2) < 0 || flagPos(2) > zLen ||
            tmpImg->operator()(flagPos(0), flagPos(1), flagPos(2)) < 1)
            it->SetVisibility(false);
    }
}

void NeuroGLWidget::Annotation_Slot()
{
    if (POINTMODE == mode_ &&pickedVertexID_.first != -1 && pickedVertexID_.second != -1){// canPickPoint_ && 
        /*cut line into 2 lines, the second is the child of the first. search the children connect to original line, and identify which connect to which.*/
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        size_t idsNum = corNeuronTree->second->m_curveList[pickedVertexID_.first].size();
        if (idsNum <= 2llu || pickedVertexID_.second < 1 || pickedVertexID_.second >= int(idsNum) - 1){
            printf("cut line is too short, or the cut dot is the end node.\n");
            return;
        }
        auto &pos = corNeuronTree->second->m_curveList[pickedVertexID_.first][pickedVertexID_.second];
        /*search if repeat*/
        bool repeat = false;
        double dist = 0.0;
        for (auto &it : annotationFlagList_) {
            dist = (it.first - pos.block(0, 0, 3, 1)).norm();
            if (dist < 1.0) {
                repeat = true;
                break;
            }
        }
        if (repeat) {
            QMessageBox::information(this, "Invalid Operation", "The position where you choose has been annotated.");
            return;
        }
        //make flag
        vtkSmartPointer<vtkPNGReader> reader =vtkSmartPointer<vtkPNGReader>::New();
        reader->SetFileName("Resource/anno.png");//QStringLiteral(":/new/tracer/Resource/open.png"
        reader->Update();
        vtkSmartPointer<vtkTexture> tex = vtkSmartPointer<vtkTexture>::New();
        tex->SetInputData(reader->GetOutput()); tex->InterpolateOn();
        tex->Update();
        vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
        planeSource->Update();
        vtkSmartPointer<vtkTextureMapToPlane> texturePlane = vtkSmartPointer<vtkTextureMapToPlane>::New();
        texturePlane->SetInputConnection(planeSource->GetOutputPort());
        vtkSmartPointer<vtkPolyDataMapper> planeSourceMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        planeSourceMapper->SetInputConnection(texturePlane->GetOutputPort());
        vtkSmartPointer<vtkFollower> actor = vtkSmartPointer<vtkFollower>::New();
        actor->SetMapper(planeSourceMapper);
        actor->SetTexture(tex);
        actor->SetScale(3);
        actor->SetCamera(camera3d_);
        actor->SetPosition(pos(0)*paramPack->xRes_, pos(1)*paramPack->yRes_, pos(2)*paramPack->zRes_);
        actorRenderer_->AddActor(actor); //update();
        GetRenderWindow()->Render();
        //save
        annotationFlagList_.emplace_back(pos.block(0, 0, 3, 1), actor);
        paramPack->activeTree->m_suspiciousPositions_.push_back(pos.block(0, 0, 3, 1));
        paramPack->activeTree->m_suspiciousCheckState_.push_back(0);
        emit NeedUpdateCheckWidget();
    }
}

void NeuroGLWidget::Unannotation_Slot()
{
    if (mode_ == POINTMODE) {
        if (annotationFlagList_.empty()) return;
        if (pickedVertexID_.first < 0 || pickedVertexID_.second < 0) return;
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        size_t idsNum = corNeuronTree->second->m_curveList[pickedVertexID_.first].size();
        if (idsNum <= 2llu || pickedVertexID_.second < 1 || pickedVertexID_.second >= int(idsNum) - 1){
            printf("cut line is too short, or the cut dot is the end node.\n");
            return;
        }
        auto &pos = corNeuronTree->second->m_curveList[pickedVertexID_.first][pickedVertexID_.second];
        /*search the chose point*/
        double dist = 0.0;
        for (size_t i = 0; i < annotationFlagList_.size(); ++i) {
            dist = (annotationFlagList_[i].first - pos.block(0, 0, 3, 1)).norm();
            if (dist < 5.0) {
                actorRenderer_->RemoveActor(annotationFlagList_[i].second);
                annotationFlagList_[i].first.setZero();
            }
        }
        annotationFlagList_.erase(std::remove_if(annotationFlagList_.begin(), annotationFlagList_.end(), [](annotationFlagClass &arg){return arg.first(0) < 0.1; }),
            annotationFlagList_.end()); 
        //update();
        GetRenderWindow()->Render();
        auto &msp = paramPack->activeTree->m_suspiciousPositions_;
        auto iter = std::find_if(msp.begin(), msp.end(), [&](Vec3d& arg){  return (arg - pos.block(0,0,3,1)).norm() < 5.0; });
        if (iter != msp.end()){
            size_t off = std::distance(msp.begin(), iter);
            msp.erase(iter);
            auto &msc = paramPack->activeTree->m_suspiciousCheckState_;
            msc.erase(msc.begin() + off);
        }
        emit NeedUpdateCheckWidget();
    }
}

void NeuroGLWidget::ShowLevelBranchID()
{
    int level, branchId; NEURITETYPE type;
    if (!GetCurveLevelAndBranch(level, branchId, type)){
        printf("error in ShowLevelBranchID\n");
        return;
    }
    if (branchId == 0) {
        QMessageBox::warning(this, "Big bug", "Debug ShowLevelBranchID");
        system("pause");
        return;
    }
    QMessageBox::information(this, "Level and Branch ID", tr("Current curve level id: %1\n branch id: %2\n Type: %3").arg(
        level).arg(branchId).arg(
        type == NEURITETYPE::UNKNOWN ? "Unknown" : (
        type == NEURITETYPE::SOMA ? "Soma" : (
        type == NEURITETYPE::AXON ? "Axon" : (
        type == NEURITETYPE::DENDRITE ? "basal Dendrite" : (
        type == NEURITETYPE::APICAL ? "apical Dendrite" : "Special"))))));
}

void  NeuroGLWidget::BuildVertexActor(const VectorVec5d& tmpPtSet, vtkSmartPointer<vtkActor> actor)
{
    VTK_CREATE(vtkPolyData, vertsPoly);
    VTK_CREATE(vtkPolyDataMapper, vertsMapper);
    VTK_CREATE(vtkCellArray, vertsCells);
    VTK_CREATE(vtkPoints, vertsPoints);
    VTK_CREATE(vtkUnsignedCharArray, vertScalars);
    vertScalars->SetNumberOfComponents(3);
    int pointNum = int(tmpPtSet.size());
    vertsCells->InsertNextCell(pointNum);
    for (vtkIdType i = 0; i < pointNum; ++i){
        const Vec5d &kPoint = tmpPtSet[i];
        vertsPoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
        vertsCells->InsertCellPoint(i);
        //vertScalars->InsertNextTupleValue(DEACTIVECOLOR);
        vertScalars->InsertNextTypedTuple(DEACTIVECOLOR);
    }
    vertsPoly->SetPoints(vertsPoints);
    vertsPoly->SetVerts(vertsCells);
    vertsPoly->GetPointData()->SetScalars(vertScalars);
    vertsMapper->SetInputData(vertsPoly);
    vertsMapper->SetScalarModeToUsePointData();
    //VTK_CREATE(vtkActor, vertsActor);
    actor->SetMapper(vertsMapper);
    actor->GetProperty()->SetPointSize(3.0);
    actor->Modified();
    //return vertsActor;
}

bool NeuroGLWidget::Draw3DPoint(double pickPoint[3])
{
    vtkSmartPointer<vtkWorldPointPicker> picker =   vtkSmartPointer<vtkWorldPointPicker>::New();
    picker->Pick(eventPos_[0], eventPos_[1], 0.0, actorRenderer_);
    double globalPos[3];
    picker->GetPickPosition(globalPos);
    globalPos[0] /= paramPack->xRes_;
    globalPos[1] /= paramPack->yRes_;
    globalPos[2] /= paramPack->zRes_;
    globalPos[0] -= paramPack->xMin_;
    globalPos[1] -= paramPack->yMin_;
    globalPos[2] -= paramPack->zMin_;
    double slope[3];
    actorRenderer_->GetActiveCamera()->GetDirectionOfProjection(slope);
    slope[0] /= paramPack->xRes_;
    slope[1] /= paramPack->yRes_;
    slope[2] /= paramPack->zRes_;
    Vec3d tmp; tmp << slope[0], slope[1], slope[2]; tmp.normalize();
    slope[0] = tmp(0); slope[1] = tmp(1); slope[2] = tmp(2);
    //6 face
    //x = 0 y = 0 z = 0 x = xsize y = ysize z = zsize
    // n1 = [1 0 0] n2 = [010] n3 = [001]
    // p1 = p2 = p3 = [0,0,0] p4 = p5 = p6 = [xsize ysize zsize]
    //line  [x y z]' = [x y z] + [slope[0] slope[1] slope[2]] * t
    //t = (n*p - n*[x y z]) / (n * slope)
    //intersect point = [x y z] + t * slope
    //if 2 points is in the range of volume, then got it
    int vX, vY, vZ;
    if (paramPack->mostdLevel_ == 1) {
        vX = volumeX_; vY = volumeY_; vZ = volumeZ_;
    }
    else if (paramPack->mostdLevel_ > 1 ) {
        int s = int(std::pow(2, paramPack->mostdLevel_ - 1));
        vX = volumeX_ * s;
        vY = volumeY_ * s;
        vZ = volumeZ_ * s;
    }
    else{
        QMessageBox::warning(this, "Draw 3D line", "mostdLevel is error. ");
        return false;
    }
    Vec3d pRay(globalPos[0], globalPos[1], globalPos[2]);
    Vec3d raySlope(slope[0], slope[1], slope[2]);
    VectorVec3d n(6);
    n[0] = Vec3d(1, 0, 0);
    n[1] = Vec3d(0, 1, 0);
    n[2] = Vec3d(0, 0, 1);
    n[3] = -n[0];
    n[4] = -n[1];
    n[5] = -n[2];
    VectorVec3d p(6);
    p[0] = Vec3d(0, 0, 0);
    p[1] = Vec3d(0, 0, 0);
    p[2] = Vec3d(0, 0, 0);
    p[3] = Vec3d(vX, vY, vZ);
    p[4] = p[3];
    p[5] = p[3];
    VectorVec3d intersectPts;
    Vec3d intersectPt;
    for (size_t i = 0; i < 6; ++i){
        double t = n[i].dot(p[i]) - n[i].dot(pRay);
        t /= n[i].dot(raySlope);
        intersectPt = pRay + t * raySlope;
        if (intersectPt.minCoeff() < 0.0 || intersectPt(0) > vX || intersectPt(1) > vY
            || intersectPt(2) > vZ){
            continue;
        }
        intersectPts.push_back(intersectPt);
    }
    if (intersectPts.size() != 2) {
        pickLineNum = 0;
        return false;
    }
    if (raySlope.dot(intersectPts[0] - intersectPts[1]) > 0){
        std::reverse(intersectPts.begin(), intersectPts.end());
    }
    //need 2 lines
    if (pickLineNum == 0) {
        ++pickLineNum;
        pick3DLineList_.clear();
        std::copy(intersectPts.begin(), intersectPts.end(), std::back_inserter(pick3DLineList_));
        //draw line
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> conn = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
        conn->InsertNextCell(2);
        pts->InsertNextPoint((intersectPts[0](0) + paramPack->xMin_)*paramPack->xRes_, (intersectPts[0](1) + paramPack->yMin_)*paramPack->yRes_
            , (intersectPts[0](2) + paramPack->zMin_)*paramPack->zRes_);
        conn->InsertCellPoint(0);
        pts->InsertNextPoint((intersectPts[1](0) + paramPack->xMin_)*paramPack->xRes_, (intersectPts[1](1) + paramPack->yMin_)*paramPack->yRes_
            , (intersectPts[1](2) + paramPack->zMin_)*paramPack->zRes_);
        conn->InsertCellPoint(1);
        poly->SetPoints(pts);
        poly->SetLines(conn);
        poly->SetVerts(conn);
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(poly);
        if(!pickLineActor_) pickLineActor_ = vtkSmartPointer<vtkActor>::New();
        pickLineActor_->SetMapper(mapper);
        pickLineActor_->GetProperty()->SetColor(255, 0, 0);
        pickLineActor_->GetProperty()->SetOpacity(1.0);
        pickLineActor_->GetProperty()->SetPointSize(5);
        actorRenderer_->AddActor(pickLineActor_);
        return false;
    }
    else if (pickLineNum == 1){
        pickLineNum = 0;
        //calculate intersect point as pick point
        Vec3d p1 = pick3DLineList_[1] - pick3DLineList_[0]; p1.normalize();//v3-v1=p1
        Vec3d p2 = intersectPts[1] - intersectPts[0]; p2.normalize();//v4-v2=p2
        //pickLineList[0] ~ v1 intersectPts[0] ~v2
        double c = p1.dot(p2);
        //AX=B
        Mat2d A; A(0, 0) = c; A(0, 1) = -1; A(1, 0) = 1; A(1, 1) = -c; 
        if (A.determinant() < 0.0001) {
            actorRenderer_->RemoveActor(pickLineActor_);
            return false;
        }
        Vec2d B; B << (pick3DLineList_[0] - intersectPts[0]).dot(p1), (pick3DLineList_[0] - intersectPts[0]).dot(p2);
        Mat2d Av = A.inverse();
        Vec2d t = Av*B;//t=[t2 t1]
        Vec3d resPt = intersectPts[0] + t(0) * p2;
        pickPoint[0] = resPt(0) + paramPack->xMin_;
        pickPoint[1] = resPt(1) + paramPack->yMin_;
        pickPoint[2] = resPt(2) + paramPack->zMin_;
        pickPoint[0] *= paramPack->xRes_;
        pickPoint[1] *= paramPack->yRes_;
        pickPoint[2] *= paramPack->zRes_;
        actorRenderer_->RemoveActor(pickLineActor_);
        return true;
    }
    return false;
}

bool NeuroGLWidget::GetCurveLevelAndBranch(int &level, int &branchNum, NEURITETYPE &type)
{
    level = branchNum = 0;
    if (!paramPack || !paramPack->activeTree) {
        NG_ERROR_MESSAGE("WTF");
        return false;
    }
    auto &topo = paramPack->activeTree->m_Topo;
    auto it = std::find(topo.begin(), topo.end(), LineID(pickedLineList_[0]));
    type = it->type;
    level = topo.depth(it);
    branchNum = 1;
    bool flag = false;
    for (auto iter = topo.begin() + 1; iter != topo.end(); ++iter) {
        if (topo.depth(iter) == level) {
            if (iter->id != it->id) ++branchNum;
            else {
                flag = true;
                break;
            }
        }
    }
    return flag;
}

void NeuroGLWidget::ToggleTargetTreeVisible(int arg)
{
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpPop, paramPack->separateTree);
    auto curTree = tmpPop->m_pop[arg];
    auto &curIt = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.second.get() == &(*curTree); });
    if (curIt->first->GetVisibility() == 0) {
        curIt->first->SetVisibility(1);
        actorRenderer_->AddActor(curIt->first);
        //if (mode_ == TREEMODE) cellPicker_->AddPickList(curIt->first);
        cellPicker_->AddPickList(curIt->first);
    }
    else{
        actorRenderer_->RemoveActor(curIt->first);
        curIt->first->SetVisibility(0);
        //if (mode_ == TREEMODE) cellPicker_->DeletePickList(curIt->first);
        cellPicker_->DeletePickList(curIt->first);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::FocusOnPoint(Vec5d &pointArg)
{
    //std::cout << "dddd" << std::endl;
    if (!camera3d_) camera3d_ = vtkSmartPointer<vtkCamera>::New();
    camera3d_->ParallelProjectionOn();
    actorRenderer_->SetActiveCamera(camera3d_);
    volumeRenderer_->SetActiveCamera(camera3d_);
    actorRenderer_->ResetCamera();
    double oldfocus[3];  camera3d_->GetFocalPoint(oldfocus);
    double oldpos[3];  camera3d_->GetPosition(oldpos);
    double focalPoint[3] = { pointArg[0] * paramPack->xRes_, pointArg[1] * paramPack->yRes_, pointArg[2] * paramPack->zRes_ };
    double posOff[3]; posOff[0] = focalPoint[0] - oldfocus[0]; posOff[1] = focalPoint[1] - oldfocus[1]; posOff[2] = focalPoint[2] - oldfocus[2];
    double cameraPos[3] = { oldpos[0] + posOff[0], oldpos[1] + posOff[1], oldpos[2] + posOff[2] };
    double paraScale = 80;
    camera3d_->SetFocalPoint(focalPoint);
    camera3d_->SetPosition(cameraPos);
    camera3d_->SetParallelScale(paraScale);
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::ToggleTraverseMode(bool arg)
{
    if (arg) {
        if (!traverseFlagActor_) traverseFlagActor_ = vtkSmartPointer<vtkActor>::New();
        actorRenderer_->AddActor(traverseFlagActor_);
        cellPicker_->InitializePickList();
        cellPicker_->AddPickList(activeTreeActor_);
        mode_ = TRAVERSEMODE;
    }
    else{
        actorRenderer_->RemoveActor(traverseFlagActor_);
        traverseFlagActor_ = vtkSmartPointer<vtkActor>::New();//delete
        cellPicker_->InitializePickList();
        if (traverseFlagActor_) {
            traverseFlagActor_->SetVisibility(0);
            actorRenderer_->RemoveActor(traverseFlagActor_);
        }
        if (pickedVertexID_.first != -1){
            if (pickedVertexID_.second != -1) {
                ChangePickedVertColor(pickedVertexID_.second, false);
            }
            auto tmp = pickedVertexID_.first;
            ShowVertex(tmp, false);
            pointPicker_->InitializePickList();
            ChangePickedCellColor(tmp, false);
            pickedVertexID_.first = -1;
            pickedVertexID_.second = -1;
        }
        mode_ = NORMALMODE;
    }
}

void NeuroGLWidget::ToggleDiffArrows(bool arg)
{
    if (!paramPack->activeTree) return;
    if (!arg) {
        for (auto &it : diffArrowList_) 
            actorRenderer_->RemoveActor(it);
        diffArrowList_.clear();
    }
    else{
        if (!diffArrowList_.empty()) {
            for (auto &it : diffArrowList_) actorRenderer_->RemoveActor(it);
            diffArrowList_.clear();
        }
        auto &curSusPos = paramPack->activeTree->m_suspiciousPositions_;
        auto &curSusState = paramPack->activeTree->m_suspiciousCheckState_;
        for (size_t k = 0; k < curSusPos.size(); ++k){
            if (curSusState[k] != 0) continue;//do not show checked points and tolerable points
            auto &it = curSusPos[k];
            diffArrowList_.push_back(CreateArrowActor(it.block(0, 0, 3, 1) - Vec3d(15, 15, 15), it.block(0, 0, 3, 1), true));
            actorRenderer_->AddActor(diffArrowList_.back());
        }
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::LightDiffArrow(int arg)
{
    if (!paramPack->activeTree) return;
    if (!diffArrowList_.empty() && arg < diffArrowList_.size()) {
        for (auto &it : diffArrowList_)
            actorRenderer_->RemoveActor(it);
        auto curActor = diffArrowList_[arg];
        actorRenderer_->AddActor(curActor);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::UpdateTraverseFlagActor()
{
    if (mode_ != TRAVERSEMODE) return;
    if (!paramPack->activeTree) {
        NG_ERROR_MESSAGE("");
        return;
    }
    //invisible other tree
    //auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.second.get() == paramPack->activeTree.get(); });
    //actorRenderer_->AddActor(corNeuronTree->first);
    //rebuild traverse flag tree
    VTK_CREATE(vtkPolyData, treePoly);
    VTK_CREATE(vtkPoints, treePoints);
    VTK_CREATE(vtkCellArray, treeCells);
    //VTK_CREATE(vtkPolyDataMapper, treeMapper);
    VTK_DEFINE(vtkPolyDataMapper, treeMapper) = vtkPolyDataMapper::SafeDownCast(traverseFlagActor_->GetMapper());
    if (!treeMapper) {
        treeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        traverseFlagActor_->SetMapper(treeMapper);
    }
    VTK_CREATE(vtkUnsignedCharArray, treeColorArray);
    treeColorArray->SetNumberOfComponents(3);
    const std::vector<Line5d>& curTree = paramPack->activeTree->m_curveList;
    int numpoint = 0;
    //count the total number of traverse points
    for (size_t i = 0; i < curTree.size(); ++i) {
        auto &curLine = curTree[i];
        for (auto &it : curLine) {
        if (it(4) > 0.5) ++numpoint;
        }
    }
    if (numpoint == 0) return;
    //
    treeCells->InsertNextCell(numpoint);
    numpoint = 0;
    for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
        auto &curLine = curTree[j];
        for (VectorVec5d::size_type k = 0; k < curLine.size(); ++k){
            const Vec5d &kPoint = curTree[j][k];
            if (kPoint(4) > 0.5) {
                treePoints->InsertNextPoint(kPoint(0) * paramPack->xRes_, kPoint(1) * paramPack->yRes_, kPoint(2) * paramPack->zRes_);
                treeCells->InsertCellPoint(numpoint);
                //treeColorArray->InsertNextTupleValue(DEACTIVECOLOR);//DEACTIVECOLOR
                treeColorArray->InsertNextTypedTuple(DEACTIVECOLOR);//DEACTIVECOLOR
                ++numpoint;
            }
        }
    }
    if (numpoint == 0) return;
    treeColorArray->Modified();
    treePoly->SetPoints(treePoints);
    treePoly->SetVerts(treeCells);
    treePoly->GetCellData()->SetScalars(treeColorArray);
    treePoly->BuildCells();
    treePoly->Modified();
    treeMapper->SetInputData(treePoly);
    treeMapper->Update();
    treeMapper->Modified();
    traverseFlagActor_->Modified();
    traverseFlagActor_->GetProperty()->SetPointSize(5.0);
    GetRenderWindow()->Render();
    
}

void NeuroGLWidget::StartTraverseFromHere_Slot()
{
    if (mode_ != TRAVERSEMODE) return;
    if (pickedVertexID_.first < 0 || pickedVertexID_.second <0) return;
    emit SetTraversePos_Signal(pickedVertexID_.first, pickedVertexID_.second);
}

/*
khtao 2017.7.31 change

*/

void NeuroGLWidget::SmoothCurrentCurve_Slot()
{
    //smooth function
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            /*if (pickedIdNum == paramPack->activeTree->m_curveList.size()) {
            QMessageBox::information(this, "Cannot Delete Branch", "You will delete a whole tree. Please turn to tree mode to delete whole tree.");
            return;
            }*/
            //Smooth  curve id = pickedLineList_[0]
            auto &curTopo = paramPack->activeTree->m_Topo;
            auto it = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedLineList_[0]));
            auto curBranch = curTopo.subtree(it, curTopo.next_sibling(it));
            std::vector<size_t> curBranchID;
            std::for_each(curBranch.begin(), curBranch.end(), [&](LineID &arg){curBranchID.push_back(arg.id); });
            if (!paramPack->UserRecord.empty()){
                paramPack->activeTree->m_ChangedBranch.clear();
                paramPack->UserRecord.clear();
                paramPack->activeTree->m_DelBranch.clear();
            }
            std::vector<std::pair<size_t, Line5d> > ChangedBranch;
            for each(size_t i in curBranchID){
                ChangedBranch.push_back(std::make_pair(i, paramPack->activeTree->m_curveList[i]));
            }
            paramPack->activeTree->m_ChangedBranch.push_back(ChangedBranch);
            paramPack->UserRecord.push_back(Record::SmoothCurve);
            pUnDel->setEnabled(true);
            NGCorrectTrace mycorrector = CorrectTrace::New();
            mycorrector->SetParam(paramPack);
            mycorrector->SetLinePoints(curBranchID);
            time_t begin = time(0);
            if (!mycorrector->Update()){
                QMessageBox::warning(this, "error!", "maybe the connect point is not in the current range!");
            }
            time_t end = time(0);
            std::cout << "SmoothCurrentCurve used :" << end - begin << "s" << std::endl;
            //
            pickedLineList_.clear();
            while (!oldColor_.empty()) {
                delete[] oldColor_.back();
                oldColor_.back() = NULL;
                oldColor_.pop_back();
            }
            RebuildActiveTreeActor();
            ExitBranchMode();
            EnterBranchMode();
        }
        //update();
        GetRenderWindow()->Render();
    }
}

/*
khtao 2017.7.31 change
*/

void NeuroGLWidget::SaveImageAndSwcForTest_Slot(){
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){

            auto &curTopo = paramPack->activeTree->m_Topo;
            auto it = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedLineList_[0]));
            auto curBranch = curTopo.subtree(it, curTopo.next_sibling(it));
            std::vector<size_t> curBranchID;
            std::for_each(curBranch.begin(), curBranch.end(), [&](LineID &arg){curBranchID.push_back(arg.id); });
            QString saveName = QFileDialog::getSaveFileName(this, "save current image and current lines", paramPack->defaultDir);
            if (!saveName.isEmpty()){
                NGCorrectTrace mycorrector = CorrectTrace::New();
                mycorrector->SetParam(paramPack);
                mycorrector->SetLinePoints(curBranchID);
                mycorrector->SaveImageAndSwcForTest(saveName.toStdString());
                pickedLineList_.clear();
            }
        }
    }
}

void NeuroGLWidget::ToggleShowTraverseFlag_Slot()
{
    if (!traverseFlagActor_) return;
    if (traverseFlagActor_->GetVisibility() != 0) {
        traverseFlagActor_->SetVisibility(0);
        actorRenderer_->RemoveActor(traverseFlagActor_);
    }
    else{
        traverseFlagActor_->SetVisibility(1);
        actorRenderer_->AddActor(traverseFlagActor_);
    }
    //update();
    GetRenderWindow()->Render();
}

void NeuroGLWidget::InsertPointBefore_Slot()
{
    if (POINTMODE == mode_ &&pickedVertexID_.first != -1 && pickedVertexID_.second != -1){// canPickPoint_ && 
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        size_t idsNum = corNeuronTree->second->m_curveList[pickedVertexID_.first].size();
        if (idsNum <= 2llu || pickedVertexID_.second < 1 || pickedVertexID_.second > int(idsNum) - 1){
            printf("line is too short, or the selected dot is the end node.\n");
            return;
        }
        auto &allCurve = corNeuronTree->second->m_curveList;
        Vec3d pt1 = allCurve[pickedVertexID_.first][pickedVertexID_.second - 1].block(0,0,3,1);
        Vec3d pt2 = allCurve[pickedVertexID_.first][pickedVertexID_.second ].block(0, 0, 3, 1);
        Vec3d tmpPt = (pt1 + pt2) / 2.0;
        Vec5d interplPt = allCurve[pickedVertexID_.first][pickedVertexID_.second];
        interplPt(0) = tmpPt(0);
        interplPt(1) = tmpPt(1);
        interplPt(2) = tmpPt(2);
        auto iter = allCurve[pickedVertexID_.first].begin() + pickedVertexID_.second;
        allCurve[pickedVertexID_.first].insert(iter, interplPt);
        auto pickedVertexIdCp = pickedVertexID_;
        RebuildActiveTreeActor();
        ExitPointMode();
        EnterPointMode();
        pickedVertexID_ = pickedVertexIdCp;

        ChangePickedCellColor(pickedVertexID_.first, true);
        ShowVertex(pickedVertexID_.first, true);
        ChangePickedVertColor(pickedVertexID_.second, true);
    }
}

void NeuroGLWidget::RemovePoint_Slot()
{
    if (POINTMODE == mode_ &&pickedVertexID_.first != -1 && pickedVertexID_.second != -1){// canPickPoint_ && 
        auto corNeuronTree = std::find_if(neuronActorList_.begin(), neuronActorList_.end(), [&](const std::pair< vtkSmartPointer<vtkActor>, NeuronTreePointer>& arg){ return arg.first.GetPointer() == activeTreeActor_; });
        size_t idsNum = corNeuronTree->second->m_curveList[pickedVertexID_.first].size();
        if (idsNum <= 2llu || pickedVertexID_.second < 1 || pickedVertexID_.second > int(idsNum) - 1){
            printf("line is too short, or the selected dot is the end node.\n");
            return;
        }
        auto &allCurve = corNeuronTree->second->m_curveList;
        auto iter = allCurve[pickedVertexID_.first].begin() + pickedVertexID_.second;
        allCurve[pickedVertexID_.first].erase(iter);
        auto pickedVertexIdCp = pickedVertexID_;
        RebuildActiveTreeActor();
        ExitPointMode();
        EnterPointMode();
        pickedVertexID_ = pickedVertexIdCp;

        ChangePickedCellColor(pickedVertexID_.first, true);
        ShowVertex(pickedVertexID_.first, true);
        ChangePickedVertColor(pickedVertexID_.second, true);
    }
}

bool NeuroGLWidget::ToggleMovePoint_Slot(bool arg)
{
    if (mode_ != POINTMODE) return false;
    if (arg) {
        if (pickedVertexID_.first < 0 || pickedVertexID_.second < 0) {
            QMessageBox::warning(this, "no selected point", "Please select point first.");
            pStartMovePt_->setChecked(false);
            return false;
        }
        if (!interactorObjectStyle_) {
            interactorObjectStyle_ = vtkSmartPointer<NGInteractorObjectStyle>::New();
        }
        auto &curTopo = paramPack->activeTree->m_Topo;
        auto &allLine = paramPack->activeTree->m_curveList;
        auto &pickPos = allLine[pickedVertexID_.first][pickedVertexID_.second];
        auto curIt = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedVertexID_.first));
        //if is tree root curve
        bool flag = false;
        std::vector<size_t> leve1BranchIDList;
        interactorObjectStyle_->isTreeRoot_ = false;
        interactorObjectStyle_->leve1BranchIDList_.clear();
        if (pickedVertexID_.second == 0) {
            int level = 1;
           
            for (auto iter = curTopo.begin() + 1; iter != curTopo.end(); ++iter) {
                if (curTopo.depth(iter) == level) {
                    leve1BranchIDList.push_back(iter->id);
                }
            }
            for (auto &it : leve1BranchIDList) {
                if (it == pickedVertexID_.first) {
                    flag = true;
                    break;
                }
            }
            
            if (flag) {//move tree root
                interactorObjectStyle_->isTreeRoot_ = true;
                interactorObjectStyle_->leve1BranchIDList_.swap(leve1BranchIDList);
            }
            else{
                QMessageBox::warning(this, "Please do not select the head node", "Please do not select the head node");
                pStartMovePt_->setChecked(false);
                return false;
            }
        }
        else if (curTopo.number_of_children(curIt) > 0) {
            int num = curTopo.number_of_children(curIt);
            for (int i = 0; i < num; ++i) {
                size_t curChildId = curTopo.child(curIt, i)->id;
                auto &curLineHead = allLine[curChildId][0];
                if ((curLineHead.block(0,0,3,1) - pickPos.block(0,0,3,1)).norm() < 1.0  ) {
                    QMessageBox::warning(this, "Cannot move point", "This point is the head of child curves. Please Change the head of child first.");
                    pStartMovePt_->setChecked(false);
                    return false;
                }
            }
        }
        
        interactorObjectStyle_->pickedVertexID_ = pickedVertexID_;
        interactorObjectStyle_->SetParam(paramPack);
        interactorObjectStyle_->SetWidget(this);
        interactorObjectStyle_->dstActor_ = activeTreeActor_;
        actorRenderer_->AddActor(interactorObjectStyle_->MoveActor);
        GetInteractor()->SetInteractorStyle(interactorObjectStyle_);
    }
    else{
        GetInteractor()->SetInteractorStyle(interactorStyle_);
        if(interactorObjectStyle_) interactorObjectStyle_->dstActor_ = NULL;
        if(interactorObjectStyle_) actorRenderer_->RemoveActor(interactorObjectStyle_->MoveActor);
    }
    isMovePointMode_ = arg;
    //update();
    GetRenderWindow()->Render();
    return true;
}

void NeuroGLWidget::SetMovePointPos(double x, double y, double z)
{
    if (isMovePointMode_) {
        movePointActor_->SetPosition(x, y, z);
        //update();
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::ToggleReconParent_Slot(bool arg)
{
    if (pickedVertexID_.second < 1 && arg) {
        QMessageBox::warning(this, "Cannot Reconnect", "Please do not select the head node.");
        pReconnectParent_->setChecked(false);
        return;
    }
    isReconnectParent_ = arg;
}

void NeuroGLWidget::BuildSomaActor(const Cell& tmpSoma, vtkSmartPointer<vtkFollower> actor)
{
    VTK_CREATE(vtkRegularPolygonSource, polygonSource);
    polygonSource->SetCenter(0, 0, 0);
    polygonSource->SetNumberOfSides(3);
    //polygonSource->SetRadius();
    VTK_CREATE(vtkPolyDataMapper, mapper);
    mapper->SetInputConnection(polygonSource->GetOutputPort());
    actor->GetProperty()->SetRepresentationToSurface();
    actor->SetMapper(mapper);
    actor->SetPosition(tmpSoma.x * paramPack->xRes_, tmpSoma.y * paramPack->yRes_, tmpSoma.z * paramPack->zRes_);
    actor->GetProperty()->SetLineWidth(LINEWIDTH_NOR);
    actor->GetProperty()->SetColor(255, 0, 0);
    actor->Modified();
}

void NeuroGLWidget::UpdateSomaList()
{
    if (!paramPack->SomaList) return;
    if (mode_ == SOMAMODE) cellPicker_->InitializePickList();
    ClearActorCollection(somaActorCollection_);
    NG_DYNAMIC_CAST(Soma, manualSoma_, paramPack->SomaList);
    for (auto i = 0; i < manualSoma_->size(); ++i) {
        VTK_NEW(vtkFollower, choosedSomaActor_);
        BuildSomaActor(manualSoma_->GetCell(i), choosedSomaActor_);
        somaActorCollection_.push_back(choosedSomaActor_);
        actorRenderer_->AddActor(choosedSomaActor_);
        if (mode_ == SOMAMODE) cellPicker_->AddPickList(choosedSomaActor_);
    }
}

void NeuroGLWidget::RefreshTree_Slot()
{
    if (mode_ == TREEMODE) {
        if (!pickedTreeList_.empty()) {
            NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
            int changNum = 0;
            int id = -1;
            for (auto &it:pickedTreeList_) {
                changNum = KPTREE::RefreshTree(*(it.second));
                if (changNum > 0){
                    printf("tree %d has been refresh for %d places.\n", ++id, changNum);
                    if (it.second.get() == paramPack->activeTree.get()) 
                        RebuildNeuronTreeActor(it, true);
                    else RebuildNeuronTreeActor(it);
                }
            }
        }
    }
}

void NeuroGLWidget::BuildCaliber_Slot()
{
    if (mode_ == BRANCHMODE) {
        if (!pickedLineList_.empty()) {
            auto &curTopo = paramPack->activeTree->m_Topo;
            auto it = std::find(curTopo.begin() + 1, curTopo.end(), LineID(pickedLineList_[0]));
            auto curBranch = curTopo.subtree(it, curTopo.next_sibling(it));
            std::vector<size_t> curBranchID;
            std::for_each(curBranch.begin(), curBranch.end(), [&](LineID &arg){curBranchID.push_back(arg.id); });
            NGCaliberBuilder caliberBuilder = CaliberBuilder::New();
            caliberBuilder->SetParam(paramPack);
            caliberBuilder->SetLinePoints(curBranchID);
            //time_t begin = time(0);
            auto res = caliberBuilder->Update();
            if (!res->success()){
                QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
                return;
            }
            //time_t end = time(0);
            if (!localCurveCaliberActor_)
                VTK_NEW(vtkActor, localCurveCaliberActor_);
            else actorRenderer_->RemoveActor(localCurveCaliberActor_);
            BuildLocalCaliberActor(*(caliberBuilder->GetCaliberPointSet()), localCurveCaliberActor_);
            actorRenderer_->AddActor(localCurveCaliberActor_);
            //update();
            GetRenderWindow()->Render();
        }
    }
}

void NeuroGLWidget::BuildLocalCaliberActor(const CellVec3d &localCaliber, vtkSmartPointer<vtkActor> actor)
{
    VTK_CREATE(vtkPolyData, vertsPoly);
    VTK_CREATE(vtkPolyDataMapper, vertsMapper);
    VTK_CREATE(vtkCellArray, vertsCells);
    VTK_CREATE(vtkPoints, vertsPoints);
    VTK_CREATE(vtkUnsignedCharArray, vertScalars);
    vertScalars->SetNumberOfComponents(3);
    int pointNum = 0;
    for (auto &it : localCaliber) 
        pointNum += int(it.size());
    vertsCells->InsertNextCell(pointNum);
    int ind = -1;
    //uchar color[3] = { 85, 20, 200 };
    uchar *color = traceColor;
    int cid = 0;
    for (auto &it : localCaliber) {
        vtkIdType sz = vtkIdType(it.size());
        for (vtkIdType i = 0; i < sz; ++i){
            const Vec3d &kPoint = it[i];
            vertsPoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
            vertsCells->InsertCellPoint(++ind);
            //vertScalars->InsertNextTupleValue(color);
            vertScalars->InsertNextTypedTuple(color + 3 * (cid % 11)); 
        }
        ++cid;
    }
    vertsPoly->SetPoints(vertsPoints);
    vertsPoly->SetVerts(vertsCells);
    vertsPoly->GetPointData()->SetScalars(vertScalars);
    vertsMapper->SetInputData(vertsPoly);
    vertsMapper->SetScalarModeToUsePointData();
    actor->SetMapper(vertsMapper);
    actor->GetProperty()->SetPointSize(3.0);
    actor->Modified();
}

void NeuroGLWidget::Snapshot(QString& path, int vx)//int vy, 
{
    //GetRenderWindow()->SetOffScreenRendering(1);
    /*int oldsz[2] = { 0, 0 };
    oldSz = this->size();
    int *s = GetRenderWindow()->GetSize();
    oldsz[0] = s[0]; oldsz[1] = s[1];*/
    //auto oldqsz = this->size();
    
    marker_->SetEnabled(0);
    outlineActor_->SetVisibility(0);
    imageBoxActor_->SetVisibility(0);
    //QSize restoreqsz = this->size();
    //QSize newqsz(vx, vy);
    GetRenderWindow()->SetOffScreenRendering(1);
    
    //GetRenderWindow()->SetSize(vx, vy);
    //GetActiveRenderer()->ResetCamera();
    //GetRenderWindow()->Render();
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(GetRenderWindow());
    windowToImageFilter->SetMagnification(vx);
    //windowToImageFilter->SetScale(vx, vy);
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(path.toStdString().c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
    GetRenderWindow()->SetOffScreenRendering(0);
    //GetRenderWindow()->SetSize(oldsz);
    //GetActiveRenderer()->ResetCamera();
    //GetRenderWindow()->Render();
    
    marker_->SetEnabled(1);
    outlineActor_->SetVisibility(1);
    imageBoxActor_->SetVisibility(1);
    //update();
    GetRenderWindow()->Render();
    
    //this->repaint();
    //update();
    //QResizeEvent event(restoreqsz, newqsz );
    //this->resizeEvent(&event);
    //QApplication::sendEvent(this, &event);
    //GetRenderWindow()->SetSize(this->width(), this->height());
}

void NeuroGLWidget::SetBranchAsDenrite_Slot()
{
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            if (pickedIdNum == paramPack->activeTree->m_curveList.size()) {
                QMessageBox::information(this, "Cannot Delete Branch", "You will delete a whole tree. Please turn to tree mode to delete whole tree.");
                return;
            }
            //get type
            QStringList items; items << "Common Dendrite" << "Apical Dendrite" << "basel Dendrite";
            bool ok;
            QString item = QInputDialog::getItem(this, tr("QInputDialog::getItem()"), tr("Dendrite type:"), items, 0, false, &ok);
            if (!ok || item.isEmpty()) return;
            //
            auto &topo = paramPack->activeTree->m_Topo;
            auto topIt = std::find_if(topo.begin(), topo.end(), [&](const LineID &arg){ return arg.id == pickedLineList_[0]; });
            if (topIt == topo.end()) {
                QMessageBox::warning(this, "Branch set error", "Cannot pick the branch, please debug.");
                return;
            }
            auto subtree = paramPack->activeTree->m_Topo.subtree(topIt, topo.next_sibling(topIt));
            std::vector<size_t> subtreeID;
            for (auto it = subtree.begin(); it != subtree.end(); ++it)
                subtreeID.push_back(it->id);
            NEURITETYPE newType;
            if (item.compare("Common Dendrite") == 0)
                newType = NEURITETYPE::DENDRITE;
            else newType = NEURITETYPE::APICAL;
            for (auto &it : subtreeID) {
                auto node = std::find_if(topo.begin(), topo.end(), [&](const LineID &arg){ return arg.id == it; });
                node->type = newType;
            }
        }
    }
}

void NeuroGLWidget::SetBranchAsAxon_Slot()
{
    if ((BRANCHMODE == mode_)){
        size_t pickedIdNum = pickedLineList_.size();
        if (pickedIdNum != 0){
            if (pickedIdNum == paramPack->activeTree->m_curveList.size()) {
                QMessageBox::information(this, "Cannot Delete Branch", "You will delete a whole tree. Please turn to tree mode to delete whole tree.");
                return;
            }
            auto &topo = paramPack->activeTree->m_Topo;
            auto topIt = std::find_if(topo.begin(), topo.end(), [&](const LineID &arg){ return arg.id == pickedLineList_[0]; });
            if (topIt == topo.end()) {
                QMessageBox::warning(this, "Branch set error", "Cannot pick the branch, please debug.");
                return;
            }
            auto subtree = paramPack->activeTree->m_Topo.subtree(topIt, topo.next_sibling(topIt));
            std::vector<size_t> subtreeID;
            for (auto it = subtree.begin(); it != subtree.end(); ++it)
                subtreeID.push_back(it->id);
            for (auto &it : subtreeID) {
                auto node = std::find_if(topo.begin(), topo.end(), [&](const LineID &arg){ return arg.id == it; });
                node->type = NEURITETYPE::AXON;
            }
        }
    }
}

void NeuroGLWidget::HideDendrite_Slot(bool arg)
{
    HideDenAxonBase(arg, NEURITETYPE::DENDRITE);
}

void NeuroGLWidget::HideAxon_Slot(bool arg)
{
    HideDenAxonBase(arg, NEURITETYPE::AXON);
}

void NeuroGLWidget::HideDenAxonBase(bool val1, NEURITETYPE val2)
{
    if (paramPack->activeTree) {
        //VTK_DEFINE(vtkPolyData, treePoly) = vtkPolyData::SafeDownCast(arg.first->GetMapper()->GetInput());
        VTK_CREATE(vtkPolyData, treePoly);
        VTK_CREATE(vtkPoints, treePoints);
        VTK_CREATE(vtkCellArray, treeCells);
        VTK_CREATE(vtkPolyDataMapper, treeMapper);
        VTK_CREATE(vtkUnsignedCharArray, treeColorArray);
        treeColorArray->SetNumberOfComponents(3);
        //save the pick line info
        if (!pickedLineList_.empty()) pickedLineList_.clear();
        pickedVertexID_.first = 0;
        pickedVertexID_.second = -1;
        const std::vector<Line5d>& curTree = paramPack->activeTree->m_curveList;
        const auto &topo = paramPack->activeTree->m_Topo;
        auto  topoHead = topo.begin() + 1;
        VectorVec5d tmpLine;
        //tmpLine.push_back(NGUtility::MakeVec5d(0, 0, 0));
        //tmpLine.push_back(NGUtility::MakeVec5d(1, 1, 1));
        for (std::vector<VectorVec5d>::size_type j = 0; j < curTree.size(); ++j){
            VTK_CREATE(vtkPolyLine, jPolyLine);
            auto topoIt = std::find(topo.begin(), topo.end(), LineID(j));
            if (val1 && val2 == NEURITETYPE::DENDRITE && (topoIt->type == NEURITETYPE::DENDRITE || topoIt->type == NEURITETYPE::APICAL)) {
                tmpLine.clear();
                tmpLine.push_back(curTree[topoHead->id][0]); tmpLine.push_back(curTree[topoHead->id][0]);
                jPolyLine->GetPointIds()->SetNumberOfIds(int(tmpLine.size()));
                for (VectorVec5d::size_type k = 0; k < tmpLine.size(); ++k){
                    const Vec5d &kPoint = tmpLine[k];
                    vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
                    jPolyLine->GetPointIds()->SetId(k, xId);
                }
            }
            else if (val1 && val2 == NEURITETYPE::AXON &&topoIt->type == val2){
                tmpLine.clear();
                tmpLine.push_back(curTree[topoHead->id][0]); tmpLine.push_back(curTree[topoHead->id][0]);
                jPolyLine->GetPointIds()->SetNumberOfIds(int(tmpLine.size()));
                for (VectorVec5d::size_type k = 0; k < tmpLine.size(); ++k){
                    const Vec5d &kPoint = tmpLine[k];
                    vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
                    jPolyLine->GetPointIds()->SetId(k, xId);
                }
            }
            else{
                jPolyLine->GetPointIds()->SetNumberOfIds(int(curTree[j].size()));
                for (VectorVec5d::size_type k = 0; k < curTree[j].size(); ++k){
                    const Vec5d &kPoint = curTree[j][k];
                    vtkIdType xId = treePoints->InsertNextPoint(kPoint[0] * paramPack->xRes_, kPoint[1] * paramPack->yRes_, kPoint[2] * paramPack->zRes_);
                    jPolyLine->GetPointIds()->SetId(k, xId);
                }
            }
            treeCells->InsertNextCell(jPolyLine);
            treeColorArray->InsertNextTypedTuple(NORMALCELLCOLOR);//InsertNextTypedTuple
        }
        treeColorArray->Modified();
        treePoly->SetPoints(treePoints);
        treePoly->SetLines(treeCells);
        treePoly->GetCellData()->SetScalars(treeColorArray);
        treePoly->BuildCells();
        treePoly->Modified();
        treeMapper->SetInputData(treePoly);
        treeMapper->Update();
        treeMapper->Modified();
        activeTreeActor_->SetMapper(treeMapper);
        GetRenderWindow()->Render();
    }
}

void NeuroGLWidget::HideCaliber_Slot(bool arg)
{
    if (localCurveCaliberActor_) {
        localCurveCaliberActor_->SetVisibility(!arg);
    }
}

void NeuroGLWidget::renderCaliber_Slot()
{
	if (!localCurveCaliberActor_)
		VTK_NEW(vtkActor, localCurveCaliberActor_);
	else 
		actorRenderer_->RemoveActor(localCurveCaliberActor_);

	BuildLocalCaliberActor(*(paramPack->caliberBulidData), localCurveCaliberActor_);
	actorRenderer_->AddActor(localCurveCaliberActor_);
	//update();
	GetRenderWindow()->Render();
}

void NeuroGLWidget::AllCaliberBulider()
{
	if (mode_ == BRANCHMODE) {
		if (!paramPack->activeTree)
		{
			std::cout << "ERROR! No active Tree!" << std::endl;
			return;
		}

		//get all branch of activeTree
		auto &curTopo = paramPack->activeTree->m_Topo;
		auto &allCurve = paramPack->activeTree->m_curveList;
		std::vector<size_t> curBranchId;
		std::for_each(curTopo.begin() + 1, curTopo.end(), [&](const LineID &arg){curBranchId.push_back(arg.id); });
		for (auto &it :curBranchId) {
			if (allCurve[it].size() < 4) {
				it = 100000000;
			}
		}
		curBranchId.erase(std::remove(curBranchId.begin(), curBranchId.end(), 100000000), curBranchId.end());
		NGCaliberBuilder caliberBuilder = CaliberBuilder::New();
		caliberBuilder->SetParam(paramPack);
		caliberBuilder->SetLinePoints(curBranchId);
		auto res = caliberBuilder->Update();
		if (!res->success()){
			QMessageBox::warning(this, QString::fromStdString(res->ErrorClassName()), QString::fromStdString(res->ErrorInfo()));
			return;
		}

		if (!localCurveCaliberActor_)
			VTK_NEW(vtkActor, localCurveCaliberActor_);
		else 
			actorRenderer_->RemoveActor(localCurveCaliberActor_);
		BuildLocalCaliberActor(*(caliberBuilder->GetCaliberPointSet()), localCurveCaliberActor_);
		actorRenderer_->AddActor(localCurveCaliberActor_);

		GetRenderWindow()->Render();
		std::cout << "Build all branch completed!" << std::endl;
	}
}
