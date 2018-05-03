/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include "nginteractorstyletrackballcamera.h"
#include <vtkRenderWindowInteractor.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkProp3DCollection.h>
#include <vtkAssemblyNode.h>
#include <vtkAssemblyPath.h>
#include <vtkCallbackCommand.h>
#include <vtkCellPicker.h>
#include <vtkWorldPointPicker.h>
#include <vtkCommand.h>
#include <vtkMath.h>
#include <vtkObjectFactory.h>
#include <vtkOutlineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
//#include <vtkTextProperty.h>
#include <vtkEventForwarderCommand.h>
#include <vtkTDxInteractorStyleCamera.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include "neuroglwidget.h"

NGInteractorStyleTrackballCamera::NGInteractorStyleTrackballCamera()
{
    //VTK_MODULE_INIT(vtkRenderingOpenGL2); 
    //VTK_MODULE_INIT(vtkInteractionStyle); 
    //VTK_MODULE_INIT(vtkRenderingFreeType);
}

void NGInteractorStyleTrackballCamera::OnChar()
{
    vtkRenderWindowInteractor *rwi = this->Interactor;
    switch (rwi->GetKeyCode())
    {
    //case 'c':
    //case 'C':{
        /* if (glbox->IsForbidden()) break;
         if (!glbox->pChooseBranchMode->isChecked()) {
         glbox->pChooseNodeMode->setChecked(false);
         glbox->pChooseSomaMode->setChecked(false);
         glbox->pChooseTreeMode->setChecked(false);
         glbox->pDrawLineMode->setChecked(false);
         glbox->pChooseBranchMode->setChecked(true);
         }
         else glbox->pChooseBranchMode->setChecked(false);*/
       // break;
    //}
    /*case 'd':
    case 'D':{
        if (glbox->IsForbidden()) break;
        glbox->pDeleteBranch->trigger();
        break;
    }*/
    /*case 'e':
    case 'E':
        if (glbox->IsForbidden()) break;
        glbox->ToggleEnableCurveTrace();
        break;*/
    /*case 'z':
    case 'Z':{
        if (glbox->IsForbidden()) break;
        if (!glbox->pChooseNodeMode->isChecked()) {
            glbox->pChooseSomaMode->setChecked(false);
            glbox->pChooseTreeMode->setChecked(false);
            glbox->pDrawLineMode->setChecked(false);
            glbox->pChooseBranchMode->setChecked(false);
            glbox->pChooseNodeMode->setChecked(true);
        }
        else glbox->pChooseNodeMode->setChecked(false);
        break;
    }*/
    /*case 'x':
    case 'X':{
        if (glbox->IsForbidden()) break;
        glbox->pCutLine->trigger();
        break;
    }*/
    /*case 'v':
    case 'V':{
        if (glbox->IsForbidden()) break;
        if (!glbox->pDrawLineMode->isChecked()) {
            glbox->pChooseSomaMode->setChecked(false);
            glbox->pChooseTreeMode->setChecked(false);
            glbox->pChooseBranchMode->setChecked(false);
            glbox->pChooseNodeMode->setChecked(false);
            glbox->pDrawLineMode->setChecked(true);
        }
        else glbox->pDrawLineMode->setChecked(false);
        break;
    }*/
    case 'm' :
    case 'M' :
        if (this->AnimState == VTKIS_ANIM_OFF)
        {
            this->StartAnimate();
        }
        else
        {
            this->StopAnimate();
        }
        break;

    case 'o':
    case 'O':
        break;

    case 'f' :
    case 'F' :
    {
        if(this->CurrentRenderer!=0)
        {
            this->AnimState = VTKIS_ANIM_ON;
            vtkAssemblyPath *path = NULL;
            this->FindPokedRenderer(rwi->GetEventPosition()[0],
                                    rwi->GetEventPosition()[1]);
            rwi->GetPicker()->Pick(rwi->GetEventPosition()[0],
                                   rwi->GetEventPosition()[1],
                                   0.0,
                                   this->CurrentRenderer);
            vtkAbstractPropPicker *picker;
            if ((picker=vtkAbstractPropPicker::SafeDownCast(rwi->GetPicker())))
            {
                path = picker->GetPath();
            }
            if (path != NULL)
            {
                rwi->FlyTo(this->CurrentRenderer, picker->GetPickPosition());
            }
            this->AnimState = VTKIS_ANIM_OFF;
        }
        else
        {
            vtkWarningMacro(<<"no current renderer on the interactor style.");
        }
    }
        break;

    case 'r' :
    case 'R' :
        this->FindPokedRenderer(rwi->GetEventPosition()[0],
                                rwi->GetEventPosition()[1]);
        if(this->CurrentRenderer!=0)
        {
            this->CurrentRenderer->ResetCamera();
        }
        else
        {
            vtkWarningMacro(<<"no current renderer on the interactor style.");
        }
        rwi->Render();
        break;

    case 'w' :
    case 'W' :
    {
        vtkActorCollection *ac;
        vtkActor *anActor, *aPart;
        vtkAssemblyPath *path;
        this->FindPokedRenderer(rwi->GetEventPosition()[0],
                                rwi->GetEventPosition()[1]);
        if(this->CurrentRenderer!=0)
        {
            ac = this->CurrentRenderer->GetActors();
            vtkCollectionSimpleIterator ait;
            for (ac->InitTraversal(ait); (anActor = ac->GetNextActor(ait)); )
            {
                for (anActor->InitPathTraversal(); (path=anActor->GetNextPath()); )
                {
                    aPart=static_cast<vtkActor *>(path->GetLastNode()->GetViewProp());
                    aPart->GetProperty()->SetRepresentationToWireframe();
                }
            }
        }
        else
        {
            vtkWarningMacro(<<"no current renderer on the interactor style.");
        }
        rwi->Render();
    }
        break;

    case 's' :
    case 'S' :
    {
        vtkActorCollection *ac;
        vtkActor *anActor, *aPart;
        vtkAssemblyPath *path;
        this->FindPokedRenderer(rwi->GetEventPosition()[0],
                                rwi->GetEventPosition()[1]);
        if(this->CurrentRenderer!=0)
        {
            ac = this->CurrentRenderer->GetActors();
            vtkCollectionSimpleIterator ait;
            for (ac->InitTraversal(ait); (anActor = ac->GetNextActor(ait)); )
            {
                for (anActor->InitPathTraversal(); (path=anActor->GetNextPath()); )
                {
                    aPart=static_cast<vtkActor *>(path->GetLastNode()->GetViewProp());
                    aPart->GetProperty()->SetRepresentationToSurface();
                }
            }
        }
        else
        {
            vtkWarningMacro(<<"no current renderer on the interactor style.");
        }
        rwi->Render();
    }
        break;

    case '3' :
        if (rwi->GetRenderWindow()->GetStereoRender())
        {
            rwi->GetRenderWindow()->StereoRenderOff();
        }
        else
        {
            rwi->GetRenderWindow()->StereoRenderOn();
        }
        rwi->Render();
        break;

    case '+':
        {
            break;
        }

    case '-':
        {
            break;
        }

    case 'p' :
    case 'P' :
        if(this->CurrentRenderer!=0)
        {
            if (this->State == VTKIS_NONE)
            {
                //vtkAssemblyPath *path = NULL;
                int *eventPos = rwi->GetEventPosition();
                this->FindPokedRenderer(eventPos[0], eventPos[1]);
                rwi->StartPickCallback();
                rwi->EndPickCallback();
            }
        }
        else
        {
            vtkWarningMacro(<<"no current renderer on the interactor style.");
        }
        break;
    }
}


NGInteractorStyleTrackballCamera *NGInteractorStyleTrackballCamera::New()
{
    return new NGInteractorStyleTrackballCamera();
}

NGInteractorStyleTrackballCamera::~NGInteractorStyleTrackballCamera()
{
}

//void NGInteractorStyleTrackballCamera::InsertLine( double x, double y, double z)
//{
//    if(!glbox) return;
//    glbox->InsertLine(x,y,z);
//}

void NGInteractorStyleTrackballCamera::OnLeftButtonDown()
{
    if (glbox->mode_ != NeuroGLWidget::NORMALMODE){
        isLeftButtonDown_ = 1;//glbox->MouseClick();
    }
    if (!glbox->Is2DView() )  vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void NGInteractorStyleTrackballCamera::SetWidget( NeuroGLWidget* arg )
{
    glbox = arg;
}

void NGInteractorStyleTrackballCamera::OnRightButtonDown()
{
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
}

void NGInteractorStyleTrackballCamera::OnRightButtonUp()
{
   vtkInteractorStyleTrackballCamera::OnRightButtonUp();
}

void NGInteractorStyleTrackballCamera::OnMiddleButtonDown()
{
        vtkInteractorStyleTrackballCamera::OnMiddleButtonDown();
        glbox->SentChangeMoveCursorSignal_Slot(true);
}

void NGInteractorStyleTrackballCamera::OnMiddleButtonUp()
{
        vtkInteractorStyleTrackballCamera::OnMiddleButtonUp();
        //isMiddleButtonDown_ = 0;
        glbox->SentChangeMoveCursorSignal_Slot(false);
}

void NGInteractorStyleTrackballCamera::OnMouseWheelForward()
{
    if (glbox->Is2DView()){
        glbox->On2DMouseWheel(true);
    }
    else{
        vtkInteractorStyleTrackballCamera::OnMouseWheelForward();
    }
}

void NGInteractorStyleTrackballCamera::OnMouseWheelBackward()
{
    if (glbox->Is2DView()){
        glbox->On2DMouseWheel(false);
    }
    //else if (glbox->isBoxSelectMode()){//sulei

    //}
    else{
        vtkInteractorStyleTrackballCamera::OnMouseWheelBackward();
    }
}

void NGInteractorStyleTrackballCamera::OnMouseMove()
{
    if (glbox->mode_ != NeuroGLWidget::NORMALMODE){
        if(isLeftButtonDown_ == 1) isLeftButtonDown_ = 0;
        numberOfClicks = 0;
    }
   vtkInteractorStyleTrackballCamera::OnMouseMove();// if (!glbox->Is2DView())  
}

//void NGInteractorStyleTrackballCamera::OnKeyPress()
//{
//    char *keySym = GetInteractor()->GetKeySym();
//    std::string keyString(keySym);
//    glbox->On2DKeyPress(keyString);
//}

void NGInteractorStyleTrackballCamera::Rotate()
{
    tmpActorRenderer = glbox->GetActorRenderer();
    tmpActiveRenderer = glbox->GetActiveRenderer();
    if (tmpActorRenderer == NULL)
    {
        return;
    }
    if (tmpActiveRenderer == NULL) return;

    vtkRenderWindowInteractor *rwi = this->Interactor;

    int dx = rwi->GetEventPosition()[0] - rwi->GetLastEventPosition()[0];
    int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];

    int *size = tmpActiveRenderer->GetRenderWindow()->GetSize();

    double delta_elevation = -20.0 / size[1];
    double delta_azimuth = -20.0 / size[0];

    double rxf = dx * delta_azimuth * this->MotionFactor;
    double ryf = dy * delta_elevation * this->MotionFactor;

    vtkCamera *camera = tmpActiveRenderer->GetActiveCamera();
    camera->Azimuth(rxf);
    camera->Elevation(ryf);
    camera->OrthogonalizeViewUp();
    //cout << camera->GetClippingRange()[0] << " " << camera->GetClippingRange()[1] << endl;
    if (this->AutoAdjustCameraClippingRange)
    {
        tmpActorRenderer->ResetCameraClippingRange();
        //camera->SetClippingRange(0.0001, 10000.00001);
    }
    //cout << camera->GetClippingRange()[0] << " " << camera->GetClippingRange()[1] << endl;
    if (rwi->GetLightFollowCamera())
    {
        tmpActiveRenderer->UpdateLightsGeometryToFollowCamera();
    }

    rwi->Render();
}

void NGInteractorStyleTrackballCamera::OnLeftButtonUp()
{
    if (glbox->mode_ != NeuroGLWidget::NORMALMODE){
        if (glbox->mode_ == NeuroGLWidget::DRAWLINEMODE ){
            ++numberOfClicks;
            glbox->MouseClick();
        }
        else if (isLeftButtonDown_ == 1){
            glbox->MouseClick();
        }
        isLeftButtonDown_ = 0;
    }
    vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
}

NGInteractorObjectStyle::NGInteractorObjectStyle()
{
    this->isMove_ = false;
    //this->PointPicker = vtkSmartPointer<vtkPointPicker>::New();

    // Setup ghost glyph
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();
    points->InsertNextPoint(0, 0, 0);
    this->MovePolyData = vtkSmartPointer<vtkPolyData>::New();
    this->MovePolyData->SetPoints(points);
    this->MoveGlyphFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
    this->MoveGlyphFilter->SetInputData(this->MovePolyData);
    this->MoveGlyphFilter->Update();

    this->MoveMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    this->MoveMapper->SetInputConnection(this->MoveGlyphFilter->GetOutputPort());

    this->MoveActor = vtkSmartPointer<vtkActor>::New();
    this->MoveActor->SetMapper(this->MoveMapper);
    //this->MoveActor->VisibilityOff();
    this->MoveActor->GetProperty()->SetPointSize(10);
    this->MoveActor->GetProperty()->SetColor(1, 0, 0);
}

void NGInteractorObjectStyle::OnMouseMove()
{
    if (!this->isMove_)
        return;
    glbox->update();
    vtkInteractorStyleTrackballActor::OnMouseMove();
}

void NGInteractorObjectStyle::OnMiddleButtonUp()
{
    this->EndPan();

    this->isMove_ = false;
    //this->MoveActor->VisibilityOff();
    //this->GetCurrentRenderer()->RemoveActor(MoveActor);

    vtkPolyData* Data = vtkPolyData::SafeDownCast(dstActor_->GetMapper()->GetInput());
   
    //auto &allCurve = param_->activeTree->m_curveList;
    std::vector<int> changeList;
    if (isTreeRoot_) {
        for (auto &it : leve1BranchIDList_) {
            changeList.push_back(int(it));
        }
    }
    else{
        changeList.push_back(pickedVertexID_.first);
    }
    if (glbox->GetLayerDisplay().empty() || glbox->GetLayerDisplay()[0] == 0) {//not level display
        for (auto &it : changeList) {
            vtkIdList *ptIds = vtkIdList::New();
            Data->GetCellPoints(it, ptIds);
            vtkIdType gid = ptIds->GetId(pickedVertexID_.second);
            ptIds->Delete();
            Data->GetPoints()->SetPoint(gid, this->MoveActor->GetPosition());
            Vec5d &curPos = param_->activeTree->m_curveList[it][pickedVertexID_.second];
            if (vtkIdType(param_->activeTree->m_curveList[it].size()) > pickedVertexID_.second + 1) {//for tree root
                Vec5d &curPos2 = param_->activeTree->m_curveList[it][pickedVertexID_.second + 1];
                //Vec5d &curPos3 = param_->activeTree->m_curveList[it][pickedVertexID_.second + 2];
                if ((curPos2.block(0, 0, 3, 1) - curPos.block(0, 0, 3, 1)).norm() < 1.0) {
                    curPos2(0) = this->MoveActor->GetPosition()[0] / param_->xRes_;
                    curPos2(1) = this->MoveActor->GetPosition()[1] / param_->yRes_;
                    curPos2(2) = this->MoveActor->GetPosition()[2] / param_->zRes_;
                    Data->GetPoints()->SetPoint(gid + 1, this->MoveActor->GetPosition());
                }
            }
            curPos(0) = this->MoveActor->GetPosition()[0] / param_->xRes_;
            curPos(1) = this->MoveActor->GetPosition()[1] / param_->yRes_;
            curPos(2) = this->MoveActor->GetPosition()[2] / param_->zRes_;
        }
    }
    else{//just search the nearest node
        vtkIdList *ptIds = vtkIdList::New();
        Data->GetCellPoints(pickedVertexID_.first, ptIds);
        vtkIdType gid = ptIds->GetId(pickedVertexID_.second);
        ptIds->Delete();
        Data->GetPoints()->SetPoint(gid, this->MoveActor->GetPosition());
        Vec5d &curPos = param_->activeTree->m_curveList[pickedVertexID_.first][pickedVertexID_.second];
        curPos(0) = this->MoveActor->GetPosition()[0] / param_->xRes_;
        curPos(1) = this->MoveActor->GetPosition()[1] / param_->yRes_;
        curPos(2) = this->MoveActor->GetPosition()[2] / param_->zRes_;
    }
    
    Data->GetPoints()->Modified();
    Data->Modified();
    //vtkInteractorStyleTrackballActor::OnMiddleButtonUp();
    this->GetCurrentRenderer()->Render();
    this->GetCurrentRenderer()->GetRenderWindow()->Render();
    //glbox->RebuildActiveTreeActor();
    glbox->update();
}

void NGInteractorObjectStyle::OnMiddleButtonDown()
{
    if (!isTreeRoot_ && pickedVertexID_.second == 0 ) {
        return;
    }
    // Get the selected point
    int x = this->Interactor->GetEventPosition()[0];
    int y = this->Interactor->GetEventPosition()[1];
    this->FindPokedRenderer(x, y);

    this->StartPan();
    this->MoveActor->VisibilityOn();
    //this->GetCurrentRenderer()->AddActor(MoveActor);
    this->isMove_ = true;
    //this->SelectedPoint = this->PointPicker->GetPointId();

    auto &pickPos = param_->activeTree->m_curveList[pickedVertexID_.first][pickedVertexID_.second];
    //double p[3];
    //vtkPolyData* Data = vtkPolyData::SafeDownCast(dstActor_->GetMapper()->GetInput());
    //Data->GetPoint(this->pickedVertexID_.second, p);
    this->MoveActor->SetPosition(pickPos(0) * param_->xRes_, pickPos(1) * param_->yRes_, pickPos(2) * param_->zRes_);

    //this->GetCurrentRenderer()->AddActor(this->MoveActor);
    this->InteractionProp = this->MoveActor;
    glbox->update();
}

void NGInteractorObjectStyle::OnLeftButtonDown()
{
    if (glbox->IsMovePointMode()) {
        glbox->ToggleMovePoint_Slot(false);
        glbox->MouseClick();
        if (!glbox->ToggleMovePoint_Slot(true))
            return;
    }
    
    vtkInteractorStyleTrackballActor::OnLeftButtonDown();
}
