/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef NGINTERACTORSTYLETRACKBALLCAMERA_H
#define NGINTERACTORSTYLETRACKBALLCAMERA_H
#include <vtkAutoInit.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCommand.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkAreaPicker.h>
#include <vtkCamera.h>
#include <vector>
#include "../../ngtypes/tree.h"
#include "../../ngtypes/ParamPack.h"
class NeuroGLWidget;
#define VTK_DEFINE(type, var) vtkSmartPointer<type> var
#define VTK_NEW(type, var) var = vtkSmartPointer<type>::New()
#define VTK_CREATE(type, var) vtkSmartPointer<type> var = vtkSmartPointer<type>::New()
#define LINEWIDTH_NOR 2
#define POINTSIZE_NOR 3
#define LINEWIDTH_CHOOSEN 5

class NGInteractorStyleTrackballCamera : public vtkInteractorStyleTrackballCamera
{
public:
    static NGInteractorStyleTrackballCamera* New();
    virtual ~NGInteractorStyleTrackballCamera();

    void SetWidget(NeuroGLWidget* arg);
    unsigned int numberOfClicks;
protected:
    NGInteractorStyleTrackballCamera();
    virtual void OnLeftButtonDown();
    virtual void OnLeftButtonUp();
    virtual void OnChar();
    virtual void OnRightButtonDown();
    virtual void OnRightButtonUp();
    virtual void OnMiddleButtonDown();
    virtual void OnMiddleButtonUp();
    virtual void OnMouseWheelForward();
    virtual void OnMouseWheelBackward();
    virtual void OnMouseMove();
    //virtual void OnKeyPress();
    virtual void Rotate();
private:
    
    NeuroGLWidget *glbox;
    bool bHasStartPoint_, isBox_, isConjoinMode_;
    int isLeftButtonDown_ = 0;
    vtkSmartPointer<vtkRenderer> tmpActorRenderer;
    vtkSmartPointer<vtkRenderer> tmpVolumeRenderer;
    vtkSmartPointer<vtkRenderer> tmpActiveRenderer;
};

class NGInteractorObjectStyle : public vtkInteractorStyleTrackballActor
{
public:
    static NGInteractorObjectStyle* New(){ return new NGInteractorObjectStyle(); }
    vtkTypeMacro(NGInteractorObjectStyle, vtkInteractorStyleTrackballActor);

    NGInteractorObjectStyle();

    void SetParam(NGParamPack arg){ param_ = arg; }

    virtual void OnLeftButtonDown();

    virtual  void OnMouseMove();

    virtual void OnMiddleButtonUp();
    virtual void OnMiddleButtonDown();

    void SetWidget(NeuroGLWidget* arg){ glbox = arg; }

    bool isMove_;
    bool isTreeRoot_;
    std::vector<size_t> leve1BranchIDList_;


    NGParamPack param_;

    //vtkPolyData* Data;
    vtkActor* dstActor_;
    vtkPolyData* GlyphData;
    NeuroGLWidget *glbox;

    vtkSmartPointer<vtkPolyDataMapper> MoveMapper;
    vtkSmartPointer<vtkActor> MoveActor;
    vtkSmartPointer<vtkPolyData> MovePolyData;
    vtkSmartPointer<vtkVertexGlyphFilter> MoveGlyphFilter;
    std::pair<vtkIdType, vtkIdType> pickedVertexID_;
    //vtkSmartPointer<vtkPointPicker> PointPicker;

    
    //vtkIdType SelectedPoint;
};

#endif // NGINTERACTORSTYLETRACKBALLCAMERA_H
