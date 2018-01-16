#ifndef __IMAGECALLBACK_H__
#define __IMAGECALLBACK_H__
#include <vtkCellArray.h>
#include <vtkCommand.h>
#include <vtkImageViewer2.h>
#include <vtkPolyData.h>
#include <vtkPropPicker.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyle.h>
#include <vtkCornerAnnotation.h>
#include <vtkAssemblyPath.h>
#include <vtkProperty.h>
#include <vtkImageProperty.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkImageActor.h>
#include <vtkObjectFactory.h>

class ImageCallBack : public vtkCommand
{
public:
    vtkTypeMacro(ImageCallBack, vtkCommand);
    static ImageCallBack *New();
    ImageCallBack();
    void SetPicker(const vtkSmartPointer<vtkPropPicker>& picker) { this->picker = picker; }
    void SetViewer(vtkImageViewer2* viewer) { this->viewer = viewer;}
    void SetCanvasSource2DActor(vtkImageActor* actor) { this->imageCanvasSource2DActor = actor; }
    virtual void Execute(vtkObject*, unsigned long event, void *);

private:
    void updateCanvasSource2DActor(int extent[6], int x0, int x1, int y0, int y1);
    void updateCanvasSource2DActor(int x0, int x1, int y0, int y1);
    void blendCanvasSource2D();

    vtkImageViewer2* viewer;
    vtkSmartPointer<vtkPropPicker> picker;
    vtkSmartPointer<vtkImageCanvasSource2D> imageCanvasSource2D;
    vtkImageActor* imageCanvasSource2DActor;

    double pickPos1[3];
    double pickPos2[3];

    bool mouseMotion;
    bool validPick;
};


#endif // !__IMAGECALLBACK_H__
