#include <vtkImageBlend.h>
#include <vtkImageData.h>
#include <vtkInteractorStyleRubberBand2D.h>

#include "ImageCallBack.h"


vtkStandardNewMacro(ImageCallBack);

ImageCallBack::ImageCallBack(): viewer(NULL), picker(NULL), mouseMotion(false), validPick(false)
{
    pickPos1[0] = pickPos1[1] = pickPos1[2] = 0.0;
    pickPos2[0] = pickPos2[1] = pickPos2[2] = 0.0;
    imageCanvasSource2D = NULL;
}

void ImageCallBack::Execute(vtkObject *, unsigned long event, void *)
{
    /*if (event == vtkCommand::DeleteEvent)
    {
        return;
    }*/
    if (imageCanvasSource2DActor == NULL)
    {
        return;
    }

    vtkImageActor* imageActor = this->viewer->GetImageActor();
    vtkRenderWindowInteractor* interactor = this->viewer->GetRenderWindow()->GetInteractor();
    vtkInteractorStyle* style = vtkInteractorStyle::SafeDownCast(interactor->GetInteractorStyle());

    // pick at the mouse location provided by the interactor
    int* screenPos = interactor->GetEventPosition();
    vtkRenderer* renderer = this->viewer->GetRenderer();
    this->picker->Pick(screenPos[0], screenPos[1], 0.0, this->viewer->GetRenderer());

    // There could be other props assigned to this picker, so make sure we picked ths image actor.
    vtkAssemblyPath* path = this->picker->GetPath();
    bool validPick = false;
    if (path)
    {
        vtkCollectionSimpleIterator sit;
        path->InitTraversal(sit);
        vtkAssemblyNode *node;
        for (int i = 0; i < path->GetNumberOfItems() && !validPick; ++i)
        {
            node = path->GetNextNode(sit);
            if (imageActor == vtkImageActor::SafeDownCast(node->GetViewProp()))
            {
                validPick = true;
            }
        }
    }

    //// same to screen pos
    //auto tmp = dynamic_cast<vtkInteractorStyleRubberBand2D*>(style);
    //if (tmp != nullptr)
    //{
    //    int* spos = tmp->GetStartPosition();
    //    int* epos = tmp->GetEndPosition();
    //    cout << "correct!" << endl;
    //}


    if (!validPick)
    {
        // Pass the event further on
        interactor->Render();
        return;
    }

    switch (event)
    {
    case vtkCommand::RightButtonPressEvent:
        renderer->RemoveActor(imageCanvasSource2DActor);
        style->OnRightButtonDown();
        break;
    case vtkCommand::LeftButtonPressEvent:
        mouseMotion = true;
        this->picker->GetPickPosition(pickPos1);
        style->OnLeftButtonDown();
        break;
    case vtkCommand::LeftButtonReleaseEvent:
        if (mouseMotion)
        {
            auto tmp = dynamic_cast<vtkInteractorStyleRubberBand2D*>(style);
            if (screenPos[0] == tmp->GetStartPosition()[0] || screenPos[1] == tmp->GetStartPosition()[1])
            {
                break;
            }
            this->picker->GetPickPosition(pickPos2);

            int x0 = 0, x1 = 0, y0 = 0, y1 = 0;
            (pickPos1[0] > pickPos2[0]) ? (x0 = round(pickPos2[0]), x1 = round(pickPos1[0])) :
                (x0 = round(pickPos1[0]), x1 = round(pickPos2[0]));
            (pickPos1[1] > pickPos2[1]) ? (y0 = round(pickPos2[1]), y1 = round(pickPos1[1])) :
                (y0 = round(pickPos1[1]), y1 = round(pickPos2[1]));

            int* ex = this->viewer->GetInput()->GetExtent();
            
            renderer->RemoveActor(imageCanvasSource2DActor);
            int ii = renderer->HasViewProp(imageCanvasSource2DActor);
            updateCanvasSource2DActor(x0, x1, y0, y1);

            renderer->AddActor(imageCanvasSource2DActor);
            ii = renderer->HasViewProp(imageCanvasSource2DActor);
        }
        mouseMotion = false;
        style->OnLeftButtonUp();
        break;

    default:
        style->OnMouseMove();
        break;
    }

    // !Corruption: add style->OnMouseMove();
    interactor->Render();
}

void ImageCallBack::updateCanvasSource2DActor(int extent[6], int x0, int x1, int y0, int y1)
{
    imageCanvasSource2D = vtkSmartPointer<vtkImageCanvasSource2D>::New();
    imageCanvasSource2D->SetScalarTypeToUnsignedChar();
    imageCanvasSource2D->SetNumberOfScalarComponents(4);
    imageCanvasSource2D->SetExtent(extent);
    imageCanvasSource2D->SetDrawColor(0, 0, 0, 0);
    imageCanvasSource2D->FillBox(extent[0], extent[1], extent[2], extent[3]);

    imageCanvasSource2D->SetDrawColor(255, 255, 0, 255);
    imageCanvasSource2D->FillBox(x0, x1, y0, y1);
    imageCanvasSource2D->SetDrawColor(255, 0, 255, 255);
    imageCanvasSource2D->DrawSegment(x0, y0, x0, y1);
    imageCanvasSource2D->DrawSegment(x0, y1, x1, y1);
    imageCanvasSource2D->DrawSegment(x1, y1, x1, y0);
    imageCanvasSource2D->DrawSegment(x1, y0, x0, y0);
    imageCanvasSource2D->Update();

    imageCanvasSource2DActor->SetInputData(imageCanvasSource2D->GetOutput());
    imageCanvasSource2DActor->GetProperty()->SetOpacity(.5);
    imageCanvasSource2DActor->Update();
    //imageCanvasSource2DActor->Modified();
}

void ImageCallBack::updateCanvasSource2DActor(int x0, int x1, int y0, int y1)
{
    imageCanvasSource2D = vtkSmartPointer<vtkImageCanvasSource2D>::New();
    imageCanvasSource2D->SetScalarTypeToUnsignedChar();
    imageCanvasSource2D->SetNumberOfScalarComponents(4);
    imageCanvasSource2D->SetExtent(x0, x1, y0, y1, 0, 0);


    imageCanvasSource2D->SetDrawColor(255, 255, 0, 255);
    imageCanvasSource2D->FillBox(x0, x1, y0, y1);
    imageCanvasSource2D->SetDrawColor(255, 0, 255, 255);
    imageCanvasSource2D->DrawSegment(x0, y0, x0, y1);
    imageCanvasSource2D->DrawSegment(x0, y1, x1, y1);
    imageCanvasSource2D->DrawSegment(x1, y1, x1, y0);
    imageCanvasSource2D->DrawSegment(x1, y0, x0, y0);
    imageCanvasSource2D->Update();
    imageCanvasSource2DActor->SetInputData(imageCanvasSource2D->GetOutput());
    imageCanvasSource2DActor->GetProperty()->SetOpacity(.5);
    imageCanvasSource2DActor->Update();
    //imageCanvasSource2DActor->Modified();

}


void ImageCallBack::blendCanvasSource2D()
{
    static vtkSmartPointer<vtkImageData> originImage = this->viewer->GetInput();
    vtkSmartPointer<vtkImageBlend> blend = vtkSmartPointer<vtkImageBlend>::New();
    blend->SetInputData(originImage);
    blend->AddInputData(imageCanvasSource2D->GetOutput());
    blend->SetOpacity(1, .5);
    blend->Update();
    this->viewer->SetInputConnection(blend->GetOutputPort());
}
