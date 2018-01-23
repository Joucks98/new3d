#include <vtkObjectFactory.h>
#include <vtkAreaPicker.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkDataSetMapper.h>
#include <vtkPlanes.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRendererCollection.h>
#include <vtkPoints.h>
#include <vtkGenericCell.h>
#include <vtkVertex.h>

#include "InteractorStylePointCloud.h"

vtkStandardNewMacro(InteractorStylePointCloud);
InteractorStylePointCloud::InteractorStylePointCloud()
{
    m_pExtractor = vtkSmartPointer<vtkExtractPolyDataGeometry>::New();
    interactData = NULL;
}

void InteractorStylePointCloud::OnLeftButtonUp()
{
    vtkInteractorStyleRubberBandPick::OnLeftButtonUp();
    if (CurrentMode != VTKISRBP_SELECT)
    {
        return;
    }
    if (interactData != NULL)
    {
        auto picker = GetInteractor()->GetPicker();
        auto pickPts = static_cast<vtkAreaPicker*>(picker)->GetClipPoints();
        auto np = pickPts->GetNumberOfPoints();
        vtkSmartPointer<vtkPolyData> voi = GenerateVoi(pickPts);

        //auto pickList = GetInteractor()->GetPicker()->GetPickList();
        vtkPlanes* frustum = static_cast<vtkAreaPicker*>(GetInteractor()->GetPicker())->GetFrustum();
        /*vtkSmartPointer<vtkAreaPicker> arePicker = vtkSmartPointer<vtkAreaPicker>::New();
        vtkPlanes* frustum = arePicker->GetFrustum();*/
        m_pExtractor->SetImplicitFunction(frustum);
        m_pExtractor->Update();


        int npt = GetRoiPoints()->GetNumberOfPoints();
        int an = Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetNumberOfItems();
        showChosenActor();
        showActor(voi);
        an = Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetNumberOfItems();
    }
}

void InteractorStylePointCloud::SetInteractData(const vtkSmartPointer<vtkPolyData>& d)
{
    interactData = d;
#if VTK_MAJOR_VERSION <= 5
    extractGeometry->SetInput(interactData);
#else
    m_pExtractor->SetInputData(interactData);
#endif // VTK_MAJOR_VERSION <= 5

}

vtkPoints * InteractorStylePointCloud::GetRoiPoints() const
{
    return m_pExtractor->GetOutput()->GetPoints();
}

void InteractorStylePointCloud::showChosenActor()
{
    if (GetRoiPoints()->GetNumberOfPoints() == 0)
        return;
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(m_pExtractor->GetOutput());
    mapper->ScalarVisibilityOff();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    actor->GetProperty()->SetOpacity(1);
    actor->GetProperty()->SetPointSize(6);

    Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
    GetInteractor()->GetRenderWindow()->Render();
}

void InteractorStylePointCloud::showActor(const vtkSmartPointer<vtkPolyData>& polyData)
{
    if (polyData == NULL)
    {
        return;
    }
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(polyData);
    mapper->ScalarVisibilityOff();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1.0, .0, 1.0);
    actor->GetProperty()->SetOpacity(1);
    actor->GetProperty()->SetPointSize(3);

    Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(actor);
    GetInteractor()->GetRenderWindow()->Render();
}

vtkSmartPointer<vtkPolyData> InteractorStylePointCloud::GenerateVoi(const vtkSmartPointer<vtkPoints>& pts)
{
    if(pts->GetNumberOfPoints() == 0)
        return vtkSmartPointer<vtkPolyData>(); // NULL
    //vtkSmartPointer<vtkGenericCell> voiCell = vtkSmartPointer<vtkGenericCell>::New();
    //voiCell->SetCellType(VTK_VOXEL); // SetCellTypeToVoxel()
    //voiCell->InstantiateCell(voiCell->GetCellType());
    //voiCell->SetPoints(pts);

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(pts);
    //cell->InsertNextCell(voiCell);
    /*polyData->Allocate(VTK_PIXEL);
    for (size_t k = 0; k < 6; ++k)
    {
        vtkSmartPointer<vtkIdList> pixelIdList = vtkSmartPointer<vtkIdList>::New();
        pixelIdList->SetNumberOfIds(4);
        pixelIdList->SetId(0, )
    }*/
    for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); i++)
    {
        /*vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
        vertex->GetPoints()->SetPoint(0, pts->GetPoint(i));
        vertices->InsertNextCell(vertex);*/
        vtkSmartPointer<vtkIdList> pid = vtkSmartPointer<vtkIdList>::New();
        pid->SetNumberOfIds(1);
        pid->SetId(0, i);
        cells->InsertNextCell(pid);
    }
    polyData->SetVerts(cells);

    cells->Reset();
    vtkSmartPointer<vtkIdList> pid = vtkSmartPointer<vtkIdList>::New();
    pid->SetNumberOfIds(4);

    /*
    (6)+--------+(7)
      /|       /|
  (4)+--------+(5)
     | |      | |
    (2)+------|-+(3)
     |/       |/
  (0)+--------+(1)
    */
    pid->SetId(0, 0); pid->SetId(1, 2); pid->SetId(2, 3); pid->SetId(3, 1);
    cells->InsertNextCell(pid);
    pid->SetId(0, 4); pid->SetId(1, 5); pid->SetId(2, 7); pid->SetId(3, 6);
    cells->InsertNextCell(pid);

    pid->SetId(0, 0); pid->SetId(1, 4); pid->SetId(2, 6); pid->SetId(3, 2);
    cells->InsertNextCell(pid);
    pid->SetId(0, 1); pid->SetId(1, 3); pid->SetId(2, 7); pid->SetId(3, 5);
    cells->InsertNextCell(pid);

    pid->SetId(0, 0); pid->SetId(1, 1); pid->SetId(2, 5); pid->SetId(3, 4);
    cells->InsertNextCell(pid);
    pid->SetId(0, 2); pid->SetId(1, 6); pid->SetId(2, 7); pid->SetId(3, 3);
    cells->InsertNextCell(pid);
    polyData->SetPolys(cells);

    int nm = polyData->GetNumberOfCells();
    return polyData;
}


