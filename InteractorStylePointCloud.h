#ifndef __INTERACTORSTYLEPOINTCLOUD_H__
#define __INTERACTORSTYLEPOINTCLOUD_H__

#include <vtkSmartPointer.h>
#include <vtkInteractorStyleRubberBandPick.h>

class vtkPolyData;
class vtkActor;
class vtkDataSetMapper;
class vtkPoints;
class vtkExtractPolyDataGeometry;

#define VTKISRBP_ORIENT 0
#define VTKISRBP_SELECT 1

class InteractorStylePointCloud : public vtkInteractorStyleRubberBandPick
{
public:
    static InteractorStylePointCloud* New();
    vtkTypeMacro(InteractorStylePointCloud, vtkInteractorStyleRubberBandPick);
    InteractorStylePointCloud();
    virtual void OnLeftButtonUp();

    void SetInteractData(const vtkSmartPointer<vtkPolyData>& d);
    vtkPoints* GetRoiPoints() const;

private:
    void showChosenActor();
    void showActor(const vtkSmartPointer<vtkPolyData>& polyData);
    vtkSmartPointer<vtkPolyData> GenerateVoi(const vtkSmartPointer<vtkPoints>& pts);
    vtkSmartPointer<vtkExtractPolyDataGeometry> m_pExtractor;
    vtkSmartPointer<vtkPolyData> interactData;

};
#endif // !__INTERACTORSTYLEPOINTCLOUD_H__
