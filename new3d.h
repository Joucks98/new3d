#ifndef __NEW3D_H__
#define __NEW3D_H__

#include <vector>

#include <QtWidgets/qmainwindow.h>

#include <vtkSmartPointer.h>

#include "ui_new3d.h"

class vtkImageViewer2;
class vtkRenderer;
class vtkCubeAxesActor;
class vtkPolyData;
class vtkLookupTable;
class vtkActor;
class vtkImageData;
class vtkRenderWindow;
class vtkPoints;
class vtkInteractorStyle;
class vtkInteractorStyleRubberBand2D;
class vtkImageActor;
class vtkScalarBarWidget;
class vtkDoubleArray;

class string;

class InteractorStylePointCloud;



namespace DIM3
{
    class Point3d
    {
    public:
        
        Point3d(double a, double b, double c):x(a), y(b), z(c), h(c), index(-1) {}
        //Point3d(const Point3d& t): x(t.x), y(t.y), z(t.z) {}
        double x;
        double y;
        double z;
        double h;
        size_t index;  // index in image
    };
}



class NEW3D : public QMainWindow
{
    Q_OBJECT

public:
    NEW3D(QWidget * parent = Q_NULLPTR);
    //~NEW3D();

private slots:
    void on_actionOpen_triggered();
    void on_actionImport_triggered();
    void on_actionZ_triggered();
    void on_actionCacheImage_triggered();
    void on_action3D_triggered();
    void on_actionRubber_toggled();
    void on_actionRubber_hovered();
    void on_actionX_triggered();
    void on_actionY_triggered();
    void on_actionFit_Plane_toggled();
    void on_actionPick_Points_toggled();
    void on_actionCorrect_triggered();

private:
    void showPointCloud(bool updateOrNot);
    void showColorImage(const vtkSmartPointer<vtkImageData>& img, int comp);
    bool readData(const std::string& fileName);


    vtkSmartPointer<vtkPolyData> toBuildPointCloudData(const vtkSmartPointer<vtkPoints>& vtkPoints, 
                                                       const vtkSmartPointer<vtkDoubleArray>& colorField = NULL) const;
    vtkSmartPointer<vtkPolyData> toBuildPointCloudData(const std::vector<DIM3::Point3d>& point3dVec) const;
    vtkSmartPointer<vtkActor> toBuildPolyDataActor(const vtkSmartPointer<vtkPolyData>& pData);
    vtkSmartPointer<vtkImageData> toBuildImageData(std::vector<DIM3::Point3d>& point3dVec);
    vtkSmartPointer<vtkImageData> initImageData(double initZ = 0);
    void updatePointCloud(bool update);
    void updateImage(bool update);

    vtkSmartPointer<vtkPoints> GetImageRoiPointsWorldData() const;
    vtkSmartPointer<vtkActor> generateFitPlaneActor(double* n, double* o, double width, double height) const;

    void updateCubeAxesActor();
    void initScalarBar(vtkSmartPointer<vtkScalarBarWidget>& scalarBarWidget);
    void initOrientationMarker();
    int initLookupTable(vtkSmartPointer<vtkLookupTable>& lut, double backGroundColor[4] = nullptr);

    void modifyPolyDataColorField(const vtkSmartPointer<vtkDoubleArray>& colorField);

    vtkSmartPointer<vtkLookupTable> m_pLookupTable;
    vtkSmartPointer<vtkImageViewer2> m_pImageViewer;
    vtkSmartPointer<vtkInteractorStyle> m_pImageStyle;
    vtkSmartPointer<vtkInteractorStyleRubberBand2D> m_pImageStyleRubber;
    vtkSmartPointer<vtkRenderWindow> m_pPointCloudWindow;
    vtkSmartPointer<vtkRenderer> m_pRenderer;
    vtkSmartPointer<InteractorStylePointCloud> m_pPointCloudStyle;
    vtkSmartPointer<vtkImageData> m_pImage;
    vtkSmartPointer<vtkPolyData> m_pPointCloudPolyData;
    vtkSmartPointer<vtkActor> m_pPointCloudActor;
    vtkSmartPointer<vtkCubeAxesActor> m_pCubeAxesActor;
    vtkMTimeType m_roi2DMTimeCache, m_roi3DMTimeCache;
    vtkSmartPointer<vtkImageActor> m_pRoi2DActor;
    vtkSmartPointer<vtkActor> m_pRoi3DActor;
    vtkSmartPointer<vtkActor> m_pFitPlaneActor;

    vtkSmartPointer<vtkScalarBarWidget> m_pScalarBarWidget;

    std::vector<DIM3::Point3d> m_point3dVec;
    vtkSmartPointer<vtkPoints> m_pPoints; // as a vtk object same to m_points3dVec




    bool isPointVecChanged;

    Ui::NEW3DClass ui;

};
#endif // !__NEW3D_H__

