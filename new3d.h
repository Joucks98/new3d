#ifndef __NEW3D_H__
#define __NEW3D_H__

#include <vector>
#include <assert.h>

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
        double value(int id) const {
            assert(id == 0 || id == 1 || id == 2);
            if (id == 0) return x;
            else if (id == 1) return y;
            else return z;
        }
        double x;
        double y;
        double z;
        double h;
        size_t index;  // index in image
    };

    class ImageParam
    {
    public:
        ImageParam()
        {
            offset[0] = -16.0;
            offset[1] = -14.99;
            offset[2] = -5.99;
            space[0] = 0.016;
            space[1] = 0.02;
            initValue = VTK_DOUBLE_MIN;
            rcNum[0] = 2000;
            rcNum[1] = 1500;
            compNum = 1;
        }
        double offset[3];
        double space[2];
        double initValue;
        int rcNum[2];
        int compNum;
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
    void showImage(const vtkSmartPointer<vtkImageData>& img, int comp);
    bool readData(const std::string& fileName);


    vtkSmartPointer<vtkPolyData> toBuildPointCloudData(const vtkSmartPointer<vtkPoints>& vtkPoints, 
                                                       const vtkSmartPointer<vtkDoubleArray>& colorField = NULL) const;
    vtkSmartPointer<vtkPolyData> toBuildPointCloudData(const std::vector<DIM3::Point3d>& point3dVec) const;
    vtkSmartPointer<vtkActor> toBuildPolyDataActor(const vtkSmartPointer<vtkPolyData>& pData);

    vtkSmartPointer<vtkImageData> toBuildImageData(std::vector<DIM3::Point3d>& point3dVec);
    vtkSmartPointer<vtkImageData> NEW3D::toBuildHeightImageData(const vtkSmartPointer<vtkImageData>& pImg);
    
    vtkSmartPointer<vtkImageData> initImageData(double initZ, int numComp);
    vtkSmartPointer<vtkImageData> initImageData(const DIM3::ImageParam&);
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

    std::map<size_t, DIM3::Point3d> m_Index_Point3_Map;


    bool isPointVecChanged;

    Ui::NEW3DClass ui;

};
#endif // !__NEW3D_H__

