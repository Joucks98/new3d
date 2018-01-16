#include <iostream>
#include <string>
#include <sstream>

#include <qdir.h>
#include <qfiledialog.h>
#include <qobject.h>
#include <qbytearray.h>
#include <qstring.h>

#include <vtkSmartPointer.h>
#include <vtkJPEGReader.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkImageActor.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkImageMapToColors.h>
#include <vtkImageData.h>
#include <vtkImageIterator.h>
#include <vtkInformation.h>
#include <vtkImageProperty.h>
#include <vtkRendererCollection.h>
#include <vtkActorCollection.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAreaPicker.h>
#include <vtkInteractorStyleRubberBand2D.h>
#include <vtkPropPicker.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageChangeInformation.h>
#include <vtkImageExtractComponents.h>


#include "new3d.h"
#include "InteractorStylePointCloud.h"
#include "ImageCallBack.h"


using std::vector;
using std::stringstream;
//using namespace DIM3;
using DIM3::Point3d;

static vtkSmartPointer<vtkPoints> vec2vtkPoints(const vector<Point3d>& ptVec)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(ptVec.size());
    for (size_t i = 0; i < ptVec.size(); ++i)
    {
        points->SetPoint(i, ptVec[i].x, ptVec[i].y, ptVec[i].z);
    }
    return points;
}



NEW3D::NEW3D(QWidget *parent)
    : QMainWindow(parent), m_pPointCloudPolyData(NULL), m_pPoints(NULL), isPointVecChanged(false)
{
    ui.setupUi(this);
    ui.actionRubber->setCheckable(true);

    m_pLookupTable = vtkSmartPointer<vtkLookupTable>::New();
    m_pLookupTable->SetBelowRangeColor(0.0, 0.0, 0.0, 1.0);
    m_pLookupTable->SetAboveRangeColor(0.0, 0.0, 0.0, 1.0);
    m_pLookupTable->UseBelowRangeColorOn();
    m_pLookupTable->UseAboveRangeColorOn();
    m_pLookupTable->Build();

    m_pImageViewer = vtkSmartPointer<vtkImageViewer2>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> imageInterator = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    m_pImageViewer->SetupInteractor(imageInterator);

    m_pImageStyle = m_pImageViewer->GetInteractorStyle();
    m_pImageViewer->GetInteractorStyle()->SetInteractionModeToImage3D();

    m_pImageStyleRubber = vtkSmartPointer<vtkInteractorStyleRubberBand2D>::New();
    imageInterator->SetInteractorStyle(m_pImageStyle);

    m_pPointCloudWindow = vtkSmartPointer<vtkRenderWindow>::New();
    m_pRenderer = vtkSmartPointer<vtkRenderer>::New();
    m_pPointCloudActor = vtkSmartPointer<vtkActor>::New();
    m_pRenderer->AddActor(m_pPointCloudActor);
    m_pPointCloudWindow->AddRenderer(m_pRenderer);
    vtkSmartPointer<vtkRenderWindowInteractor> pointCloudInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    pointCloudInteractor->SetRenderWindow(m_pPointCloudWindow);
    m_pPointCloudStyle = vtkSmartPointer<InteractorStylePointCloud>::New();
    vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
    pointCloudInteractor->SetPicker(areaPicker);
    pointCloudInteractor->SetInteractorStyle(m_pPointCloudStyle);


    m_pROIActor = vtkSmartPointer<vtkImageActor>::New();



    ui.actionZ->setEnabled(false);
    ui.action3D->setEnabled(false);
}



void NEW3D::on_actionOpen_triggered()
{
    QString filter = "Data File(*.jpg *.jpeg *.*)";
    QDir dir;
    // ... dig into
    QString fileName = QFileDialog::getOpenFileName(this,
        QString(tr("Open File")), dir.absolutePath(), filter);
    if (fileName.isEmpty())
        return;

    QByteArray ba = fileName.toLocal8Bit();
    const char* fileNameStr = ba.data();

    vtkSmartPointer<vtkJPEGReader> reader = vtkSmartPointer<vtkJPEGReader>::New();
    reader->SetFileName(fileNameStr);
    int flag = reader->CanReadFile(fileNameStr);

    reader->Update();

    m_pImageViewer->SetInputConnection(reader->GetOutputPort());
    m_pImageViewer->UpdateDisplayExtent();

    m_pImageViewer->SetSliceOrientationToXY();
    m_pImageViewer->GetImageActor()->InterpolateOff();
    m_pImageViewer->GetRenderer()->ResetCamera();

    ui.m_qVTKViewer->SetRenderWindow(m_pImageViewer->GetRenderWindow());
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

void NEW3D::on_actionImport_triggered()
{
    QString filter = "Data File(*.txt *.*)";
    QDir dir;
    // ...dig into
    QString fileName = QFileDialog::getOpenFileName(this,
        QString(tr("Open File")), dir.absolutePath(), filter);
    if (fileName.isEmpty())
        return;
    //... dig into : check file format,give feedback.
    bool flag = readData(fileName.toStdString());
    if (flag)
    {
        m_pPoints = vec2vtkPoints(m_point3dVec);
        updatePointCloud(1);
        updateImage(1);
    }

    ui.actionZ->setEnabled(!m_point3dVec.empty());
    ui.action3D->setEnabled(!m_point3dVec.empty());
}

void NEW3D::on_actionZ_triggered()
{
    showColorImage(2); // ...dig into
    ui.actionZ->setEnabled(false);
    ui.action3D->setEnabled(true);
}

void NEW3D::on_actionCacheImage_triggered()
{
    auto img = ui.m_qVTKViewer->cachedImage();
    m_pImageViewer->SetInputData(img);
    ui.m_qVTKViewer->SetRenderWindow(m_pImageViewer->GetRenderWindow());
    m_pImageViewer->GetRenderer()->ResetCamera();
    m_pImageViewer->Render();
}

void NEW3D::on_action3D_triggered()
{
    showPointCloud(0);
    ui.action3D->setEnabled(false);
    ui.actionZ->setEnabled(true);
}

void NEW3D::on_actionRubber_toggled()
{
    bool checked = ui.actionRubber->isChecked();
    auto iteractor = m_pImageViewer->GetRenderWindow()->GetInteractor();
    if (checked)
    {
        iteractor->SetInteractorStyle(m_pImageStyleRubber);
        // Picker to pick pixels
        vtkSmartPointer<vtkPropPicker> propPicker = vtkSmartPointer<vtkPropPicker>::New();
        propPicker->PickFromListOn();
        // Give the picker a prop to pick
        vtkSmartPointer<vtkImageActor> imageActor = m_pImageViewer->GetImageActor();
        propPicker->AddPickList(imageActor);
        // disable interpolation, so we can see each pixel
        imageActor->InterpolateOff();
        vtkSmartPointer<ImageCallBack> imageCallback = vtkSmartPointer<ImageCallBack>::New();
        imageCallback->SetViewer(m_pImageViewer);
        imageCallback->SetPicker(propPicker);
        imageCallback->SetCanvasSource2DActor(m_pROIActor);
        m_pImageStyleRubber->AddObserver(vtkCommand::LeftButtonPressEvent, imageCallback);
        m_pImageStyleRubber->AddObserver(vtkCommand::LeftButtonReleaseEvent, imageCallback);
        m_pImageStyleRubber->AddObserver(vtkCommand::RightButtonPressEvent, imageCallback);
    }
    else
    {
        m_pImageViewer->GetRenderer()->RemoveActor(m_pROIActor);
        m_pImageViewer->GetRenderWindow()->RemoveAllObservers();
        iteractor->SetInteractorStyle(m_pImageStyle);
    }
    // update m_qvtkWidget window
    m_pImageViewer->GetRenderWindow()->Render();
}

void NEW3D::on_actionRubber_hovered()
{
    if (ui.actionRubber->isChecked())
    {
        /*double* origin0 = m_pImageViewer->GetInput()->GetOrigin();
        int* d0 = m_pImageViewer->GetInput()->GetDimensions();
        int* e0 = m_pImageViewer->GetInput()->GetExtent();*/

        double* origin1 = m_pROIActor->GetInput()->GetOrigin(); // default always zero point unless change image information.
        int* d1 = m_pROIActor->GetInput()->GetDimensions();
        int* e1 = m_pROIActor->GetInput()->GetExtent(); // the global coordiates
        bool vis = ui.actionRubber->isVisible();


        double* space = m_pImage->GetSpacing();
        //double space[3] = { 1, 1, 1 };
        int y0 = static_cast<int>(e1[2] / space[1]);
        int y1 = static_cast<int>(e1[3] / space[1]);
        int x0 = static_cast<int>(e1[0] / space[0]);
        int x1 = static_cast<int>(e1[1] / space[0]);

        double* ptr = (double*)m_pImage->GetScalarPointer();
        int dimX = m_pImage->GetDimensions()[0];
        for (size_t j = y0; j <= y1; ++j)
        {
            for (size_t i = x0; i <= x1; ++i)
            {
                ptr[dimX*j + i] = -15;
            }
            
        }
    }   
}

void NEW3D::on_actionX_triggered()
{
    showColorImage(0);
}


void NEW3D::showPointCloud(bool updateOrNot)
{
    if (m_pPointCloudActor == NULL)
        return;
    // update m_pPointCloudActor
    m_pRenderer->RemoveActor(m_pPointCloudActor);
    updatePointCloud(updateOrNot);
    m_pRenderer->AddActor(m_pPointCloudActor);

    m_pRenderer->ResetCamera();
    m_pRenderer->SetBackground(0, 0, 0);
    // change window to m_pPointCloudWindow
    ui.m_qVTKViewer->SetRenderWindow(m_pPointCloudWindow);

    // m_pRenderer->Render() is not enough to update window.
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

void NEW3D::showColorImage(int comp, bool updateOrNot)
{
    updateImage(updateOrNot);
    if (m_pImage == NULL)
        return;

    vtkSmartPointer<vtkImageExtractComponents> extractCompFilter =
        vtkSmartPointer<vtkImageExtractComponents>::New();
    extractCompFilter->SetInputData(m_pImage);
    extractCompFilter->SetComponents(comp);
    extractCompFilter->Update();

    vtkSmartPointer<vtkImageMapToColors> colorMap =
        vtkSmartPointer<vtkImageMapToColors>::New();
    colorMap->SetInputConnection(extractCompFilter->GetOutputPort());
    colorMap->SetLookupTable(m_pLookupTable);
    m_pLookupTable->SetRange(m_pPoints->GetBounds()[4], m_pPoints->GetBounds()[5]);
    colorMap->Update();

    vtkSmartPointer<vtkImageChangeInformation> changer =
        vtkSmartPointer<vtkImageChangeInformation>::New();
    changer->SetInputConnection(colorMap->GetOutputPort());
    changer->SetOutputOrigin(0, 0, 0);
    //changer->SetSpacingScale(1/ m_pImage->GetSpacing()[0], 1/m_pImage->GetSpacing()[1], 1);
    changer->Update();

    // clear roi actor
    m_pImageViewer->GetRenderer()->RemoveActor(m_pROIActor);
    // !error: 
    //m_pImageViewer->GetImageActor()->SetInputData(changer->GetOutput());
    m_pImageViewer->SetInputConnection(changer->GetOutputPort());
    // change window to imageviewer
    // !error: m_pImageViewer->SetRenderWindow(ui.m_qVTKViewer->GetRenderWindow());
    ui.m_qVTKViewer->SetRenderWindow(m_pImageViewer->GetRenderWindow());
    ui.m_qVTKViewer->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->ResetCamera();
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

bool NEW3D::readData(const std::string & fileName)
{
    ifstream filesteram(fileName);
    std::string line;

    m_point3dVec.clear();
    while (getline(filesteram, line))
    {
        double x, y, z;
        char c;
        stringstream linestream;
        linestream << line;
        linestream >> x >> c >> y >> c >> z >> c >> z;
        m_point3dVec.emplace_back(x, y, z);
    }
    filesteram.close();
    return true;
}

vtkSmartPointer<vtkPolyData> NEW3D::toBuildPointCloudData(const vtkSmartPointer<vtkPoints>& vtkPoints)
{
    if(vtkPoints->GetNumberOfPoints() == 0)
        return vtkSmartPointer<vtkPolyData>();
    // initialize cells: one point one cell
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vertices->Reset();
    vtkIdType num = vtkPoints->GetNumberOfPoints();
    for (vtkIdType i = 0; i < num; ++i)
    {
        vtkSmartPointer<vtkIdList> pid = vtkSmartPointer<vtkIdList>::New();
        pid->SetNumberOfIds(1);
        pid->SetId(0, i);
        vertices->InsertNextCell(pid);
    }

    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(vtkPoints);
    polyData->SetVerts(vertices);

    // set coler field
    vtkSmartPointer<vtkDoubleArray> colerField = vtkSmartPointer<vtkDoubleArray>::New();
    colerField->SetName("Color Field");
    colerField->SetNumberOfComponents(1);
    colerField->SetNumberOfTuples(num);
    for (vtkIdType i = 0; i < num; ++i)
    {
        colerField->SetValue(i, vtkPoints->GetPoint(i)[2]);
    }
    polyData->GetPointData()->AddArray(colerField);

    return polyData;
}

vtkSmartPointer<vtkActor> NEW3D::toBuildPolyDataActor(const vtkSmartPointer<vtkPolyData>& pData)
{
    if(pData == NULL)
        return vtkSmartPointer<vtkActor>();

    vtkSmartPointer<vtkPolyDataMapper> pMapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    pMapper->SetInputData(pData);
    pMapper->SetLookupTable(m_pLookupTable);
    m_pLookupTable->SetRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    //pMapper->SetScalarRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    pMapper->UseLookupTableScalarRangeOn();
    pMapper->SetScalarModeToUsePointFieldData();
    pMapper->SelectColorArray("Color Field");

    vtkSmartPointer<vtkActor> pActor = vtkSmartPointer<vtkActor>::New();
    pActor->SetMapper(pMapper);
    return pActor;
}

vtkSmartPointer<vtkImageData> NEW3D::toBuildImageData()
{
    if(m_point3dVec.empty())
        return vtkSmartPointer<vtkImageData>();

    // init imagedata parameter
    auto img = initImageData();
    double* ptr = (double*)img->GetScalarPointer();
    for (size_t i = 0; i < m_point3dVec.size(); ++i)
    {
        size_t yId = round((m_point3dVec[i].y - img->GetOrigin()[1]) / img->GetSpacing()[1]);
        size_t xId = round((m_point3dVec[i].x - img->GetOrigin()[0]) / img->GetSpacing()[0]);
        size_t index = img->GetDimensions()[0] * yId + xId;
        *(ptr + index*img->GetNumberOfScalarComponents() + 0) = m_point3dVec[i].x;
        *(ptr + index*img->GetNumberOfScalarComponents() + 1) = m_point3dVec[i].y;
        *(ptr + index*img->GetNumberOfScalarComponents() + 2) = m_point3dVec[i].z;

    }
    return img;
}

vtkSmartPointer<vtkImageData> NEW3D::initImageData(double initZ)
{
    // ... dig into
    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
    double x_offset = -16.0, y_offset = -14.99, z_offset = -5.99;
    double dist_x = 0.016, dist_y = 0.02;
    int iNum = 2000, jNum = 1500;
    int componentNum = 3;
    img->SetDimensions(iNum, jNum, 1);
    img->SetOrigin(x_offset, y_offset, z_offset);
    img->SetSpacing(dist_x, dist_y, 0);
    auto info = vtkSmartPointer<vtkInformation>::New();
    img->SetScalarType(VTK_DOUBLE, info);
    img->SetNumberOfScalarComponents(componentNum, info); // x, y, z
    img->AllocateScalars(info);

    vtkImageIterator<double> iter(img, img->GetExtent());
    while (!iter.IsAtEnd())
    {
        double* spanIter = iter.BeginSpan();
        double* spanEnd = iter.EndSpan();
        *(spanEnd - 1) = initZ;
        /*while (spanIter != spanEnd)
        {
            *spanIter = initZ;
            ++spanIter;
        }*/
        iter.NextSpan();
    }
    return img;
}

void NEW3D::updatePointCloud(bool update)
{
    if (update)
    {
        m_pPointCloudPolyData = toBuildPointCloudData(m_pPoints);
        m_pPointCloudStyle->SetInteractData(m_pPointCloudPolyData);
        m_pPointCloudActor = toBuildPolyDataActor(m_pPointCloudPolyData);
    }
}

void NEW3D::updateImage(bool update)
{
    if (update)
    {
        m_pImage = toBuildImageData();
    }
}
