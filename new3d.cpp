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
#include <vtkExtractVOI.h>
#include <vtkPlaneSource.h>
#include <vtkDataArray.h>
#include <vtkCallbackCommand.h>
#include <vtkCubeAxesActor.h>
#include <vtkTextProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkScalarBarWidget.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>


#include "new3d.h"
#include "InteractorStylePointCloud.h"
#include "ImageCallBack.h"
#include "Algorithm3D.h"

using std::vector;
using std::stringstream;
using namespace DIM3;
//using DIM3::Point3d;


// 定义回调函数。注意回调函数的签名，不能更改。
static void CallbackFunc(vtkObject* obj, unsigned long eid, void* clientdata, void *calldata)
{
    
    std::cout << " ZYH " << std::endl;
}


NEW3D::NEW3D(QWidget *parent) : QMainWindow(parent), 
                                m_pPointCloudPolyData(NULL), 
                                m_pRoi3DActor(NULL), 
                                m_pFitPlaneActor(NULL),
                                isPointVecChanged(false)
{
    ui.setupUi(this);
    ui.actionRubber->setCheckable(true);
    ui.actionFit_Plane->setCheckable(true);
    ui.actionPick_Points->setCheckable(true);

    initLookupTable(m_pLookupTable);

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
    m_pCubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
    m_pRenderer->AddActor(m_pPointCloudActor);
    m_pRenderer->AddActor(m_pCubeAxesActor);
    m_pPointCloudWindow->AddRenderer(m_pRenderer);
    vtkSmartPointer<vtkRenderWindowInteractor> pointCloudInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    pointCloudInteractor->SetRenderWindow(m_pPointCloudWindow);
    m_pPointCloudStyle = vtkSmartPointer<InteractorStylePointCloud>::New();
    vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
    pointCloudInteractor->SetPicker(areaPicker);
    pointCloudInteractor->SetInteractorStyle(m_pPointCloudStyle);

    
    m_pRoi2DActor = vtkSmartPointer<vtkImageActor>::New();

    m_pScalarBarWidget = vtkSmartPointer<vtkScalarBarWidget>::New();

    m_roi2DMTimeCache = m_pRoi2DActor->GetMTime();
    m_roi3DMTimeCache = 0;

    ui.actionZ->setEnabled(false);
    ui.menuStyle->setEnabled(false);
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
        updatePointCloud(1);
        updateImage(1);
    }

    ui.actionZ->setEnabled(!m_point3dVec.empty());
    ui.action3D->setEnabled(!m_point3dVec.empty());
}

void NEW3D::on_actionZ_triggered()
{
    showColorImage(m_pImage, 2); // ...dig into
    ui.actionZ->setEnabled(false);
    ui.menuStyle->setEnabled(true);
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
    ui.menuStyle->setEnabled(false);
}

void NEW3D::on_actionRubber_toggled()
{
    bool checked = ui.actionRubber->isChecked();
    auto iteractor = m_pImageViewer->GetRenderWindow()->GetInteractor();
    if (checked)
    {
        m_pRoi2DActor->SetVisibility(1);
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
        imageCallback->SetCanvasSource2DActor(m_pRoi2DActor);
        m_pImageStyleRubber->AddObserver(vtkCommand::LeftButtonPressEvent, imageCallback);
        m_pImageStyleRubber->AddObserver(vtkCommand::LeftButtonReleaseEvent, imageCallback);
        m_pImageStyleRubber->AddObserver(vtkCommand::RightButtonPressEvent, imageCallback);
        // 将回调函数和 vtkCallbackCommand 联系起来
        vtkSmartPointer<vtkCallbackCommand> mouseCallback = vtkSmartPointer<vtkCallbackCommand>::New();
        mouseCallback->SetCallback(CallbackFunc);
        iteractor->AddObserver(vtkCommand::RightButtonPressEvent, mouseCallback);
    }
    else
    {
        m_pRoi2DActor->SetVisibility(0);
        m_pImageViewer->GetRenderWindow()->RemoveAllObservers();
        iteractor->SetInteractorStyle(m_pImageStyle);
    }
    // update m_qvtkWidget window
    m_pImageViewer->GetRenderWindow()->Render();
}

void NEW3D::on_actionRubber_hovered()
{
}

void NEW3D::on_actionX_triggered()
{
    showColorImage(m_pImage, 0);
    ui.menuStyle->setEnabled(true);
}

void NEW3D::on_actionY_triggered()
{
    showColorImage(m_pImage, 1);
    ui.menuStyle->setEnabled(true);
}

void NEW3D::on_actionFit_Plane_toggled()
{
    if (ui.actionFit_Plane->isChecked())
    {
        if ((m_pRoi3DActor->GetMTime() > m_roi3DMTimeCache)
            || (m_pFitPlaneActor == NULL))
        {
            // update m_pFitPlaneActor
            m_pRenderer->RemoveActor(m_pFitPlaneActor);
            auto roiData = GetImageRoiPointsWorldData();
            if (roiData != NULL)
            {
                double normal[3], origin[3];
                planeFitting(filterPoints(roiData), normal, origin);
                double *bound = roiData->GetBounds();
                m_pFitPlaneActor = generateFitPlaneActor(normal, origin, bound[1] - bound[0], bound[3] - bound[2]);
                // this all 0
                /*auto oo = m_pFitPlaneActor->GetOrigin();
                auto pp = m_pFitPlaneActor->GetPosition();*/
            }
            m_pRenderer->AddActor(m_pFitPlaneActor);

            m_roi3DMTimeCache = m_pRoi3DActor->GetMTime();
        }
        if (m_pFitPlaneActor != NULL)
        {
            m_pFitPlaneActor->SetVisibility(1);
        }        
    }
    else
    {
        m_pFitPlaneActor->SetVisibility(0);
    }
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

void NEW3D::on_actionPick_Points_toggled()
{
    m_pRenderer->RemoveActor(m_pRoi3DActor);
    if (ui.actionPick_Points->isChecked())
    {
        if (m_pRoi2DActor->GetMTime() > m_roi2DMTimeCache)
        {
            // update m_pRoi3DActor
            //m_pRenderer->RemoveActor(m_pRoi3DActor);
            vtkSmartPointer<vtkPoints> roiPts = GetImageRoiPointsWorldData();
            if (roiPts != NULL)
            {
                // save roiData
                ofstream outfile("../roiData.txt", ios::trunc);
                for (int i = 0; i < roiPts->GetNumberOfPoints(); ++i)
                {
                    outfile << roiPts->GetPoint(i)[0] << "," <<
                        roiPts->GetPoint(i)[1] << "," <<
                        roiPts->GetPoint(i)[2] << endl;
                }
                outfile.close();

                // build m_pRoi3DActor
                auto polyData = toBuildPointCloudData(roiPts);
                m_pRoi3DActor = toBuildPolyDataActor(polyData);                
                // set m_pRoi3DActor size & color
                if (m_pRoi3DActor == nullptr)
                {
                    return;
                }
                m_pRoi3DActor->GetProperty()->SetPointSize(3);
                m_pRoi3DActor->GetProperty()->SetColor(1, .3, .3);
                m_pRoi3DActor->GetProperty()->SetOpacity(.3);
            }
            //m_pRenderer->AddActor(m_pRoi3DActor);

            m_roi2DMTimeCache = m_pRoi2DActor->GetMTime();
        }
        m_pRoi3DActor->SetVisibility(1);
        m_pRenderer->AddActor(m_pRoi3DActor);
    }
    else
    {
        //m_pRoi3DActor->SetVisibility(0); // will modify m_pRoi3DActor
    }
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

void NEW3D::on_actionCorrect_triggered()
{
    auto roiData = GetImageRoiPointsWorldData();
    //vtkSmartPointer<vtkDoubleArray> colorScalarArray = vtkSmartPointer<vtkDoubleArray>::New();
    if (roiData != NULL)
    {
        double normal[3], origin[3];
        if (planeFitting(filterPoints(roiData), normal, origin))
        {
            std::cerr << "fitting error!" << endl;
            return;
        }
        vtkMath::Normalize(normal);

        double* ptr = (double*)m_pImage->GetScalarPointer(); // change m_pImage
        
        int numComp = m_pImage->GetNumberOfScalarComponents();
        vector<double> tmp(numComp);        
        for (size_t i = 0; i < m_point3dVec.size(); ++i)
        {
            size_t scalarId = numComp * m_point3dVec[i].index;
            for (size_t j = 0; j < numComp; ++j)
            {
                tmp[j] = ptr[scalarId + j] - origin[j];
            }
            double height = vtkMath::Dot(normal, &tmp[0]);
            ptr[scalarId + 2] = height;
            m_point3dVec[i].h = height;
        }

        /*ofstream outfile("../height.txt",ios::trunc);
        for (int i = 0; i < m_point3dVec.size(); ++i)
        {
            outfile << m_point3dVec[i].h<< endl;
        }
        outfile.close();*/
    }
    //m_pImage->GetPointData()->SetActiveScalars("Height_Field");

    // the following line return all three components of points in m_pImage
    // auto vv = m_pImage->GetPointData()->GetScalars();
    showColorImage(m_pImage, 2);
    //ui.m_qVTKViewer->GetRenderWindow()->Render();

    
    vector<double> tmp(m_point3dVec.size());
    std::transform(m_point3dVec.begin(), m_point3dVec.end(), tmp.begin(), [](const auto& a) {return a.h; });
    modifyPolyDataColorField(vec2vtkDoubleArray(tmp));
}


void NEW3D::showPointCloud(bool updateOrNot)
{
    if (m_pPointCloudActor == NULL)
        return;
    // update m_pPointCloudActor
    m_pRenderer->RemoveActor(m_pPointCloudActor);
    updatePointCloud(updateOrNot);
    m_pRenderer->AddActor(m_pPointCloudActor);

    // set color
    m_pPointCloudActor->GetMapper()->SetLookupTable(m_pLookupTable);
    //m_pLookupTable->SetRange(m_pPoints->GetBounds()[4], m_pPoints->GetBounds()[5]);



    //m_pPointCloudActor->GetMapper()->GetInput()->
    //pMapper->SetScalarRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    m_pPointCloudActor->GetMapper()->UseLookupTableScalarRangeOn();
    m_pPointCloudActor->GetMapper()->SetScalarModeToUsePointFieldData();
    m_pPointCloudActor->GetMapper()->SelectColorArray("Color_Field");    
    
    // update cubeAxes
    updateCubeAxesActor();
    // update scalerBarWidget
    
    initScalarBar(m_pScalarBarWidget);

    // rendering
    m_pRenderer->ResetCamera();
    m_pRenderer->SetBackground(0, 0, 0);
    // change window to m_pPointCloudWindow
    ui.m_qVTKViewer->SetRenderWindow(m_pPointCloudWindow);

    // m_pRenderer->Render() is not enough to update window.
    ui.m_qVTKViewer->GetRenderWindow()->Render();
}

void NEW3D::showColorImage(const vtkSmartPointer<vtkImageData>& pImg, int comp)
{
    //updateImage(updateOrNot);
    if (pImg == NULL)
        return;

    vtkSmartPointer<vtkImageExtractComponents> extractCompFilter =
        vtkSmartPointer<vtkImageExtractComponents>::New();
    extractCompFilter->SetInputData(pImg);
    extractCompFilter->SetComponents(comp);
    extractCompFilter->Update();

    vtkSmartPointer<vtkImageMapToColors> colorMap =
        vtkSmartPointer<vtkImageMapToColors>::New();
    colorMap->SetInputConnection(extractCompFilter->GetOutputPort());
    colorMap->SetLookupTable(m_pLookupTable);
    // change range to involve all value

    vtkDataArray* tmp = extractCompFilter->GetOutput()->GetPointData()->GetScalars(); // return active scalar data
    // this only work on z component
    // vtkDataArray* tmp = img->GetPointData()->GetScalars("Height_Field");
    
    double minHeight = VTK_DOUBLE_MAX, maxHeight = VTK_DOUBLE_MIN;
    for (vtkIdType i = 0; i < tmp->GetNumberOfTuples(); ++i)
    {
        if (tmp->GetTuple1(i) > maxHeight)
        {
            maxHeight = tmp->GetTuple1(i);
        }
        if (tmp->GetTuple1(i) > VTK_DOUBLE_MIN && tmp->GetTuple1(i) < minHeight)
        {
            minHeight = tmp->GetTuple1(i);
        }
    }
    
    //m_pLookupTable->SetRange(minHeight, maxHeight);

    /* !error:
    double* bounds = extractCompFilter->GetOutput()->GetBounds();
    m_pLookupTable->SetRange(bounds[2], bounds[3]);*/

    colorMap->Update();

    vtkSmartPointer<vtkImageChangeInformation> changer =
        vtkSmartPointer<vtkImageChangeInformation>::New();
    changer->SetInputConnection(colorMap->GetOutputPort());
    changer->SetOutputOrigin(0, 0, 0);
    //changer->SetSpacingScale(1/ pImg->GetSpacing()[0], 1/pImg->GetSpacing()[1], 1);
    changer->Update();

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

vtkSmartPointer<vtkPolyData> NEW3D::toBuildPointCloudData(const vtkSmartPointer<vtkPoints>& vtkPoints, 
                                                          const vtkSmartPointer<vtkDoubleArray>& colorField) const
{
    if(vtkPoints == NULL ||
       vtkPoints->GetNumberOfPoints() == 0)
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
    if (colorField != nullptr && colorField->GetNumberOfTuples() == num)
    {
        colorField->SetName("Color_Field");
        polyData->GetPointData()->AddArray(colorField);
    }
    else
    {
        vtkSmartPointer<vtkDoubleArray> tmp = vtkSmartPointer<vtkDoubleArray>::New();
        tmp->SetName("Color_Field");
        tmp->SetNumberOfComponents(1);
        tmp->SetNumberOfTuples(num);
        for (vtkIdType i = 0; i < num; ++i)
        {
            tmp->SetValue(i, vtkPoints->GetPoint(i)[2]);
        }
        polyData->GetPointData()->AddArray(tmp);
    }   

    return polyData;
}

vtkSmartPointer<vtkPolyData> NEW3D::toBuildPointCloudData(const std::vector<DIM3::Point3d>& point3dVec) const
{
    if (point3dVec.empty())
    {
        return nullptr;
    }
    vtkSmartPointer<vtkPoints> vtkPoints = vec2vtkPoints(point3dVec);
    auto m = vtkPoints->GetData()->GetNumberOfValues();
    auto mm = vtkPoints->GetData()->GetSize();
    vector<double> tmp(point3dVec.size());
    std::transform(point3dVec.begin(), point3dVec.end(), tmp.begin(), [](const auto& a) {return a.h; });
    return toBuildPointCloudData(vtkPoints, vec2vtkDoubleArray(tmp));
}

vtkSmartPointer<vtkActor> NEW3D::toBuildPolyDataActor(const vtkSmartPointer<vtkPolyData>& pData)
{
    if(pData == NULL)
        return vtkSmartPointer<vtkActor>();

    vtkSmartPointer<vtkPolyDataMapper> pMapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();
    pMapper->SetInputData(pData);

    //pMapper->SetLookupTable(m_pLookupTable);
    //m_pLookupTable->SetRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    ////pMapper->SetScalarRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    //pMapper->UseLookupTableScalarRangeOn();
    //pMapper->SetScalarModeToUsePointFieldData();
    //pMapper->SelectColorArray("Color_Field");

    vtkSmartPointer<vtkActor> pActor = vtkSmartPointer<vtkActor>::New();
    pActor->SetMapper(pMapper);
    return pActor;
}

vtkSmartPointer<vtkImageData> NEW3D::toBuildImageData(vector<Point3d>& point3dVec)
{
    if(point3dVec.empty())
        return vtkSmartPointer<vtkImageData>();

    // init imagedata parameter
    auto img = initImageData(VTK_DOUBLE_MIN);


    vtkSmartPointer<vtkDoubleArray> heightField = vtkSmartPointer<vtkDoubleArray>::New();
    heightField->SetName("Height_Field");
    heightField->SetNumberOfComponents(1);
    heightField->SetNumberOfTuples(img->GetNumberOfPoints());
    for (vtkIdType i = 0; i < img->GetNumberOfPoints(); ++i)
    {
        heightField->SetValue(i, VTK_DOUBLE_MIN);
    }


    double* ptr = (double*)img->GetScalarPointer();
    int numComp = img->GetNumberOfScalarComponents();
    for (size_t i = 0; i < point3dVec.size(); ++i)
    {
        size_t yId = round((point3dVec[i].y - img->GetOrigin()[1]) / img->GetSpacing()[1]);
        size_t xId = round((point3dVec[i].x - img->GetOrigin()[0]) / img->GetSpacing()[0]);
        size_t index = img->GetDimensions()[0] * yId + xId;
        *(ptr + index*numComp + 0) = point3dVec[i].x;
        *(ptr + index*numComp + 1) = point3dVec[i].y;
        *(ptr + index*numComp + 2) = point3dVec[i].z;
        //*(ptr + index*img->GetNumberOfScalarComponents() + 3) = 1;

        heightField->SetValue(index, point3dVec[i].h);
        // double* tmp = img->GetPoint(index);  // this point has the same (x,y), but not the z .
        
        point3dVec[i].index = index;
    }

    
    img->GetPointData()->AddArray(heightField);


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
        //*(spanEnd - 1) = initZ;
        while (spanIter != spanEnd)
        {
            *spanIter = initZ;
            ++spanIter;
        }
        iter.NextSpan();
    }
    return img;
}

void NEW3D::updatePointCloud(bool update)
{
    if (update)
    {
        m_pPointCloudPolyData = toBuildPointCloudData(m_point3dVec);
        m_pPointCloudStyle->SetInteractData(m_pPointCloudPolyData);
        m_pPointCloudActor = toBuildPolyDataActor(m_pPointCloudPolyData);
    }
}

void NEW3D::updateImage(bool update)
{
    if (update)
    {
        m_pImage = toBuildImageData(m_point3dVec);
    }
}

vtkSmartPointer<vtkPoints> NEW3D::GetImageRoiPointsWorldData() const
{
    vtkSmartPointer<vtkPoints> result = NULL;
    if (/*!ui.actionRubber->isChecked() ||*/ 
        m_pImage == NULL || 
        m_pRoi2DActor->GetInput()->GetNumberOfPoints() == 0
        ) 
    {
        return result;
    }

    vtkImageData* roiImg = m_pRoi2DActor->GetInput();    
    int* extent = roiImg->GetExtent(); // here extent keeps roiImg global coordiates

    double* space = m_pImage->GetSpacing();
    int minIdx = static_cast<int>(extent[0] / space[0]);
    int maxIdx = static_cast<int>(extent[1] / space[0]);
    int minIdy = static_cast<int>(extent[2] / space[1]);
    int maxIdy = static_cast<int>(extent[3] / space[1]);    

    /* //function: GetScalarPointerForExtent isn't same to extractVOI
    int dim[6] = { minIdx, maxIdx, minIdy, maxIdy, 0, 0 };
    auto pts = static_cast<double*>(m_pImage->GetScalarPointerForExtent(dim));*/

    vtkSmartPointer<vtkExtractVOI> extractVOI =
        vtkSmartPointer<vtkExtractVOI>::New();
    extractVOI->SetInputData(m_pImage);
    extractVOI->SetVOI(minIdx, maxIdx, minIdy, maxIdy, 0, 0);
    extractVOI->Update();

    vtkImageData* extractImg = extractVOI->GetOutput();
    int num = extractImg->GetNumberOfPoints();
    if (num == 0)
    {
        return result;
    }
    result = vtkSmartPointer<vtkPoints>::New();
    // You cannot determine the points you picked have value or not.
    //result->SetNumberOfPoints(num);

    auto pts2 = static_cast<double*>(extractImg->GetScalarPointer());
    int numComp = extractImg->GetNumberOfScalarComponents();
    for (size_t i = 0; i < num; ++i)
    {
        //double* tmp = extractImg->GetPoint(i);
        if (pts2[numComp * i+2] > VTK_DOUBLE_MIN)
            result->InsertNextPoint(&pts2[numComp*i]);
        //result->SetPoint(i, &pts2[numComp * i]);
    }
        
    /*
    double* ptr = (double*)m_pImage->GetScalarPointer();
    double* pts3 = new double[comp*(maxIdy - minIdy + 1)*(maxIdx - minIdx + 1)];
    int k = 0;
    for (size_t i = minIdy; i <= maxIdy; ++i)
    {
        for (size_t j = minIdx; j <= maxIdx; ++j)
        {
            size_t index = m_pImage->GetDimensions()[0] * i + j;
            double x = ptr[index*m_pImage->GetNumberOfScalarComponents() + 0];
            double y = ptr[index*m_pImage->GetNumberOfScalarComponents() + 1];
            double z = ptr[index*m_pImage->GetNumberOfScalarComponents() + 2];
            pts3[k++] = x;
            pts3[k++] = y;
            pts3[k++] = z;
        }
    }
    delete[]pts3;*/
    return result;
}


vtkSmartPointer<vtkActor> NEW3D::generateFitPlaneActor(double * n, double * o, double width, double height) const
{
    //auto ptsData = GetImageRoiPointsWorldData();
    //if (ptsData == NULL)
    //    return vtkSmartPointer<vtkActor>();

    //double normal[3], origin[3];
    //planeFitting(filterPoints(ptsData), normal, origin);

    if (n == nullptr || o == nullptr 
        || width < 0 || height < 0 )
    {
        return vtkSmartPointer<vtkActor>();
    }

    // build plane actor
    vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
    // set plane size
    planeSource->SetOrigin(0, 0, 0);
    planeSource->SetPoint1(0, height, 0);
    planeSource->SetPoint2(width, 0, 0);
    // set origin and normal
    planeSource->SetCenter(o);
    planeSource->SetNormal(n);
    planeSource->Update();
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(planeSource->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    actor->GetProperty()->SetPointSize(1);
    actor->GetProperty()->SetColor(0, 1, 1);
    actor->GetProperty()->SetOpacity(0.6);
    return actor;
}

void NEW3D::updateCubeAxesActor()
{
    double *range = m_pPointCloudActor->GetBounds();//= { -1, 5, -1, 5, 0, 1 };
    m_pCubeAxesActor->SetBounds(range);
    m_pCubeAxesActor->SetAxisOrigin(range[0], range[2], range[4]);
    m_pCubeAxesActor->SetCamera(m_pRenderer->GetActiveCamera());
    m_pCubeAxesActor->GetTitleTextProperty(0)->SetColor(1.0, 0.0, 0.0);
    m_pCubeAxesActor->GetLabelTextProperty(0)->SetColor(1.0, 0.0, 0.0);
    m_pCubeAxesActor->GetTitleTextProperty(1)->SetColor(0.0, 1.0, 0.0);
    m_pCubeAxesActor->GetLabelTextProperty(1)->SetColor(0.0, 1.0, 0.0);
    m_pCubeAxesActor->GetTitleTextProperty(2)->SetColor(0.0, 0.0, 1.0);
    m_pCubeAxesActor->GetLabelTextProperty(2)->SetColor(0.0, 0.0, 1.0);

    //m_pCubeAxesActor->DrawXGridlinesOn();
    //m_pCubeAxesActor->DrawYGridlinesOn();
    //m_pCubeAxesActor->DrawZGridlinesOn();
    m_pCubeAxesActor->DrawXInnerGridlinesOn();
    m_pCubeAxesActor->DrawYInnerGridlinesOn();

    //m_pCubeAxesActor->SetGridLineLocation(2);
    m_pCubeAxesActor->XAxisMinorTickVisibilityOn();
    m_pCubeAxesActor->YAxisMinorTickVisibilityOn();
    m_pCubeAxesActor->ZAxisMinorTickVisibilityOn();

    //m_pCubeAxesActor->SetFlyModeToFurthestTriad();
    //m_pCubeAxesActor->SetFlyModeToOuterEdges();
    m_pCubeAxesActor->SetFlyModeToStaticTriad();

    //m_pCubeAxesActor->StickyAxesOn();
    //m_pCubeAxesActor->CenterStickyAxesOff();
    
}

void NEW3D::initScalarBar(vtkSmartPointer<vtkScalarBarWidget>& scalarBarWidget)
{
    // show color bar
    vtkSmartPointer<vtkScalarBarActor> scalarBarActor = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBarActor->SetOrientationToHorizontal();
    //scalarBarActor->SetOrientationToVertical();
    scalarBarActor->SetLookupTable(m_pLookupTable);    
    scalarBarWidget->SetInteractor(m_pPointCloudWindow->GetInteractor());
    scalarBarWidget->SetScalarBarActor(scalarBarActor);
    scalarBarWidget->On();
}

void NEW3D::initOrientationMarker()
{
    ///////////////////////////// vtkOrientationMarkerWidget
    vtkSmartPointer<vtkAxesActor> iconActor = vtkSmartPointer<vtkAxesActor>::New();
    vtkSmartPointer<vtkOrientationMarkerWidget> orientationWidget =
        vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    orientationWidget->SetOutlineColor(0.9300, 0.5700, 0.1300);
    orientationWidget->SetInteractor(m_pPointCloudWindow->GetInteractor());
    orientationWidget->SetOrientationMarker(iconActor);
    orientationWidget->SetViewport(0.0, 0.0, 0.2, 0.2);
    orientationWidget->SetEnabled(1);
    orientationWidget->InteractiveOn();
}

int NEW3D::initLookupTable(vtkSmartPointer<vtkLookupTable>& lut, double backGroundColor[4])
{
    lut = vtkSmartPointer<vtkLookupTable>::New();
    double color[4] = { 0.0, 0.0, 0.0, 1.0 };
    if (backGroundColor != nullptr)
    {
        std::copy(backGroundColor, backGroundColor + 4, color);
    }    
    lut->SetBelowRangeColor(color);
    lut->SetAboveRangeColor(color);
    lut->SetUseBelowRangeColor(1);
    lut->SetUseAboveRangeColor(1); 
    lut->SetRange(-5, 3);
    lut->Build();

    return 0;
}

void NEW3D::modifyPolyDataColorField(const vtkSmartPointer<vtkDoubleArray>& colorField)
{
    if (colorField == nullptr || 
        m_pPointCloudPolyData->GetNumberOfPoints() != colorField->GetNumberOfTuples())
    {
        return;
    }    
    auto colors = m_pPointCloudPolyData->GetPointData()->GetScalars("Color_Field");
    if (colors != nullptr)
    {
        colors->Resize(colorField->GetNumberOfTuples());
        for (vtkIdType i = 0; i < colors->GetNumberOfTuples(); ++i)
        {
            colors->SetTuple1(i, colorField->GetTuple1(i));
        }
    }    
}
