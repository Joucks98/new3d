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
#include <vtkTable.h>
#include <vtkPCAStatistics.h>
#include <vtkDataArray.h>

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



NEW3D::NEW3D(QWidget *parent) : QMainWindow(parent), 
                                m_pPointCloudPolyData(NULL), 
                                m_pRoi3DActor(NULL), 
                                m_pFitPlaneActor(NULL),
                                m_pPoints(NULL), 
                                isPointVecChanged(false)
{
    ui.setupUi(this);
    ui.actionRubber->setCheckable(true);
    ui.actionFit_Plane->setCheckable(true);

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

    
    m_pRoi2DActor = vtkSmartPointer<vtkImageActor>::New();

    m_roiMTimeCache = m_pRoi2DActor->GetMTime();

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
    }
    else
    {
        //m_pImageViewer->GetRenderer()->RemoveActor(m_pRoi2DActor);
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
    showColorImage(0);
    ui.menuStyle->setEnabled(true);
}

void NEW3D::on_actionY_triggered()
{
    showColorImage(1);
    ui.menuStyle->setEnabled(true);
}

void NEW3D::on_actionFit_Plane_toggled()
{
    if (ui.actionFit_Plane->isChecked())
    {
        if ((m_pRoi2DActor->GetMTime() > m_roiMTimeCache) || (m_pFitPlaneActor == NULL))
        {
            m_pRenderer->RemoveActor(m_pFitPlaneActor);
            m_pFitPlaneActor = generateFitPlaneActor();
            m_pRenderer->AddActor(m_pFitPlaneActor);
        }
        m_pFitPlaneActor->SetVisibility(1);
    }
    else
    {
        m_pFitPlaneActor->SetVisibility(0);
    }
    ui.m_qVTKViewer->GetRenderWindow()->Render();
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
    m_pLookupTable->SetRange(m_pPoints->GetBounds()[4], m_pPoints->GetBounds()[5]);
    //pMapper->SetScalarRange(pData->GetPoints()->GetBounds()[4], pData->GetPoints()->GetBounds()[5]);
    m_pPointCloudActor->GetMapper()->UseLookupTableScalarRangeOn();
    m_pPointCloudActor->GetMapper()->SetScalarModeToUsePointFieldData();
    m_pPointCloudActor->GetMapper()->SelectColorArray("Color_Field");

    // update m_pRoi3DActor
    if (m_pRoi2DActor->GetMTime() > m_roiMTimeCache)
    {
        m_pRenderer->RemoveActor(m_pRoi3DActor);
        vtkSmartPointer<vtkPoints> roiPts = GetImageRoiPointsWorldData();
        if (roiPts != NULL)
        {
            auto polyData = toBuildPointCloudData(roiPts);
            m_pRoi3DActor = toBuildPolyDataActor(polyData);
            m_pRenderer->AddActor(m_pRoi3DActor);
            // set roi Actor size & color
            m_pRoi3DActor->GetProperty()->SetPointSize(3);
            m_pRoi3DActor->GetProperty()->SetColor(1, .3, .3);
            m_pRoi3DActor->GetProperty()->SetOpacity(.3);
        }
        m_roiMTimeCache = m_pRoi2DActor->GetMTime();
    }
    
    /*if (m_pFitPlaneActor != NULL)
        m_pRenderer->AddActor(m_pFitPlaneActor);*/

    // rendering
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
    // change range to involve all value
    m_pLookupTable->SetRange(m_pPoints->GetBounds()[2*comp], m_pPoints->GetBounds()[2*comp+1]);
    colorMap->Update();

    vtkSmartPointer<vtkImageChangeInformation> changer =
        vtkSmartPointer<vtkImageChangeInformation>::New();
    changer->SetInputConnection(colorMap->GetOutputPort());
    changer->SetOutputOrigin(0, 0, 0);
    //changer->SetSpacingScale(1/ m_pImage->GetSpacing()[0], 1/m_pImage->GetSpacing()[1], 1);
    changer->Update();

    // update roi actor
    if (m_pRoi2DActor->GetMTime() > m_roiMTimeCache)
    {
        m_pImageViewer->GetRenderer()->RemoveActor(m_pRoi2DActor);
    }

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
    if(vtkPoints == NULL || vtkPoints->GetNumberOfPoints() == 0)
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
    colerField->SetName("Color_Field");
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

vtkSmartPointer<vtkImageData> NEW3D::toBuildImageData()
{
    if(m_point3dVec.empty())
        return vtkSmartPointer<vtkImageData>();

    // init imagedata parameter
    auto img = initImageData(VTK_DOUBLE_MIN);
    double* ptr = (double*)img->GetScalarPointer();
    for (size_t i = 0; i < m_point3dVec.size(); ++i)
    {
        size_t yId = round((m_point3dVec[i].y - img->GetOrigin()[1]) / img->GetSpacing()[1]);
        size_t xId = round((m_point3dVec[i].x - img->GetOrigin()[0]) / img->GetSpacing()[0]);
        size_t index = img->GetDimensions()[0] * yId + xId;
        *(ptr + index*img->GetNumberOfScalarComponents() + 0) = m_point3dVec[i].x;
        *(ptr + index*img->GetNumberOfScalarComponents() + 1) = m_point3dVec[i].y;
        *(ptr + index*img->GetNumberOfScalarComponents() + 2) = m_point3dVec[i].z;
        //*(ptr + index*img->GetNumberOfScalarComponents() + 3) = 1;


        double* tmp = img->GetPoint(index);
        
    }

    /*vtkSmartPointer<vtkIntArray> indicateField = vtkSmartPointer<vtkIntArray>::New();
    indicateField->SetName("Indicate_Field");
    indicateField->SetNumberOfComponents(1);
    indicateField->SetNumberOfTuples(img->GetNumberOfPoints());
    
    for (vtkIdType i = 0; i < img->GetNumberOfPoints(); ++i)
    {
        indicateField->SetValue(i, vtkPoints->GetPoint(i)[2]);
    }
    polyData->GetPointData()->AddArray(indicateField);*/


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

vtkSmartPointer<vtkPoints> NEW3D::GetImageRoiPointsWorldData() const
{
    vtkSmartPointer<vtkPoints> result = NULL;
    if (!ui.actionRubber->isChecked() || 
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

    /* //function of the following lines isn't same to extractVOI
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
namespace {
    int planeFitting(vtkPoints* pts, double normal[3], double origin[3], int iterNum = 3)
    {
        if (pts == NULL || pts->GetNumberOfPoints() == 0)
        {
            return 1; // error occur
        }

        vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
        vtkIdType num = pts->GetNumberOfPoints();
        double alpha = sqrt(num-1);
        
        //int nComp = pts->GetData()->GetNumberOfComponents();
        //char name[] = { 'x', '\0' };
        //for (size_t i = 0; i < nComp; ++i)
        //{
        //    vtkSmartPointer<vtkDoubleArray> array =
        //        vtkSmartPointer<vtkDoubleArray>::New();
        //    name[0] += i;
        //    array->SetName(name);
        //    array->SetNumberOfTuples(num);
        //    array->SetNumberOfComponents(1);
        //    datasetTable->AddColumn(array);

        //    origin[i] = 0;
        //    for (vtkIdType k = 0; k < num; ++k)
        //    {
        //        //double* tmp = pts->GetPoint(i); same to  double* tmp1 = pts->GetData()->GetTuple(i);
        //        double tmp = pts->GetData()->GetComponent(k, i);
        //        array->SetValue(k, tmp*alpha);       
        //        origin[i] += array->GetComponent(k, 0); // same to : tmp*alpha
        //    }
        //    origin[i] /= num;
        //}





        //-----------------------------------------------------
        // These would be all of your "x" values.
        vtkSmartPointer<vtkDoubleArray> xArray =
            vtkSmartPointer<vtkDoubleArray>::New();
        xArray->SetNumberOfComponents(1);
        xArray->SetName("x");

        // These would be all of your "y" values.
        vtkSmartPointer<vtkDoubleArray> yArray =
            vtkSmartPointer<vtkDoubleArray>::New();
        yArray->SetNumberOfComponents(1);
        yArray->SetName("y");

        // These would be all of your "z" values.
        vtkSmartPointer<vtkDoubleArray> zArray =
            vtkSmartPointer<vtkDoubleArray>::New();
        zArray->SetNumberOfComponents(1);
        zArray->SetName("z");

        for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i)
        {
            double* tmp = pts->GetPoint(i);
            xArray->InsertNextValue(tmp[0]*alpha);
            yArray->InsertNextValue(tmp[1]*alpha);
            zArray->InsertNextValue(tmp[2]*alpha);
        }

        /*vtkSmartPointer<vtkTable> datasetTable =
            vtkSmartPointer<vtkTable>::New();*/
        datasetTable->AddColumn(xArray);
        datasetTable->AddColumn(yArray);
        datasetTable->AddColumn(zArray);

        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;
        for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i)
        {
            double* tmp = pts->GetPoint(i);
            origin[0] += tmp[0];
            origin[1] += tmp[1];
            origin[2] += tmp[2];
        }
        origin[0] /= pts->GetNumberOfPoints();
        origin[1] /= pts->GetNumberOfPoints();
        origin[2] /= pts->GetNumberOfPoints();

        //-----------------------------------------------------







        vtkSmartPointer<vtkPCAStatistics> pcaStatistics =
            vtkSmartPointer<vtkPCAStatistics>::New();
        pcaStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
        pcaStatistics->SetColumnStatus("x", 1);
        pcaStatistics->SetColumnStatus("y", 1);
        pcaStatistics->SetColumnStatus("z", 1);
        pcaStatistics->RequestSelectedColumns();
        pcaStatistics->SetDeriveOption(true);
        //pcaStatistics->SetMedianAbsoluteDeviation(true);
        pcaStatistics->Update();

        ///////// Eigenvalues ////////////
        vtkSmartPointer<vtkDoubleArray> eigenvalues =
            vtkSmartPointer<vtkDoubleArray>::New();
        pcaStatistics->GetEigenvalues(eigenvalues);
        double ev[3] = { eigenvalues->GetValue(0), eigenvalues->GetValue(1), eigenvalues->GetValue(2) };
        // min value of eigenvalues
        double d = eigenvalues->GetValue(2);
        ///////// Eigenvectors ////////////
        /*vtkSmartPointer<vtkDoubleArray> eigenvectors =
            vtkSmartPointer<vtkDoubleArray>::New();
        pcaStatistics->GetEigenvectors(eigenvectors);*/
        vtkSmartPointer<vtkDoubleArray> eigenVector = vtkSmartPointer<vtkDoubleArray>::New();
        pcaStatistics->GetEigenvector(2, eigenVector);
        for (vtkIdType i = 0; i < eigenVector->GetNumberOfTuples(); ++i)
        {
            eigenVector->GetTuple(i, &normal[i]);
        }
        vtkMath::Normalize(normal);
        return 0;
    }
}

vtkSmartPointer<vtkActor> NEW3D::generateFitPlaneActor() const
{
    auto ptsData = GetImageRoiPointsWorldData();
    if (ptsData == NULL)
        return vtkSmartPointer<vtkActor>();
    double normal[3], origin[3];
    planeFitting(ptsData, normal, origin);
    // build plane actor
    vtkSmartPointer<vtkPlaneSource> planeSource = vtkSmartPointer<vtkPlaneSource>::New();
    // set plane size
    planeSource->SetOrigin(0, 0, 0);
    planeSource->SetPoint1(0, ptsData->GetBounds()[3]- ptsData->GetBounds()[2], 0);
    planeSource->SetPoint2(ptsData->GetBounds()[1] - ptsData->GetBounds()[0], 0, 0);
    // set origin and normal
    planeSource->SetCenter(origin);
    planeSource->SetNormal(normal);
    //planeSource->SetResolution(100, 100);
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
