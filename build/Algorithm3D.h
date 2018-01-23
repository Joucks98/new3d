#ifndef __ALGORITHM3D_H_
#define __ALGORITHM3D_H_
#include <vtkPoints.h>
#include <vtkPCAStatistics.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPlane.h>
#include <vtkKMeansStatistics.h>
#include <vtkIdList.h>

#include <numeric>
#include <algorithm>


namespace DIM3 {
    vtkSmartPointer<vtkPoints> vec2vtkPoints(const std::vector<Point3d>& ptVec)
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(ptVec.size());
        for (size_t i = 0; i < ptVec.size(); ++i)
        {
            points->SetPoint(i, ptVec[i].x, ptVec[i].y, ptVec[i].z);
        }
        return points;
    }

    // dim: number of components of tuples in return vector
    void vtkPoints2Vec(vtkPoints* vtkPts, std::vector<double>* ptVec, int begin = 0, int dim = 3)
    {
        ptVec->clear();
        if (vtkPts == NULL || ptVec == NULL)
            return;
        ptVec->reserve(dim * vtkPts->GetNumberOfPoints());
        for (vtkIdType i = 0; i < vtkPts->GetNumberOfPoints(); ++i)
        {
            double* tuple = vtkPts->GetPoint(i);
            std::copy(tuple+begin, tuple + dim, std::back_inserter(*ptVec));
        }
    }

    double calcStandardVariance(const std::vector<double>& data)
    {
        if (data.size() < 2)
            return -1;
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
        double sigma = 0;
        for (size_t i = 0; i < data.size(); ++i)
        {
            sigma += (data[i] - mean)*(data[i] - mean);
        }
        sigma /= (data.size() - 1);
        sigma = std::sqrt(sigma);
        return sigma;
    }
    vtkSmartPointer<vtkPoints> filterPoints(vtkPoints* pts, int numberOfClusters = 2)
    {
        std::vector<double> zVec;
        vtkPoints2Vec(pts, &zVec, 2);
        // get rid of outliers
        double sigma = calcStandardVariance(zVec);

        std::sort(zVec.begin(), zVec.end());
        double midZValue = zVec[(size_t)(zVec.size()*.5)]; // middle value

        vtkSmartPointer<vtkIdList> pickedIdList0 = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> pickedIdList1 = vtkSmartPointer<vtkIdList>::New();
        for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i)
        {
            double* tmp = pts->GetPoint(i);
            if (std::abs(tmp[2] - midZValue) > 2 * sigma)
            {
                continue;
            }
            if (tmp[2] < midZValue)
            {
                pickedIdList0->InsertNextId(i);
            }
            else
            {
                pickedIdList1->InsertNextId(i);
            }                
        }

        vtkSmartPointer<vtkPoints> result = vtkSmartPointer<vtkPoints>::New();
        if (pickedIdList0->GetNumberOfIds() > pickedIdList1->GetNumberOfIds())
        {
            pts->GetPoints(pickedIdList0, result);
        }
        else
        {
            pts->GetPoints(pickedIdList1, result);
        }

        return result;

        // Get the points into the format needed for KMeans
        vtkSmartPointer<vtkTable> inputData =
            vtkSmartPointer<vtkTable>::New();

        for (int c = 0; c < 3; ++c)
        {
            /*std::stringstream colName;
            colName << "coord " << c;*/
            vtkSmartPointer<vtkDoubleArray> doubleArray =
                vtkSmartPointer<vtkDoubleArray>::New();
            //doubleArray->SetName(colName.str().c_str());
            doubleArray->SetNumberOfComponents(1);
            doubleArray->SetNumberOfTuples(pts->GetNumberOfPoints());


            for (int r = 0; r < pts->GetNumberOfPoints(); ++r)
            {
                double p[3];
                pts->GetPoint(r, p);

                doubleArray->SetValue(r, p[c]);
            }

            inputData->AddColumn(doubleArray);
        }

        // kmean clustering
        vtkSmartPointer<vtkKMeansStatistics> kMeansStatistics =
            vtkSmartPointer<vtkKMeansStatistics>::New();

#if VTK_MAJOR_VERSION <= 5
        kMeansStatistics->SetInput(vtkStatisticsAlgorithm::INPUT_DATA, datasetTable);
#else
        kMeansStatistics->SetInputData(vtkStatisticsAlgorithm::INPUT_DATA, inputData);
#endif
        kMeansStatistics->SetColumnStatus(inputData->GetColumnName(0), 1);
        kMeansStatistics->SetColumnStatus(inputData->GetColumnName(1), 1);
        kMeansStatistics->SetColumnStatus(inputData->GetColumnName(2), 1);
        kMeansStatistics->RequestSelectedColumns();
        kMeansStatistics->SetAssessOption(true);
        kMeansStatistics->SetDefaultNumberOfClusters(numberOfClusters);
        kMeansStatistics->Update();

        // Analyze clusters
        /*vtkSmartPointer<vtkIntArray> clusterArray =
            vtkSmartPointer<vtkIntArray>::New();
        clusterArray->SetNumberOfComponents(1);
        clusterArray->SetName("ClusterId");*/

        vtkSmartPointer<vtkIdList> cluster0IdList = vtkSmartPointer<vtkIdList>::New();
        vtkSmartPointer<vtkIdList> cluster1IdList = vtkSmartPointer<vtkIdList>::New();
        for (vtkIdType r = 0; r < kMeansStatistics->GetOutput()->GetNumberOfRows(); ++r)
        {
            vtkVariant v = kMeansStatistics->GetOutput()->GetValue(r, kMeansStatistics->GetOutput()->GetNumberOfColumns() - 1);
            //clusterArray->InsertNextValue(v.ToInt());
            std::cout << "Point " << r << " is in cluster " << v.ToInt() << std::endl;
            if (v.ToInt() == 0)
            {
                cluster0IdList->InsertNextId(r);
            }
            else
            {
                cluster1IdList->InsertNextId(r);
            }
        }

        vtkSmartPointer<vtkPoints> reducedPts = vtkSmartPointer<vtkPoints>::New();
        if (cluster0IdList->GetNumberOfIds() > pts->GetNumberOfPoints() * .5)
        {
            pts->GetPoints(cluster0IdList, reducedPts);
        }
        else
        {
            pts->GetPoints(cluster1IdList, reducedPts);
        }

        return reducedPts;
    }

    int planeFitting(vtkPoints* pts, double normal[3], double origin[3], int iterNum = 100)
    {
        if (pts == NULL || pts->GetNumberOfPoints() == 0)
        {
            return 1; // error occur
        }

        vtkSmartPointer<vtkTable> datasetTable = vtkSmartPointer<vtkTable>::New();
        vtkIdType num = pts->GetNumberOfPoints();
        double alpha = std::sqrt(num - 1);

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
            xArray->InsertNextValue(tmp[0] * alpha);
            yArray->InsertNextValue(tmp[1] * alpha);
            zArray->InsertNextValue(tmp[2] * alpha);
        }

        /*vtkSmartPointer<vtkTable> datasetTable =
        vtkSmartPointer<vtkTable>::New();*/
        datasetTable->AddColumn(xArray);
        datasetTable->AddColumn(yArray);
        datasetTable->AddColumn(zArray);


        // calc origin

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
        // Display the results
        //pcaStatistics->GetOutput()->Dump();

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

        double tmpNormal[3];
        for (vtkIdType i = 0; i < eigenVector->GetNumberOfTuples(); ++i)
        {
            eigenVector->GetTuple(i, &tmpNormal[i]);
        }
        double len = vtkMath::Normalize(tmpNormal);

        double meanDist = len*d / num;


        vtkSmartPointer<vtkPlane> tmpPlane = vtkSmartPointer<vtkPlane>::New();
        tmpPlane->SetOrigin(origin);
        tmpPlane->SetNormal(tmpNormal);
        std::vector<double> distVec(num);
        for (vtkIdType i = 0; i < num; ++i)
        {
            distVec[i] = tmpPlane->DistanceToPlane(pts->GetPoint(i));
        }

        double meanDist1 = std::accumulate(distVec.begin(), distVec.end(), 0.0) / num;
        double sigma = calcStandardVariance(distVec);
        vtkSmartPointer<vtkIdList> chosenIdList = vtkSmartPointer<vtkIdList>::New();
        int k = 0;
        for (vtkIdType i = 0; i < num; ++i)
        {
            if (std::abs(distVec[i] - meanDist1) < 2 * sigma)
            {
                chosenIdList->InsertId(k++, i);
            }
        }        
        vtkSmartPointer<vtkPoints> chosenPts = vtkSmartPointer<vtkPoints>::New();
        pts->GetPoints(chosenIdList, chosenPts);

        // if there still left interNum-1 times to execute 
        // and the result tmpNormal is not convergent, recursive. 
        if ((std::abs(vtkMath::Dot(tmpNormal, normal) - 1) > 1e-12)
            && (iterNum-1) > 0)
        {
            planeFitting(chosenPts, tmpNormal, origin, (iterNum-1));
        }
        normal[0] = tmpNormal[0];
        normal[1] = tmpNormal[1];
        normal[2] = tmpNormal[2];

        return 0;
    }
}




#endif // !__ALGORITHM3D_H_

