#ifndef __ALGORITHM3D_H__
#define __ALGORITHM3D_H__
#include <vtkPoints.h>
#include <vtkPCAStatistics.h>
#include <vtkTable.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkPlane.h>
#include <vtkKMeansStatistics.h>
#include <vtkIdList.h>

#include <numeric> // accumulate
#include <algorithm>
//#include "D:\\Program Files\\eigen\\Eigen\\Core"
#include "D:\\Program Files\\eigen\\Eigen\\Eigen"

using namespace Eigen;
using std::vector;

namespace DIM3 {
    vtkSmartPointer<vtkDoubleArray> vec2vtkDoubleArray(const std::vector<double> vec)
    {
        vtkSmartPointer<vtkDoubleArray> result = vtkSmartPointer<vtkDoubleArray>::New();
        result->SetNumberOfTuples((vtkIdType)vec.size());
        result->SetNumberOfComponents(1);
        for (size_t i = 0; i < vec.size(); ++i)
        {
            result->SetValue(i, vec[i]);
        }
        return result;
    }

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

    void vtkPoints2Point3dVec(vtkPoints* vtkPts, std::vector<Point3d>* ptVec)
    {
        ptVec->clear();
        if (vtkPts == NULL || ptVec == NULL)
            return;
        ptVec->reserve(vtkPts->GetNumberOfPoints());
        for (vtkIdType i = 0; i < vtkPts->GetNumberOfPoints(); ++i)
        {
            double* tuple = vtkPts->GetPoint(i);
            ptVec->emplace_back(tuple[0], tuple[1], tuple[2]);
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
        if (pts == nullptr || pts->GetNumberOfPoints() == 0)
        {
            return nullptr;
        }
        
        std::vector<Point3d> originalPointsVec;
        vtkPoints2Point3dVec(pts, &originalPointsVec);
        auto num = originalPointsVec.size();
        std::sort(originalPointsVec.begin(), originalPointsVec.end(), [](const Point3d& a, const Point3d& b) {return a.z < b.z; });
        auto midIter = originalPointsVec.begin() + static_cast<size_t>(num*.5);
        auto lowIter = originalPointsVec.begin() + static_cast<size_t>(.025*num);
        auto highIter = originalPointsVec.begin() + static_cast<size_t>(.975*num);
        std::vector<Point3d> resultPoints;
        if (midIter->z < .5*(lowIter->z + highIter->z))
        {
            resultPoints.assign(lowIter, midIter);
        }
        else
        {
            resultPoints.assign(midIter, highIter);
        }
        return vec2vtkPoints(resultPoints);


        //std::vector<double> zVec;
        //vtkPoints2Vec(pts, &zVec, 2);
        //// get rid of outliers
        //double sigma = calcStandardVariance(zVec);
        //std::sort(zVec.begin(), zVec.end());
        //double midZValue = zVec.size() > 0 ? zVec[(size_t)(zVec.size()*.5)]: 0.0; // middle value
        //vtkSmartPointer<vtkIdList> pickedIdList0 = vtkSmartPointer<vtkIdList>::New();
        //vtkSmartPointer<vtkIdList> pickedIdList1 = vtkSmartPointer<vtkIdList>::New();
        //for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i)
        //{
        //    double* tmp = pts->GetPoint(i);
        //    if (std::abs(tmp[2] - midZValue) > 2 * sigma)
        //    {
        //        continue;
        //    }
        //    if (tmp[2] < midZValue)
        //    {
        //        pickedIdList0->InsertNextId(i);
        //    }
        //    else
        //    {
        //        pickedIdList1->InsertNextId(i);
        //    }                
        //}
        //vtkSmartPointer<vtkPoints> result = vtkSmartPointer<vtkPoints>::New();
        //if (pickedIdList0->GetNumberOfIds() > pickedIdList1->GetNumberOfIds())
        //{
        //    pts->GetPoints(pickedIdList0, result);
        //}
        //else
        //{
        //    pts->GetPoints(pickedIdList1, result);
        //}
        //return result;

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

    bool modeFilter(const std::vector<Point3d> & pointsVec, std::vector<double>* pickRange)
    {
        assert(pickRange != nullptr);
        pickRange->clear();
        auto num = pointsVec.size();
        if (pointsVec.empty())
        {
            return false;
        }
        /*VectorXi indexes = VectorXi::LinSpaced(num, 0, num - 1);
        std::sort(indexes.data(), indexes.data() + indexes.size(),
        [&pointsVec](int a, int b) {return pointsVec[a].z < pointsVec[b].z; });*/
        std::vector<int> indexes(pointsVec.size());
        std::generate(indexes.begin(), indexes.end(), [n = 0]() mutable {return n++; });
        std::sort(indexes.data(), indexes.data() + indexes.size(),
            [&pointsVec](int a, int b) {return pointsVec[a].z < pointsVec[b].z; });
        /*auto mm = MY1.array().sum() / num;
        auto absDer = (MY1.array() - mm).abs();
        auto r = 2*sqrt(absDer.square().sum() / (num - 1));
        W1 = ((MY1.array() - MY1(indexes(num/2),0)).array() > r).select(0, W1);*/


        auto mm = num >> 1;
        double medianZ = pointsVec[indexes[mm]].z;
        auto meanZ = std::accumulate(pointsVec.begin(), pointsVec.end(), 0.0,
            [](double re, auto& a) { return re += a.z; }) / num;  // use mean avoid outliers.
                                                                  //int tm = static_cast<int>(num / 40);
        if (medianZ < meanZ)
        {
            int low = 0;//tm;
            int high = mm;
            pickRange->push_back(pointsVec[indexes[low]].z);
            pickRange->push_back(2 * pointsVec[indexes[high]].z - pointsVec[indexes[low]].z);
        }
        else
        {
            int low = mm; /*mm+1: may exceed range*/
            int high = num - 1;// -tm;
            pickRange->push_back(2 * pointsVec[indexes[low]].z - pointsVec[indexes[high]].z);
            pickRange->push_back(pointsVec[indexes[high]].z);
        }        
        return true;
    }


    int planeFitting2(vtkPoints* pts, double normal[3], double origin[3], int iterNum = 100)
    {
        if (pts == NULL || pts->GetNumberOfPoints() == 0)
        {
            return 1; // error occur
        }
        std::vector<Point3d> pointsVec;
        vtkPoints2Point3dVec(pts, &pointsVec);
        std::vector<double> pickRange;
        for (size_t i = 0; i < 3; ++i)
        {
            auto tmp = pointsVec.size();
            std::vector<double> tmpRange;
            if (modeFilter(pointsVec, &tmpRange))
            {
                pickRange = tmpRange;
                pointsVec.erase(std::remove_if(pointsVec.begin(), pointsVec.end(), [&tmpRange](auto& p) {
                    return p.z < tmpRange[0] || p.z > tmpRange[1]; }), pointsVec.end());
                double percent = pointsVec.size() * 1.0 / tmp;
                if (percent > .9)
                {
                    break;
                }
            }
        }

        auto num = pointsVec.size();

        MatrixXf MX1(num, 3);
        MatrixXf MY1(num, 1);
        for (size_t i = 0; i < num; ++i)
        {
            MX1(i, 0) = pointsVec[i].x;  MX1(i, 1) = pointsVec[i].y;  MX1(i, 2) = 1;
            MY1(i, 0) = pointsVec[i].z;
        }

        /*MatrixXf W1 = MatrixXf::Zero(num, 1);
        for (size_t i = 0; i < num; ++i)
        {
            if (pointsVec[i].z > pickRange[0] && pointsVec[i].z < pickRange[1])
            {
                W1(i) = 1;
            }
        }
        MatrixXf W(num, 3);
        W << W1, W1, W1;*/


        MatrixXf tmp1 = MX1.transpose()*MX1;//W.cwiseProduct(MX1);
        Vector3f P2 = tmp1.inverse()*(MX1.transpose()*MY1);//W1.cwiseProduct(MY1));
        normal[0] = P2(0); normal[1] = P2(1); normal[2] = -1;
        vtkMath::Normalize(normal);
        origin[0] = std::accumulate(pointsVec.begin(), pointsVec.end(), 0.0, [](double r, const Point3d& a) {return r += a.x; })/num;
        origin[1] = std::accumulate(pointsVec.begin(), pointsVec.end(), 0.0, [](double r, const Point3d& a) {return r += a.y; })/num;
        origin[2] = P2(0)*origin[0] + P2(1)*origin[1] + P2[2];
        return 0;
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


    void heightNormalization(const vector<double>& input, vector<double>* output)
    {
        assert(output != nullptr);
        output->clear();
        output->resize(input.size());
        auto mm = std::minmax_element(input.begin(), input.end());
        double range = *mm.second - *mm.first;
        double scale = 255 / range;
        std::transform(input.begin(), input.end(), output->begin(),
            [scale, &mm](double h) {return floor((h - *mm.first)*scale); });
    }
    void histogram256(const vector<int>& normalizedInput, vector<int>* countOutput)
    {
        assert(countOutput != nullptr);
        countOutput->clear();
        countOutput->resize(256);
        for (auto i:normalizedInput)
        {
            assert(i <= 255);
            (*countOutput)[i] += 1;
        }
    }
    int otsuThreshold(const vector<int>& hist)
    {
        vector<double> w(256), u(256);
        int num = std::accumulate(hist.begin(), hist.end(), 0);
        w[0] = 1.0 * hist[0] / num;
        u[0] = w[0];
        for (int i = 0; i < 256; ++i)
        {
            auto pi = 1.0* hist[i] / num;
            w[i] = w[i - 1] + pi;
            u[i] = u[i - 1] + (i + 1)*pi;
        }
        vector<double> val(255);//255 separation points
        for (int i = 1; i < 256; ++i)
        {
            auto tmp = u[255] * w[i - 1] - u[i - 1];
            val[i - 1] = (tmp*tmp) / (w[i - 1] * (1 - w[i - 1]));
        }
        return std::max_element(val.begin(), val.end()) - val.begin() + 1;
    }
    
    void getMaskArr(const double* data, int widthStep, int width, int height,
        int maskWidth, int maskHeight, int r, int c, vector<double>* re)
    {
        assert(re != nullptr);
        re->reserve(maskHeight*maskWidth);
        int wr = maskWidth >> 1;
        int hr = maskHeight >> 1;
        for (size_t v = -wr; v <= wr; ++v)
        {
            for (size_t u = -hr; u <= hr; ++u)
            {
                int row = r + u;
                row = std::max(0, row);
                row = std::min(row, height - 1);
                int col = c + v;
                col = std::max(0, col);
                col = std::min(col, width - 1);
                int id = row*widthStep + col;
                re->emplace_back(data[id]);
            }
        }
    }
    //void medianOnline()
    //{
    //    vector<double> maskArr;
    //    getMaskArr(zData, widthStep, width, height, maskWidth, maskHeight, r, wRadius, &maskArr);
    //    std::sort(maskArr.begin(), maskArr.end());
    //    dstData[idxBase + wRadius] = maskArr[maskArr.size() >> 1];
    //    for (int c = wRadius+1; c < width - wRadius; ++c)
    //    {
    //        auto remvoePtr = zData + (r - hRadius)*widthStep + c - wRadius - 1;
    //        auto addPtr = removePtr + maskWidth;
    //        for (int k = 0; k < maskHeight; ++k)
    //        {
    //            // the last element equal to remove item
    //            auto rId = std::upper_bound(maskArr.begin(), maskArr.end(), *removePtr) - maskArr.begin() - 1;
    //            // rId must be exist and non-negative
    //            if (rId == -1) rId = 0;
    //            // the insert index for new element
    //            auto iId = std::upper_bound(maskArr.begin(), maskArr.end(), *addPtr) - maskArr.begin() - 1;
    //            if (iId != maskArr.size() - 1) // when no element greater than new input
    //            {
    //                // move the item greater than new input
    //                // when iId == -1, move from 0 to rId
    //                while (rId > iId +1)
    //                {
    //                    maskArr[rId] = maskArr[rId - 1];
    //                    --rId;
    //                }
    //                while (rId < iId)
    //                {
    //                    maskArr[rId] = maskArr[rId + 1];
    //                    ++rId;
    //                }
    //                maskArr[rId] = *addPtr;
    //            }
    //            else
    //            {
    //                while (rId < maskArr.size() - 1)
    //                {
    //                    maskArr[rId] = maskArr[rId + 1];
    //                    ++rId;
    //                }
    //                maskArr[rId] = *addPtr;
    //            }
    //            removePtr += widthStep;
    //            addPtr += widthStep;                                
    //        }
    //        dstData[idxBase + c] = maskArr[maskArr.size() >> 1];
    //    }// end c
    //}
}




#endif // !__ALGORITHM3D_H__

