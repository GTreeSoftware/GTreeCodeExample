/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/
#include <cstdlib>
#include <stdio.h>
#include <algorithm>
#include <numeric>
#include <deque>
#include <iterator>
#include <utility>
#include <Eigen/Core>
#include <Eigen/LU>
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
using namespace google;
#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <ctime>
#else
#include <math.h>
#include <cmath>
#include <sys/time.h>
#endif
#include "sparsetracefilter.h"
#include "../../../ngtypes/tree.h"
#include "../../../ngtypes/soma.h"
#include "../../volumealgo.h"
#include "../../contourutil.h"
#include "../traceutil.h"
#include "../../NGUtility.h"
#include "../Function/IO/imagewriter.h"
#include "../CaliberBuilder.h"
const double EXPO = 0.000001;

SparseTraceFilter::SparseTraceFilter()
{
    className_ = std::string("SparseTraceFilter");
    globalID_ = 1;
    initialDirection_ << 0,0,0;
    m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));
}

SparseTraceFilter::~SparseTraceFilter()
{

}

ProcStatPointer SparseTraceFilter::Update()
{
    if (!m_Input || !m_Back || !m_Bin){
        printf("error occured in %s\n", className_.c_str());
        //LOG(ERROR) << "error occured in " << className_;
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if(m_Input->GetProcessObject()){//|| !m_Soma
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    //2015-6-8

    origImgPointer =
            std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Input);//2015-6-8
    backImgPointer =
            std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Back);//2015-6-8
    binImgPointer =
            std::dynamic_pointer_cast<const Volume<NGCHAR> >(m_Bin);
    std::shared_ptr<const Soma > tmpSoma =
            std::dynamic_pointer_cast<const Soma>(m_Soma);
    std::shared_ptr<TreeCurve > tmpTree =
            std::dynamic_pointer_cast<TreeCurve>(m_Source);

    if(!origImgPointer || !backImgPointer || !binImgPointer || !tmpSoma) {
        printf("error backnoise image!\n");
        LOG(ERROR)<<"error backnoise image!";
        printf("error occured in %s\n", className_.c_str());
        LOG(ERROR)<<"error occured in " <<className_;
        MAKEPROCESSSTATUS(resSta, false, className_, "input data is wrong.");
        return resSta;
    }

    Vec3d soma;
    soma << tmpSoma->GetCell(0).x, tmpSoma->GetCell(0).y, tmpSoma->GetCell(0).z;
    std::vector<VectorVec3d> initialDendrites;
    std::vector<std::vector<VectorVec3d> > treeDataset;
    TraceDendritesConnectedSoma(soma, initialDendrites );
    treeDataset.push_back(initialDendrites);

    size_t traceDeep = 10;//2015-6-8
    for(size_t i = 0; i < traceDeep; ++i){
        //resultDendrites.clear();
        //RecursiveTrace(*tmpOrig, *tmpBack, traceLabelMatrix, tmpDendrites, resultDendrites);
        //2015-6-8
        if(!treeDataset.back().empty()){
            RecursiveTrace(treeDataset);
        }else{
            break;
        }
    }

    //rebuilt into one tree
    std::vector<VectorVec5d> tmpCurves;
    for(size_t i = 0; i < treeDataset.size(); ++i){
        for(size_t j = 0; j < treeDataset[i].size(); ++j){
            tmpCurves.push_back(VectorVec5d());
            VectorVec5d tmpCurve;
            Vec5d tmp;
            tmpCurve.clear();
            for(size_t ij = 0; ij < treeDataset[i][j].size(); ++ij){
                tmp << treeDataset[i][j][ij](0), treeDataset[i][j][ij](1),
                        treeDataset[i][j][ij](2), 1,1;
                tmpCurve.push_back(tmp);
            }
            tmpCurves.back().swap(tmpCurve);
        }
    }
    tmpTree->SetCurve(tmpCurves);

    //destroy
    indexImg_.SetSize(0,0,0);
    traceLabelMatrix_.SetSize(0,0,0);

    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

ConstIDataPointer SparseTraceFilter::GetOutput()
{
    if(!m_Source)
#ifdef _WIN32
        m_Source = std::tr1::shared_ptr<TreeCurve>(new TreeCurve(this));
#else
        m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));
#endif   
    return m_Source;
}

IDataPointer SparseTraceFilter::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData(m_Source);
    m_Source.reset();
    return tData;
}

void SparseTraceFilter::SetInputBack(ConstIDataPointer p)
{
    m_Back = p;
}

void SparseTraceFilter::SetInputBin(ConstIDataPointer p)
{
    m_Bin = p;
}

void SparseTraceFilter::SetSoma(ConstIDataPointer p)
{
    m_Soma = p;
}


//2015-6-8 make use of c++ class function, origImg, backImg, indexImg, traceLabelMatrix are all global varient
void SparseTraceFilter::TraceDendritesConnectedSoma(const Vec3d &soma,
                                                    std::vector<VectorVec3d> &initialDendrites)
{
    const Volume<unsigned short> &origImg = *origImgPointer;
    //const Volume<unsigned short> &backImg = *backImgPointer;

    //initial
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();
    //2015-6-8
    indexImg_.SetSize(nxx, nyy,nzz);//global varient
    indexImg_.SetZero();
    traceLabelMatrix_.SetSize(nxx, nyy, nzz);//global varient
    traceLabelMatrix_.SetZero();
    //get contour and label
    VectorVec3i innerSomaPts;
    ReconstructSomaShapeQuanReviV2(soma, innerSomaPts);//2015-6-8
    LOG(INFO)<<"for";
    for(size_t i = 0; i < innerSomaPts.size(); ++i){
        Vec3i &curPt = innerSomaPts[i];
        if(curPt(0) > -1 && curPt(0) < nxx && curPt(1) > -1 && curPt(1) < nyy
                && curPt(2) > -1 && curPt(2) < nzz){
            traceLabelMatrix_(curPt(0), curPt(1), curPt(2)) = 1;
            indexImg_(curPt(0), curPt(1), curPt(2)) = 1;//2015-6-8
        }
    }
    //2015-6-8 inflat 9 level
    std::vector<VectorVec3i> growPointSet;
    growPointSet.resize(9);
    CVolume growLabelMatrix;//2015-6-8 warning!! below there is a Volume<int> growLabelMatrix
    growLabelMatrix.SetSize(nxx, nyy, nzz);
    for(int i = 0 ; i < nxx; ++i)//2015-6-8
        for(int j = 0 ; j < nyy; ++j)
            for(int ij = 0 ; ij < nzz; ++ij)
                growLabelMatrix(i,j,ij) = traceLabelMatrix_(i,j,ij);
    //growLabelMatrix.QuickCopy(traceLabelMatrix);
    VectorVec3i curPoints(innerSomaPts);
    //size_t loop= isOutOfSomaImg_ ? 6 : 9;
    for(size_t i = 0; i < 9; ++ i){//2015-6-8
        VectorVec3i growPoints;
        RegionInflationModify(curPoints, growLabelMatrix, 35, growPoints);//2016-3-13
        
        curPoints.swap(growPoints);
        if(!curPoints.empty()){
            growPointSet[i] = curPoints;
        } else{
            break;
        }
    }
    /*
     * clustre the grow point of each set. growPointSet has 9 sets.
     * To adapt to C++, I set growLabelMatrix as global varient
    */
    growLabelMatrix.SetZero();
    //all points in each level
    std::vector<VectorVec3i> clustrePtSet;
    //label for clustres in each level
    std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> > clustreLabel;
    ClusterGrowPtSet(growPointSet, growLabelMatrix, clustrePtSet, clustreLabel);
    growLabelMatrix.SetSize(0,0,0);//2015-6-8 destroy growLabelMatrix
    //find the connection between i-level sets and (i + 1)-level sets.
    //there may be some sets which terminate prematurely.
    std::vector<std::vector<int> > subRegionCon;
    RegionConnections(clustrePtSet, clustreLabel, subRegionCon);//2015-6-8 need new label matrix
    /*from end nodes to top nodes. each row in pathMatrix represent
     *the path of 1-th to 9-th node.
     */
    MatXi pathMatrix;
    ConnectPathInTree(subRegionCon, pathMatrix);//2015-6-8

    /*total 9 level. just compare 4-6th level to remove the overlap curve.
     * why not compare 7-9th? 1-6th is long enough. overlap 4-6th level is
     * enough long.
     */
    std::vector<int> redundancyLabel;
    DeleteRedundancyPath(pathMatrix, redundancyLabel);//2015-6-8
    /*
     * according to region connection, get the initial curves from soma.
     * clustreLabel includes the center of each clustre. Connect the center point
     * and get the curves from soma.
    */
    std::vector<VectorVec4d> connectPoints;
    std::vector<VectorVec4d> shortConnectPoints;//it is not use
    ConnectGrowClustreFromSoma(subRegionCon, clustreLabel, connectPoints, shortConnectPoints);//2015-6-8
    LOG(INFO)<<"cool";
    /*
     * just add soma position to initial curves.
    */
    //std::vector<VectorVec3d> initialDendrites; //it is global varient
    //AddSomaToInitDendrites(connectPoints, shortConnectPoints, soma, initialDendrites);
    /*
     * trace 5 steps
    */
    //IdentifyingRayburst(origImg, initialDendrites);
    //size_t initialDendritesNum = initialDendrites.size();
    /*
     * calculate main direction of curve
    */
    //2015-6-8
    //int remainPathNum = std::accumulate(redundancyLabel.begin(), redundancyLabel.end(),0);
    //std::vector<VectorVec3d> initialDendrites;
    VectorVec3d directionSet;
    //directionSet.resize(remainPathNum);
    for(int i = 0; i < pathMatrix.rows(); ++i){
        if(redundancyLabel[i]==1){
            VectorVec3d curDendrite;
            NGUtility::GetReversePartVectorVec3d(connectPoints[i], 0, 8, curDendrite);
            VectorVec3d tmp;
            tmp.push_back(soma);
            std::copy(curDendrite.begin(), curDendrite.begin() + 7, back_inserter(tmp));
            initialDendrites.push_back(VectorVec3d());
            initialDendrites.back().swap(tmp);

            Vec3d curDir;
            CurveDirection(curDendrite, curDir);
            tmp.clear();
            std::copy(curDendrite.begin()+3, curDendrite.end(), back_inserter(tmp));
            CurveDirection(tmp, curDir);
            directionSet.push_back(curDir);
        }
    }
    /*
     * trace more 2015-6-8
    */
    globalID_ = 1;
    std::vector<VectorVec3d> initialDendritesCp;//(initialDendrites);
    VectorVec3d directionSetCp;
    if (initialDirection_.norm() > 0.1) {
        for (size_t i = 0; i < initialDendrites.size(); ++i) {
            if (directionSet[i].dot(initialDirection_) > 0) {
                initialDendritesCp.push_back(initialDendrites[i]);
                directionSetCp.push_back(directionSet[i]);
            }
        }
    } else {//all direction
        initialDendritesCp = initialDendrites;
        directionSetCp = directionSet;
    }
    LOG(INFO)<<"cool";
    ContrainedPCATraceSomaBifur(initialDendritesCp, directionSetCp, initialDendrites);
}

/*
2015-6-8 origImg, backImg, traceLabelMatrix, indexImg is global varient
*/
void SparseTraceFilter::RecursiveTrace(std::vector<std::vector<VectorVec3d> > &treeDataSet)
{
    /*initialize
    * const Volume<unsigned short> &origImg = *origImgPointer;
     *const Volume<unsigned short> &backImg = *backImgPointer;
     *std::vector<VectorVec3d> resultDendrites;
    */
    const std::vector<VectorVec3d> &tmpDendrites = treeDataSet.back();
    size_t tmpDendritesNum = tmpDendrites.size();
    std::vector<VectorVec3d> sortClustre;
    sortClustre.reserve(tmpDendritesNum);
    for(size_t i = 0; i < tmpDendritesNum; ++i){
        if(!tmpDendrites[i].empty())
            sortClustre.push_back(tmpDendrites[i]);//cannot swap , cuz tmpDendrites is constant.
    }
    /*
     * detect bifurcation
    */
    VectorVec3d bifurcationDirections;
    std::vector<VectorVec3d> bifurcationClustre;
    CVolume traceLabelMatrixCp;
    traceLabelMatrixCp.QuickCopy(traceLabelMatrix_);
    //2015-6-8
    VectorVec3i boundaryLabel;
    DetectBifurcation(traceLabelMatrixCp, sortClustre, bifurcationDirections, bifurcationClustre, boundaryLabel);

    /*
    * remove cross curves
    */


    /*
     * according to bifurcation set , traverse to trace the dendrites.
    */
    if(!bifurcationDirections.empty()){
        ContrainedPCATraceBifur(bifurcationClustre, bifurcationDirections, treeDataSet);
    }else{
        treeDataSet.push_back(std::vector<VectorVec3d>());
    }

    //2015-6-8
    if(!boundaryLabel.empty()){
        for(size_t i = 0; i < boundaryLabel.size(); ++i){
            traceLabelMatrix_(boundaryLabel[i](0),boundaryLabel[i](1),boundaryLabel[i](2)) = 1;
        }
    }
}

//2015-6-8
void SparseTraceFilter::ReconstructSomaShapeQuanReviV2(const Vec3d &initialPoint,
                                                     VectorVec3i &innerSomaPts)
{
    //2015-6-8
    const Volume<unsigned short> &origImg = *origImgPointer;
    const Volume<unsigned short> &backImg = *backImgPointer;
    //
    innerSomaPts.clear();
    SVolume locOrigImg;
    SVolume locBackImg;
    Vec3d locPoint;
    const int Theta = 40;
    const int Phi = 20;
    ExtractLocalDomainV2(initialPoint, origImg, locOrigImg, locPoint);//2015-6-8
    ExtractLocalDomainV2(initialPoint, backImg, locBackImg, locPoint);//2015-6-8
    double slice = 0.5;
    const int blocksum = 81;//2016-3-13

    std::vector<double> raySlice;
    for(int i = 0; i < blocksum; ++i){
        raySlice.push_back(slice * double(i));
    }

    /**/
    Volume<double> sphereRayWet;
    sphereRayWet.SetSize(int(raySlice.size()), Theta, Phi);
    Volume<double> sphereBackWet;
    sphereBackWet.SetSize(int(raySlice.size()), Theta, Phi);

    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;

    double rayNodeDense(0.0), rayNodeWet(0.0);
    double x,y,z;

    for (std::vector<double>::size_type k = 0; k < raySlice.size(); ++k){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                x = raySlice[k] * std::sin(b * (double)j * PI_180) * std::cos(a * (double)i * PI_180) + locPoint(0);
                y = raySlice[k] * std::sin(b * (double)j * PI_180) * std::sin(a * (double)i * PI_180) + locPoint(1);
                z = raySlice[k] * std::cos(b * (double)j * PI_180) + locPoint(2);
                
                rayNodeDense = rayNodeWet = 0.0;
                ContourUtil::CalculateSphereOneNode(locOrigImg, 1.0, x, y, z, rayNodeDense, rayNodeWet);
                sphereRayWet(int(k), i-1, j-1) = (double)(rayNodeDense / (rayNodeWet + 0.0001));

                //2015-6-8
                rayNodeDense = rayNodeWet = 0.0;
                ContourUtil::CalculateSphereOneNode(locBackImg, 1.0, x, y, z, rayNodeDense, rayNodeWet);
                sphereBackWet(int(k), i-1, j-1) = (double)(rayNodeDense / (rayNodeWet + 0.0001));
            }
        }
    }//for

    std::vector<std::vector<double> > rayLimit;
    GetRayLimitV2(sphereRayWet, sphereBackWet, 8, rayLimit);

    Volume<double> smoothRay;

    TraceUtil::GetGradientVectorFlowForTrace(sphereRayWet, smoothRay);

    double gradThrev(0.0);//max value of 
    GetAreaMaxValue(smoothRay, 0,  smoothRay.x() - 1, 0, smoothRay.y() - 1, 0, smoothRay.z() - 1, gradThrev);
    gradThrev *= 0.1;
    MatXd resultRayLength;//2015-2-7
    resultRayLength = 5.0 * MatXd::Ones(Theta + 2, Phi + 2);//ray length
    for(size_t i = 0; i < rayLimit.size();++i){
        for(size_t j = 0; j < rayLimit[0].size(); ++j){
            resultRayLength(i+1,j+1)=rayLimit[i][j];
        }
    }
    MatXd AAss = resultRayLength;//2015-6-8

    std::vector<double> raySliceIndex;
    for(int i = 0; i < blocksum; ++i){
        raySliceIndex.push_back(1.0 + double(i));
    }

    double currentRayLength(0);
    std::vector<double> currentRayGrad(smoothRay.x());
    std::vector<double> rayGradWet(smoothRay.x());

    double pullParm[2];
    pullParm[0] = 0.8;
    pullParm[1] = 1.0 - pullParm[0];

    int repeat = 20;//2015-6-8
    double rayGradWetSum;
    LOG(INFO)<<repeat<<" "<<Theta<<" "<<Phi<<" "<<smoothRay.x()<<" "<<rayGradWet.size();
    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                currentRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    currentRayGrad[ij] = smoothRay(ij, i - 1, j - 1);
                    rayGradWet[ij] = currentRayGrad[ij]
                            * std::exp( -0.05 *std::pow( raySliceIndex[ij] - currentRayLength, 2.0));
                }
                //rayGradWetSum = std::accumulate(rayGradWet.begin(), rayGradWet.end(), 0.0);
                rayGradWetSum = 0.0;
                for(std::vector<double>::size_type k = 0; k < rayGradWet.size(); ++k){
                    rayGradWetSum += rayGradWet[k];
                }
                //if (std::abs(rayGradWetSum) > 0.001){
                resultRayLength(i,j) += pullParm[0] * ( ( (std::inner_product(raySliceIndex.begin(), raySliceIndex.end(), rayGradWet.begin(), 0.0)) /
                                                          (rayGradWetSum ) ) - resultRayLength(i,j) )//2015-6-8+0.001
                        - pullParm[1] * (4.0 * resultRayLength(i,j)
                                         - resultRayLength(i-1,j) - resultRayLength(i+1,j) - resultRayLength(i,j-1) - resultRayLength(i,j+1));
                //}

            }
            //resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength = resultRayLength.cwiseMax(AAss);//2015-6-8
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for
    LOG(INFO)<<"for";
    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                currentRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    currentRayGrad[ij] = smoothRay(ij, i - 1, j - 1);
                    rayGradWet[ij] = currentRayGrad[ij] * std::exp( -0.05 * std::pow( raySliceIndex[ij] - currentRayLength, 2.0));
                }

                rayGradWetSum = 0.0;
                for(std::vector<double>::size_type k = 0; k < rayGradWet.size(); ++k){
                    rayGradWetSum += rayGradWet[k];
                }
                if ( smoothRay( (std::min)(NGUtility::Round(currentRayLength + 0.5) -1, (int)raySlice.size() - 1), i - 1, j - 1 ) >  gradThrev ){
                    resultRayLength(i,j) += pullParm[0] * ( ( (inner_product(raySliceIndex.begin(), raySliceIndex.end(), rayGradWet.begin(), 0.0)) /
                                                              (rayGradWetSum ) ) - resultRayLength(i,j) )//2015-6-8+0.001
                            - pullParm[1] * (4.0 * resultRayLength(i,j) - resultRayLength(i-1,j) - resultRayLength(i+1,j) - resultRayLength(i,j-1) - resultRayLength(i,j+1));
                }
                else{
                    resultRayLength(i,j) += (-pullParm[1]) * (4.0 * resultRayLength(i,j) - resultRayLength(i-1,j) - resultRayLength(i+1,j) - resultRayLength(i,j-1) - resultRayLength(i,j+1));
                }
            }
            //resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength = resultRayLength.cwiseMax(AAss);//2015-6-8
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for

    /*
     * 2015-2-5 new inner point collect method
    */
    Vec3d curPt, curNormPt;
    double curPtLength;
    double altitudeAngle, azimuthAngle;
    int altitudeIndicator, azimuthIndicator;
    double ray0, ray1, ray2, ray3, ray4, rayMean;
    Vec3i soma;
    soma << NGUtility::Round(initialPoint(0)), NGUtility::Round(initialPoint(1)), NGUtility::Round(initialPoint(2));
    //2015-6-8
    for(int i = -40; i <= 40; ++i){
        for(int j = -40; j <= 40; ++j){
            for(int ij = -40; ij <= 40; ++ij){
                curPt << i,j,ij;
                curPtLength = curPt.norm();
                if(curPtLength > 0.01){
                    curNormPt = curPt.normalized();
                    CalVectorAngle(curNormPt, altitudeAngle, azimuthAngle);
                    altitudeIndicator=std::min(NGUtility::Round(altitudeAngle*180.0/M_PI/b+0.5),Phi);//+1;
                    azimuthIndicator=std::min(NGUtility::Round(azimuthAngle*180.0/M_PI/a+0.5),Theta);//+1;//C++ minus 1
                    ray0 = resultRayLength(azimuthIndicator,altitudeIndicator)*0.5;
                    ray1 = resultRayLength(azimuthIndicator+1,altitudeIndicator)*0.5;
                    ray2 = resultRayLength(azimuthIndicator-1,altitudeIndicator)*0.5;
                    ray3 = resultRayLength(azimuthIndicator,altitudeIndicator+1)*0.5;
                    ray4 = resultRayLength(azimuthIndicator,altitudeIndicator-1)*0.5;
                    rayMean = 0.2 * (ray0 + ray1 + ray2 + ray3 + ray4);
                    if(curPtLength < rayMean + 0.5)
                        innerSomaPts.push_back(Vec3i(i,j,ij) + soma);
                }else{
                    innerSomaPts.push_back(Vec3i(i,j,ij) + soma);//2015-6-8
                }
            }
        }
    }
    {
        VectorVec3i tmp(innerSomaPts);
        innerSomaPts.swap(tmp);
    }
}

void SparseTraceFilter::RegionInflationModify(const VectorVec3i &curPoints,
                                              CVolume &growLabelMatrix, double threv,
                                              VectorVec3i &growPoints)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int xBdy = origImg.x() - 1;
    int yBdy = origImg.y() - 1;
    int zBdy = origImg.z() - 1;

    size_t nxx = curPoints.size();
    growPoints.clear();
    /*
     * area grow
    */
    bool binImgLabel;
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VectorVec3i tmpGrowPts;
    for(size_t i = 0; i < nxx; ++i){

        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPoints[i](0), curPoints[i](1), curPoints[i](2),
                    2,2,1, 0, xBdy, 0, yBdy, 0, zBdy);//2015-2-8
                tmpGrowPts.clear();
        tmpGrowPts.reserve(75);
        //if(!isOutOfSomaImg_){//in soma
        for(int ii = xMin; ii <= xMax; ++ii){
            for(int jj = yMin; jj <= yMax; ++jj){
                for(int ij = zMin; ij <= zMax; ++ij){
                    binImgLabel = double(origImg(ii,jj,ij)) >
                            double(backImg(ii,jj,ij)) + threv * std::sqrt(double(backImg(ii,jj,ij)));
                    if(0 == growLabelMatrix(ii,jj,ij) && (binImgLabel))
                        tmpGrowPts.push_back(Vec3i(ii,jj,ij));
                }
            }
        }//for
        //}else{//2015-6-16
        //    for(int ii = xMin; ii <= xMax; ++ii){
        //        for(int jj = yMin; jj <= yMax; ++jj){
        //            for(int ij = zMin; ij <= zMax; ++ij){
        //                binImgLabel = double(origImg(ii,jj,ij)) >
        //                        double(backImg(ii,jj,ij)) + threv * std::sqrt(double(backImg(ii,jj,ij)));
        //                bool lishiweiLabel = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + diffuseValue_;//2016-3-13 diffuseValue_
        //                if(0 == growLabelMatrix(ii,jj,ij) && (binImgLabel || lishiweiLabel))
        //                    tmpGrowPts.push_back(Vec3i(ii,jj,ij));
        //            }
        //        }
        //    }//for
        //}
        if(tmpGrowPts.size() > 2){
            std::copy(tmpGrowPts.begin(), tmpGrowPts.end(), std::back_inserter(growPoints));
            for(size_t j = 0; j < tmpGrowPts.size(); ++j)
                growLabelMatrix(tmpGrowPts[j](0), tmpGrowPts[j](1), tmpGrowPts[j](2)) = 1;
        }
    }
}

void SparseTraceFilter::RegionInflationModifyV2(const VectorVec3i &curPoints,
                                              CVolume &growLabelMatrix, double threv,
                                              VectorVec3i &growPoints)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int xBdy = origImg.x() - 1;
    int yBdy = origImg.y() - 1;
    int zBdy = origImg.z() - 1;

    size_t nxx = curPoints.size();
    growPoints.clear();
    /*
     * area grow
    */
    bool flag1;
    bool flag2;
    bool flag3;
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VectorVec3i tmpGrowPts;
    for(size_t i = 0; i < nxx; ++i){

        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPoints[i](0), curPoints[i](1), curPoints[i](2),
                    2,2,1, 0, xBdy, 0, yBdy, 0, zBdy);//2015-2-8
                tmpGrowPts.clear();
        tmpGrowPts.reserve(75);
        for(int ii = xMin; ii <= xMax; ++ii){
            for(int jj = yMin; jj <= yMax; ++jj){
                for(int ij = zMin; ij <= zMax; ++ij){
                    flag1 = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + threv * std::sqrt(double(backImg(ii,jj,ij)));
                    flag2 = double(origImg(ii,jj,ij)) > 1.5*double(backImg(ii,jj,ij));
                    flag3 = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + paramPack->diffuseValue_;
                    if(0 == growLabelMatrix(ii,jj,ij) && (flag1 || flag2 || flag3))
                        tmpGrowPts.push_back(Vec3i(ii,jj,ij));
                }
            }
        }//for
        if(tmpGrowPts.size() > 2){
            std::copy(tmpGrowPts.begin(), tmpGrowPts.end(), std::back_inserter(growPoints));
            for(size_t j = 0; j < tmpGrowPts.size(); ++j)
                growLabelMatrix(tmpGrowPts[j](0), tmpGrowPts[j](1), tmpGrowPts[j](2)) = 1;
        }
    }
}

void SparseTraceFilter::ClusterGrowPtSet(const std::vector<VectorVec3i> &growPointSet, CVolume &growLabelMatrix,
                                         std::vector<VectorVec3i> &clustrePtSet,
                                         std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> > &clustreLabel)
{
    //initialize
    size_t nxx = growPointSet.size();
    clustreLabel.clear();
    clustreLabel.resize(nxx);
    clustrePtSet = growPointSet;
    /*
     * clustre each level
    */
    VectorVec3i inflatPtSet;
    VectorVec3i curConnectPt;
    Vec3d curCenter, curStd;
    for(size_t j = 0; j < nxx; ++j){
        inflatPtSet=clustrePtSet[j];
        //growLabelMatrix is initialized in function ClustersExtractionPoints
        VectorVec3i subClustrePtSet;
        std::vector<int> subClustrePtSetNum;
        ClustersExtractionPoints(inflatPtSet, growLabelMatrix, subClustrePtSet, subClustrePtSetNum);
        std::vector<int> subClustrePtSetCumNum(subClustrePtSetNum.size()+1);
        subClustrePtSetCumNum[0] = 0;
        std::partial_sum(subClustrePtSetNum.begin(), subClustrePtSetNum.end(), subClustrePtSetCumNum.begin()+1);
        size_t subClustrePtSetNumLength = subClustrePtSetNum.size();
        std::vector<SparseTraceFilter::CLUSTRELABEL> curPtSetFeature(subClustrePtSetNumLength);
        for(size_t i = 0; i < curPtSetFeature.size(); ++i)
            curPtSetFeature[i].index = subClustrePtSetNum[i];
        for(size_t i = 0; i < subClustrePtSetNumLength; ++i){
            curConnectPt.clear();
            std::copy(subClustrePtSet.begin() + subClustrePtSetCumNum[i],
                      subClustrePtSet.begin() + subClustrePtSetCumNum[i + 1],
                    std::back_inserter(curConnectPt));
            //the same as isempty function
            if(subClustrePtSetCumNum[i + 1] - 1 > subClustrePtSetCumNum[i]){
                curCenter.setZero();
                // get center
                for(size_t ij = 0; ij < curConnectPt.size(); ++ij){
                    curCenter += Vec3d(curConnectPt[ij](0), curConnectPt[ij](1), curConnectPt[ij](2));
                }
                curCenter /= double(curConnectPt.size());
                // get standard deviation
                curStd.setZero();
                for(size_t ij = 0; ij < curConnectPt.size(); ++ij){
                    curStd += (Vec3d(curConnectPt[ij](0), curConnectPt[ij](1), curConnectPt[ij](2))
                                     - curCenter).array().square().matrix();
                }
                curStd = curStd.array() / double(curConnectPt.size() - 1);
                curStd = curStd.array().sqrt();
                // save to feature
                curPtSetFeature[i].center = curCenter;
                curPtSetFeature[i].standardDeviation = curStd;
            } else{
                curPtSetFeature[i].center = Vec3d(curConnectPt[0](0), curConnectPt[0](1),
                        curConnectPt[0](2));
                curPtSetFeature[i].standardDeviation.setZero();
            }
        }
        clustrePtSet[j].swap(subClustrePtSet);
        clustreLabel[j].swap(curPtSetFeature);
    }
}

void SparseTraceFilter::RegionConnections(const std::vector<VectorVec3i> &clustrePtSet,
                                          const std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> > &clustreLabel,
                                          std::vector<std::vector<int> > &clustreConInfo)
{
    //initial
    size_t clustreNum = clustreLabel.size();
    clustreNum = std::min(clustreNum, size_t(10));
    const SVolume &origImg = *origImgPointer;
    Volume<int> growLabelMatrix;//2015-6-8 create a new growLabelMatrix of Volume<int>
    growLabelMatrix.SetSize(origImg.x(), origImg.y(), origImg.z());
    growLabelMatrix.SetZero();
    //
    VectorVec3i curConPtSet;
    int curConPtSetFlag;
    for(size_t i = 0; i < clustreNum; ++i){
        const VectorVec3i& curClustrePtSet = clustrePtSet[i];
        const std::vector<SparseTraceFilter::CLUSTRELABEL>& curClustreFeature = clustreLabel[i];
        std::vector<int> tmpList;
        for(size_t j = 0; j < curClustreFeature.size(); ++j)
            tmpList.push_back(curClustreFeature[j].index);
        std::vector<int> curClustreConNum(curClustreFeature.size()+1);
        curClustreConNum[0] = 0;
        std::partial_sum(tmpList.begin(), tmpList.end(), curClustreConNum.begin()+1);
        for(size_t j = 0; j < curClustreConNum.size() - 1; ++j){
            curConPtSet.clear();
            std::copy(curClustrePtSet.begin() + curClustreConNum[j],
                      curClustrePtSet.begin() + curClustreConNum[j + 1],
                    std::back_inserter(curConPtSet));
            //just generate a unique global ID. C++ must plus 1 to avoid 0 ID.
            curConPtSetFlag = int(i+1) * 10000 + int(j)+1;//2015-6-8
            for(size_t mm = 0; mm < curConPtSet.size(); ++mm)
                growLabelMatrix(curConPtSet[mm](0), curConPtSet[mm](1), curConPtSet[mm](2)) = curConPtSetFlag;
        }
    }
    //check collide area ID
    int connectionFlag = 0;
    clustreConInfo.resize(clustreNum - 1);
    std::vector<int> curClustreConInfo;
    LOG(INFO)<<clustreNum;
    for(size_t i = 1; i < clustreNum; ++i){
        connectionFlag = 0;
        const VectorVec3i& curClustrePtSet = clustrePtSet[i];
        const std::vector<SparseTraceFilter::CLUSTRELABEL>& curClustreFeature = clustreLabel[i];
        std::vector<int> tmpList;
        for(size_t j = 0; j < curClustreFeature.size(); ++j)
            tmpList.push_back(curClustreFeature[j].index);
        std::vector<int> curClustreConNum(curClustreFeature.size()+1);
        curClustreConNum[0] = 0;
        std::partial_sum(tmpList.begin(), tmpList.end(), curClustreConNum.begin()+1);
        curClustreConInfo.clear();
        curClustreConInfo.resize(curClustreConNum.size() - 1, 0);
        for(size_t j = 0; j < curClustreConNum.size() - 1; ++j){
            curConPtSet.clear();
            std::copy(curClustrePtSet.begin() + curClustreConNum[j],
                      curClustrePtSet.begin() + curClustreConNum[j + 1],
                    std::back_inserter(curConPtSet));
            // check region connection info
            // region id has plused one
            RegionConnectionsSub(curConPtSet, growLabelMatrix, int(i), connectionFlag);
            curClustreConInfo[j] = connectionFlag;
//            if(connectionFlag < 0)
//                printf("nia\n");
        }
//        if(curClustreConInfo.size() > 1 && i-1 == 4  )
//            if(curClustreConInfo[1] < 0)
//                printf("nia\n");

        clustreConInfo[i - 1].swap(curClustreConInfo);
    }
}

//2015-6-8
void SparseTraceFilter::ConnectGrowClustreFromSoma(const std::vector<std::vector<int> > &subRegionCon,
                                                   const std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> > &clustreLabel,
                                                   std::vector<VectorVec4d> &connectPoints,
                                                   std::vector<VectorVec4d> &shortConnectPoints)
{
    size_t nxx = clustreLabel.size();
    /*
     * clustreLabel is consist of index, three dimension coordinate of center,
     * three dimension coordinate of standard deviation, of all 7 elements.
    */
    std::vector<SparseTraceFilter::CLUSTRELABEL> terminalPtFeature = clustreLabel.back();
    size_t endPtNum = terminalPtFeature.size();
    //nxx is the depth level of clutre grow sets.
    /*
     * endPtNum is the num of 6-th level clustre grow sets.
     * connectPoints is consist of three dimension coordinate of center and a flag value of
     * current center point.
    */
    connectPoints.clear();
    connectPoints.resize(endPtNum);
    for(size_t i = 0; i < endPtNum; ++i){
        connectPoints[i].resize(nxx);
    }
    //traverse  all end Points and connect the points.
    size_t curRecursiveIndex;
    for(size_t j = 0; j < endPtNum; ++j){
        const Vec3d& curCenterPt = terminalPtFeature[j].center;
        //initialize the end node.
        connectPoints[j][0] << curCenterPt(0), curCenterPt(1), curCenterPt(2), 1.0;
        curRecursiveIndex = j;
        //warning!!! C++ index
        for(size_t i = 0; i < nxx - 1; ++i){
            curRecursiveIndex = subRegionCon[nxx - i - 2][curRecursiveIndex] - 1;
            connectPoints[j][i +1].block(0,0,3,1) = clustreLabel[nxx - i -2][curRecursiveIndex].center;
            const Vec3d &maxCoordinate = clustreLabel[nxx - i -2][curRecursiveIndex].standardDeviation;
            if(maxCoordinate.maxCoeff() < 15.1)//2015-6-8 3.1
                connectPoints[j][i +1](3) = 1.0;
            else
                connectPoints[j][i +1](3) = 0.0;
        }
    }
    shortConnectPoints.clear();//2015-6-8
    /*
     * nxx -1 = 5, the size of subRegionCon is 5 while the size of clustreLabel is 6.
    */
//    std::vector<int> shortTraversalLabel(clustreLabel[nxx -2].size());
//    for(size_t i = 0; i < subRegionCon[nxx - 2].size(); ++i){
//        shortTraversalLabel[subRegionCon[nxx -2][i] - 1] = 1;
//    }
//    //the second traversal begins in 5-th level.
//    //find the 5-th level end points which have no connect points in 6-th level.
//    endPtNum = 0;
//    for(size_t i = 0; i < shortTraversalLabel.size(); ++i){
//        if(shortTraversalLabel[i] == 0) ++endPtNum;
//    }
//    shortConnectPoints.clear();
//    shortConnectPoints.reserve(endPtNum);
////    for(size_t i = 0; i < endPtNum; ++i){
////        shortConnectPoints[i].reserve(nxx - 1);
////    }
//    terminalPtFeature = clustreLabel[nxx-2];
//    for(size_t j = 0; j < terminalPtFeature.size(); ++j){
//        if(0 == shortTraversalLabel[j]){
//            shortConnectPoints.push_back(VectorVec4d());
//            const Vec3d& curCenterPt = terminalPtFeature[j].center;
//            shortConnectPoints.back().push_back(Vec4d(curCenterPt(0), curCenterPt(1), curCenterPt(2), 1.0));
//            curRecursiveIndex = j;
//            for(size_t i = 0; i < nxx -2; ++i){//because from 5-th level
//                curRecursiveIndex = subRegionCon[nxx - i - 3][curRecursiveIndex] - 1;
//                shortConnectPoints.back().push_back(Vec4d(0,0,0,0));
//                shortConnectPoints.back().back().block(0,0,3,1) = clustreLabel[nxx - i -3][curRecursiveIndex].center;
//                const Vec3d &maxCoordinate = clustreLabel[nxx - i -3][curRecursiveIndex].standardDeviation;
//                if(maxCoordinate.maxCoeff() < 3.1)
//                    shortConnectPoints.back().back()(3) = 1.0;
//                else
//                    shortConnectPoints.back().back()(3) = 0.0;
//            }
//        }
//    }
}

void SparseTraceFilter::AddSomaToInitDendrites(const std::vector<VectorVec4d> &connectPoints,
                                               const std::vector<VectorVec4d> &shortConnectPoints,
                                               const Vec3d &soma, std::vector<VectorVec3d> &initDendrites)
{
    //initialize
    initDendrites.clear();
    initDendrites.resize(connectPoints.size() + shortConnectPoints.size());
    //connect connectPoints
    //VectorVec4d curDendrites;
    for(size_t i = 0; i < connectPoints.size(); ++i){
        const VectorVec4d &curDendrites = connectPoints[i];
        int idexx = -1;
        for(size_t j = 0; j < connectPoints[i].size(); ++j){
            if(connectPoints[i][j](3) < 0.5){//==0
                idexx = int(j);
                break;
            }
        }
        if(-1 == idexx){//all curve node is usable
            initDendrites[i].push_back(soma);
            size_t len = curDendrites.size();
            for(size_t ij = 0; ij < len; ++ij){
                initDendrites[i].push_back(connectPoints[i][len - ij - 1].block(0,0,3,1));
            }
        } else{
            initDendrites[i].push_back(soma);
            //int len = int(curDendrites.size()) - 1;
            for(int ij = idexx - 1; ij > -1; --ij){//cuz idexx maybe -1, donot size_t(idexx)
                initDendrites[i].push_back(connectPoints[i][ij].block(0,0,3,1));
            }
        }
    }
    size_t curNum = connectPoints.size();
    //connect shortConnectPoints
    for(size_t i = 0; i < shortConnectPoints.size(); ++i){
        const VectorVec4d &curDendrites = shortConnectPoints[i];
        int idexx = -1;
        for(size_t j = 0; j < shortConnectPoints[i].size(); ++j){
            if(shortConnectPoints[i][j](3) < 0.5){//==0
                idexx = int(j);
                break;
            }
        }
        if(-1 == idexx){//all curve node is usable
            initDendrites[i + curNum].push_back(soma);
            size_t len = curDendrites.size();
            for(size_t ij = 0; ij < len; ++ij){
                initDendrites[i + curNum].push_back(curDendrites[len - ij - 1].block(0,0,3,1));
            }
        } else{
            initDendrites[i + curNum].push_back(soma);
            //int len = int(curDendrites.size()) - 1;
            for(int ij = idexx - 1; ij > -1; --ij){
                initDendrites[i + curNum].push_back(curDendrites[ij].block(0,0,3,1));
            }
        }
    }
}

void SparseTraceFilter::IdentifyingRayburst(const Volume<unsigned short> &origImg, std::vector<VectorVec3d> &initialDendrites)
{
    size_t initialDendritesNum = initialDendrites.size();
    LOG(INFO)<<initialDendritesNum;
    for(size_t i = 0; i < initialDendritesNum; ++i){
        const VectorVec3d& curInitDendrite = initialDendrites[i];
        size_t curInitDendriteNum = curInitDendrite.size();
        VectorVec3d forwardArea;
        RayburstShapeTrackV2(curInitDendrite[curInitDendriteNum - 2], curInitDendrite.back(),
                origImg, 5, forwardArea);
        if(!forwardArea.empty()){
            std::copy(forwardArea.begin(), forwardArea.end(), std::back_inserter(initialDendrites[i]));
        }
    }
}

void SparseTraceFilter::CurveDirection(const VectorVec3d &curDendrite, Vec3d &curDir)
{
    size_t nxx = curDendrite.size();
    curDir << 0.0,0.0,0.0;
    if(nxx > 2){
        for(size_t i = 0; i < nxx - 1; ++i){
            curDir += curDendrite[i + 1] - curDendrite[0];
        }
        curDir /= std::max(curDir.norm(), 0.001);
    }
}

//2015-6-8 threValue,origImg, backImg, globalID is global varient
void SparseTraceFilter::ContrainedPCATraceSomaBifur(const std::vector<VectorVec3d> &initialDendritesCp,
                                                const VectorVec3d &directionSet,
                                                std::vector<VectorVec3d> &resultDendrites)
{
    //initial 2015-6-8
    //const SVolume &origImg = *origImgPointer;
    //const SVolume &backImg = *backImgPointer;
    //int collideId;
    size_t initialDendritesCpNum = initialDendritesCp.size();
    resultDendrites.clear();
    resultDendrites.resize(initialDendritesCpNum);
    VectorVec5d traceDendrite;
    //
    LOG(INFO)<<initialDendritesCpNum;
    for(size_t i = 0; i < initialDendritesCpNum; ++i){
        const VectorVec3d& curDendrite = initialDendritesCp[i];
        const Vec3d& seedPoint = curDendrite.back();
        const Vec3d& curDirection = directionSet[i];
        if(curDirection.array().abs().sum() > 0.0){
            traceDendrite.clear();
            //2017-12-5
            //MeanShiftGrowTraceTotal(seedPoint, curDirection, origImg, backImg, traceDendrite, collideId);
            //2015-6-8
            TraceCurveForwardFromSomaV2(seedPoint, curDirection, traceDendrite);
            if(!traceDendrite.empty()){
                resultDendrites[i].resize(curDendrite.size() + traceDendrite.size());
                std::copy(curDendrite.begin(), curDendrite.end(), resultDendrites[i].begin());
                //std::copy(traceDendrite.begin(), traceDendrite.end(), resultDendrites[i].begin() + curDendrite.size());
                for(size_t j = 0; j < traceDendrite.size(); ++j){
                    resultDendrites[i][j + curDendrite.size()] = traceDendrite[j].block(0,0,3,1);
                }
                ++globalID_;
                //2017-12-4
                CubeLabelL1Opt(resultDendrites[i]);
                //2016-3-13
                //CubeLabelSoma(resultDendrites[i]);
                //CubeLabel(traceLabelMatrix, origImg, backImg, resultDendrites[i], globalID, indexImg);
                //2015-5-7
                //VectorVec3d tmpDend;
                //ThreeDimdCurveResample(resultDendrites[i], 2.0, tmpDend);
                //resultDendrites[i].swap(tmpDend);
            } else{//2015-5-7
                resultDendrites[i].clear();
                ++globalID_;
                //2017-12-4
                //CubeLabelL1Opt(curDendrite);
                CubeLabelSoma(curDendrite);//2016-3-13
                //CubeLabel(traceLabelMatrix, origImg, backImg, curDendrite, globalID, indexImg);
            }
        }
    }
}

//2015-6-8 threValue,origImg, backImg, globalID is global variant
void SparseTraceFilter::ContrainedPCATraceBifur(const std::vector<VectorVec3d> &initialDendritesCp,
                                                const VectorVec3d &directionSet,
                                                std::vector<std::vector<VectorVec3d> > &treeDataSet)
{
    //initial 2015-6-8
    size_t initialDendritesCpNum = initialDendritesCp.size();
    treeDataSet.push_back(std::vector<VectorVec3d>());
    VectorVec3d resultDendrites;
    resultDendrites.clear();
    //resultDendrites.resize(initialDendritesCpNum);
    VectorVec5d traceDendrite;
    //
    //LOG(INFO)<<initialDendritesCpNum;
    //int collideIDnouse;
    for(size_t i = 0; i < initialDendritesCpNum; ++i){
        const VectorVec3d& curDendrite = initialDendritesCp[i];
        const Vec3d& seedPoint = curDendrite.back();
        const Vec3d& curDirection = directionSet[i];
        if(curDirection.array().abs().sum() > 0.0){
            traceDendrite.clear();
            //2017-12-5
            //MeanShiftGrowTraceTotal(seedPoint, curDirection, *origImgPointer, *backImgPointer, traceDendrite, collideIDnouse);
            //2015-6-8
            TraceCurveForwardFromSomaV2(seedPoint, curDirection, traceDendrite);
            if(!traceDendrite.empty()){
                //2015-6-8
                resultDendrites.resize(curDendrite.size() + traceDendrite.size());
                std::copy(curDendrite.begin(), curDendrite.end(), resultDendrites.begin());
                for(size_t j = 0; j < traceDendrite.size(); ++j){
                    resultDendrites[j + curDendrite.size()] = traceDendrite[j].block(0,0,3,1);
                }
                //treeDataSet.back().push_back(VectorVec3d());
                //treeDataSet.back().back().swap(resultDendrites);
                //2015-6-8 add reverse connect parent dendrites
                VectorVec3d reverseHead;
                NGUtility::GetReversePartVectorVec3d(resultDendrites, 0, std::min(5, int(resultDendrites.size() - 1)), reverseHead);
                Vec3d reverseDirection;
                CalcPtCurveDirection(reverseHead, reverseDirection);
                Vec3d seedPointNew = curDendrite[0];
                VectorVec5d reverseCurve;
                int collideID;
                //change to mean shift trace method
                {
                    /* IDataPointer tmpIndex = std::make_shared<SVolume>();
                     NG_CREATE_DYNAMIC_CAST(SVolume, tmp, tmpIndex);
                     tmp->SetSize(indexImg_.x(), indexImg_.y(), indexImg_.z());
                     for (int x = 0; x < indexImg_.x(); ++x) {
                     for (int y = 0; y < indexImg_.y(); ++y) {
                     for (int z = 0; z < indexImg_.z(); ++z) {
                     if (indexImg_(x,y,z) != 0) {
                     tmp->operator()(x, y, z) = 255;
                     }
                     }
                     }
                     }
                     NGImageWriter writer = ImageWriter::New();
                     writer->SetOutputFileName("F:/wocao.tif");
                     writer->SetInput(tmpIndex);
                     writer->Update();
                     system("pause");*/
                }
                //TraceCurveForwardFromSomaV3(seedPointNew, reverseDirection, reverseCurve, collideID, 15);
                MeanShiftGrowTraceTotal(seedPointNew, reverseDirection, *origImgPointer, *backImgPointer, reverseCurve, collideID);

                if(collideID < 1){
                    //collide nothing
                    /*VectorVec3d tmp;
                    NGUtility::GetReversePartVectorVec3d(reverseCurve, 0, int(reverseCurve.size()) - 1, tmp);
                    tmp.resize(reverseCurve.size() + resultDendrites.size());
                    std::copy(resultDendrites.begin(), resultDendrites.end(), tmp.begin() + reverseCurve.size());
                    treeDataSet.back().push_back(VectorVec3d());
                    treeDataSet.back().back().swap(tmp);*/
                    continue;
                } else if(collideID == 1){
                    //collide soma
                    Vec5d tmpp;
                    tmpp.block(0,0,3,1) = treeDataSet[0][0][0];
                    tmpp(3)=1;tmpp(4)=1;
                    reverseCurve.push_back(tmpp);//add soma position
                    VectorVec3d tmp;
                    NGUtility::GetReversePartVectorVec3d(reverseCurve, 0, int(reverseCurve.size()) - 1, tmp);
                    std::copy(resultDendrites.begin(), resultDendrites.end(), back_inserter(tmp));
                    treeDataSet.back().push_back(VectorVec3d());
                    treeDataSet.back().back().swap(tmp);
                }
                else if (collideID > 1000000){
                    //collide to the previous part of curve.
                    continue;
                }
                else{
                    int cellLevel, cellInd;
                    Ind2CellInd(treeDataSet, collideID-1, cellLevel, cellInd);//soma is excluded.
                    double minDist = 10000;
                    int minIdx = -1;
                    double tmpDist;
                    Vec3d recurback;
                    if (reverseCurve.empty()) {
                        recurback = traceDendrite.front().block(0,0,3,1);
                    }else
                        recurback = reverseCurve.back().block(0, 0, 3, 1);
                    for(size_t k = 0; k < treeDataSet[cellLevel][cellInd].size(); ++k){
                        tmpDist = 0.0;
                        tmpDist += paramPack->xRes_ * std::abs((treeDataSet[cellLevel][cellInd][k](0) - recurback(0)));
                        tmpDist += paramPack->yRes_ * std::abs((treeDataSet[cellLevel][cellInd][k](1) - recurback(1)));
                        tmpDist += paramPack->zRes_ * std::abs((treeDataSet[cellLevel][cellInd][k](2) - recurback(2)));
                        if(tmpDist < minDist){
                            minIdx = int(k);
                            minDist = tmpDist;
                        }
                    }
                    if(minDist < 5.0){//wrong cubelabel when weak signal
                        //if (minIdx < 2) minIdx = 0;
                        /*if (minIdx > treeDataSet[cellLevel][cellInd].size() - 3){ //minIdx = treeDataSet[cellLevel][cellInd].size() - 1;
                            treeDataSet[cellLevel][cellInd].erase(treeDataSet[cellLevel][cellInd].begin() + minIdx + 1, treeDataSet[cellLevel][cellInd].end());
                        }*/
                        Vec5d tmpp;
                        tmpp.block(0,0,3,1) = treeDataSet[cellLevel][cellInd][minIdx];
                        tmpp(3)=1;tmpp(4)=1;
                        reverseCurve.push_back(tmpp);//add soma position
                        VectorVec3d tmp;
                        NGUtility::GetReversePartVectorVec3d(reverseCurve, 0, int(reverseCurve.size()) - 1, tmp);
                        tmp.resize(reverseCurve.size() + resultDendrites.size());
                        std::copy(resultDendrites.begin(), resultDendrites.end(), tmp.begin() + reverseCurve.size());
						//std::copy(tmp.begin(), tmp.end(), std::back_inserter(treeDataSet[cellLevel][cellInd]));
                        treeDataSet.back().push_back(VectorVec3d());
                        treeDataSet.back().back().swap(tmp);
                    }else
                        continue;
                }

                ++globalID_;
                CubeLabelL1Opt(treeDataSet.back().back());
                //2015-6-8
                //if(!treeDataSet.back().empty())
                    //CubeLabel(treeDataSet.back().back());
            } else{//2015-5-7
            }
        }
    }
}

void SparseTraceFilter::CalVectorAngle(const Vec3d &curNormPt, double &altitudeAngle, double &azimuthAngle)
{
    altitudeAngle = std::acos(curNormPt(2));
    if(std::abs(curNormPt(0)) < 0.01){
        if(curNormPt(1) > 0.0)
            azimuthAngle = M_PI_2;
        else
            azimuthAngle = M_PI_2 + M_PI;
    }
    else{
        azimuthAngle=std::atan2(curNormPt(1), curNormPt(0));
        if(azimuthAngle <= 0.0)
            azimuthAngle += 2 * M_PI;
    }
}

void SparseTraceFilter::ClustersExtractionPoints(const VectorVec3i &inflatPtSet, CVolume &growLabelMatrix,
                                                 VectorVec3i &subClustrePtSet, std::vector<int> &subClustrePtSetNum)
{
    growLabelMatrix.SetZero();
    size_t inflatPtSetNum = inflatPtSet.size();
    for(size_t i = 0; i < inflatPtSetNum; ++i){
        growLabelMatrix(inflatPtSet[i](0), inflatPtSet[i](1), inflatPtSet[i](2)) = 1;
    }
    subClustrePtSet.clear();
    subClustrePtSetNum.clear();
    ClustersExtraction(growLabelMatrix, inflatPtSet, subClustrePtSet, subClustrePtSetNum);
}

void SparseTraceFilter::RegionConnectionsSub(const VectorVec3i &curConPtSet, Volume<int> &growLabelMatrix,
                                             int threvValue, int &connectionFlag)
{
    size_t curConPtSetNum = curConPtSet.size();
    std::vector<int> indexWindowVal;
    indexWindowVal.reserve(2000);
    int nx = growLabelMatrix.x() - 1;
    int ny = growLabelMatrix.y() - 1;
    int nz = growLabelMatrix.z() - 1;//warning!!!
    int xMin, xMax, yMin, yMax, zMin, zMax;
    int pixelValue;
    int headFlag;
    //int kk = -1;
    for(size_t i = 0; i < curConPtSetNum; ++i){
        const Vec3i& curPt = curConPtSet[i];
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPt(0), curPt(1), curPt(2), 2,2,2,
                    0, nx, 0, ny, 0, nz);
        for(int i1 = xMin; i1 <= xMax; ++i1){
            for(int j1 = yMin; j1 <= yMax; ++j1){
                for(int k1 = zMin; k1 <= zMax; ++k1){
                   pixelValue = growLabelMatrix(i1,j1,k1);
                   headFlag = pixelValue / 10000;//2015-6-8
                   if( headFlag == threvValue){
                       //++kk;
                       indexWindowVal.push_back(pixelValue - 10000 * headFlag);
                       if(indexWindowVal.size() == 2000)
                           indexWindowVal.pop_back();
                   }
                }
            }
        }
    }
    //find the largest collide area
    int maxVal, minVal;
    if(!indexWindowVal.empty()){
        std::vector<int> indexWindowValSort;
        for(size_t i = 0; i < indexWindowVal.size(); ++i)
            indexWindowValSort.push_back(indexWindowVal[i]);
        std::sort(indexWindowValSort.begin(), indexWindowValSort.end());
        maxVal = indexWindowValSort.back();
        minVal = indexWindowValSort.front();
        std::vector<int> indexHist;
        indexHist.resize(maxVal - minVal + 1, 0);
        for(size_t j = 0; j < indexWindowValSort.size(); ++j){
            ++indexHist[indexWindowValSort[j] - minVal];
        }
        connectionFlag = int(std::max_element(indexHist.begin(), indexHist.end()) - indexHist.begin()) + minVal;
    }else{
        connectionFlag = 0;
    }
}

void SparseTraceFilter::RayburstShapeTrackV2(const Vec3d &end2ndPt, const Vec3d &endPt,
                                             const Volume<unsigned short> &origImg, int curveLen,
                                             VectorVec3d &forwardArea)
{
    forwardArea.clear();
    //forwardArea.resize(curveLen);
    Vec3d initDir = endPt - end2ndPt;
    Vec3d tmpVec;
    VectorVec3i rayNode;
    double intersectAngle;
    Vec3d initPtCp;
    Vec3d recursiveEndPt = endPt;
    for(int j = 0; j < curveLen; ++j){
        rayNode.clear();
        RayBurstShapeV2(recursiveEndPt, origImg, rayNode);
        size_t rayNodeNum = rayNode.size();
        if(rayNodeNum > 2){
            std::vector<int> validFlag(rayNodeNum, 0);
            for(size_t i = 0; i < rayNodeNum; ++i){
                tmpVec << rayNode[i](0) - recursiveEndPt(0), rayNode[i](1) - recursiveEndPt(1),
                          rayNode[i](2) - recursiveEndPt(2);
                intersectAngle = initDir.dot(tmpVec);
                if(intersectAngle > 0.0)
                    validFlag[i] = 1;
            }
            std::vector<size_t> validNodeList;
            for(size_t i = 0; i < validFlag.size(); ++i){
                if(validFlag[i] == 1)
                    validNodeList.push_back(i);
            }
            if(validNodeList.size() > 1){
                VectorVec3d tmpRayNode;
                Vec3d tmpVec3d;
                for(size_t k = 0; k < validNodeList.size(); ++k){
                    tmpVec3d << rayNode[validNodeList[k]](0), rayNode[validNodeList[k]](1),
                            rayNode[validNodeList[k]](2);
                    tmpRayNode.push_back(tmpVec3d);
                }
                std::vector<double> rayNodeWet;
                WeighRayValue(tmpRayNode, origImg, rayNodeWet);
                double sumRayNodeWet = 0.0;
                for(size_t k = 0; k < rayNodeWet.size(); ++k){
                    sumRayNodeWet += rayNodeWet[k];
                }
                for(size_t k = 0; k < rayNodeWet.size(); ++k){
                    rayNodeWet[k] /= sumRayNodeWet;
                }
                initPtCp.setZero();
                for(size_t k = 0; k < rayNodeWet.size(); ++k){
                    tmpVec3d << rayNode[validNodeList[k]](0), rayNode[validNodeList[k]](1),
                            rayNode[validNodeList[k]](2);
                    initPtCp += rayNodeWet[k] * tmpVec3d;
                }
            } else
                break;
            forwardArea.push_back(initPtCp);
            recursiveEndPt = forwardArea[j];
        }else
            break;
    }
}

void SparseTraceFilter::RayBurstShapeV2(const Vec3d &initPt, const Volume<unsigned short> &origImg, VectorVec3i &rayNode)
{
    rayNode.clear();

    double slice = 0.3;//(double)(minLen) / 82.0;
    const int blocksum = 26;

    std::vector<double> lineSegment;
    for(int i = 0; i < blocksum; ++i){
        lineSegment.push_back(double(i) * slice);
    }

    const int Theta = 20;
    const int Phi = 10;

    const SVolume& locOrigImg = origImg;//2015-6-8
    const Vec3d& locSoma = initPt;

    DVolume sphereRayWet;//(lineSegment.size(), Theta, Phi);
    sphereRayWet.SetSize(int(lineSegment.size()), Theta, Phi);
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;

    double segmentDense(0.0), segmentWet(0.0);
    double x,y,z;

    size_t numLineSegment;
    numLineSegment = lineSegment.size();
    for (size_t k = 0; k < numLineSegment; ++k){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                x = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::cos(a * (double)i * PI_180) + locSoma(0);
                y = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::sin(a * (double)i * PI_180) + locSoma(1);
                z = lineSegment[k] * std::cos(b * (double)j * PI_180) + locSoma(2);
                segmentDense = segmentWet = 0.0;
                ContourUtil::CalculateSphereOneNode(locOrigImg, 1.0, x, y, z, segmentDense, segmentWet);
                sphereRayWet(k, i-1, j-1) = segmentDense / (segmentWet + 0.0001);
            }
        }
    }//for
    //
    int lenR0_1 = int(lineSegment.size()) - 1;
    std::vector<double> outerShell;
    for (int j = 0; j < Phi; ++j){
        for (int i = 0; i < Theta; ++i){
            outerShell.push_back(sphereRayWet(lenR0_1, i, j));
        }
    }

    std::vector<double> boundaryBack;
    ContourUtil::GetBoundaryBack(outerShell, 4, boundaryBack);

    //three_vs = mean(Lssx1)+3.5*std(Lssx1)
    double constrictionThrev(0.0);
    double constrictionThrevMean =
        std::accumulate(boundaryBack.begin(), boundaryBack.end(), 0.0)
            / double(boundaryBack.size());
    //standard variation---------------------------------------
    double constrictionThrevMeanSqrtSum(0);
    int numBoundaryBack = int(boundaryBack.size());
    for (int i = 0; i< numBoundaryBack; ++i)
    {
        constrictionThrevMeanSqrtSum +=
            (boundaryBack[i] - constrictionThrevMean) * (boundaryBack[i] - constrictionThrevMean);
    }
    double stdThrev = std::sqrt(constrictionThrevMeanSqrtSum / double(boundaryBack.size() - 1));
    constrictionThrev = 3.5 * std::sqrt(constrictionThrevMeanSqrtSum / (boundaryBack.size() - 1));
    constrictionThrev += constrictionThrevMean;
    //boundary value
    for (int i = 0; i < Phi; ++i){
        for (int j = 0; j < Theta; ++j)
            sphereRayWet(lenR0_1, j, i) = boundaryBack[i * Theta + j];
    }

    std::vector<std::vector<double> > rayLimit;
    ContourUtil::GetRayLimit(sphereRayWet, constrictionThrev, rayLimit);
    //2015-2-7
    std::vector<std::vector<double> > rayLimit1;
    ContourUtil::GetRayLimit(sphereRayWet, constrictionThrev + 4.5*stdThrev, rayLimit1);
    Volume<double> smoothRay;//2015-2-7
    TraceUtil::GetGradientVectorFlowForTrace(sphereRayWet, smoothRay);

    std::vector<double> lineSegLength;//
    //generate_n(back_inserter(lineSegLength), blocksum, GenArray<double>(1.0, 1.0));
    for(int i = 0; i < blocksum; ++i){
        lineSegLength.push_back(double(i + 1));
    }

    double curRayLength(0);
    std::vector<double> curSmoothRay( smoothRay.x() );
    std::vector<double> distWet( smoothRay.x() );

    double reduceSmooth[2];
    reduceSmooth[0] = 0.95;//2015-2-7
    reduceSmooth[1] = 1.0 - reduceSmooth[0];

    MatXd resultRayLength = 2.0*MatXd::Ones(Theta+2, Phi+2);// 2015-2-7
    for(size_t i = 0; i < rayLimit.size();++i){
        for(size_t j = 0; j < rayLimit[0].size(); ++j){
            resultRayLength(i+1,j+1) = std::max(rayLimit[i][j] - 2.0, 2.0);
        }
    }

    int repeat = 50;//2015-2-7
    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                curRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    curSmoothRay[ij] = smoothRay(ij, i - 1, j - 1);
                    distWet[ij] = curSmoothRay[ij] * std::exp( -0.05 *std::pow( lineSegLength[ij] - curRayLength, 2.0));
                }

                resultRayLength(i,j) += reduceSmooth[0] * ( ( (std::inner_product(lineSegLength.begin(), lineSegLength.end(),
                    distWet.begin(), 0.0)) /
                    (std::accumulate(distWet.begin(), distWet.end(), 0.0) + 0.001 ) ) - resultRayLength(i,j) )
                    - reduceSmooth[1] * (4.0 * resultRayLength(i,j)
                    - resultRayLength(i-1,j)
                    - resultRayLength(i+1,j)
                    - resultRayLength(i,j-1)
                    - resultRayLength(i,j+1));

                resultRayLength(i,j) = std::min(resultRayLength(i,j), rayLimit[i - 1][ j - 1]);
                resultRayLength(i,j) = std::max(resultRayLength(i,j), rayLimit1[i - 1][ j - 1]);
            }
            resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for

    rayNode.reserve(Theta * Phi * lineSegment.size());
    double xx(0.0), yy(0.0), zz(0.0);
    Vec3i tmp;
    for(int j = 1; j < Phi + 1; ++j){
        for(int i = 1; i < Theta + 1; ++i){
            if(resultRayLength(i, j) > 1.0){
                for(double k = 1.0; k <= 0.3 * resultRayLength(i,j); k += 0.5){
                    xx = (k * std::sin(b * (double)(j) * PI_180) * std::cos(a * (double)(i) * PI_180) + initPt(0));
                    yy = (k * std::sin(b * (double)(j) * PI_180) * std::sin(a * (double)(i) * PI_180) + initPt(1));
                    zz = (k * std::cos(b * (double)(j) * PI_180) + initPt(2));
                    tmp(0) = NGUtility::Round(xx);
                    tmp(1) = NGUtility::Round(yy);
                    tmp(2) = NGUtility::Round(zz);
                    rayNode.push_back(tmp);
                }
            }
        }
    }
#ifdef _WIN32
    std::sort(rayNode.begin(), rayNode.end(), Vec3i_less());
#else
    std::sort(rayNode.begin(), rayNode.end(), [](const Vec3i& lhs, const Vec3i &rhs){
        if (lhs(0) != rhs(0))  return lhs(0) < rhs(0);
        else if (lhs(1) != rhs(1))  return lhs(1) < rhs(1);
        else if (lhs(2) != rhs(2))  return lhs(2) < rhs(2);
        return false;
    });
#endif
    
    rayNode.erase(std::unique(rayNode.begin(), rayNode.end()), rayNode.end());

    int nxxp = int(rayNode.size());
    VectorVec3i rayNodeCp;
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();
    for(int i = 0; i < nxxp; ++i){
        const Vec3i& curRayNode = rayNode[i];
        if(curRayNode.minCoeff() > -1.0 && curRayNode(0) < nxx
                && curRayNode(1) < nyy  && curRayNode(2) < nzz ){
            rayNodeCp.push_back(curRayNode);
        }
    }
    rayNode.swap(rayNodeCp);
}

//2015-6-8 threValue,origImg, backImg, traceLabelMatrix is global variant
void SparseTraceFilter::TraceCurveForwardFromSomaV2(const Vec3d &seedPoint, const Vec3d &curDirection,
                                                    VectorVec5d &resultSeedCurve)
{
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    resultSeedCurve.clear();
    Vec5d tmp1;
    tmp1 << seedPoint(0), seedPoint(1), seedPoint(2), 0, 0;
    resultSeedCurve.push_back(tmp1);

    //Vec2i collideFlag;
    //collideFlag.setZero();
    int isNextNodeValid(1);
    int i = 0;

    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();

    Vec3d nextDendDir(curDirection);
    //Vec3d fixInitDendDir(curDirection);

    while(isNextNodeValid == 1 && i < 3000){
        int xMin, xMax, yMin, yMax, zMin, zMax;
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax,
            NGUtility::Round(resultSeedCurve[i](0) ), NGUtility::Round(resultSeedCurve[i](1)), NGUtility::Round(resultSeedCurve[i](2)), 10, 10, 7,
            0, nxx - 1, 0, nyy - 1, 0, nzz - 1);

        Volume<unsigned short> locOrigImg;
        ExtractArea(origImg, xMin, xMax, yMin, yMax, zMin, zMax, locOrigImg);

        if((xMax - xMin) * (yMax - yMin) * (zMax - zMin) <= 0){
            isNextNodeValid = 0;
            break;
        }
        Vec3d lefttopAxis(xMin, yMin, zMin);
        //2014-4-22
        int rdrvx = std::min(std::max(NGUtility::Round(resultSeedCurve[i](0) - 1.0 ), 0), nxx - 1);
        int rdrvy = std::min(std::max(NGUtility::Round(resultSeedCurve[i](1) - 1.0 ), 0), nyy - 1);
        int rdrvz = std::min(std::max(NGUtility::Round(resultSeedCurve[i](2) - 1.0 ), 0), nzz - 1);//2015-6-8
        double threv = backImg(rdrvx, rdrvy, rdrvz);
        //double threv = backImg(rdrvx, rdrvy, NGUtility::Round(double(zMax - zMin)/2.0 + zMin));
        Vec3d tmpCurrentCurveNode(resultSeedCurve[i](0) - xMin, resultSeedCurve[i](1) - yMin, resultSeedCurve[i](2) - zMin);

        Vec3d nextCurveNode;
        nextCurveNode.setZero();
        Vec3d nextDenDirCopy;
        nextDenDirCopy.setZero();
        Vec3d tmpx1;
        //2015-6-8 threValue is global varient, ttk_Label is used here.
        if(i==0)
            TraceNextCurveNodeV2(locOrigImg, tmpCurrentCurveNode, threv, nextDendDir,0,nextCurveNode, tmpx1,
                             isNextNodeValid);
        else
            TraceNextCurveNodeV2(locOrigImg, tmpCurrentCurveNode, threv, nextDendDir,1,nextCurveNode, tmpx1,
                                 isNextNodeValid);
        nextDendDir = tmpx1;
        nextDenDirCopy = nextDendDir;

        if ( nextCurveNode(0) != 0.0 && nextCurveNode(1) != 0.0 && nextCurveNode(2) != 0.0 )
            nextCurveNode += lefttopAxis;

        if (isNextNodeValid == 1){
            // in fact it is not nextCurveNode but corrected current node according to forward direction.
            resultSeedCurve[i](0) = nextCurveNode(0);
            resultSeedCurve[i](1) = nextCurveNode(1);
            resultSeedCurve[i](2) = nextCurveNode(2);

            double radius(0.0);
            double rav(0.0);

            CalcParmOfCurveNode(origImg, backImg, nextCurveNode, radius, rav);
            resultSeedCurve[i](3) = radius;
            //if (radius > 1.0){//2015-2-7
            ++i;//warning!!
            Vec3d rayNode2 = nextCurveNode + 2.0 * nextDenDirCopy;
            Vec3d rayNode1 = nextCurveNode + nextDenDirCopy;
            VectorVec3d tmpRayNode;
            tmpRayNode.push_back(rayNode1);
            tmpRayNode.push_back(rayNode2);
            std::vector<double> rayValue;
            WeighRayValue(tmpRayNode, origImg, rayValue);

            if ( rayValue[0] > rayValue[1] + 3.5 * std::sqrt(threv) ){
                Vec3d tmpVec = nextCurveNode + 2.0 * nextDenDirCopy;
                Vec5d tmp; tmp.setZero();
                tmp(0) = tmpVec(0);
                tmp(1) = tmpVec(1);
                tmp(2) = tmpVec(2);
                tmp(4) = rayValue[0];
                resultSeedCurve.push_back(tmp);
            }
            else{
                Vec3d tmpVec = nextCurveNode + 2.0 * nextDenDirCopy;
                Vec5d tmp; tmp.setZero();
                tmp(0) = tmpVec(0);
                tmp(1) = tmpVec(1);
                tmp(2) = tmpVec(2);
                tmp(4) = rayValue[1];
                resultSeedCurve.push_back(tmp);
            }//if
            //else isNextNodeValid = 0;
        }

        Vec3d tmpCurCurveNode(resultSeedCurve[i](0) * paramPack->xScale_, 
            resultSeedCurve[i](1) * paramPack->yScale_, resultSeedCurve[i](2) * paramPack->zScale_);
        if (tmpCurCurveNode.minCoeff() < 2.0 || tmpCurCurveNode(0) > nxx * paramPack->xScale_ - 2
            || tmpCurCurveNode(1) > nyy * paramPack->yScale_ - 2
            || tmpCurCurveNode(2) > nzz * paramPack->zScale_ - 2) {//2017-3-27
            isNextNodeValid = 0;
        }

        if(i < 4)//2015-6-8
            nextDendDir.normalize();

        if (i > 3){//matlab is 4
            //2015-2-7
            Vec3d nearNodeDist;
            nearNodeDist.setZero();
            for (int iij = 1; iij < 4; ++iij){
                Vec3d tmp11(resultSeedCurve[i - iij](0), resultSeedCurve[i - iij](1), resultSeedCurve[i - iij](2));
                Vec3d tmp22(resultSeedCurve[i - iij + 1](0), resultSeedCurve[i - iij + 1](1), resultSeedCurve[i - iij + 1](2));
                nearNodeDist(iij - 1) = (tmp22 - tmp11).norm();
            }

            VectorVec3d dataPrincipald;
            Vec3d tmpDataPrin;
            for ( int n = std::max(i - 10, 0) ; n <= i - 2; ++n ){
                tmpDataPrin(0) = resultSeedCurve[n](0);
                tmpDataPrin(1) = resultSeedCurve[n](1);
                tmpDataPrin(2) = resultSeedCurve[n](2);
                dataPrincipald.push_back(tmpDataPrin);
            }

            Vec3d sx1;//, sx2, sx3;
            CalcPtCurveDirection(dataPrincipald, sx1);
            nextDendDir = sx1;
            //x2 = sx2;
            //x3 = sx3;

            Vec3d forwardAngle(sx1);
            Vec3d dataLLi, dataLLi_1;
            dataLLi(0) = resultSeedCurve[i](0);
            dataLLi(1) = resultSeedCurve[i](1);
            dataLLi(2) = resultSeedCurve[i](2);

            dataLLi_1(0) = resultSeedCurve[i - 1](0);
            dataLLi_1(1) = resultSeedCurve[i - 1](1);
            dataLLi_1(2) = resultSeedCurve[i - 1](2);

            Vec3d backwardAngle = dataLLi - dataLLi_1;

            double crossAngle = forwardAngle.transpose() * backwardAngle;
            crossAngle /=  forwardAngle.norm() * backwardAngle.norm();

            if (nearNodeDist.minCoeff() < 0.05 || crossAngle < 0.1 )//2015-2-7
                isNextNodeValid = 0;
            if(isNextNodeValid == 1)
                IsCurvesCollideV2(resultSeedCurve[i-1], isNextNodeValid);//2015-6-8
        }
    }
    //2015-6-8
    if(i<4)
        resultSeedCurve.clear();
}

//2015-6-8 threValue,origImg, backImg, traceLabelMatrix is global variant
void SparseTraceFilter::TraceCurveForwardFromSomaV3(const Vec3d &seedPoint, const Vec3d &curDirection,
                                                    VectorVec5d &resultSeedCurve, int& collideID, int maxStep)
{
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    collideID = 0;//2015-6-8
    resultSeedCurve.clear();
    Vec5d tmp1;
    tmp1 << seedPoint(0), seedPoint(1), seedPoint(2), 0, 0;
    resultSeedCurve.push_back(tmp1);

    Vec2i collideFlag;
    collideFlag.setZero();
    int isNextNodeValid(1);
    int i = 0;

    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();

    Vec3d nextDendDir(curDirection);
    //Vec3d fixInitDendDir(curDirection);

    while (isNextNodeValid == 1 && i < maxStep){
        int xMin, xMax, yMin, yMax, zMin, zMax;
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax,
            NGUtility::Round(resultSeedCurve[i](0) ), NGUtility::Round(resultSeedCurve[i](1)), NGUtility::Round(resultSeedCurve[i](2)), 10, 10, 7,
            0, nxx - 1, 0, nyy - 1, 0, nzz - 1);

        Volume<unsigned short> locOrigImg;
        ExtractArea(origImg, xMin, xMax, yMin, yMax, zMin, zMax, locOrigImg);

        if((xMax - xMin) * (yMax - yMin) * (zMax - zMin) == 0){
            isNextNodeValid = 0;
            break;
        }
        Vec3d lefttopAxis(xMin, yMin, zMin);
        //2014-4-22
        int rdrvx = std::min(std::max(NGUtility::Round(resultSeedCurve[i](0) - 1.0 ), 0), nxx - 1);
        int rdrvy = std::min(std::max(NGUtility::Round(resultSeedCurve[i](1) - 1.0 ), 0), nyy - 1);
        int rdrvz = std::min(std::max(NGUtility::Round(resultSeedCurve[i](2) - 1.0 ), 0), nzz - 1);//2015-6-8
        double threv = backImg(rdrvx, rdrvy, rdrvz);
        Vec3d tmpCurrentCurveNode(resultSeedCurve[i](0) - xMin, resultSeedCurve[i](1) - yMin, resultSeedCurve[i](2) - zMin);

        Vec3d nextCurveNode;
        nextCurveNode.setZero();
        Vec3d nextDenDirCopy;
        nextDenDirCopy.setZero();
        Vec3d tmpx1;
        //2015-6-8 threValue is global variant, ttk_Label is used here.
        if(i==0)
            TraceNextCurveNodeV2(locOrigImg, tmpCurrentCurveNode, threv, nextDendDir,0,nextCurveNode, tmpx1,
                             isNextNodeValid);
        else
            TraceNextCurveNodeV2(locOrigImg, tmpCurrentCurveNode, threv, nextDendDir,1,nextCurveNode, tmpx1,
                                 isNextNodeValid);
        nextDendDir = tmpx1;
        nextDenDirCopy = nextDendDir;

        if ( nextCurveNode(0) != 0.0 && nextCurveNode(1) != 0.0 && nextCurveNode(2) != 0.0 )
            nextCurveNode += lefttopAxis;

        if (isNextNodeValid == 1){
            // in fact it is not nextCurveNode but corrected current node according to forward direction.
            resultSeedCurve[i](0) = nextCurveNode(0);
            resultSeedCurve[i](1) = nextCurveNode(1);
            resultSeedCurve[i](2) = nextCurveNode(2);

            double radius(0.0);
            double rav(0.0);

            CalcParmOfCurveNode(origImg, backImg, nextCurveNode, radius, rav);
            resultSeedCurve[i](3) = radius;
            //if (radius > 1.0){//2015-2-7
            ++i;//warning!!
            Vec3d rayNode2 = nextCurveNode + 2.0 * nextDenDirCopy;
            Vec3d rayNode1 = nextCurveNode + nextDenDirCopy;
            VectorVec3d tmpRayNode;
            tmpRayNode.push_back(rayNode1);
            tmpRayNode.push_back(rayNode2);
            std::vector<double> rayValue;
            WeighRayValue(tmpRayNode, origImg, rayValue);

            if ( rayValue[0] > rayValue[1] + 3.5 * std::sqrt(threv) ){
                Vec3d tmpVec = nextCurveNode + 2.0 * nextDenDirCopy;//2015-6-8
                Vec5d tmp; tmp.setZero();
                tmp(0) = tmpVec(0);
                tmp(1) = tmpVec(1);
                tmp(2) = tmpVec(2);
                tmp(4) = rayValue[0];
                resultSeedCurve.push_back(tmp);
            }
            else{
                Vec3d tmpVec = nextCurveNode + 2.0 * nextDenDirCopy;
                Vec5d tmp; tmp.setZero();
                tmp(0) = tmpVec(0);
                tmp(1) = tmpVec(1);
                tmp(2) = tmpVec(2);
                tmp(4) = rayValue[1];
                resultSeedCurve.push_back(tmp);
            }//if
            //else isNextNodeValid = 0;
        }

        Vec3d tmpCurCurveNode(resultSeedCurve[i](0), resultSeedCurve[i](1), resultSeedCurve[i](2));
        if ( tmpCurCurveNode.minCoeff() < 0.0 || tmpCurCurveNode(0) > nxx-1
            || tmpCurCurveNode(1) > nyy-1
            || tmpCurCurveNode(2) > nzz-1  ) {//2015-5-7
            isNextNodeValid = 0;
        }

        if(i < 4)//2015-6-8
            nextDendDir.normalize();

        if (i > 3){//matlab is 4
            //2015-2-7
            Vec3d nearNodeDist;
            nearNodeDist.setZero();
            for (int iij = 1; iij < 4; ++iij){
                Vec3d tmp11(resultSeedCurve[i - iij](0), resultSeedCurve[i - iij](1), resultSeedCurve[i - iij](2));
                Vec3d tmp22(resultSeedCurve[i - iij + 1](0), resultSeedCurve[i - iij + 1](1), resultSeedCurve[i - iij + 1](2));
                nearNodeDist(iij - 1) = (tmp22 - tmp11).norm();
            }

            VectorVec3d dataPrincipald;
            Vec3d tmpDataPrin;
            for ( int n = std::max(i - 10, 0) ; n <= i - 2; ++n ){
                tmpDataPrin(0) = resultSeedCurve[n](0);
                tmpDataPrin(1) = resultSeedCurve[n](1);
                tmpDataPrin(2) = resultSeedCurve[n](2);
                dataPrincipald.push_back(tmpDataPrin);
            }

            Vec3d sx1;//, sx2, sx3;
            CalcPtCurveDirection(dataPrincipald, sx1);
            nextDendDir = sx1;
            //x2 = sx2;
            //x3 = sx3;

            Vec3d forwardAngle(sx1);
            Vec3d dataLLi, dataLLi_1;
            dataLLi(0) = resultSeedCurve[i](0);
            dataLLi(1) = resultSeedCurve[i](1);
            dataLLi(2) = resultSeedCurve[i](2);

            dataLLi_1(0) = resultSeedCurve[i - 1](0);
            dataLLi_1(1) = resultSeedCurve[i - 1](1);
            dataLLi_1(2) = resultSeedCurve[i - 1](2);

            Vec3d backwardAngle = dataLLi - dataLLi_1;

            double crossAngle = forwardAngle.transpose() * backwardAngle;
            crossAngle /=  forwardAngle.norm() * backwardAngle.norm();

            if (nearNodeDist.minCoeff() < 0.05 || crossAngle < 0.1 )//2015-2-7
                isNextNodeValid = 0;
            if(isNextNodeValid == 1)
                IsCurvesCollideV3(resultSeedCurve[i-1], isNextNodeValid, collideID);//2015-6-8
        }
    }
    //2015-6-8
    if(i<4)
        resultSeedCurve.clear();
}

//2015-6-8 add indexImg
void SparseTraceFilter::CubeLabel(const VectorVec3d &curDendrite)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();

    VectorVec5d rawSomaCurve;
    rawSomaCurve.clear();
    Vec5d tmpVec5d;
    for(size_t i = 0; i < curDendrite.size(); ++i){
        tmpVec5d << curDendrite[i](0), curDendrite[i](1), curDendrite[i](2), 0, 0;
        rawSomaCurve.push_back(tmpVec5d);
    }
    VectorVec3d tmpPtCurve;
    Vec3d tmpVec3d;
    for (size_t i = 0; i < rawSomaCurve.size(); ++i){
        tmpVec3d << rawSomaCurve[i](0), rawSomaCurve[i](1), rawSomaCurve[i](2);
        tmpPtCurve.push_back(tmpVec3d);
    }
    //
    if(!curDendrite.empty()){
        std::vector<double> R;//here is list
        std::vector<double> V;
        VectorVec5d resultCurveCopy;
        CalcParmOfCurveNodeList(origImg, backImg, tmpPtCurve, R, V );
        for (size_t i = 0; i < rawSomaCurve.size(); ++i){
            rawSomaCurve[i](3) = R[i];
            rawSomaCurve[i](4) = V[i];
        }
        Vec5d half;
        //
        for (size_t i = 0; i < rawSomaCurve.size() - 1; ++i){
            resultCurveCopy.push_back(rawSomaCurve[i]);
            half = 0.5 * (rawSomaCurve[i] + rawSomaCurve[i + 1]);
            resultCurveCopy.push_back(half);
        }
        resultCurveCopy.push_back(rawSomaCurve[rawSomaCurve.size() - 1]);

        std::vector<double> three;
        for (VectorVec5d::size_type i = 0; i < resultCurveCopy.size(); ++i){
            three.push_back((double)std::min(4, NGUtility::Round(resultCurveCopy[i](3) + 3.5) ) );//2017-11
        }
        int id1(0), id2(0), id3(0);
        int xMin, xMax, yMin, yMax, zMin, zMax;
        double threBack;
        for (int ik = 1; ik < 2*int(rawSomaCurve.size())-2; ++ik){
            id1 = NGUtility::Round(resultCurveCopy[ik](0) );
            id2 = NGUtility::Round(resultCurveCopy[ik](1) );
            id3 = NGUtility::Round(resultCurveCopy[ik](2) );

            xMin = std::max(0, id1 - std::max(0, (int)three[ik]));
            xMax = std::min(nxx - 1, id1 + std::max(0, (int)three[ik]));
            yMin = std::max(0, id2 - std::max(0, (int)three[ik]));
            yMax = std::min(nyy - 1, id2 + std::max(0, (int)three[ik]));
            zMin = std::max(0, id3 - std::max(0, (int)three[ik]));
            zMax = std::min(nzz - 1, id3 + std::max(0, (int)three[ik]));

            //double threBack(0);
            threBack = std::max(double(origImg(id1, id2, id3) + backImg(id1, id2, id3)) / 2, 0.85*double(origImg(id1, id2, id3)));
            threBack = std::max(threBack, backImg(id1, id2, id3) + paramPack->axonDiffuseValue_);
            for (int ii = xMin; ii <= xMax; ++ii){
                for (int jj = yMin; jj <= yMax; ++jj){
                    for (int kk = zMin; kk <= zMax; ++kk){
                        //2015-6-8
                        threBack = double(backImg(ii, jj, kk)) + 4.0 * std::sqrt(double(backImg(ii, jj, kk)));//2016-3-18
                        threBack = std::min(threBack, 1.5 * double(backImg(ii, jj, kk)) );//2016-3-18
                        threBack = std::max(300.0, threBack);//2016-3-18
                        if (traceLabelMatrix_(ii, jj, kk) == 0 && double(origImg(ii, jj, kk)) > threBack){
                            traceLabelMatrix_(ii, jj, kk) = 1;
                            indexImg_(ii,jj,kk) = globalID_;//2015-6-8
                        }
                    }
                }
            }//for
        }//for
        //printf("%d\n",int(test.size()));
    }
}

//2015-6-8 here the traceLabelMatrix is the copy of global traceLabelMatrix
void SparseTraceFilter::DetectBifurcation(CVolume &traceLabelMatrix, const std::vector<VectorVec3d> &curveCluster,
                                          VectorVec3d &bifurcationDirections, std::vector<VectorVec3d> &bifurcationClustre,
                                          VectorVec3i& boundaryLabel)
{
    //traceLabelMatrix is traceLabelMatrixCp;
    std::vector<VectorVec3i> growPointSet;
    ExtractCubeBoundary(curveCluster, traceLabelMatrix, growPointSet, boundaryLabel);
    bool isContinue = true;
    for(size_t i = 0; i < growPointSet.size(); ++i){
        if(growPointSet[i].empty())
            isContinue = false;
    }
    if(isContinue){
        std::vector<VectorVec3i> clustrePtSet;
        std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> > clustreLabel;
        ClusterGrowPtSet(growPointSet, traceLabelMatrix, clustrePtSet, clustreLabel);
        size_t flagNum = 0;
        for(size_t i =0; i < clustreLabel.size(); ++i){
            flagNum = std::max(flagNum, clustreLabel[i].size());
        }
        if(flagNum > 10000){//2015-6-8
            bifurcationDirections.clear();
            bifurcationDirections.clear();
        } else{
            std::vector<std::vector<int> > subRegionCon;
            LOG(INFO)<<flagNum;
            RegionConnections(clustrePtSet, clustreLabel, subRegionCon);//2015-6-8

            MatXi pathMatrix;
            ConnectPathInTree(subRegionCon, pathMatrix);//2015-6-8

            /*total 9 level. just compare 4-6th level to remove the overlap curve.
             * why not compare 7-9th? 1-6th is long enough. overlap 4-6th level is
             * enough long.
             */
            std::vector<int> redundancyLabel;
            DeleteRedundancyPathV2(pathMatrix, redundancyLabel);//2015-6-8

            std::vector<VectorVec4d> connectPoints;
            std::vector<VectorVec4d> shortConnectPoints;
            ConnectGrowClustreFromSoma(subRegionCon, clustreLabel, connectPoints, shortConnectPoints);
            //size_t connectPointsNum = connectPoints.size();
            //size_t shortConnectPointsNum = shortConnectPoints.size();

            //2015-6-8
            //int remainPathNum = std::accumulate(redundancyLabel.begin(), redundancyLabel.end(),0);
            //bifurcationDirections.resize(remainPathNum);
            //bifurcationClustre.resize(remainPathNum);
            for(int i = 0; i < pathMatrix.rows(); ++i){
                if(redundancyLabel[i]==1){
                    VectorVec3d curDendrite;
                    NGUtility::GetReversePartVectorVec3d(connectPoints[i], 0, 5, curDendrite);//2015-6-8
                    VectorVec3d tmp;
                    //tmp.push_back(soma);
                    std::copy(curDendrite.begin(), curDendrite.begin()+5, back_inserter(tmp));
                    bifurcationClustre.push_back(VectorVec3d());
                    bifurcationClustre.back().swap(tmp);

                    Vec3d curDir;
                    //VectorVec3d tmp;
                    tmp.clear();
                    std::copy(curDendrite.begin()+1, curDendrite.end(), back_inserter(tmp));//2015-6-8
                    CurveDirection(tmp, curDir);
                    bifurcationDirections.push_back(curDir);
                }
            }
        }
    }else{
        bifurcationClustre.clear();
        bifurcationDirections.clear();
    }
}

//2015-6-8
void SparseTraceFilter::ExtractCubeBoundary(const std::vector<VectorVec3d> &curveCluster, CVolume &traceLabelMatrix,
                                            std::vector<VectorVec3i> &growPointSet, VectorVec3i& boundaryLabel)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int nxx = origImg.x() - 1;
    int nyy = origImg.y() - 1;
    int nzz = origImg.z() - 1;//warning!!
    boundaryLabel.clear();
    for(size_t i0 = 0; i0 < curveCluster.size(); ++i0){
        VectorVec3d curvePoints;
        ThreeDimdCurveResample(curveCluster[i0], 2.0, curvePoints);
        std::vector<double> curvePointsRadius(curvePoints.size(), 0.0);
        for(size_t i1 = 0; i1 < curvePoints.size(); ++ i1){
            double R,W;
            CalcParmOfCurveNode(origImg, backImg, curvePoints[i1], R, W);
            curvePointsRadius[i1] = R;
        }
        //
        Vec3d curPoint;
        double pointRadius;
        int xMin, xMax, yMin, yMax, zMin, zMax;
        bool curFlag = false;
        for(size_t i2 = 0; i2 < curvePoints.size(); ++ i2){
            curPoint << NGUtility::Round(curvePoints[i2](0)), NGUtility::Round(curvePoints[i2](1)), NGUtility::Round(curvePoints[i2](2));
            pointRadius = std::max(4.0, double(NGUtility::Round(curvePointsRadius[i2] + 2.0)) + 1.0);//2016-3-13
            Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPoint[0], curPoint[1], curPoint[2],
                    pointRadius, pointRadius, pointRadius, 0, nxx, 0, nyy, 0, nzz);

            for(int i_x  = xMin; i_x <= xMax; ++i_x){
                for(int i_y = yMin; i_y <= yMax; ++i_y){
                    for(int i_z = zMin; i_z <= zMax; ++i_z){
                        //2015-6-8
                        if(origImg(i_x, i_y, i_z) - backImg(i_x, i_y, i_z) > 4.0 * std::sqrt(double(backImg(i_x, i_y, i_z))))
                            curFlag = true;
                        else
                            curFlag = false;
                        bool curFlag1 = origImg(i_x, i_y, i_z) > 1.5*double(backImg(i_x, i_y, i_z));
                        bool curFlag2 = origImg(i_x, i_y, i_z) > double(backImg(i_x, i_y, i_z)) + paramPack->diffuseValue_;
                        curFlag = curFlag || curFlag1 || curFlag2;
                        if(curFlag && 0 == traceLabelMatrix(i_x, i_y, i_z)){
                            traceLabelMatrix(i_x, i_y, i_z) = 1;//just fill traceLabelMatrix, not indexImage_;
                            boundaryLabel.push_back(Vec3i(i_x, i_y, i_z));//2015-6-8
                        }
                    }
                }
            }
        }
    }
    //
    {
        VectorVec3i tmp(boundaryLabel);
        boundaryLabel.swap(tmp);
    }
    VectorVec3i boundaryLabelCp(boundaryLabel);
    //
    growPointSet.clear();
    growPointSet.resize(6);
    for(size_t i = 0; i < 6; ++i){
        VectorVec3i curPoints;
        RegionInflationModifyV2(boundaryLabelCp, traceLabelMatrix, 15, curPoints);//2016-3-13
        if(!curPoints.empty()){
            growPointSet[i] = curPoints;
            boundaryLabelCp.swap(curPoints);
        } else{
            break;
        }
    }
}

void SparseTraceFilter::ThreeDimdCurveResample(const VectorVec3d &curvePoints, double incValue, VectorVec3d &resampleCurve)
{
    if(curvePoints.empty()){
        printf("error occured in ThreeDimdCurveResample.\n");
        LOG(ERROR) << "error occured in ThreeDimdCurveResample.";
    }
    //
    size_t curvePointsNum = curvePoints.size();
    resampleCurve.clear();
    resampleCurve.reserve(5 * curvePointsNum);
    resampleCurve.push_back(curvePoints[0]);
    std::vector<size_t> resampleIndex(5 * curvePointsNum, 0);
    resampleIndex.reserve(5 * curvePointsNum);
    resampleIndex[0] = 0;//matlab is 0, while C++ is 1
    for(size_t i = 0; i < 5* curvePointsNum - 1; ++i){
        if(resampleIndex[i] < curvePointsNum - 1){
            VectorVec3d adjacentPoints;
            adjacentPoints.assign(curvePoints.begin() + resampleIndex[i] + 1,
                                       curvePoints.begin() + std::min(resampleIndex[i] + 4, curvePointsNum));

            if(adjacentPoints.size() > 1){
                std::vector<double> vectorDist(adjacentPoints.size(), 0.0);
                for(size_t j = 0; j < adjacentPoints.size(); ++j){
                    vectorDist[j] = (adjacentPoints[j] - resampleCurve[i]).norm();
                }
                //
                std::vector<size_t> resampleFlag;
                for(size_t k = 0; k < vectorDist.size(); ++k){
                    if(vectorDist[k] > incValue) resampleFlag.push_back(k);
                }
                //
                if(resampleFlag.empty()){
                    size_t sz = std::min(resampleIndex[i] + 3, curvePointsNum - 1);//C++
                    resampleIndex[i+1] = sz;
                    resampleCurve.push_back(curvePoints[sz]);
                } else{
                    if(resampleFlag[0] == 0){
                        resampleIndex[i+1] = resampleIndex[i];
                        Vec3d curDir = (curvePoints[resampleIndex[i] + 1] - resampleCurve[i]) / vectorDist[0];
                        resampleCurve.push_back(resampleCurve[i] + 1.5 * curDir);
                    } else{
                        Vec3d curDir = adjacentPoints[resampleFlag[0]] - adjacentPoints[resampleFlag[0] - 1];
                        curDir.normalize();
                        double tmp = (adjacentPoints[resampleFlag[0] -1 ] - resampleCurve[i]).norm();
                        tmp = tmp * tmp;
                        double incStep = incValue * incValue - tmp;
                        resampleCurve.push_back(adjacentPoints[resampleFlag[0] - 1] + std::sqrt(incStep) * curDir);
                        //resampleIndex.push_back(resampleFlag[0] + resampleIndex[i] - 1);
                        resampleIndex[i+1] = resampleFlag[0] + resampleIndex[i];//2015-5-7
                    }
                }
            }
            if(adjacentPoints.size() == 1){
                resampleIndex[i+1] = curvePointsNum - 1;
                resampleCurve.push_back(adjacentPoints[0]);
            }
        }else{
            break;
        }
    }
    {
        VectorVec3d tmp(resampleCurve);
        resampleCurve.swap(tmp);
        resampleCurve.pop_back();//2015-6-8
    }
}

void SparseTraceFilter::WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,
                                std::vector<double> &rayNodeWet)
{
    typedef double spheredist;
    int nxss = int(rayNode.size());
    rayNodeWet.clear();

    typedef double spheredist;
    spheredist x,y,z;
    spheredist segmentDense, segmentWet;//, w1;
    for (int i = 0; i < nxss; ++i){
        x = rayNode[i](0);
        y = rayNode[i](1);
        z = rayNode[i](2);
        segmentDense = segmentWet = 0.0;
        ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
        rayNodeWet.push_back(segmentDense / (segmentWet + 0.0001));
    }
}

void SparseTraceFilter::CalcParmOfCurveNodeList(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
                                          const VectorVec3d &curNode,
                                          std::vector<double> &radius, std::vector<double> &rav)
{
    radius.clear();
    rav.clear();

    VectorVec3d::size_type nxx = curNode.size();
    double tmpR(0.0);
    double tmpV(0.0);
    for (VectorVec3d::size_type i = 0; i < nxx; ++i){
        Vec3d tmpdata1;
        tmpdata1(0) = NGUtility::Round(curNode[i](0));//int(curNode[i](0) + 0.5);
        tmpdata1(1) = NGUtility::Round(curNode[i](1));
        tmpdata1(2) = NGUtility::Round(curNode[i](2));
        CalcParmOfCurveNode(origImg, backImg, tmpdata1, tmpR, tmpV);
        radius.push_back(tmpR);
        rav.push_back(tmpV);
    }
}

void SparseTraceFilter::CalcParmOfCurveNode(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
                                      const Vec3d &curNode, double &radius, double &wet)
{
    int vx = backImg.x();
    int vy = backImg.y();
    int vz = backImg.z();

    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, NGUtility::Round(curNode(0) ), NGUtility::Round(curNode(1) ), NGUtility::Round(curNode(2) ),
        5, 5, 5, 0, vx - 1, 0, vy - 1, 0, vz - 1);

    int xLen = xMax - xMin + 1;
    int yLen = yMax - yMin + 1;
    int zLen = zMax - zMin + 1;

    Vec3d center(curNode(0) - xMin, curNode(1) - yMin, curNode(2) - zMin);//SL
    Volume<double> procImg;
    procImg.SetSize(xLen, yLen, zLen);//MM1

    for (int i = xMin; i <= xMax; ++i){
        for (int j = yMin; j <= yMax; ++j){
            for (int ij = zMin; ij <= zMax; ++ij){
                procImg(i-xMin, j-yMin, ij-zMin) = double(origImg(i, j, ij)) - double(backImg(i, j, ij))
                    - 3.0 * std::sqrt(double(backImg(i, j, ij)));
            }
        }
    }//for

    radius = 0.05;//r
    wet = 0.0;//v
    MatXd distWet(2, 4); distWet.setZero();
    double vv0 = 0.0;
    double vv1 = 0.0;

    for (double i = 0; i < xLen; ++i){
        for (double j = 0; j < yLen; ++j){
            for (double ij = 0; ij < zLen; ++ij){
                Vec3d dist(i - center(0), j - center(1), ij - center(2));
                double distNorm = dist.norm();
                if (distNorm <= 1.0){ //
                    distWet(0, 0) += 1.0;
                    if (procImg(i, j, ij) > 0.0){
                        vv0 += procImg(i, j, ij);
                        distWet(1, 0) += 1.0;
                    }
                }//if

                if (distNorm <= 2.0){ //
                    distWet(0, 1) += 1.0;
                    if (procImg(i, j, ij) > 0){
                        vv1 += procImg(i, j, ij);
                        distWet(1, 1) += 1.0;
                    }
                }

                if (distNorm <= 3.0){ //
                    distWet(0, 2) += 1.0;
                    if (procImg(i, j, ij) > 0)
                        distWet(1, 2) += 1.0;
                }

                if (distNorm <= 4.0){ //
                    distWet(0, 3) += 1.0;
                    if (procImg(i, j, ij) > 0)
                        distWet(1, 3) += 1.0;
                }
            }//for
        }
    }//for

    Vec4d procDistWet = distWet.row(1).array() / ( distWet.row(0).array() + 0.0001);
    for (int i = 1; i <= 4; ++i){
        if (procDistWet(4-i) > 0.5){
            radius = 5.0 - (double)i;
            break;
        }
    }

    if (radius > 2.1){//warning !!
        radius = radius / std::pow(distWet(0, int(radius + 0.5) - 1) / distWet(1, int(radius + 0.5) - 1), 1.0/3.0);
        wet = vv1 / distWet(1, 1);
    }

    if (radius < 2.1){
        if (distWet(1,1) > 0.0){
            radius = 2.0 / std::pow(distWet(0, 1) / distWet(1, 1), 1.0/3.0);
            wet = vv1 / distWet(1, 1);
        }
        if (std::abs(distWet(1, 1) - 0.0) < EXPO && distWet(1, 0) > 0.0 ){
            radius = 1.0 / std::pow(distWet(0, 0) / distWet(1, 0), 1.0/3.0);
            wet = vv0 / distWet(1, 0);
        }
    }
}

void SparseTraceFilter::TraceNextCurveNodeV2(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode,
                                         const double &threv, const Vec3d &initDir, bool flag,
                                         Vec3d &nextCurveNode, Vec3d &nextDenDir,  int &isNextNodeValid)
{
    nextCurveNode.setZero();
    isNextNodeValid = 1;
    nextDenDir = initDir;

    VectorVec3d neighborPtSet;//
    std::vector<double> neighborWet;//W

    Vec3d x10;//x10
    x10.setZero();
    Vec3d x11;//x11
    x11.setZero();
    //2015-6-8 add flag
    CalcNeighborSignalV2(locOrigImg, curSeedNode, nextDenDir, threv, flag, neighborPtSet, neighborWet, x10, x11);

    size_t neighborWetNum = neighborPtSet.size();
    /*W1 = sort(W, 'desend')*/
    std::vector<double> sortedNeighborWet;//W1
    sortedNeighborWet.assign(neighborWet.begin(), neighborWet.end());
    std::sort(sortedNeighborWet.begin(), sortedNeighborWet.end(), std::greater<double>());

    size_t dsw = (std::min)((size_t)10, sortedNeighborWet.size());//2015-7-7
    //size_t dsw = (std::min)((size_t)11, sortedNeighborWet.size());//2015-6-15 TDI079
    /**/
    double thrdk(0.0);
    if (dsw > 5){//2015-7-7
    //if(dsw > 10){//2015-6-15 TDI079
        thrdk = 0.0;
        for(int i = 0; i < int(dsw); ++i){
            thrdk += sortedNeighborWet[i];
        }
        thrdk /= double(dsw);
    }
    else thrdk = 0.0;

    Mat3d sigmaH;//sigmaH
    sigmaH.setZero();
    Vec3d kk;//kk
    kk.setZero();

    double P = 0.0;
    //2015-6-8
    if (neighborWetNum > std::max(0.1 * threv, 10.0) &&
        thrdk > std::max(std::min(std::min(0.6 * threv, 4.5 * std::sqrt(threv)), paramPack->diffuseValue_), paramPack->traceValue_)){
    /*if (neighborWetNum > std::max(0.1 * threv, 10.0) &&
        thrdk > std::max(std::min(std::min(60 * threv, 20.0 * std::sqrt(threv)), paramPack->diffuseValue_), paramPack->traceValue_)){*/

        //2015-6-8 thrdk>max(min([0.6*threv,4.5*sqrt(threv),threValue]),1)
        CalcConstraintPCA(neighborPtSet, curSeedNode, 0.150 * Eigen::Matrix3d::Identity(), neighborWet, P, sigmaH, kk);
        if (P>0){
            if ( std::abs(sigmaH.determinant() / ((sigmaH.inverse()).determinant()) + 0.0001) < 10.0 ) {
                CalcPCADirections(NGUtility::sign(sigmaH.determinant()) * sigmaH.inverse(), initDir.normalized(),
                              x11, 0.8, nextDenDir);//2015-6-8
            }
        } else{
            isNextNodeValid=0;
        }
        Mat3d tmpResultDirections;
        Vec3d x2, x3;
        CalcOrthoBasis(nextDenDir, x2, x3);
        tmpResultDirections.col(0) = x3;
        tmpResultDirections.col(1) = x2;
        tmpResultDirections.col(2) = nextDenDir;
        //vmmk = nextDenDir;
        nextCurveNode = - (tmpResultDirections.col(0) * tmpResultDirections.col(0).transpose() + tmpResultDirections.col(1) * tmpResultDirections.col(1).transpose()) * kk + curSeedNode;
    }
    else isNextNodeValid = 0;
}

//2015-6-8 add flag
void SparseTraceFilter::CalcNeighborSignalV2(const Volume<unsigned short> &origImg,
                                            const Vec3d &curSeedNode, const Vec3d &initVec, const double threv, bool flag,
                                            VectorVec3d &neighborPtSet, std::vector<double> &neighborWet,
                                            Vec3d &firDir, Vec3d &secDir)
{
    neighborPtSet.clear();
    neighborWet.clear();
    firDir.setZero();
    secDir.setZero();
    int vx = origImg.x();
    int vy = origImg.y();
    int vz = origImg.z();

    /**/
    int xMin = std::max(NGUtility::Round(curSeedNode(0) - 7.0 ), 0);
    int xMax = std::min(NGUtility::Round(curSeedNode(0) + 7.0), vx - 1);
    int yMin = std::max(NGUtility::Round(curSeedNode(1) - 7.0 ), 0);
    int yMax = std::min(NGUtility::Round(curSeedNode(1) + 7.0 ), vy - 1);
    int zMin = 0;//std::max(NGUtility::Round(curSeedNode(2) - 6.0 ), 0);//2015-6-8
    int zMax = vz-1;//std::min(NGUtility::Round(curSeedNode(2) + 6.0 ), vz - 1);

    Vec3d curCenter(curSeedNode(0) - (double)xMin,
        curSeedNode(1) - (double)yMin, curSeedNode(2) - (double)zMin);//ML1
    /**/
    Volume<unsigned short> locOrigImg;
    //Change : 2014-3-19-10-16
    locOrigImg.SetSize(xMax - xMin + 1, yMax - yMin + 1, zMax - zMin + 1);
    //const Volume<unsigned short> &pptr = origImg;//
    for (int i = xMin; i <= xMax; ++i){
        for (int j = yMin; j <= yMax; ++j){
            for (int ij = zMin; ij <= zMax; ++ij){
                //int nima = origImg(i, j, ij);
                locOrigImg(i - xMin, j - yMin, ij - zMin) = origImg(i, j, ij);
            }
        }
    }
    /**/
    int sx = xMax - xMin + 1;
    int sy = yMax - yMin + 1;
    int sz = zMax - zMin + 1;
    //Volumn subVol(xMax - xMin + 1, yMax - yMin + 1, vol.z);//nx1,ny1,nz1 zMax - zMin + 1
    /**/
    int minLen = (sx) < (sy) ?
        ((sx) < (sz) ? (sx) : (sz)) :
        ((sy) < (sz) ? (sy) : (sz));

    double ii(0.0);
    double jj(0.0);
    double angleDiff(0.0);

    //double angleThrev1=-0.9;//2015-2-7
    //double angleThrev2= 0.85;
    //2015-6-8
    double angleThrev1;
    double angleThrev2;
    if(flag){//==1
        angleThrev1=-0.95;
        angleThrev2=0.95;
    }else{//==0
        angleThrev1=-0.9;
        angleThrev2=0.85;
    }

    if (minLen > 3 && vz > 3){
        /*AAj = zeros(4,216)*/
        VectorVec4d backwardNodeSet;//AAj
        for (int i = 0; i < 36; ++i){
            for ( int j = 0; j < 18; ++j){
                ii = (double)i * M_PI / 18.0;
                jj = (double)j * M_PI / 18.0;
                Vec3d polarPosition((std::sin(jj)) * (std::cos(ii)), (std::sin(ii)) * (std::sin(jj)), (std::cos(jj)));
                angleDiff = polarPosition.dot(initVec);/*jd*/
                if (angleDiff < angleThrev1) //if jd < -0.85
                {
                    VectorVec3d rayNode;
                    rayNode.push_back(curCenter + polarPosition);
                    rayNode.push_back(curCenter + 2 * polarPosition);
                    rayNode.push_back(curCenter + 3 * polarPosition);
                    rayNode.push_back(curCenter + 4 * polarPosition);
                    rayNode.push_back(curCenter + 5 * polarPosition);
                    std::vector<double> rayNodeWet;
                    WeighRayValue(rayNode, locOrigImg, rayNodeWet);
                    double meanRayNodeWet = std::accumulate(rayNodeWet.begin(), rayNodeWet.begin() + 3, 0.0);
                    meanRayNodeWet /= 3.0;//akmmv.size();
                    backwardNodeSet.push_back(Vec4d(polarPosition(0), polarPosition(1), polarPosition(2), meanRayNodeWet));
                }//if
            }
        }//for
        /**/
        if (backwardNodeSet.size() < 2) firDir =  - initVec;
        else{
#ifdef _WIN32
            VectorVec4d::iterator maxItem = std::max_element(backwardNodeSet.begin(), backwardNodeSet.end(),
                Vec4d_3th_less());//Node_4TH_MAX
#else
            VectorVec4d::iterator maxItem = std::max_element(backwardNodeSet.begin(), backwardNodeSet.end(),
                [](const Vec4d& lhs, const Vec4d& rhs){
                    return lhs(3) < rhs(3);
            });//Node_4TH_MAX
#endif
            
            Vec3d maxVec((*maxItem)(0), (*maxItem)(1), (*maxItem)(2));//xss
            firDir = maxVec.normalized();//x10 =
        }
    }
    else firDir = - initVec;

    if (minLen > 3 && vz > 3){
        /*AAj2 = zeros(4,216)*/
        VectorVec4d forwardNodeSet;//AAj2
        for (int i = 0; i < 36; ++i){
            for ( int j = 0; j < 18; ++j){
                ii = (double)i * M_PI / 18.0;
                jj = (double)j * M_PI / 18.0;
                Vec3d ds1((std::sin(jj)) * (std::cos(ii)), (std::sin(ii)) * (std::sin(jj)) , (std::cos(jj)));//ds1
                angleDiff = ds1.dot(initVec);/*jd*/
                if (angleDiff > angleThrev2){ //if jd > -0.85
                    VectorVec3d rayNode;
                    rayNode.push_back(curCenter + ds1);
                    rayNode.push_back(curCenter + 2 * ds1);
                    rayNode.push_back(curCenter + 3 * ds1);
                    rayNode.push_back(curCenter + 4 * ds1);
                    rayNode.push_back(curCenter + 5 * ds1);
                    rayNode.push_back(curCenter + 6 * ds1);
                    std::vector<double> rayNodeWet;
                    WeighRayValue(rayNode, locOrigImg, rayNodeWet);
                    double meanRayNodeWet = accumulate(rayNodeWet.begin(), rayNodeWet.begin()+3, 0.0);
                    meanRayNodeWet /= 3.0;//akmmv.size();
                    forwardNodeSet.push_back(Vec4d(ds1[0], ds1[1], ds1[2], meanRayNodeWet));
                }
            }
        }//for
        /*exception*/
        if (forwardNodeSet.size() < 2) secDir = initVec;
        else{
#ifdef _WIN32
            VectorVec4d::iterator maxItem = max_element(forwardNodeSet.begin(), forwardNodeSet.end(), Vec4d_3th_less());//Node_4TH_MAX
#else
            VectorVec4d::iterator maxItem = max_element(forwardNodeSet.begin(), forwardNodeSet.end(), [](const Vec4d& lhs, const Vec4d& rhs){
                return lhs(3) < rhs(3);
            });//Node_4TH_MAX
#endif
            
            Vec3d maxVec((*maxItem).x(), (*maxItem).y(), (*maxItem).z());//xss
            secDir = maxVec.normalized();
        }
    }
    else secDir = initVec;

    /*dataSS = zeros(4,8*125);//*/
    VectorVec4d rayArea;//dataSS,dataS
    for (int i = 0; i < 3; ++i){//2014-4-21
        //Vec3d data10(curSeedNode(0), curSeedNode(1), curSeedNode(2));
        //Vec3d ray = i * firDir + data10; //data10 + i * x10
        Vec3d rayNode = i * firDir + curSeedNode;
        int minX = std::max(std::min(NGUtility::Round(rayNode(0) - 4.0), vx - 1), 0);//2016-3-17
        int maxX = std::min(std::max(NGUtility::Round(rayNode[0] + 4.0), 0), vx - 1);
        int minY = std::max(std::min(NGUtility::Round(rayNode[1] - 4.0), vy - 1), 0);
        int maxY = std::min(std::max(NGUtility::Round(rayNode[1] + 4.0), 0), vy - 1);
        int minZ = std::max(std::min(NGUtility::Round(rayNode[2] - 4.0), vz - 1), 0);
        int maxZ = std::min(std::max(NGUtility::Round(rayNode[2] + 4.0), 0), vz - 1);
        //Volumn localVol(maxX - minX + 1, maxY - minY + 1, maxZ - minZ + 1);
        int lx = maxX - minX + 1;
        int ly = maxY - minY + 1;
        int lz = maxZ - minZ + 1;
        /**/
        Vec4d tmpNode;
        for (int ix = 0; ix < lx; ++ix){
            for (int iy = 0; iy < ly; ++iy){
                for (int iz = 0; iz < lz; ++iz){
                    int newX = minX + ix;
                    int newY = minY + iy;
                    int newZ = minZ + iz;
                    tmpNode << newX, newY , newZ , std::max(double(origImg(newX, newY, newZ)) - threv, 0.0) ;
                    rayArea.push_back(tmpNode);
                }
            }
        }
    }
    /**/
    for (int i = 1; i < 3; ++i){
        //Vec3d data10(curSeedNode(0), curSeedNode(1), curSeedNode(2));
        //Vec3d ray = i * secDir + data10; //data10 + i * x11
        Vec3d rayNode = i * secDir + curSeedNode;
        int minX = std::max(std::min(int(rayNode[0] - 4.0 + 0.5), vx - 1), 0);//2016-3-17
        int maxX = std::min(std::max(int(rayNode[0] + 4.0 + 0.5), 0), vx - 1);
        int minY = std::max(std::min(int(rayNode[1] - 4.0 + 0.5), vy - 1), 0);
        int maxY = std::min(std::max(int(rayNode[1] + 4.0 + 0.5), 0), vy - 1);
        int minZ = std::max(std::min(int(rayNode[2] - 4.0 + 0.5), vz - 1), 0);
        int maxZ = std::min(std::max(int(rayNode[2] + 4.0 + 0.5), 0), vz - 1);
        //Volumn localVol(maxX - minX + 1, maxY - minY + 1, maxZ - minZ + 1);
        int lx = maxX - minX + 1;
        int ly = maxY - minY + 1;
        int lz = maxZ - minZ + 1;
        /**/
        Vec4d tmpNode;
        for (int ix = 0; ix < lx; ++ix){
            for (int iy = 0; iy < ly; ++iy){
                for (int iz = 0; iz < lz; ++iz){
                    int newX  = minX + ix;
                    int newY = minY + iy;
                    int newZ = minZ + iz;
                    tmpNode << newX, newY , newZ , (std::max)((double)origImg(newX, newY, newZ) - threv, 0.0) ;
                    rayArea.push_back(tmpNode);
                }
            }
        }
    }

    if (!rayArea.empty()){
        /*rayArea*/
        //sort from x-y-z-v
#ifdef _WIN32
        std::sort(rayArea.begin(), rayArea.end(), Vec4d_0123less());//Node_v_x_y_z_MAX
#else
        std::sort(rayArea.begin(), rayArea.end(), [](const Vec4d& lhs, const Vec4d& rhs){
            if(lhs(0) != rhs(0)) return lhs(0) < rhs(0);
            else if(lhs(1) != rhs(1)) return lhs(1) < rhs(1);
            else if(lhs(2) != rhs(2)) return lhs(2) < rhs(2);
            return lhs(3) < rhs(3);
        });//Node_v_x_y_z_MAX
#endif
        
        rayArea.erase(std::unique(rayArea.begin(), rayArea.end()), rayArea.end());
        Vec3d tmp;
        for (VectorVec4d::size_type i = 0; i < rayArea.size(); ++i){
            neighborWet.push_back(rayArea[i](3));
            //tmp(0) = (rayArea[i](0), rayArea[i](1), rayArea[i](2));
            tmp << rayArea[i](0), rayArea[i](1), rayArea[i](2);
            neighborPtSet.push_back(tmp);
        }
    }
    else{
        neighborWet.clear();
        neighborPtSet.clear();
    }
}

//2015-6-8
void SparseTraceFilter::CalcConstraintPCA(const VectorVec3d &neighborPtSet, const Vec3d &curSeedNode, const Mat3d &convMat,
                           const std::vector<double> &neighborWet,
                           double &P, Mat3d &sigmaH, Vec3d &mdata3)
{
    Vec3d g;//g
    Mat3d H;//H HH
    g.setZero();
    H.setZero();

    VectorVec3d::size_type length = neighborPtSet.size();//nx
    Vec3d X(curSeedNode(0), curSeedNode(1), curSeedNode(2)); //X
    std::vector<double> C;//C
    VectorVec3d U;//U

    Mat3d MMs;//MMs
    MMs.setZero();
    Vec3d MMs1;//MMs1
    MMs1.setZero();
    Vec3d dd(0.0, 0.0, 0.0);//dd
    //vec3d tmpCenter;//data1
    for (VectorVec3d::size_type i = 0; i < length; ++i){
        dd = X - neighborPtSet[i];//dd=X-data1(:,i)
        //dd(2) *= 2.0;//2015-6-8
        if (dd.norm() < 14.5){//2015-6-8
            C.push_back(neighborWet[i] * std::exp(-0.5 * dd.transpose() * convMat * dd));
        }
        else{
            C.push_back(0.0);//C(:,i)
        }

        U.push_back(convMat * dd);//U(:,i) = data2 *dd
        g += C[i] * U[i];//g
        H += C[i] * (U[i] * U[i].transpose() - convMat);//HH
        MMs += C[i] * convMat;//MMs
        MMs1 += (C[i] * convMat) * dd;//MMs1
    }

    //P = std::accumulate(C.begin(), C.end(), 0.0);//P = sum(C)
    P = 0.0;
    for(std::vector<double>::size_type i = 0; i < C.size(); ++i){
        P += C[i];
    }
    g *= -1.0;

    /*sigmaH = -H./P + (P^(-2)) * g * g'*/
    sigmaH = - H / P + std::pow(P, -2.0) * (g * g.transpose());
    if (std::abs(MMs.determinant()) > 0.1){
        Mat3d tmpcorrection = MMs.inverse();
        mdata3 = tmpcorrection * MMs1;//inv(MMs)*MMs1
    } else{
        mdata3 = MMs1;
    }
}

void SparseTraceFilter::CalcPCADirections(const Mat3d &sigmaH, const Vec3d &initVec, const Vec3d &T2, const double threv,
                                    Vec3d &vec1)
{
    Vec3d vec0 = initVec;//x0 = T1
    for (int i = 1; i < 16; ++i){
        CalcPCAMainDirection(vec0, sigmaH, initVec, T2, threv, vec1);//gproPCAA
        vec0 = vec1;
        double angleDiff = vec0.transpose() * initVec;
        if (angleDiff < 0.7) break;
    }
    //directionc(vec1, vec2, vec3);//
}

void SparseTraceFilter::CalcPCAMainDirection(const Vec3d &x0, const Mat3d &sigmaH,
                                       const Vec3d &T1, const Vec3d &T2, const double threv, Vec3d &x1)
{
    double init = 0.5;//a
    Vec3d g = sigmaH * x0;//g
    Mat3d A;
    A.col(0) = x0; A.col(1) = -T1; A.col(2) = -T2;

    Mat3d H;
    H.setZero();

    H(0, 0) = x0.transpose() * x0 - 1.0;//x0' * x0 - 1
    H(1, 1) = - T1.transpose() * x0 + threv;
    H(2, 2) = - T2.transpose() * x0 + threv;

    Mat3d ss1 = (A.transpose() * A - H);
    Mat3d ss = ss1.inverse();
    Mat3d B = ss * A.transpose();
    Vec3d U = B * g;
    Mat3d P = Mat3d::Identity() - (A * ss) * A.transpose();
    if ( ((P * g).cwiseAbs()).sum() <0.001 &&
        U.minCoeff()  > -0.0001){
        x1 = x0;
    }
    else{
        Vec3d v; v.setZero();
        for (int n = 0; n < 3; ++n){//2 previously
            if (U(n) < -0.0001) {
                v(n) = U(n);
            }
            else{
                v(n) = -H(n, n);
            }
        }//for

        Vec3d S = P * g + B.transpose() * v;
        /*p1 = (1-a) * g' * S / (abs(U' * ones(2,1)) + 1)*/
        double tp1 = (1.0 - init) * g.transpose() * S;
        double tp2 = std::abs(U.transpose() * Vec3d::Ones()) + 1.0;
        double p1 = tp1 / tp2;
        Vec3d d = S - p1 * (B.transpose() * Vec3d::Ones());
        d /= d.cwiseAbs().maxCoeff();

        int i = 1;
        for (i = 1; i <= 20; ++i){
            x1 = x0 + pow(0.5, (double)i) * d;
            x1.normalize();
            double tLs1 = x1.transpose() * sigmaH * x1;
            double tLs2 =  x0.transpose() * sigmaH * x0;
            double Ls = tLs1 - tLs2;
            double tdkkss1 = -T1.transpose() * x1;
            tdkkss1 /= x1.norm();
            double dkkss1 = tdkkss1 + threv + 0.001;
            double dkkss2 = pow(0.5, double(i+1)) * g.transpose() * d;
            double tdkkss3 = -T2.transpose() * x1;
            tdkkss3 /= x1.norm();
            double dkkss3 = tdkkss3 + threv + 0.001;
            if (Ls > dkkss2 && x1.norm() < 1.0001 && x1.norm() > 0.999 && dkkss1 < 0.0 &&
                 dkkss3 < 0.0) {
                break;
            }
        }
        if ( i == 20)  x1 = x0;
    }//else
}

//2015-6-8 traceLabelMatrix is global varient
void SparseTraceFilter::IsCurvesCollideV2(//const Volume<int> &traceLabelMatrix,
                                  const Vec5d &seed, int& isNextNodeValid)
{
    int nxx = traceLabelMatrix_.x() - 1;
    int nyy = traceLabelMatrix_.y() - 1;
    int nzz = traceLabelMatrix_.z() - 1;//warning!!
    int id1 = NGUtility::Round(seed(0));
    int id2 = NGUtility::Round(seed(1));
    int id3 = NGUtility::Round(seed(2));
    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, id1, id2, id3, 2,2,2,0, nxx, 0,nyy,0, nzz);
    std::vector<int> indexWindowVal;
    int kk11 = 0;
    for (int i = xMin; i <= xMax; ++i){
        for (int j = yMin; j <= yMax; ++j){
            for (int k = zMin; k <= zMax; ++k){
                if ( 0 != traceLabelMatrix_(i, j, k) ) indexWindowVal.push_back( traceLabelMatrix_(i, j, k) );
                if( traceLabelMatrix_(i,j,k) < 0) ++kk11;
            }
        }
    }//for
    //2014-4-24
    if ( indexWindowVal.size() >  31 || kk11 > 60 ){
        isNextNodeValid = 0;
    } else{
        isNextNodeValid = 1;
    }
}

//2015-6-8 traceLabelMatrix is global varient, return collide index
void SparseTraceFilter::IsCurvesCollideV3(const Vec5d &seed, int& isNextNodeValid, int& collideID)
{
    int nxx = indexImg_.x() - 1;
    int nyy = indexImg_.y() - 1;
    int nzz = indexImg_.z() - 1;//warning!!
    int id1 = NGUtility::Round(seed(0));
    int id2 = NGUtility::Round(seed(1));
    int id3 = NGUtility::Round(seed(2));
    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, id1, id2, id3, 6,6,6,0, nxx, 0,nyy,0, nzz);//2015-9-19
    std::vector<int> indexWindowVal;
    int kk11 = 0;
    for (int i = xMin; i <= xMax; ++i){
        for (int j = yMin; j <= yMax; ++j){
            for (int k = zMin; k <= zMax; ++k){
                if ( 0 != indexImg_(i, j, k) ) indexWindowVal.push_back( indexImg_(i, j, k) );
                if( indexImg_(i,j,k) < 0) ++kk11;//it is impossible
            }
        }
    }//for
    //2015-6-8
    collideID = 0;
    if ( indexWindowVal.size() >  18 || kk11 > 60 ){//2015-9-19
        isNextNodeValid = 0;
        int maxVal = *(std::max_element(indexWindowVal.begin(), indexWindowVal.end()));
        std::vector<int> countArray(maxVal,0);
        for(size_t i = 0; i < indexWindowVal.size(); ++i){
            if(indexWindowVal[i] > 0)//no value is less than 0
                ++countArray[indexWindowVal[i] - 1 ];
        }
        collideID = *(std::max_element(indexWindowVal.begin(), indexWindowVal.end()));
    } else{
        isNextNodeValid = 1;
    }
}

void SparseTraceFilter::CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir)
{
    int nx = int(ptCurve.size());
    Vec3d tmpFir;
    tmpFir.setZero();
    double tmp(0.0);
    Vec3d tmpVec;
    for (int i = 0; i < nx - 1; ++i){
        tmpVec = ptCurve[i + 1] - ptCurve[i];
        tmp = tmpVec.norm();
        tmpFir += tmp * tmpVec;
    }
    fir = tmpFir.normalized();
    //CalcOrthoBasis(fir, x2, x3);
}

void SparseTraceFilter::CalcOrthoBasis(const Vec3d &vec1, Vec3d &vec2, Vec3d &vec3)
{
    vec2 = vec1;
    vec3.setZero();
    std::vector<std::pair<int, double> > tmp;
    tmp.push_back(std::pair<int, double>(0, std::abs(vec1(0))));
    tmp.push_back(std::pair<int, double>(1, std::abs(vec1(1))));
    tmp.push_back(std::pair<int, double>(2, std::abs(vec1(2))));
    /*[idxv,idexx]=sort(abs(vec1))*/
#ifdef _WIN32
    std::sort(tmp.begin(), tmp.end(), Pair_less());
#else
    std::sort(tmp.begin(), tmp.end(),[](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){
        return lhs.second < rhs.second;
    });
#endif
    
    double cdd = (tmp[0].second * tmp[0].second + tmp[1].second * tmp[1].second) / tmp[2].second;
    vec2(tmp[2].first) = -(double)NGUtility::sign(vec1[tmp[2].first]) * cdd;

    if (std::abs(cdd - 0.0) > EXPO)
        vec2.normalize();//cdd != 0.0
    else vec2[tmp[1].first] = 1.0;

    /*3th primary element*/
    vec3(tmp[0].first) = -1.0;
    Mat2d denMat;
    denMat << vec1(tmp[1].first), vec1(tmp[2].first), vec2(tmp[1].first), vec2(tmp[2].first);

    //Matrix2f num1;
    Mat2d num1;
    num1 << vec1(tmp[0].first), vec1(tmp[2].first), vec2(tmp[0].first), vec2(tmp[2].first);

    //Matrix2f num2;
    Mat2d num2;
    num2 << vec1(tmp[1].first), vec1(tmp[0].first), vec2(tmp[1].first), vec2(tmp[0].first);

    double den = denMat.determinant();
    vec3(tmp[1].first) = num1.determinant() / den;
    vec3(tmp[2].first) = num2.determinant() / den;
    vec3.normalize();
}

void SparseTraceFilter::ClustersExtraction(CVolume& growLabelMatrix, const VectorVec3i& inflatPtSet,
                                           VectorVec3i& subClustrePtSet, std::vector<int>& subClustrePtSetNum)
{
    subClustrePtSet.clear();
    subClustrePtSetNum.clear();
    for(size_t i = 0; i < inflatPtSet.size(); ++i){
        if(growLabelMatrix(inflatPtSet[i](0), inflatPtSet[i](1), inflatPtSet[i](2)) == 1){
            VectorVec3i extractedPtSet;
            ExtractSubRegionOP(inflatPtSet[i], extractedPtSet, growLabelMatrix);
            size_t sz = subClustrePtSet.size();
            subClustrePtSet.resize(sz + extractedPtSet.size());
            std::copy(extractedPtSet.begin(), extractedPtSet.end(), subClustrePtSet.begin() + sz);
            subClustrePtSetNum.push_back(int(extractedPtSet.size()));
        }
    }
}

void SparseTraceFilter::ExtractSubRegionOP(const Vec3i &currentPoint, VectorVec3i &extractedPtSet, CVolume &growLabelMatrix)
{
    int i = 0;
    extractedPtSet.clear();//
    VectorVec3i tempPos1;		//
    VectorVec3i tempPos2;		//
    tempPos2.push_back(currentPoint);			//
    ///--------------------------------
    growLabelMatrix( currentPoint(0) , currentPoint(1) ,  currentPoint(2) ) = 0;
    extractedPtSet.push_back(currentPoint);//
    tempPos1.push_back(currentPoint);//
    ///---------------begin-------------------------
    while(!tempPos2.empty() && i < 1280 ){
        ++i;
        ///////////////growing area////////////////////////
        ExtractSubRegionAboveThrevModify( tempPos1, tempPos2, growLabelMatrix);
        //
        VectorVec3i::size_type sz = extractedPtSet.size();
        extractedPtSet.resize(sz + tempPos2.size());
        std::copy(tempPos2.begin(), tempPos2.end(), extractedPtSet.begin() + sz);
        tempPos1 = tempPos2;
    }
}

void SparseTraceFilter::ExtractSubRegionAboveThrevModify(const VectorVec3i &currentPtSet, VectorVec3i &extractedPtSet,
                                                         CVolume &indexImg)
{
    int i,j,ij;
    int n, t_Sum;;
    int xMin,xMax;
    int yMin,yMax;
    int zMin,zMax;
    int nx	= indexImg.x();
    int ny	= indexImg.y();
    int nz	= indexImg.z();
    int nxx	= int(currentPtSet.size());

    Vec3i pt;
    pt.setZero();
    extractedPtSet.clear();//must clear
    ///-------------------------------------
    for (n = 0 ; n < nxx; ++n){
        const Vec3i& it = currentPtSet[n];
        int id1 = it(0);
        int id2 = it(1);
        int id3 = it(2);

        if (id1 > 1 && id1 < nx-2 &&id2 > 1 && id2 < ny-2 && id3 > 1 && id3 < nz-2){
            xMin = id1 -1; xMax = id1 +1;
            yMin = id2 -1; yMax = id2 +1;
            zMin = id3 -1; zMax = id3 +1;

            //
            t_Sum = 0;
            for (i = xMin; i <= xMax; ++i)
                for (ij = zMin; ij <= zMax; ++ij)
                    for (j = yMin; j <= yMax; ++j){
                        t_Sum += indexImg(i , j , ij );
                    }
            //t_Sum /= 255;indexImg value is 0-1, not 0-255

            //extract
            if (t_Sum > 0){//2015-2-7
                for (i = xMin; i <= xMax; ++i)
                    for (j = yMin; j <= yMax; ++j)
                        for (ij = zMin; ij <= zMax; ++ij)
                            if (indexImg(i ,j , ij ) > 0){//-2015-2-9
                                indexImg(i ,j , ij ) = 0;
                                pt << i, j, ij;
                                extractedPtSet.push_back(pt);
                            }
                        }
        }
    }//for
}

//2015-6-8
void SparseTraceFilter::GetRayLimitV2(const Volume<double> &sphereRayWet, const Volume<double> &sphereBackWet,
                                      const double constrictionThrev, std::vector<std::vector<double> > &rayLimit)
{
    int nx = sphereRayWet.x();
    int ny = sphereRayWet.y();
    int nz = sphereRayWet.z();
    //std::vector<double> sphereRay;
    //std::vector<double> sphereBackRay;
    std::vector<double> tmpRay;
    std::vector<double> tmp;
    int arra(0);
    double value;

    for (int i = 0; i < ny; ++i){
        tmp.clear();
        for (int j = 0; j < nz; ++j){
            for (int ij = 0; ij < nx; ++ij){
                //sphereRay.push_back(sphereRayWet(ij, i, j));
                //sphereBackRay.push_back(sphereBackWet(ij, i, j));
                value = sphereRayWet(ij, i, j) - sphereBackWet(ij, i, j) -
                        constrictionThrev * std::sqrt(sphereBackWet(ij, i, j));
                tmpRay.push_back(value);
            }

            ContourUtil::CalculateOneRayLimit(tmpRay, 0, arra);
            tmp.push_back(arra);
            tmpRay.clear();
            //sphereBackRay.clear();
        }
        rayLimit.push_back(tmp);
    }
}

//2015-6-8
void SparseTraceFilter::ConnectPathInTree(const std::vector<std::vector<int> > &subRegionCon, MatXi &pathMatrix)
{
    std::vector<std::vector<int> > subRegionConCp(subRegionCon);
    int nxx = int(subRegionConCp.size());
    int pathNum = int(subRegionConCp.back().size());
    std::reverse(subRegionConCp.begin(), subRegionConCp.end());
    pathMatrix = MatXi(pathNum, nxx +1);
    pathMatrix.setZero();
    for(int i = 0; i < pathNum; ++i){
        pathMatrix(i,0) = i+1;//warning!
        for(int j = 1; j < nxx + 1; ++j){
            const std::vector<int>& curData = subRegionConCp[j-1];
            pathMatrix(i,j) = curData[pathMatrix(i,j-1)-1];//C++
        }
    }
    MatXi pathMatrixCp = pathMatrix.rowwise().reverse();//resolve eigen alias
    pathMatrix = pathMatrixCp;
}

//2015-6-8
void SparseTraceFilter::DeleteRedundancyPath(const MatXi &pathMatrix, std::vector<int> &redundancyLabel)
{
    int pathNum = pathMatrix.rows();
    redundancyLabel.resize(pathNum, 1);
    for(int i = 0; i < pathNum; ++i){
        if(redundancyLabel[i] == 1){//warning!
            Eigen::RowVector3i referenceData = pathMatrix.block(i,3,1,3);
            for(int j = i+1; j < pathNum; ++j){
                Eigen::RowVector3i curData = pathMatrix.block(j,3,1,3);
                if( (referenceData-curData).array().abs().sum() == 0)
                    redundancyLabel[j]=0;
            }
        }
    }
}

//2015-6-8
void SparseTraceFilter::DeleteRedundancyPathV2(const MatXi &pathMatrix, std::vector<int> &redundancyLabel)
{
    int pathNum = pathMatrix.rows();
    redundancyLabel.resize(pathNum, 1);
    for(int i = 0; i < pathNum; ++i){
        if(redundancyLabel[i] == 1){//warning!
            Eigen::RowVector3i referenceData = pathMatrix.block(i,1,1,3);//different from DeleteRedundancyPath
            for(int j = i+1; j < pathNum; ++j){
                Eigen::RowVector3i curData = pathMatrix.block(j,1,1,3);
                if( (referenceData-curData).array().abs().sum() == 0)
                    redundancyLabel[j]=0;
            }
        }
    }
}

void SparseTraceFilter::Ind2CellInd(std::vector<std::vector<VectorVec3d> >& treeDataSet, int ind, int &cellLevel, int &cellInd)
{
    cellLevel = cellInd = 0;
//    int c_CellInd = cellInd - 1;//for C++ index, the index should minus 1
//    if(c_CellInd < 0){
//        printf("error in Ind2CellInd.\n");
//        return;
//    }
    if(ind < 0){
        printf("error in Ind2CellInd.\n");
        LOG(ERROR) <<"error in Ind2CellInd.";
        return;
    }

    std::vector<int> numList(1,0);
    for(size_t i = 0; i < treeDataSet.size(); ++i){
        numList.push_back(int(treeDataSet[i].size()));
    }
    std::vector<int> cumList(numList.size());
    std::partial_sum(numList.begin(), numList.end(), cumList.begin());

    for(size_t i = 0; i < cumList.size(); ++i){
        if(ind <= cumList[i]){//cumList[0] = 0, so i should be greater than 0
            cellLevel = int(i)-1;
            cellInd = ind - cumList[i-1] -1;//for C++, index - 1
            break;
        }
    }
}

void SparseTraceFilter::ExtractLocalDomainV2(const Vec3d &initPoint, const SVolume &origImg,
                                       SVolume &locOrigImg, Vec3d &locPoint)
{
    //Volumn vol(v.Width(), v.Height(), v.Depth());
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    int minX = std::max(NGUtility::Round(initPoint(0) - 40.0 ), 0);
    int maxX = std::min(NGUtility::Round(initPoint(0) + 40.0 ), nx - 1);
    int minY = std::max(NGUtility::Round(initPoint(1) - 40.0 ), 0);
    int maxY = std::min(NGUtility::Round(initPoint(1) + 40.0 ), ny - 1);
    int minZ = std::max(NGUtility::Round(initPoint(2) - 40.0 ), 0);//2015-6-18
    int maxZ = std::min(NGUtility::Round(initPoint(2) + 40.0 ), nz - 1);

    locPoint = Vec3d(initPoint(0) - minX, initPoint(1) - minY, initPoint(2) - minZ);

    //Volumn SubVol(maxX - minX + 1, maxY - minY + 1, maxZ - minZ + 1);
    int sx = maxX - minX + 1;
    int sy = maxY - minY + 1;
    int sz = maxZ - minZ + 1;

    locOrigImg.SetSize(sx, sy, sz);
    for (int i = minX; i <= maxX; ++i){
        for (int j = minY; j <= maxY; ++j){
            for (int ij = minZ; ij <= maxZ; ++ij)
                locOrigImg(i - minX, j - minY, ij - minZ) = origImg(i, j, ij);
        }
    }//for
}

void SparseTraceFilter::SetInitialDirection( const Vec3d& arg)
{
    initialDirection_ = arg;
}

//2016-3-13
void SparseTraceFilter::CubeLabelSoma( const VectorVec3d &curDendrite )
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();

    VectorVec5d rawSomaCurve;
    rawSomaCurve.clear();
    Vec5d tmpVec5d;
    for(size_t i = 0; i < curDendrite.size(); ++i){
        tmpVec5d << curDendrite[i](0), curDendrite[i](1), curDendrite[i](2), 0, 0;
        rawSomaCurve.push_back(tmpVec5d);
    }
    VectorVec3d tmpPtCurve;
    Vec3d tmpVec3d;
    for (size_t i = 0; i < rawSomaCurve.size(); ++i){
        tmpVec3d << rawSomaCurve[i](0), rawSomaCurve[i](1), rawSomaCurve[i](2);
        tmpPtCurve.push_back(tmpVec3d);
    }
    //
    if(!curDendrite.empty()){
        std::vector<double> R;//here is list
        std::vector<double> V;
        VectorVec5d resultCurveCopy;
        CalcParmOfCurveNodeList(origImg, backImg, tmpPtCurve, R, V );
        for (size_t i = 0; i < rawSomaCurve.size(); ++i){
            rawSomaCurve[i](3) = R[i];
            rawSomaCurve[i](4) = V[i];
        }
        Vec5d half;
        //
        for (size_t i = 0; i < rawSomaCurve.size() - 1; ++i){
            resultCurveCopy.push_back(rawSomaCurve[i]);
            half = 0.5 * (rawSomaCurve[i] + rawSomaCurve[i + 1]);
            resultCurveCopy.push_back(half);
        }
        resultCurveCopy.push_back(rawSomaCurve[rawSomaCurve.size() - 1]);

        std::vector<double> three;
        for (VectorVec5d::size_type i = 0; i < resultCurveCopy.size(); ++i){
            three.push_back((double)std::min(4, NGUtility::Round(resultCurveCopy[i](3) + 3.5) ) );//2016-3-13
        }
        int id1(0), id2(0), id3(0);
        int xMin, xMax, yMin, yMax, zMin, zMax;
        //VectorVec3i test;//test
        double threBack(0);
        for (int ik = 1; ik < 2*int(rawSomaCurve.size())-2; ++ik){
            id1 = NGUtility::Round(resultCurveCopy[ik](0) );
            id2 = NGUtility::Round(resultCurveCopy[ik](1) );
            id3 = NGUtility::Round(resultCurveCopy[ik](2) );

            xMin = std::max(0, id1 - std::max(0, (int)three[ik]));
            xMax = std::min(nxx - 1, id1 + std::max(0, (int)three[ik]));
            yMin = std::max(0, id2 - std::max(0, (int)three[ik]));
            yMax = std::min(nyy - 1, id2 + std::max(0, (int)three[ik]));
            zMin = std::max(0, id3 - std::max(0, (int)three[ik]));
            zMax = std::min(nzz - 1, id3 + std::max(0, (int)three[ik]));

            //double threBack(0);
            threBack = std::max(double(origImg(id1, id2, id3) + backImg(id1, id2, id3)) / 2, 0.85*double(origImg(id1, id2, id3)));
            threBack = std::max(threBack, backImg(id1, id2, id3) + paramPack->axonDiffuseValue_);
            for (int ii = xMin; ii <= xMax; ++ii){
                for (int jj = yMin; jj <= yMax; ++jj){
                    for (int kk = zMin; kk <= zMax; ++kk){
                        //2015-6-8
                        //threBack = double(backImg(ii, jj, kk)) + 4.0 * std::sqrt(double(backImg(ii, jj, kk)));
                        //threBack = std::min(threBack, 1.5 * double(backImg(ii, jj, kk)) );
                        //threBack = std::max(300.0, threBack);//2016-3-18
                        if (traceLabelMatrix_(ii, jj, kk) == 0 && double(origImg(ii, jj, kk)) > threBack){
                            traceLabelMatrix_(ii, jj, kk) = 1;
                            indexImg_(ii,jj,kk) = globalID_;//2015-6-8
                            //test.push_back(Vec3i(ii,jj,kk));
                        }
                    }
                }
            }//for
        }//for
    }
}

void SparseTraceFilter::BuildCurTreeLevel(const std::vector<std::vector<VectorVec3d> >& curTree)
{
    if (curTree.empty()) {
        printf("no curve!\n");
        return;
    }
    treeLevel_.clear();
    treeLevel_.resize(curTree.size());
    size_t level = 0;
    size_t id = 0;
    for (size_t i = 0; i < curTree.size();++i) {
         for (size_t j = 0; j < curTree[i].size(); ++j) {
            treeLevel_[i].push_back(id++);
         }
         ++level;
    }
}

bool SparseTraceFilter::IsNearBoundary(const Vec3d& arg, int x, int y, int z)
{
    return arg(0) <= paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(0) > double(x) - paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(1) <= paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(1) > double(y) - paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(2) <= 3.0
        || arg(2) > double(z - 3);
}


void SparseTraceFilter::MeanShiftGrowTraceTotal(Vec3d initPt, Vec3d initDir, const SVolume& origImg, const SVolume& backImg,
    VectorVec5d& resultCurve, int& collideIndex)
{
    if (IsNearBoundary(initPt, origImg.x(), origImg.y(), origImg.z())) return;
    Vec5d tmpPt; tmpPt << initPt(0), initPt(1), initPt(2), 0, 0;
    resultCurve.push_back(tmpPt);
    bool traceFlag = true;
    int i(0);
    //Vec3d tmpDir = initDir;
    collideIndex = 0;
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    Vec3d nextNode, nextDir, curPt, tmpDir; nextDir.setZero();
    VecXd tmpValueList, signalIntensity, backIntensity;
    MatXd curFeature, tmp;
    VectorVec3d tmpCurve;
    double internalDist;
    while (traceFlag) {
        //LOG(INFO)<<"MeanShiftGrowTrace";
        MeanShiftGrowTrace(initPt, initDir, origImg, backImg, nextNode, nextDir);
        if ((nextDir - Vec3d(0, 0, 0)).norm() - 0.0 > 0.001) {
            ++i;
            tmpPt.block(0, 0, 3, 1) = nextNode;
            resultCurve.push_back(tmpPt);
            tmpDir = initDir;
            initDir = nextDir;
            initPt = nextNode;
            if (initPt.minCoeff() < 1 || initPt(0) > nx - 2 || initPt(1) > ny - 2 || initPt(2) > nz - 2) {//
                traceFlag = false;
            }
            if (i > 1) {
                VectorVec3d tmpResultCurve;
                NGUtility::GetPartVectorVec3d(resultCurve, i - 2, i, tmpResultCurve);
                NGUtility::WeighRayValue(tmpResultCurve, origImg, signalIntensity);
                NGUtility::WeighRayValue(tmpResultCurve, backImg, backIntensity);

                tmpValueList = signalIntensity - backIntensity - backIntensity.cwiseSqrt() *  3.0;
                if (tmpValueList.maxCoeff() < 0) {
                    //printf("Signal weak\n");
                    int nxx = origImg.x();
                    int nyy = origImg.y();
                    int nzz = origImg.z();
                    int rdrvx = std::min(std::max(NGUtility::Round(initPt(0) - 1.0), 0), nxx - 1);
                    int rdrvy = std::min(std::max(NGUtility::Round(initPt(1) - 1.0), 0), nyy - 1);
                    int rdrvz = std::min(std::max(NGUtility::Round(initPt(2) - 1.0), 0), nzz - 1);//2015-6-8
                    double threv = backImg(rdrvx, rdrvy, rdrvz);
                    Vec3d nextCurveNode, tmpx1;
                    int isNextNodeValid;
                    //TraceNextCurveNodeV2(origImg, initPt, threv, nextDir, 0, nextCurveNode, tmpx1, isNextNodeValid);
                    SignalIdentityMeanshift(origImg, initPt, threv, nextDir, 0, isNextNodeValid);
                    if (isNextNodeValid != 0) 
                        traceFlag = true;
                    else{
                        traceFlag = 0;
                        i -= 3;
                    }
                }

            }
            if (i > 1) {
                internalDist = (resultCurve[i].block(0, 0, 3, 1) - resultCurve[i - 1].block(0, 0, 3, 1)).norm();
                if (internalDist < 0.02) {
                    printf("SVM 2\n");
                    traceFlag = false;
                }
            }

            if (i > maxSteps_ - 1) {
                traceFlag = false;
            }

            if (tmpDir.transpose() * initDir < -0.5) {
                printf("svm 3\n");
                traceFlag = false;
            }

            if (i > 1) {
                curPt = resultCurve[i].block(0, 0, 3, 1);
                if (traceFlag) {
                    DetectCollisionByMeanShift(curPt, indexImg_, origImg, traceFlag, collideIndex);
                }
            }
        }
        else{
            printf("svm 4\n");
            traceFlag = 0;
        }
    }
    {
        int n = i + 1;
        VectorVec5d tmpResultCurve(n);
        std::copy(resultCurve.begin(), resultCurve.begin() + n, tmpResultCurve.begin());
        resultCurve.swap(tmpResultCurve);
    }
}

void SparseTraceFilter::MeanShiftGrowTrace(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg,
    Vec3d& nextNode, Vec3d& nextDir)
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    VectorMatXd samplePlaneSet;
    std::vector<VectorVec4d> extractPtSet;
    VectorVec3d centerPoints, newPositions;
    std::vector<double> centerPtValue;
    ExtractSampleConePtSet(initPt, initDir, origImg, backImg, samplePlaneSet, extractPtSet, centerPoints);
    WeighRayValue(centerPoints, origImg, centerPtValue);
    bool flag(true);
    double centerPtMin(100000000.0), centerPtXMax(0.0), centerPtYMax(0.0), centerPtZMax(0.0);
    for (size_t i = 0; i < centerPoints.size(); ++i) {
        centerPtMin = std::min(centerPtMin, centerPoints[i].minCoeff());
        centerPtXMax = std::max(centerPtXMax, centerPoints[i](0));
        centerPtYMax = std::max(centerPtYMax, centerPoints[i](1));
        centerPtZMax = std::max(centerPtZMax, centerPoints[i](2));
    }

    if (centerPtMin < 1.0 || int(centerPtXMax) > nx - 2 || int(centerPtYMax) > ny - 2 || int(centerPtZMax) > nz - 2)
        flag = false;

    if (std::accumulate(centerPtValue.begin(), centerPtValue.end(), 0.0) / double(centerPtValue.size()) > 20.0 && flag) {
        CalcPositionByMeanshift(extractPtSet, centerPoints, newPositions);
        nextNode = newPositions[1];
        Vec3d nextDir1 = newPositions[1] - newPositions[0]; nextDir1.normalize();
        //NGUtility::Principald(newPositions, nextDir1);
        nextDir = 0.1 * initDir + 0.9 * nextDir1;
        nextDir /= nextDir.norm() + 0.0001;
    }
    else{
        nextDir.setZero();
        nextNode.setZero();
    }
}

void SparseTraceFilter::ExtractSampleConePtSet(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg,
    VectorMatXd& samplePlaneSet, std::vector<VectorVec4d> &extractPtSet, VectorVec3d& centerPoints)
{
    Vec3d x2, x3;
    CalcOrthoBasis(initDir, x2, x3);
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    int cpNum = 5;
    samplePlaneSet.resize(cpNum-1);
    extractPtSet.resize(cpNum-1);
    centerPoints.resize(cpNum);
    centerPoints[0] = initPt + initDir;
    Vec3d curPt;
    int rx, ry, rz;
    MatXd perpendicularPlane;
    bool collideFlag = false;
    for (int jj = 1; jj < cpNum; ++jj) {
        collideFlag = false;
        int r = jj * 2;
        perpendicularPlane.setZero(3 + 2 * r, 3 + 2 * r);
        Vec3d centerPt = initPt + (1 + 0.5*double(jj)) * initDir;
        VectorVec4d curExtractPtSet;
        centerPoints[jj] = centerPt;
        for (int i = -r - 1; i <= r + 1; ++i) {
            for (int j = -r - 1; j <= r + 1; ++j) {
                curPt = centerPt + double(i) * x2 + double(j) * x3;
                curPt(0) = std::max(std::min(curPt(0), double(nx - 1)), 0.0);
                curPt(1) = std::max(std::min(curPt(1), double(ny - 1)), 0.0);
                curPt(2) = std::max(std::min(curPt(2), double(nz - 1)), 0.0);
                rx = NGUtility::Round(curPt(0));
                ry = NGUtility::Round(curPt(1));
                rz = NGUtility::Round(curPt(2));
                if (indexImg_(rx, ry, rz) != 0){//region growing function fill the trace label matrix
                    collideFlag = true;
                    continue;
                }
                double signalIntensity = std::max(0.0, NGUtility::WeighRayValue(curPt, origImg) - NGUtility::WeighRayValue(curPt, backImg));
                curExtractPtSet.push_back(Vec4d(curPt(0), curPt(1), curPt(2), signalIntensity));
                perpendicularPlane(i + r + 1, j + r + 1) = signalIntensity;
            }
        }
        if (collideFlag) {
            int num = 0; double sum = 0.0;
            for (int i = 0; i < perpendicularPlane.rows(); ++i) {
                for (int j = 0; j < perpendicularPlane.cols(); ++j){
                    if (perpendicularPlane(i,j) < 0.1) 
                        continue;
                    ++num;
                    sum += perpendicularPlane(i, j);
                }
            }
            double mean = sum / double(num);
            for (int i = 0; i < perpendicularPlane.rows(); ++i) {
                for (int j = 0; j < perpendicularPlane.cols(); ++j){
                    if (perpendicularPlane(i, j) < 0.1)
                        perpendicularPlane(i, j) = mean;
                }
            }
        }
        extractPtSet[jj - 1].swap(curExtractPtSet);
        samplePlaneSet[jj - 1].swap(perpendicularPlane);
    }
}

//centerPoints is passed by value
void SparseTraceFilter::CalcPositionByMeanshift(const std::vector<VectorVec4d> &extractPtSet, VectorVec3d centerPoints,//centerPoints will be modified
    VectorVec3d& newPositions)
{
    size_t centerPointsNum = centerPoints.size();
    Vec3d pos0, pos1, pos2, pos3;
    Vec3d pos0bak, pos1bak, pos2bak, pos3bak;
    pos0bak.setZero();
    pos1bak.setZero();
    pos2bak.setZero();
    pos3bak.setZero();
    for (size_t kkj = 0u; kkj < 20lu; ++kkj) {
        pos0bak = pos0;
        pos1bak = pos1;
        pos2bak = pos2;
        pos3bak = pos3;

        CalcOnePosByMeanShift(extractPtSet[0], centerPoints[1], pos0);
        pos0 = 0.8 * pos0 + 0.1 * (centerPoints[0] + centerPoints[2]);

        CalcOnePosByMeanShift(extractPtSet[1], centerPoints[2], pos1);
        pos1 = 0.8 * pos1 + 0.1 * (centerPoints[1] + centerPoints[3]);

        CalcOnePosByMeanShift(extractPtSet[2], centerPoints[3], pos2);
        pos2 = 0.8 * pos2 + 0.1 * (centerPoints[2] + centerPoints[4]);

        CalcOnePosByMeanShift(extractPtSet[3], centerPoints[4], pos3);

        centerPoints[1] = pos0;
        centerPoints[2] = pos1;
        centerPoints[3] = pos2;
        centerPoints[4] = pos3;

        if (kkj > 10u){
            if ((pos0bak - pos0).norm() < 0.001 && (pos1bak - pos1).norm() < 0.001 &&
                (pos2bak - pos2).norm() < 0.001 && (pos3bak - pos3).norm() < 0.001) {
                break;
            }
        }
    }
    newPositions.resize(centerPointsNum - 1);
    std::copy(centerPoints.begin() + 1, centerPoints.begin() + centerPointsNum, newPositions.begin());
}

void SparseTraceFilter::CalcOnePosByMeanShift(const VectorVec4d& curExtractPtSet, const Vec3d& centerPoint, Vec3d& resultPos)
{
    size_t ptSetNum = curExtractPtSet.size();
    resultPos.setZero();
    double den(0.0);
    double curDist, curPixel, gaussianDist;
    double wet;
    for (size_t i = 0; i < ptSetNum; ++i) {
        curDist = (curExtractPtSet[i].block(0, 0, 3, 1) - centerPoint).norm();
        gaussianDist = std::exp(-std::pow(curDist, 2.0));
        curPixel = curExtractPtSet[i](3) / 100.0;
        wet = curPixel * gaussianDist;
        resultPos += wet * curExtractPtSet[i].block(0, 0, 3, 1);
        den += wet;
    }
    if (den > 0.1) {
        resultPos = resultPos / den;
    }
    else {
        resultPos = centerPoint;
    }
}

void SparseTraceFilter::DetectCollisionByMeanShift(const Vec3d &curPt, const IVolume& indexImg, const SVolume& origImg, 
    bool& traceFlag, int& collideIndex)
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    collideIndex = 0;
    traceFlag = true;
    int id1 = NGUtility::Round(curPt(0));
    int id2 = NGUtility::Round(curPt(1));
    int id3 = NGUtility::Round(curPt(2));
    int r = 1;
    int xMin = std::max(0, id1 - r);
    int xMax = std::min(nx - 1, id1+r);
    int yMin = std::max(0, id2-r);
    int yMax = std::min(ny - 1, id2+r);
    int zMin = std::max(0, id3-r);
    int zMax = std::min(nz - 1, id3+r);
    VectorVec5d collideNodeSet;
    Vec5d tmp5d;
    int curIndex;
    double dist;
    for (int ii = xMin; ii <= xMax; ++ii) {
        for (int jj = yMin; jj <= yMax; ++jj) {
            for (int ij = zMin; ij <= zMax; ++ij) {
                curIndex = indexImg(ii, jj, ij);
                if (curIndex != 0 ) {
                    tmp5d << double(ii), double(jj), double(ij), 0.0, 0.0;
                    dist = (tmp5d.block(0, 0, 3, 1) - curPt).norm();
                    tmp5d(3) = dist;
                    tmp5d(4) = curIndex;
                    collideNodeSet.push_back(tmp5d);
                }
            }
        }
    }
    if (collideNodeSet.size() > 0) {
        double threvDist = 0;
        double tmpIndex = 0;
        /*if (collideNodeSet.size() > 3) {
            std::partial_sort(collideNodeSet.begin(), collideNodeSet.begin() + 4, collideNodeSet.end(), Vec5d_3th_less());
            threvDist = collideNodeSet[3](3);
            tmpIndex = collideNodeSet[3](4);
            }
            else{
            VectorVec5d::iterator it = std::min_element(collideNodeSet.begin(), collideNodeSet.end(), Vec5d_3th_less());
            threvDist = (*it)(3);
            tmpIndex = (*it)(4);
            }*/
        VectorVec5d::iterator it = std::min_element(collideNodeSet.begin(), collideNodeSet.end(), Vec5d_3th_less());
        threvDist = (*it)(3);
        tmpIndex = (*it)(4);
        if (threvDist < 11.5) {
            //printf("svm 5\n");
            traceFlag = false;
            collideIndex = tmpIndex;//int((*it)(4));
        }
    }
}

void SparseTraceFilter::CalcParmOfCurveNodeListSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
    const VectorVec3d &curNode,
    std::vector<double> &radius, std::vector<double> &rav)
{
    radius.clear();
    rav.clear();

    VectorVec3d::size_type nxx = curNode.size();
    double tmpR(0.0);
    double tmpV(0.0);
    for (VectorVec3d::size_type i = 0; i < nxx; ++i){
        Vec3d tmpdata1;
        tmpdata1(0) = NGUtility::Round(curNode[i](0));//int(curNode[i](0) + 0.5);
        tmpdata1(1) = NGUtility::Round(curNode[i](1));
        tmpdata1(2) = NGUtility::Round(curNode[i](2));
        //LOG(INFO)<<"CalcParmOfCurveNodeSVM";
        CalcParmOfCurveNodeSVM(origImg, backImg, tmpdata1, tmpR, tmpV);
        radius.push_back(tmpR);
        rav.push_back(tmpV);
    }
}

void SparseTraceFilter::CalcParmOfCurveNodeSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
    const Vec3d &curNode, double &radius, double &wet)
{
    int vx = backImg.x();
    int vy = backImg.y();
    int vz = backImg.z();

    int xMin, xMax, yMin, yMax, zMin, zMax;
    //2016-4-6
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, NGUtility::Round(curNode(0)), NGUtility::Round(curNode(1)), NGUtility::Round(curNode(2)),
        5, 5, 5, 0, vx - 1, 0, vy - 1, 0, vz - 1);

    int xLen = xMax - xMin + 1;
    int yLen = yMax - yMin + 1;
    int zLen = zMax - zMin + 1;

    Vec3d center(curNode(0) - xMin, curNode(1) - yMin, curNode(2) - zMin);//SL
    Volume<double> procImg;
    procImg.SetSize(xLen, yLen, zLen);//MM1

    for (int i = xMin; i <= xMax; ++i){
        for (int j = yMin; j <= yMax; ++j){
            for (int ij = zMin; ij <= zMax; ++ij){
                procImg(i - xMin, j - yMin, ij - zMin) = double(origImg(i, j, ij)) - double(backImg(i, j, ij))
                    - 3.0 * std::sqrt(double(backImg(i, j, ij)));
            }
        }
    }//for

    radius = 0.05;//r
    wet = 0.0;//v
    MatXd distWet(2, 4); distWet.setZero();
    double vv0 = 0.0;
    double vv1 = 0.0;

    //2015-8-13
    Vec4d level; level << 1.0, 2.0, 3.0, 4.0;

    for (double i = 0; i < xLen; ++i){
        for (double j = 0; j < yLen; ++j){
            for (double ij = 0; ij < zLen; ++ij){
                Vec3d dist(i - center(0), j - center(1), ij - center(2));
                double distNorm = dist.norm();
                if (distNorm <= level(0)){ //距离中心为1
                    distWet(0, 0) += 1.0;
                    if (procImg(i, j, ij) > 0.0){
                        vv0 += procImg(i, j, ij);
                        distWet(1, 0) += 1.0;
                    }
                }//if

                if (distNorm <= level(1)){ //距离中心为2
                    distWet(0, 1) += 1.0;
                    if (procImg(i, j, ij) > 0){
                        vv1 += procImg(i, j, ij);
                        distWet(1, 1) += 1.0;
                    }
                }

                if (distNorm <= level(2)){ //距离中心为3
                    distWet(0, 2) += 1.0;
                    if (procImg(i, j, ij) > 0)
                        distWet(1, 2) += 1.0;
                }

                if (distNorm <= level(3)){ //距离中心为4
                    distWet(0, 3) += 1.0;
                    if (procImg(i, j, ij) > 0)
                        distWet(1, 3) += 1.0;
                }
            }//for
        }
    }//for

    //2015-8-13
    Vec4d procDistWet = distWet.row(1).array() / (distWet.row(0).array() + 0.0001);
    for (int i = 1; i <= 4; ++i){
        if (procDistWet(4 - i) > 0.5){
            radius = 5.0 - (double)i;
            break;
        }
    }

    if (radius > 2.1){//warning!!
        radius = radius / std::pow(distWet(0, int(radius + 0.5) - 1) / distWet(1, int(radius + 0.5) - 1), 1.0 / 3.0);
        wet = vv1 / distWet(1, 1);
    }

    if (radius < 2.1){
        if (distWet(1, 1) > 0.0){
            radius = 2.0 / std::pow(distWet(0, 1) / distWet(1, 1), 1.0 / 3.0);
            wet = vv1 / distWet(1, 1);
        }
        if (std::abs(distWet(1, 1) - 0.0) < EXPO && distWet(1, 0) > 0.0){
            radius = 1.0 / std::pow(distWet(0, 0) / distWet(1, 0), 1.0 / 3.0);
            wet = vv0 / distWet(1, 0);
        }
    }
}

void SparseTraceFilter::CubeLabelL1Opt(const VectorVec3d &curDendrite)
{
    //initial
    //const SVolume &origImg = *origImgPointer;
    //int nxx = origImg.x();
    //int nyy = origImg.y();
    //int nzz = origImg.z();
    //
    std::vector<Line3d> tmp;
    tmp.push_back(curDendrite);
    auto &tmp1 = tmp[0];
    for (auto &it : tmp1) {
        it(0) += paramPack->xMin_;
        it(1) += paramPack->yMin_;
        it(2) += paramPack->zMin_;
    }
    NGCaliberBuilder caliberBuilder = CaliberBuilder::New();
    caliberBuilder->SetParam(paramPack);
    caliberBuilder->SetLinePoints(tmp);
    auto res = caliberBuilder->Update();
    if (!res->success()){
        return;
    }
    auto resContourPt = caliberBuilder->GetCaliberPointSet();
    for (auto &it : *resContourPt) {
        for (auto &iter :it) {
            iter(0) -= paramPack->xMin_;
            iter(1) -= paramPack->yMin_;
            iter(2) -= paramPack->zMin_;
        }
    }
    /*NGUtility::WriteVectorVec3d("F:/nima.swc", "dot", resContourPt[0]);
    NGImageWriter writer = ImageWriter::New();
    writer->SetInput(paramPack->OrigImage);
    writer->SetOutputFileName("F:/nima.tif");
    writer->Update();
    system("pause");*/
    int rx, ry, rz;
    for (auto &it : *resContourPt) {
        for (auto &iter : it) {
            rx = NGUtility::Round(iter(0));
            ry = NGUtility::Round(iter(1));
            rz = NGUtility::Round(iter(2));
            if (traceLabelMatrix_(rx, ry, rz) == 0 || indexImg_(rx,ry,rz) == 0){//because the extract cube label do not fill the indeximg_
                traceLabelMatrix_(rx, ry, rz) = 1;
                indexImg_(rx, ry, rz) = globalID_;//2015-6-8
                //test.push_back(Vec3i(ii,jj,kk));
            }
        }
    }
}

void SparseTraceFilter::SignalIdentityMeanshift(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode, const double &threv, const Vec3d &initDir, bool flag, int &isNextNodeValid)
{
    //nextCurveNode.setZero();
    isNextNodeValid = 1;
    Vec3d nextDenDir = initDir;

    VectorVec3d neighborPtSet;//
    std::vector<double> neighborWet;//W

    Vec3d x10;//x10
    x10.setZero();
    Vec3d x11;//x11
    x11.setZero();
    //2015-6-8 add flag
    CalcNeighborSignalV2(locOrigImg, curSeedNode, nextDenDir, threv, flag, neighborPtSet, neighborWet, x10, x11);

    size_t neighborWetNum = neighborPtSet.size();
    /*W1 = sort(W, 'descend)*/
    std::vector<double> sortedNeighborWet;//W1
    sortedNeighborWet.assign(neighborWet.begin(), neighborWet.end());
    std::sort(sortedNeighborWet.begin(), sortedNeighborWet.end(), std::greater<double>());

    size_t dsw = (std::min)((size_t)10, sortedNeighborWet.size());//2015-7-7
    //size_t dsw = (std::min)((size_t)11, sortedNeighborWet.size());//2015-6-15 TDI079
    /**/
    double thrdk(0.0);
    if (dsw > 5){//2015-7-7
        //if(dsw > 10){//2015-6-15 TDI079
        thrdk = 0.0;
        for (int i = 0; i < int(dsw); ++i){
            thrdk += sortedNeighborWet[i];
        }
        thrdk /= double(dsw);
    }
    else thrdk = 0.0;

    Mat3d sigmaH;//sigmaH
    sigmaH.setZero();
    Vec3d kk;//kk
    kk.setZero();

    double P = 0.0;
    //2015-6-8
    if (neighborWetNum > std::max(0.1 * threv, 10.0) &&
        thrdk > std::max(std::min(std::min(0.6 * threv, 4.5 * std::sqrt(threv)), paramPack->diffuseValue_), paramPack->traceValue_)){
        /*if (neighborWetNum > std::max(0.1 * threv, 10.0) &&
        thrdk > std::max(std::min(std::min(60 * threv, 20.0 * std::sqrt(threv)), paramPack->diffuseValue_), paramPack->traceValue_)){*/

        //2015-6-8 thrdk>max(min([0.6*threv,4.5*sqrt(threv),threValue]),1)
        CalcConstraintPCA(neighborPtSet, curSeedNode, 0.150 * Eigen::Matrix3d::Identity(), neighborWet, P, sigmaH, kk);
        if (P>0){
            isNextNodeValid = 1;
        }
        else{
            isNextNodeValid = 0;
        }
    }
    else isNextNodeValid = 0;
}
