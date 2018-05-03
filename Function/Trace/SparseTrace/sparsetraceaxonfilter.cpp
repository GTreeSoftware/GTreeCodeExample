/* SparseAxonTracer
*	function: trace single neuron
*	Author : zhouhang, lishiwei
*	2015-10-28
*/
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <deque>
#include <iterator>
#include <utility>
#include <Eigen/Core>
#include <Eigen/LU>
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
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
using namespace google;
#include "sparsetraceaxonfilter.h"
#include "../../../ngtypes/tree.h"
#include "../../../ngtypes/soma.h"
#include "../../volumealgo.h"
#include "../../contourutil.h"
#include "../../Function/IO/imagewriter.h"
#include "CrossBranchFilter.h"
#include "CrossFiberFilter.h"
#include "../../NGUtility.h"

const double EXPO = 0.000001;

SparseTraceAxonFilter::SparseTraceAxonFilter()
{
    className_ = std::string("SparseTraceAxonFilter");
    //diffuseValue_ = 6;//2015-6-15 TDI079
    globalID_ = 1;
    //boundaryDistanceThreshold_ = 5;
    initialDirection_ << 0,0,0;
    m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));
}

SparseTraceAxonFilter::~SparseTraceAxonFilter()
{

}

ProcStatPointer SparseTraceAxonFilter::Update()
{
    if(!m_Input || !m_Back || !m_Bin ){
        printf("error occurred in %s\n", className_.c_str());
        LOG(ERROR)<<"error occurred in "<<className_;
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
    //2016-4-23
    if(!m_Source){
        m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));
    }
    //2015-6-8
    origImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Input);//2015-6-8
    backImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Back);//2015-6-8
    binImgPointer =  std::dynamic_pointer_cast<const Volume<NGCHAR> >(m_Bin);
    std::shared_ptr<const Soma > tmpSoma =  std::dynamic_pointer_cast<const Soma>(m_Soma);
    std::shared_ptr<TreeCurve > tmpTree =  std::dynamic_pointer_cast<TreeCurve>(m_Source);

    if(!origImgPointer || !backImgPointer || !binImgPointer || !tmpSoma) {
        printf("error back noise image!\n");
        //LOG(ERROR)<<"error back noise image!";
        printf("error occurred in %s\n", className_.c_str());
        //LOG(ERROR)<<"error occurred in "<<className_;
        MAKEPROCESSSTATUS(resSta, false, className_, "input data is wrong.");
        return resSta;
    }

    const Volume<unsigned short> &origImg = *origImgPointer;
    //initial
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();
    //2015-6-8
    indexImg_.SetSize(nxx, nyy, nzz);//global variant
    indexImg_.SetZero();
    traceLabelMatrix_.SetSize(nxx, nyy, nzz);//global variant
    traceLabelMatrix_.SetZero();
    //Vec3d soma_;
    soma_ << tmpSoma->GetCell(0).x, tmpSoma->GetCell(0).y, tmpSoma->GetCell(0).z;
    if (paramPack->activeTree) {
        NG_CREATE_DYNAMIC_CAST(NeuronTree, allTree, paramPack->activeTree);
        if (allTree) {
            //if (paramPack->isLocalTrace_) {
            // fill current image to avoid trace repeatedly. Because the image is not scaled to low anisotropy, so the soma contour reconstruction is invalid.
            //Vec3i curPt;
            //int rx, ry, rz;
            //if (allTree) {//has previous tree
            //    auto &allCurve = allTree->m_curveList;
            //    for (auto &curve : allCurve) {
            //        for (auto &pt : curve) {
            //            rx = NGUtility::Round(pt(0));
            //            ry = NGUtility::Round(pt(1));
            //            rz = NGUtility::Round(pt(2));
            //            if (rx < paramPack->xMin_ || rx > paramPack->xMax_ || ry < paramPack->yMin_ || ry > paramPack->yMax_
            //                || rz < paramPack->zMin_ || rz > paramPack->zMax_) {
            //                continue;
            //            }
            //            curPt << std::max(std::min(NGUtility::Round(pt(0)) - paramPack->xMin_, nxx - 1), 0),
            //                std::max(std::min(NGUtility::Round(pt(1)) - paramPack->yMin_, nyy - 1), 0),
            //                std::max(std::min(NGUtility::Round(pt(2)) - paramPack->zMin_, nzz - 1), 0);
            //            traceLabelMatrix_(curPt(0), curPt(1), curPt(2)) = 1;
            //            indexImg_(curPt(0), curPt(1), curPt(2)) = std::numeric_limits<int>::max();
            //        }
            //    }
            //}
            //}
            //else{
            Vec3i curPt;
            int rx, ry, rz;
            //if (allTree) {//has previous tree
                auto &allCurve = allTree->m_curveList;
                int r, xMin, xMax, yMin, yMax, zMin, zMax;
                size_t sz = 0;
                for (auto &curve : allCurve) {
                    sz = curve.size();
                    for (auto k = 0; k < sz; k += 2llu) {//not all pt need to fill
                        auto &pt = curve[k];
                        //for (auto &pt : curve) {
                        rx = NGUtility::Round(pt(0));
                        ry = NGUtility::Round(pt(1));
                        rz = NGUtility::Round(pt(2));
                        if (rx < paramPack->xMin_ || rx > paramPack->xMax_ || ry < paramPack->yMin_ || ry > paramPack->yMax_
                            || rz < paramPack->zMin_ || rz > paramPack->zMax_) {
                            continue;
                        }
                        if (std::abs(soma_(0) - rx) < 4.0 && std::abs(soma_(1) - ry) < 4.0 && std::abs(soma_(2) - rz) < 4.0) {
                            continue;
                        }
                        curPt << std::max(std::min(NGUtility::Round(pt(0)) - paramPack->xMin_, nxx - 1), 0),
                            std::max(std::min(NGUtility::Round(pt(1)) - paramPack->yMin_, nyy - 1), 0),
                            std::max(std::min(NGUtility::Round(pt(2)) - paramPack->zMin_, nzz - 1), 0);
                        r = std::min(4, NGUtility::Round(NGUtility::Round(pt(3)) + 1.5));
                        xMin = std::max(0, curPt(0) - r);
                        xMax = std::min(nxx - 1, curPt(0) + r);
                        yMin = std::max(0, curPt(1) - r);
                        yMax = std::min(nyy - 1, curPt(1) + r);
                        zMin = std::max(0, curPt(2) - r);
                        zMax = std::min(nzz - 1, curPt(2) + r);
                        for (int ii = xMin; ii <= xMax; ++ii){
                            for (int jj = yMin; jj <= yMax; ++jj){
                                for (int kk = zMin; kk <= zMax; ++kk){
                                    traceLabelMatrix_(ii, jj, kk) = 1;
                                    indexImg_(ii, jj, kk) = std::numeric_limits<int>::max();
                                }
                            }
                        }
                    }
                }
            //}
        }
    }

    globalID_ = 1;
    std::vector<VectorVec3d> initialDendrites;
    std::vector<std::vector<VectorVec3d> > treeDataset;
    TraceDendritesConnectedAxon(soma_, initialDendrites );
    treeDataset.push_back(initialDendrites);

    size_t traceDeep = paramPack->enableSVM_? 0:2;//2016-8-2 many many bugs in RecursiveTrace
    for(size_t i = 0; i < traceDeep; ++i){
        if(!treeDataset.back().empty()){
            RecursiveTrace(treeDataset);
        }else  break;
    }
    //clear empty item
    if (treeDataset.back().empty()) 
        treeDataset.pop_back();

    //2016-4-19-------------------------//
    //get rawDendList
    Vec5d tmpVec5d; tmpVec5d.setZero();
    VectorVec5d tmpVectorVec5d;
    std::vector<VectorVec5d> rawDendList;
    for (size_t i = 0; i < treeDataset.size(); ++i) {
        for (size_t j = 0; j < treeDataset[i].size(); ++j) {
            tmpVectorVec5d.clear();
            for (size_t ij = 0; ij < treeDataset[i][j].size(); ++ij) {
                tmpVec5d.block(0,0,3,1) = treeDataset[i][j][ij];
                tmpVectorVec5d.push_back(tmpVec5d);
            }
            if (tmpVectorVec5d.empty()) continue;
            rawDendList.push_back(VectorVec5d());
            rawDendList.back().swap(tmpVectorVec5d);
        }
    }
    //get rawDendConInfo
    VectorMat2i rawDendConInfo;
    Mat2i tmpRawDendConInfo; tmpRawDendConInfo.setZero();
    for (size_t i = 0; i < rawDendList.size(); ++i) {
        const VectorVec5d& curCurve = rawDendList[i];
        for (size_t j = 0; j < rawDendList.size(); ++j) {
            if (j != i) {
                const VectorVec5d conCurve = rawDendList[j];
                for (size_t ij = 0; ij < conCurve.size(); ++ij) {
                    if ( (curCurve[0].block(0,0,3,1) - conCurve[ij].block(0,0,3,1)).norm() < 0.1 ) {
                        tmpRawDendConInfo(0,0) = int(j) + 1;
                    }
                    if ( (curCurve.back().block(0,0,3,1) - conCurve[ij].block(0,0,3,1)).norm() < 0.1 ) {
                        tmpRawDendConInfo(0,1) = int(j) + 1;
                    }
                }
            }
        }
        rawDendConInfo.push_back(tmpRawDendConInfo);
    }
    size_t num = 0;
    for (size_t i = 0; i < rawDendList.size(); ++i) {
        num += rawDendList[i].size();
    }
    printf("before SVM, there are %u points.\n", num);
    //SVM
    if(!rawDendList.empty() && paramPack->enableSVM_){
        if (!svmFilter) svmFilter = SVMTraceFilter::New(indexImg_, traceLabelMatrix_);
        svmFilter->SetInput(m_Input);
        svmFilter->SetParam(paramPack);
        svmFilter->SetInputBack(m_Back);
        svmFilter->SetInputBin(m_Bin);
        svmFilter->SetParam(paramPack);
        svmFilter->SetInputRawCurves(rawDendList);
        svmFilter->SetInputRawConInfo(rawDendConInfo);
        svmFilter->SetTraceInitPoint(false);
        svmFilter->SetSoma(soma_);
        svmFilter->SetInitialDirection(initialDirection_);
        svmFilter->SetCrossFiberFilter(cff);
        auto sres = svmFilter->Update();
        if (!sres->success()) 
            return sres;
        svmFilter->SwapRawCurves(rawDendList);
    }
    else if ( paramPack->enableSVM_){
        if (!svmFilter) svmFilter = SVMTraceFilter::New(indexImg_, traceLabelMatrix_);
        rawDendList.push_back(VectorVec5d());
        rawDendList.back().push_back(NGUtility::MakeVec5d(soma_));
        rawDendList.back().push_back(NGUtility::MakeVec5d(soma_ + initialDirection_ * 2.0));
        Mat2i tmp; tmp.setZero();
        rawDendConInfo.push_back(tmp);
        svmFilter->SetInput(m_Input);
        svmFilter->SetInputBack(m_Back);
        svmFilter->SetInputBin(m_Bin);
        svmFilter->SetParam(paramPack);
        svmFilter->SetInputRawCurves(rawDendList);
        svmFilter->SetInputRawConInfo(rawDendConInfo);
        svmFilter->SetTraceInitPoint(true);
        svmFilter->SetSoma(rawDendList.back().back().block(0,0,3,1));// !!!
        svmFilter->SetInitialDirection(initialDirection_);// !!!
        svmFilter->SetCrossFiberFilter(cff);
        auto sres = svmFilter->Update();
        if (!sres->success())
            return sres;
        svmFilter->SwapRawCurves(rawDendList);
    }
    rawDendList.erase(std::remove_if(rawDendList.begin(), rawDendList.end(), [](const Line5d& arg){return arg.size() < 3llu; }), rawDendList.end());
    //cross branch
    {
        NGCrossFiberFilter crossFilter = CrossFiberFilter::New();
        crossFilter->SetParam(paramPack);
        crossFilter->SetSoma(soma_);
        crossFilter->SetTraceLine(rawDendList);
        auto cres = crossFilter->Update();
        if (!cres->success())
            return cres;
    }
    
    tmpTree->SetCurve(rawDendList);//rawDendInfo is useless
    //destroy
    indexImg_.SetSize(0,0,0);
    traceLabelMatrix_.SetSize(0,0,0);
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}
  
ConstIDataPointer SparseTraceAxonFilter::GetOutput()
{
    if (!m_Source)   NG_SMART_POINTER_NEW(TreeCurve, m_Source, this);
    return m_Source;
}

IDataPointer SparseTraceAxonFilter::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData(m_Source);
    m_Source.reset();
    return tData;
}


//2015-6-8 make use of c++ class function, origImg, backImg, indexImg, traceLabelMatrix are all global varient
void SparseTraceAxonFilter::TraceDendritesConnectedAxon(const Vec3d& soma, std::vector<VectorVec3d> &initialDendrites)
{
    std::vector<VectorVec3d> tmpDendrites(1);
    tmpDendrites[0].push_back(soma);
    VectorVec3d directionSet;
    directionSet.push_back(initialDirection_);
    ContrainedPCATraceSomaBifur(tmpDendrites, directionSet, initialDendrites);
}

void SparseTraceAxonFilter::TraceNextCurveNodeV2(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode,
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
    /*if (neighborWetNum > std::max(0.1 * threv, 10.0) &&
        thrdk > std::max(std::min(std::min(0.6 * threv, 2.5 * std::sqrt(threv)), paramPack->axonDiffuseValue_*0.8), paramPack->axonTraceValue_)){//2017-3-27*/
    if (neighborWetNum > std::max(0.1 * threv, 20.0) &&
        thrdk > std::max(std::min(std::min(0.6 * threv, 2.5 * std::sqrt(threv)), paramPack->axonDiffuseValue_*0.8), paramPack->axonTraceValue_)){//2017-3-27
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

void SparseTraceAxonFilter::ExtractCubeBoundary( const std::vector<VectorVec3d>& curveCluster, CVolume& traceLabelMatrix,
                                                std::vector<VectorVec3i>& growPointSet,VectorVec3i& boundaryLabel )
{
    //initial
    //std::cout << "ExtractCubeBoundary bifurcation value_: " << paramPack->bifurcationValue_ << std::endl;
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
            pointRadius = std::max(2.0, double(NGUtility::Round(curvePointsRadius[i2] + 3.5)) + 1.0);//2016-3-18
            Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, int(curPoint[0]), int(curPoint[1]), int(curPoint[2]),
                int(pointRadius), int(pointRadius), int(pointRadius), 0, nxx, 0, nyy, 0, nzz);

            for(int i_x  = xMin; i_x <= xMax; ++i_x){
                for(int i_y = yMin; i_y <= yMax; ++i_y){
                    for(int i_z = zMin; i_z <= zMax; ++i_z){
                        //2015-6-8
                        if(origImg(i_x, i_y, i_z) - backImg(i_x, i_y, i_z) > paramPack->bifurcationValue_ * std::sqrt(double(backImg(i_x, i_y, i_z))))
                            curFlag = true;
                        else
                            curFlag = false;
						bool curFlag1 = origImg(i_x, i_y, i_z) > paramPack->axonBinaryThreshold_*double(backImg(i_x, i_y, i_z));
                        bool curFlag2 = origImg(i_x, i_y, i_z) > double(backImg(i_x, i_y, i_z)) + paramPack->axonDiffuseValue_;
                        if (paramPack->strongSignalMode_) curFlag = curFlag && curFlag1 && curFlag2;
                        else curFlag = curFlag || curFlag1 || curFlag2;
                        if(curFlag && 0 == traceLabelMatrix(i_x, i_y, i_z)){
                            traceLabelMatrix(i_x, i_y, i_z) = 1;
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
        RegionInflationModifyV2(boundaryLabelCp, traceLabelMatrix, curPoints);//2016-3-18
        if(!curPoints.empty()){
            growPointSet[i] = curPoints;
            boundaryLabelCp.swap(curPoints);
        } else{
            break;
        }
    }
}

void SparseTraceAxonFilter::CubeLabel( const VectorVec3d &curDendrite )
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
            //three.push_back((double)std::min(4, NGUtility::Round(resultCurveCopy[i](3) + 2.5) ) );//2016-3-13
            three.push_back((double)std::max(3, NGUtility::Round(resultCurveCopy[i](3) + 0.5) ) );//2016-3-18
        }
        int id1(0), id2(0), id3(0);
        int xMin, xMax, yMin, yMax, zMin, zMax;
        double threBack(0);
        //VectorVec3i test;//test
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

            for (int ii = xMin; ii <= xMax; ++ii){
                for (int jj = yMin; jj <= yMax; ++jj){
                    for (int kk = zMin; kk <= zMax; ++kk){
                        threBack = std::max(double(origImg(id1, id2, id3) + backImg(id1, id2, id3)) / 2, 0.85*double(origImg(id1, id2, id3)));
                        threBack = std::max(threBack, backImg(id1, id2, id3) + paramPack->axonDiffuseValue_);//TODO: 2018-04-26
                        //2015-6-8
                        //threBack = double(backImg(ii, jj, kk)) + 4.0 * std::sqrt(double(backImg(ii, jj, kk)));//2016-3-18
                        //threBack = std::min(threBack, 1.5 * double(backImg(ii, jj, kk)) );//2016-3-18
                        //threBack = std::max(300.0, threBack);//2016-3-18
                        if (traceLabelMatrix_(ii, jj, kk) == 0 && (double(origImg(ii, jj, kk)) > threBack)){
                            // || double(origImg(ii, jj, kk)) > double(backImg(ii,jj,kk))+paramPack->axonDiffuseValue_
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

void SparseTraceAxonFilter::RecursiveTrace( std::vector<std::vector<VectorVec3d> > &treeDataSet )
{
    //initialize
    const std::vector<VectorVec3d> &tmpDendrites = treeDataSet.back();
    size_t tmpDendritesNum = tmpDendrites.size();
    std::vector<VectorVec3d> sortClustre;
    sortClustre.reserve(tmpDendritesNum);
    for(size_t i = 0; i < tmpDendritesNum; ++i){
        if(!tmpDendrites[i].empty())
            sortClustre.push_back(tmpDendrites[i]);//cannot swap , because tmpDendrites is constant.
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
    if (bifurcationDirections.size() > 50) {
        NG_ERROR_MESSAGE("wrong bifurcation value.");
        bifurcationDirections.clear();
        bifurcationClustre.clear();
        boundaryLabel.clear();
    }
    //remove reverse trace point and direction
    for (auto k = 0; k < bifurcationDirections.size(); ++k) {
        if ((bifurcationClustre[k][0] - soma_ ).norm() < 8.0 ) {
            bifurcationClustre[k].clear();
            bifurcationDirections[k].setZero();
        }
        if ((bifurcationClustre[k][0] - soma_).norm() < 20.0 && bifurcationDirections[k].dot(initialDirection_) < -0.75) {
            bifurcationClustre[k].clear();
            bifurcationDirections[k].setZero();
        }
    }
    bifurcationClustre.erase(std::remove_if(bifurcationClustre.begin(), bifurcationClustre.end(), 
        [](const VectorVec3d &arg){return arg.empty(); }), bifurcationClustre.end());
    bifurcationDirections.erase(std::remove_if(bifurcationDirections.begin(), bifurcationDirections.end(),
        [](const Vec3d& arg){return arg.norm() < 0.1; }), bifurcationDirections.end());
    
    /*
    * filter bi fur
    */
    {
        //Vec3d soma;
        //std::tr1::shared_ptr<const Soma > tmpSoma =  std::tr1::dynamic_pointer_cast<const Soma>(m_Soma);
        //soma << tmpSoma->GetCell(0).x, tmpSoma->GetCell(0).y, tmpSoma->GetCell(0).z;
        /*
        * remove cross curves
        */
        //if(!cff) cff = CrossBranchFilter::New();
        //cff->SetParam(paramPack);
        //cff->SetData(tmpDendrites, soma, initialDirection_, bifurcationClustre, bifurcationDirections, boundaryLabel);
        //cff->Update();
    }
    
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

void SparseTraceAxonFilter::Train( std::vector<VectorVec5d>& arg)
{
    if (!svmFilter) {
        return;
    }
    svmFilter->Train(arg);
}

void SparseTraceAxonFilter::DetectBifurcation( CVolume &traceLabelMatrix, const std::vector<VectorVec3d> &curveCluster, VectorVec3d &bifurcationDirections, std::vector<VectorVec3d> &bifurcationClustre, VectorVec3i& boundaryLabel )
{
    //traceLabelMatrix is traceLabelMatrixCp;
    std::vector<VectorVec3i> growPointSet;
    ExtractCubeBoundary(curveCluster, traceLabelMatrix, growPointSet, boundaryLabel);
    {
        /*NGImageWriter writer = ImageWriter::New();
        writer->SetInput(paramPack->OrigImage);
        writer->SetOutputFileName("F:/wocaonima.tif");
        writer->Update();
        VectorVec3i tmp;
        for (auto &it : growPointSet) {
        std::copy(it.begin(), it.end(), std::back_inserter(tmp));
        }
        NGUtility::WriteVectorVec3i("F:/wocao.swc", "dot", tmp);
        system("pause");*/
    }

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
    {
        /*VectorVec3d tmp;trace
        for (auto &it : bifurcationClustre) {
        std::copy(it.begin(), it.end(), std::back_inserter(tmp));
        }
        NGUtility::WriteVectorVec3d("E:/Didem_BigNeuron/wocao.swc", "dot", tmp);
        system("pause"); */
    }
}

void SparseTraceAxonFilter::RegionInflationModifyV2( const VectorVec3i &curPoints, CVolume &growLabelMatrix, VectorVec3i &growPoints )
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
    bool finalFlag;
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
                    flag1 = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + paramPack->bifurcationValue_ * std::sqrt(double(backImg(ii,jj,ij)));
					flag2 = double(origImg(ii, jj, ij)) > paramPack->axonBinaryThreshold_*double(backImg(ii, jj, ij));
                    flag3 = double(origImg(ii, jj, ij)) > double(backImg(ii, jj, ij)) + paramPack->axonDiffuseValue_;
                    if(paramPack->strongSignalMode_) finalFlag = flag1 && flag2 && flag3;
                    else finalFlag = flag1 || flag2 || flag3;
                    if(0 == growLabelMatrix(ii,jj,ij) && finalFlag)
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



/*
* C++11
*/

