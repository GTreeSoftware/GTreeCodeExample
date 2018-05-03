#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <deque>
#include <utility>
#include <Eigen/LU>
#include <iostream>
#include <fstream>
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
#include "ngtypes/tree.h"
#include "ngtypes/soma.h"
#include "ngtypes/volume.h"
#include "Function/volumealgo.h"
#include "Function/contourutil.h"
#include "Function/Trace/traceutil.h"
#include "../../Function/NGUtility.h"
#include "../../Function/binaryfilter.h"
const double EXPO = 0.000001;
#ifdef _WIN32
#include <ctime>
#define M_PI       3.14159265358979323846
#else
#include <sys/time.h>
#endif
#include "../../ngtypes/tree.h"
#include "../../ngtypes/soma.h"
#include "../../ngtypes/volume.h"
#include "../../volumealgo.h"
#include "../../contourutil.h"
#include "../traceutil.h"
#include "SVMTraceFilter.h"
#include "CrossBranchFilter.h"
#include "Function/IO/imagewriter.h"

SVMTraceFilter::SVMTraceFilter(Volume<int>& p, Volume<unsigned char>& t):indexImg_(p), traceLabelMatrix_(t)
{
    className_ = std::string("SVMTraceFilter");
    foreFeatureVector_old_.resize(0,0);
    backFeatureVector_old_.resize(0,0);
    isTraceInitialPoint_ = false;
    libsvmclassifier.param_.svm_type = C_SVC;
    libsvmclassifier.param_.kernel_type = LINEAR;
    libsvmclassifier.param_.degree = 0;
    libsvmclassifier.param_.gamma = 0;	// 1/num_features
    libsvmclassifier.param_.coef0 = 0;
    libsvmclassifier.param_.nu = 0;
    libsvmclassifier.param_.cache_size = 100;
    libsvmclassifier.param_.C = 0.5;
    libsvmclassifier.param_.eps = 1e-3;
    libsvmclassifier.param_.p = 0.1;
    libsvmclassifier.param_.shrinking = 1;
    libsvmclassifier.param_.probability = 1;
    libsvmclassifier.param_.nr_weight = 0;
    libsvmclassifier.param_.weight_label = NULL;
    libsvmclassifier.param_.weight = NULL;
}

SVMTraceFilter::~SVMTraceFilter(void)
{
}

void SVMTraceFilter::SetInputRawCurves( const std::vector<VectorVec5d>& p )
{
    rawDendList_ = p;
}

void SVMTraceFilter::SetInputRawConInfo( const VectorMat2i& p)
{
    rawDendConInfo_ = p;
}

ProcStatPointer SVMTraceFilter::Update()
{
    if(!m_Input || !m_Back || !m_Bin || !paramPack){
        printf("error occurred in %s\n", className_.c_str());
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

    origImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Input);
    backImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short> >(m_Back);
    binImgPointer =  std::dynamic_pointer_cast<const Volume<NGCHAR> >(m_Bin);

    if(!origImgPointer || !backImgPointer || !binImgPointer ) {
        printf("error occurred in %s\n", className_.c_str());
        //LOG(ERROR) <<" error occurred in " <<className_;
        MAKEPROCESSSTATUS(resSta, false, className_, "input data is wrong.");
        return resSta;
    }
    //2016-3-22
    VectorVec5d allEndNodes;
    VectorVec3d allEndNodesDirection;
    if (isTraceInitialPoint_) {//2016-7-21
        allEndNodes.push_back(NGUtility::MakeVec5d(soma_, 1.0, 1.0));
        allEndNodesDirection.push_back(initDir_);
    }else {
        //LOG(INFO)<<"TileTreeRemoveExtraCurveEnd";
        TileTreeRemoveExtraCurveEnd(*origImgPointer, *backImgPointer, rawDendConInfo_, rawDendList_ );
        //2016-5-23 remove previous block branches
        /*VectorVec5d allEndNodesCp;
        VectorVec3d allEndNodesDirectionCp;
        for (size_t i = 0; i < allEndNodes.size(); ++ i) {
        Vec3d tmpVec = allEndNodes[i].block(0,0,3,1) - soma_;
        if (tmpVec.norm() < 5.0 || (tmpVec.norm() < 10.0 && (tmpVec.normalized()).dot(initDir_) < -0.5)) continue;
        allEndNodesCp.push_back(allEndNodes[i]);
        allEndNodesDirectionCp.push_back(allEndNodesDirection[i]);
        }
        allEndNodes.swap(allEndNodesCp);
        allEndNodesDirection.swap(allEndNodesDirectionCp);*/
    }
    //Guo Rui
    int nx = origImgPointer->x();
    int ny = origImgPointer->y();
    int nz = origImgPointer->z();
    indexImg_.SetZero();
    traceLabelMatrix_.SetZero();
    size_t curvesNum = rawDendList_.size();
    int xMin, xMax, yMin, yMax, zMin, zMax;
    //LOG(INFO)<<"for";
    for (size_t i = 0;i < curvesNum;i++) {
        const VectorVec5d &curCurve = rawDendList_[i];
        //vector<double> curRadius = curCurve(3);
        std::vector<int> curRadius;
        for (size_t k = 0; k < curCurve.size(); ++k) 
            curRadius.push_back( std::min(std::max(NGUtility::Round(curCurve[k](3)+0.25),1),2) );
        int currNumPoints = int(curCurve.size());

        for (int ik = 1;ik < currNumPoints-2;ik++) {
            int id1 = NGUtility::Round(curCurve[ik](0));
            int id2 = NGUtility::Round(curCurve[ik](1));
            int id3 = NGUtility::Round(curCurve[ik](2));
            int offset = std::max(curRadius[ik],0);
            Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, id1, id2, id3, offset, offset, offset, 0, nx -1, 0, ny -1, 0, nz - 1);
            for (int x = xMin;x <= xMax;++x) {
                for (int xx = yMin;xx <= yMax;++xx) {
                    for (int xxx = zMin;xxx <= zMax;++xxx) {
                        if (indexImg_(x,xx,xxx) == 0) {
                            indexImg_(x,xx,xxx) = int(i)+1;
                            traceLabelMatrix_(x,xx,xxx) = 1;
                        }
                    }
                }
            }
        }
    }
    //VecXd  w(16); double b;//double tmpw;
    //clock_t debugbeg = clock();
    if(!GetSVMClassifier( rawDendList_, *origImgPointer)){//, paramPack->w, paramPack->b
        printf("SVM failed and skip.\n");
        LOG(ERROR) <<"SVM failed and skip .";
        //return false;
        MAKEPROCESSSTATUS(resSta, false, className_, "SVM failed and skip.");
        return resSta;
    }
    //clock_t debugend = clock();
    //printf("%d ms eclipsed in SVM Train. \n", int(debugend - debugbeg));
    allEndNodes.clear();
    allEndNodesDirection.clear();
    TileTreeRemoveExtraCurveEndBySVM(*origImgPointer, rawDendConInfo_, rawDendList_, allEndNodes, allEndNodesDirection);
    //if delete soma, then add soma
    {
        bool addFlag = true;
        size_t minID; bool headFlag; double minDist = 100000.0, tmpDist1, tmpDist2;
        for (size_t k = 0; k < rawDendList_.size(); ++k) {
            auto &it = rawDendList_[k];
            tmpDist1 = (it[0].block(0, 0, 3, 1) - soma_).norm();
            tmpDist2 = (it.back().block(0, 0, 3, 1) - soma_).norm();
            if (tmpDist1 < 1.0 || tmpDist2 < 1.0) {
                addFlag = false;
                break;
            }
            if (tmpDist1  < minDist) {
                minID = k; minDist = tmpDist1; headFlag = true;
            }
            if (tmpDist2 < minDist) {
                minID = k; minDist = tmpDist2; headFlag = false;
            }
        }
        if (addFlag) {
            if (headFlag) {
                rawDendList_[minID].insert(rawDendList_[minID].begin(), NGUtility::MakeVec5d(soma_) );
            }
            else rawDendList_[minID].push_back( NGUtility::MakeVec5d(soma_));
        }
    }
    {
        /*VectorVec3d tmp;
        TraceUtil::GetPartVectorVec3d(rawDendList_[0], 0, rawDendList_[0].size() - 1, tmp);
        MatXd featureMat;  GetSVMSampleFeatures(*origImgPointer, tmp, featureMat);
        FILE *fp = fopen("F:/Linhuimin/hehe.txt", "w");
        for (int k = 0; k < featureMat.rows(); ++k) {
        for (int m = 0; m < featureMat.cols(); ++m) {
        fprintf(fp, "%lf ", featureMat(k, m));
        }
        fprintf(fp,"\n");
        }fclose(fp);
        system("pause");*/
    }
    {//2016-5-23 remove previous block branches
        VectorVec5d allEndNodesCp;
        VectorVec3d allEndNodesDirectionCp;
        for (size_t i = 0; i < allEndNodes.size(); ++i) {
            Vec3d tmpVec = allEndNodes[i].block(0, 0, 3, 1) - soma_;
            if (tmpVec.norm() < 5.0 || (tmpVec.norm() < 10.0 && (tmpVec.normalized()).dot(initDir_) < -0.5)) continue;
            allEndNodesCp.push_back(allEndNodes[i]);
            allEndNodesDirectionCp.push_back(allEndNodesDirection[i]);
        }
        allEndNodes.swap(allEndNodesCp);
        allEndNodesDirection.swap(allEndNodesDirectionCp);
    }
    
    paramPack->bifurcationValue_ = std::max(1.0, CalcBifurcationThrev(rawDendList_));
    printf("bifBinThrev:%lf\n", paramPack->bifurcationValue_);
    //trace
    for (size_t ii = 0; ii < allEndNodes.size(); ++ii) {
        int curEndNodeIndex = int(allEndNodes[ii](4));
        VectorVec5d resultCurve;
        int collideIndex;
        MeanShiftGrowTraceTotal(allEndNodes[ii].block(0,0,3,1), allEndNodesDirection[ii], *origImgPointer, *backImgPointer, curEndNodeIndex, 
            resultCurve, collideIndex);//paramPack->w, paramPack->b, 
        if (!resultCurve.empty()) {
            if ( int(allEndNodes[ii](3)) == 0 ) {
                VectorVec5d curCurve;// = rawDendList[curEndNodeIndex];
                std::copy(resultCurve.rbegin(), resultCurve.rend(), std::back_inserter(curCurve));
                std::copy(rawDendList_[curEndNodeIndex-1].begin(),rawDendList_[curEndNodeIndex-1].end(),std::back_inserter(curCurve));
                rawDendList_[curEndNodeIndex-1].swap(curCurve);
                Mat2i& tmpMat = rawDendConInfo_[curEndNodeIndex-1];
                tmpMat(0,0) = collideIndex;
            }
            if ( int(allEndNodes[ii](3)) == 1 ) {
                VectorVec5d curCurve;// = rawDendList[curEndNodeIndex];
                std::copy(rawDendList_[curEndNodeIndex-1].begin(),rawDendList_[curEndNodeIndex-1].end(),std::back_inserter(curCurve));
                std::copy(resultCurve.begin(), resultCurve.end(), std::back_inserter(curCurve));
                rawDendList_[curEndNodeIndex-1].swap(curCurve);
                Mat2i& tmpMat = rawDendConInfo_[curEndNodeIndex-1];
                tmpMat(0,1) = collideIndex;
            }
        }
    }
    std::vector<VectorVec5d> newTracedCurves(rawDendList_);
    VectorVec3i CubeLabelPts;
    globalID_ = 0;
    //find new bifurcation and trace
    int iterNum = 0;
    int maxIterNum = paramPack->isLocalTrace_ ? 3 : 1;
    while (!newTracedCurves.empty() && iterNum < maxIterNum) {
        ++iterNum;
        std::vector<VectorVec5d> newTracedCurvesCp;
        for (size_t ij = 0; ij < newTracedCurves.size(); ++ij) {
            if (!newTracedCurves[ij].empty()) {
                ++globalID_;
                VectorVec3d tmpTracedCurve;
                NGUtility::GetPartVectorVec3d(newTracedCurves[ij], 0, int(newTracedCurves[ij].size()) - 1, tmpTracedCurve);
                CubeLabel(tmpTracedCurve, CubeLabelPts);
                VectorVec3d bifurcationDirections;
                std::vector<VectorVec3d> bifurcationClustre;
                DetectBifurcation(CubeLabelPts, bifurcationDirections, bifurcationClustre);
                if (bifurcationClustre.size() > 10) continue;
                /*
                * remove cross curves
                */
                {
                    if(!cff) cff = CrossBranchFilter::New();
                    cff->SetParam(paramPack);
                    std::vector<VectorVec3d> tmpP; tmpP.push_back(tmpTracedCurve);
                    cff->SetData(tmpP, soma_, initDir_, bifurcationClustre, bifurcationDirections);
                    auto cfres = cff->Update();
                    if (!cfres->success())
                        return cfres;
                    /* std::vector<VectorVec3d> bifurcationClustreCp;
                     VectorVec3d bifurcationDirectionsCp;
                     std::set<size_t> removeInd;
                     std::vector<VectorVec3d> tmpP; tmpP.push_back(tmpTracedCurve);
                     RemoveCrossBifurcation(tmpP, bifurcationClustre, removeInd);
                     for (size_t i = 0; i < bifurcationClustre.size(); ++ i) {
                     Vec3d tmpVec = bifurcationClustre[i].front().block(0,0,3,1) - soma_;
                     if (tmpVec.norm() < 5.0 || (tmpVec.norm() < 10.0 &&
                     (tmpVec.normalized()).dot(initDir_) < -0.5)) continue;
                     if (removeInd.count(i) != 0) continue;
                     bifurcationClustreCp.push_back(bifurcationClustre[i]);
                     bifurcationDirectionsCp.push_back(bifurcationDirections[i]);
                     }
                     bifurcationClustre.swap(bifurcationClustreCp);
                     bifurcationDirections.swap(bifurcationDirectionsCp);*/
                }
                
                Mat2i tmpMat2i;
                for (size_t ii = 0; ii < bifurcationClustre.size(); ++ii) {
                    VectorVec5d resultCurve;
                    VectorVec5d resultCurveModify;
                    int collideIndex;
                    int id = int(rawDendList_.size()) + 1;
                    MeanShiftGrowTraceTotal(bifurcationClustre[ii][3].block(0,0,3,1), bifurcationDirections[ii], *origImgPointer, *backImgPointer, id,
                        resultCurve, collideIndex);//paramPack->w, paramPack->b, 
                    if(resultCurve.size() < 4) continue;
                    if (!resultCurve.empty()) {
                        //find nearest point in parent curve
                        //fill in bifurcationClustre 0~2 points
                        rawDendList_.push_back(VectorVec5d());
                        size_t headInd;
                        double dist = 100000.0;
                        Vec3d tmpFirPt= bifurcationClustre[ii][0].block(0,0,3,1);
                        for (size_t k = 0; k < newTracedCurves[ij].size(); ++k) {
                            if ( (tmpFirPt - newTracedCurves[ij][k].block(0,0,3,1)).norm() < dist  ) {
                                headInd = k;
                                dist = (tmpFirPt - newTracedCurves[ij][k].block(0,0,3,1)).norm();
                            }
                        }
                        rawDendList_.back().push_back(newTracedCurves[ij][headInd]);
                        rawDendList_.back().push_back(
                            NGUtility::MakeVec5d(bifurcationClustre[ii][0](0), bifurcationClustre[ii][0](1), bifurcationClustre[ii][0](2), 0, 0));
                        rawDendList_.back().push_back(
                            NGUtility::MakeVec5d(bifurcationClustre[ii][1](0), bifurcationClustre[ii][1](1), bifurcationClustre[ii][1](2), 0, 0));
                        std::copy(resultCurve.begin(), resultCurve.end(), std::back_inserter(rawDendList_.back()));
                        tmpMat2i << globalID_, 0, 0, 0;
                        rawDendConInfo_.push_back(tmpMat2i);
                        newTracedCurvesCp.push_back(rawDendList_.back());
                    }
                    
                }
            }
        }
        newTracedCurves.swap(newTracedCurvesCp);
    }
    paramPack->bifurcationValue_ = std::max(1.0, CalcBifurcationThrev(rawDendList_));
    printf("svm update bifBinThrev:%lf\n", paramPack->bifurcationValue_);
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

void SVMTraceFilter::GetSVMSample( const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, 
                                  MatXd &foreSamples, MatXd& backSamples )
{
    size_t rawDendListNum = rawDendList.size();
    VectorVec3d forePtSet;
    int kk = 0,nx,ny,nz,foreSamplesNum,backSamplesNum;
    //system("pause");
    //time_t beg = time(0);
    for (size_t i = 0; i < rawDendListNum; ++ i) {
        const VectorVec5d &curForeCurve = rawDendList[i];
        if (curForeCurve.size() < 7) continue;
        VectorVec3d tmpCurForeCurve;
        NGUtility::GetPartVectorVec3dStep(curForeCurve, 4, int(curForeCurve.size()) - 6, 1, tmpCurForeCurve);
        std::copy(tmpCurForeCurve.begin(), tmpCurForeCurve.end(), std::back_inserter(forePtSet));
        kk += int(tmpCurForeCurve.size());
        if( kk > 1000000-1000) break;
    }
    if (forePtSet.empty() && foreFeatureVector_old_.cols() == 0) {
        printf("error in GetSVMSample, the forePtSet is empty.\n");
        LOG(ERROR)<<"error in GetSVMSample, the forePtSet is empty.";
        return;
    }
    
    int maxNum = 500;
    if ( int(forePtSet.size()) < maxNum)
        GetSVMSampleFeatures(origImg,forePtSet,foreSamples);
    else{
        std::random_shuffle(forePtSet.begin(),forePtSet.end());
        VectorVec3d tmpForePtSet(maxNum);
        std::copy(forePtSet.begin(), forePtSet.begin() + maxNum, tmpForePtSet.begin());
        GetSVMSampleFeatures(origImg,tmpForePtSet, foreSamples);
    }   
    nx = origImg.x();
    ny = origImg.y();
    nz = origImg.z();
    int featureNum = foreSamples.rows();
    VectorVec3d backPtSet;
    srand(time(NULL));
    for(int i = 1; i <= foreSamples.cols(); ++i){
        double x = (double(rand())/double(RAND_MAX)*(nx-2)+1);
        double y = (double(rand())/double(RAND_MAX)*(ny-2)+1);
        double z = (double(rand())/double(RAND_MAX)*(nz-2)+1);
        backPtSet.push_back(Vec3d(x,y,z));//round([rand*(nx-2)+1,rand*(ny-2)+1,rand*(nz-2)+1]')
        //backPtSet.push_back(Vec3d(Round(double(i)/double(maxNum)*(nx-2)),Round(double(i)/double(maxNum)*(ny-2)),Round(double(i)/double(maxNum)*(nz-2))));
    }
    //LOG(INFO)<<"GetSVMSampleFeatures";
    GetSVMSampleFeatures(origImg,backPtSet,backSamples); 
    MatXd tmpAllSamples(backSamples.rows(), backSamples.cols() + foreSamples.cols());
    //tmpAllSamples << backSamples, foreSamples;
    tmpAllSamples.block(0,0,featureNum, backSamples.cols()) = backSamples;
    tmpAllSamples.block(0,backSamples.cols(),featureNum, foreSamples.cols()) = foreSamples;
    MatXd convMatrix = tmpAllSamples.transpose() * backSamples;//convMatrix=([backSamples,foreSamples]'*backSamples);
    foreSamplesNum = int(foreSamples.cols());
    backSamplesNum = int(backSamples.cols());
    std::vector<int> removeFlagList(backSamplesNum, 0);
    for (size_t i = 0;i < backSamplesNum; ++i){
        double maxv1 = (convMatrix.block(0,i,backSamplesNum, 1)).maxCoeff();
        double maxv2 = (convMatrix.block(backSamplesNum, i , foreSamplesNum, 1)).maxCoeff();
        if (maxv1 < maxv2) removeFlagList[i] = 1;
    }
    //backSamples=backSamples(:,removeFlagList==0);
    int newBackSampleNum = int(removeFlagList.size()) - std::accumulate(removeFlagList.begin(), removeFlagList.end(), 0);
    MatXd cpBackSamples(backSamples.rows(), newBackSampleNum);
    int idx = 0;
    for (size_t i = 0;i < backSamplesNum; ++i){
        if (removeFlagList[i] == 0) {
            cpBackSamples.col(idx++) = backSamples.col(i);
        }
    }
    backSamples.swap(cpBackSamples);
    int backSamplesNumFinal = backSamples.cols();
    if (foreFeatureVector_old_.size()>0) {
        if (foreSamples.cols() < maxNum) {
            int lostNum = maxNum - foreSamples.cols();
            if (lostNum > foreFeatureVector_old_.cols()) {
                MatXd tmpForeFeature(featureNum, foreSamples.cols() + foreFeatureVector_old_.cols());
                tmpForeFeature << foreSamples, foreFeatureVector_old_;
                foreSamples.swap(tmpForeFeature);
                MatXd tmpBackFeature(featureNum, backSamples.cols() + backFeatureVector_old_.cols());
                tmpBackFeature << backSamples, backFeatureVector_old_;
                backSamples.swap(tmpBackFeature);
            } else{
                MatXd tmpForeFeature(featureNum, maxNum);
                MatXd tmpBackFeature(featureNum, maxNum);
                tmpForeFeature.block(0,0,featureNum,foreSamplesNum) = foreSamples;
                tmpBackFeature.block(0,0,featureNum,backSamplesNumFinal) = backSamples;
                for (int i = foreSamplesNum; i < maxNum; ++i) {
                    int randomNum = std::max(0, std::min(int(double(rand()) / double(RAND_MAX) * foreFeatureVector_old_.cols()), int(foreFeatureVector_old_.cols())-1));
                    tmpForeFeature.col(i) = foreFeatureVector_old_.col( randomNum );
                }//2016-7-9
                for (int i = backSamplesNumFinal; i < maxNum; ++i) {
                    int randomNum = std::max(0, std::min(int(double(rand()) / double(RAND_MAX) * backFeatureVector_old_.cols()), int(backFeatureVector_old_.cols())-1));
                    tmpBackFeature.col(i) = backFeatureVector_old_.col(randomNum);
                }
                foreSamples.swap(tmpForeFeature);
                backSamples.swap(tmpBackFeature);
            }
        }
    }
    //save current sample
    foreFeatureVector_old_ = foreSamples;
    backFeatureVector_old_ = backSamples;

    /*{
        FILE *fp = fopen("F:/Linhuimin/tmpForePtSet.swc", "w");
        int ind = 0;
        for (auto &it : forePtSet) {
        fprintf(fp, "%d 1 %lf %lf %lf 1 -1\n", ++ind, it(0), it(1), it(2));
        }
        fclose(fp);
        }

        {
        FILE *fp = fopen("F:/Linhuimin/tmpBackPtSet.swc", "w");
        int ind = 0;
        for (auto &it : backPtSet) {
        fprintf(fp, "%d 1 %lf %lf %lf 1 -1\n", ++ind, it(0), it(1), it(2));
        }
        fclose(fp);
        }
        system("pause");*/
}



void SVMTraceFilter::GetSVMSampleBySoma(const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, MatXd &foreSamples, MatXd& backSamples)
{
    size_t rawDendListNum = rawDendList.size();
    VectorVec4d forePtSet;
    int kk = 0, nx, ny, nz, foreSamplesNum, backSamplesNum;
    //system("pause");
    //time_t beg = time(0);
    for (size_t i = 0; i < rawDendListNum; ++i) {
        const VectorVec5d &curForeCurve = rawDendList[i];
        if (curForeCurve.size() < 7) continue;
        VectorVec3d tmpCurForeCurve;
        //TraceUtil::GetPartVectorVec3dStep(curForeCurve, 4, curForeCurve.size() - 6, 1, tmpCurForeCurve);
        NGUtility::GetPartVectorVec3dStep(curForeCurve, 4, int(curForeCurve.size()) - 6, 2, tmpCurForeCurve);//2017-5-19
        forePtSet.reserve(tmpCurForeCurve.size());
        for (auto &it : tmpCurForeCurve) {
            forePtSet.push_back(Vec4d(it(0), it(1), it(2), WeighRayValue(Vec3d(it(0), it(1), it(2)), origImg)));
        }
        kk += int(tmpCurForeCurve.size());
        if (kk > 1000000 - 1000) break;
    }

    if (forePtSet.empty()) {
        printf("error in GetSVMSample, the forePtSet is empty.\n");
        return;
    }

    size_t maxNum = 500;
    if (int(forePtSet.size()) < maxNum){
        printf("svm no 500\n");
        VectorVec3d tmpForePtSet(forePtSet.size());
        for (size_t i = 0; i < tmpForePtSet.size(); ++i) tmpForePtSet[i] = forePtSet[i].block(0, 0, 3, 1);
        GetSVMSampleFeatures(origImg, tmpForePtSet, foreSamples);
    }
    else{
        //std::random_shuffle(forePtSet.begin(), forePtSet.end());
        VectorVec3d tmpForePtSet(maxNum);
        std::sort(forePtSet.begin(), forePtSet.end(), [&](const Vec4d& lhs, const Vec4d& rhs){return lhs(3) < rhs(3); });
        size_t beg = forePtSet.size() / 2lu - 250;
        size_t end = beg + 500;
        for (size_t i = beg; i <end; ++i) tmpForePtSet[i - beg] = forePtSet[i].block(0, 0, 3, 1);
        GetSVMSampleFeatures(origImg, tmpForePtSet, foreSamples);
    }

    nx = origImg.x();
    ny = origImg.y();
    nz = origImg.z();
    int featureNum = foreSamples.rows();
    VectorVec3d backPtSet;
    srand(time(NULL));
    for (int i = 1; i <= foreSamples.cols(); ++i){
        double x = (double(rand()) / double(RAND_MAX)*(nx - 2) + 1);
        double y = (double(rand()) / double(RAND_MAX)*(ny - 2) + 1);
        double z = (double(rand()) / double(RAND_MAX)*(nz - 2) + 1);
        backPtSet.push_back(Vec3d(x, y, z));//round([rand*(nx-2)+1,rand*(ny-2)+1,rand*(nz-2)+1]')
        //backPtSet.push_back(Vec3d(Round(double(i)/double(maxNum)*(nx-2)),Round(double(i)/double(maxNum)*(ny-2)),Round(double(i)/double(maxNum)*(nz-2))));
    }

    //LOG(INFO)<<"GetSVMSampleFeatures";
    GetSVMSampleFeatures(origImg, backPtSet, backSamples);
    MatXd tmpAllSamples(backSamples.rows(), backSamples.cols() + foreSamples.cols());
    //tmpAllSamples << backSamples, foreSamples;
    tmpAllSamples.block(0, 0, featureNum, backSamples.cols()) = backSamples;
    tmpAllSamples.block(0, backSamples.cols(), featureNum, foreSamples.cols()) = foreSamples;
    MatXd convMatrix = tmpAllSamples.transpose() * backSamples;//convMatrix=([backSamples,foreSamples]'*backSamples);
    foreSamplesNum = int(foreSamples.cols());
    backSamplesNum = int(backSamples.cols());
    std::vector<int> removeFlagList(backSamplesNum, 0);
    for (size_t i = 0; i < backSamplesNum; ++i){
        double maxv1 = (convMatrix.block(0, i, backSamplesNum, 1)).maxCoeff();
        double maxv2 = (convMatrix.block(backSamplesNum, i, foreSamplesNum, 1)).maxCoeff();
        if (maxv1 < maxv2) removeFlagList[i] = 1;
    }
    //backSamples=backSamples(:,removeFlagList==0);
    int newBackSampleNum = int(removeFlagList.size()) - std::accumulate(removeFlagList.begin(), removeFlagList.end(), 0);
    MatXd cpBackSamples(backSamples.rows(), newBackSampleNum);
    int idx = 0;
    for (size_t i = 0; i < backSamplesNum; ++i){
        if (removeFlagList[i] == 0) {
            cpBackSamples.col(idx++) = backSamples.col(i);
        }
    }
    backSamples.swap(cpBackSamples);
}


void SVMTraceFilter::GetSVMSampleFeatures( const SVolume &origImg, VectorVec3d &samplePtSet, MatXd &featureVector )
{
    //LOG(INFO)<<"GetSVMSampleFeatures start";
    size_t samplePtSetNum = samplePtSet.size();
    int allPtNum,validSignalSum;
    std::vector<double> rayNodeWet, extractThrevList;
    int radiuList = 9;
    VecXd signalVec(9);
    double curSeedValue;
    featureVector = Eigen::MatrixXd::Zero(9,samplePtSetNum);
    MatXd samplePtSetFeature;
    //time_t beg = time(0);
    double ratioVal(0.0);
    for (size_t ii = 0;ii < samplePtSetNum; ++ii){
        const Vec3d &curSeedPt = samplePtSet[ii];
        VectorVec3d tmp;
        tmp.push_back(curSeedPt);
        WeighRayValue(tmp,origImg,rayNodeWet);
        curSeedValue = rayNodeWet[0];    
        extractThrevList.clear();
        VectorDecorator<double>(extractThrevList) << 1.0*curSeedValue,0.975*curSeedValue,0.95*curSeedValue,0.925*curSeedValue,
            0.9*curSeedValue, 0.875*curSeedValue, 0.85*curSeedValue, 0.825*curSeedValue, 0.8*curSeedValue;
        //size_t ny = radiuList.size();
        size_t nx = extractThrevList.size();
        samplePtSetFeature.setZero(nx,1);//ny
        Vec3i roundCurSeedPt(NGUtility::Round(curSeedPt(0)), NGUtility::Round(curSeedPt(1)), NGUtility::Round(curSeedPt(2)));
        //LOG(INFO)<<"GetSVMSampleFeatureByRadiusThrev";
        for (size_t i =0;i < nx;++i){
            //for (size_t j = 0;j < ny;++j){
            GetSVMSampleFeatureByRadiusThrev(origImg, roundCurSeedPt ,extractThrevList[i], radiuList,allPtNum,validSignalSum);
            ratioVal = double(validSignalSum) / double(allPtNum);
            //samplePtSetFeature(i,0) = ratioVal;
            if (ratioVal > 0.7 ) {
                for (size_t j =i;j < nx;++j) samplePtSetFeature(j,0) = 1.0;
                break;
            }
            else samplePtSetFeature(i,0) = ratioVal;
            //}
        }
        //LOG(INFO)<<"GetSVMSampleFeatureByRadiusThrev complete";
        /*signalVec << samplePtSetFeature.row(0).transpose(), samplePtSetFeature.row(1).transpose(),
            samplePtSetFeature.row(2).transpose(), samplePtSetFeature.row(3).transpose();*/
        signalVec << samplePtSetFeature.col(0);
        featureVector.col(ii) = signalVec;
    }
    //LOG(INFO)<<"GetSVMSampleFeatures complete";
    //time_t end = time(0);
    //time_t eclapse = end - beg;
}

void SVMTraceFilter::GetSVMSampleFeatureByRadiusThrev( const SVolume &origImg, Vec3i seedPt, double extractThrev, 
                                                      int radius, int& allPtNum, int& validSignalSum )
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();

    seedPt(0) = std::min(std::max(seedPt(0),0), nx - 1);
    seedPt(1) = std::min(std::max(seedPt(1),0), ny - 1);
    seedPt(2) = std::min(std::max(seedPt(2),0), nz - 1);
    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, seedPt(0), seedPt(1), seedPt(2), radius, radius, radius, 0, nx-1, 0, ny - 1, 0, nz -1);

    int xLen = xMax - xMin + 1, yLen = yMax - yMin + 1, zLen = zMax - zMin + 1;

    //SVolume subOrigImg;
    ExtractArea(origImg, xMin, xMax, yMin, yMax, zMin, zMax, subOrigImg_);
    allPtNum = xLen * yLen * zLen;
    VectorVec3i growPtSet;
    VectorVec3i incPtSet;
    incPtSet.push_back(Vec3i(seedPt(0) - xMin, seedPt(1) - yMin, seedPt(2) - zMin));

    //CVolume signalMatrix;
    signalMatrix_.SetSize(xLen,yLen,zLen);
    signalMatrix_(incPtSet[0](0), incPtSet[0](1), incPtSet[0](2)) = 1;
    validSignalSum = 1;
    for (int i = 1; i <= 8;i++) {//for safe use int
        SignalPointsIdentifyAreaGrowing(growPtSet,incPtSet, subOrigImg_,signalMatrix_,extractThrev);
        //printf("%lf\n",extractThrev);
        //system("pause");
        //SignalPointsIdentifyAreaGrowing(growPtSet,origImg,signalMatrix,extractThrev);
        if (growPtSet.empty()) {
            --i;
            break;
        }
    }
    GetAreaSum(signalMatrix_, 0, signalMatrix_.x() - 1, 0, signalMatrix_.y() - 1, 0, signalMatrix_.z() - 1, validSignalSum);
}

void SVMTraceFilter::SignalPointsIdentifyAreaGrowing(VectorVec3i& growPtSet, VectorVec3i& incPtSet, const SVolume& subOrigImg, CVolume& signalMatrix, double threv, int offset, int ptNumThrev)
{
    int nx = subOrigImg.x();
    int ny = subOrigImg.y();
    int nz = subOrigImg.z();
    size_t nxx = incPtSet.size();
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VectorVec3i resultGrowPtSet;//TODO:debug
    VectorVec3i tmpGrowPtSet;
    Vec3i tmp;
    for (size_t i=0; i < nxx; i++) {
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, incPtSet[i](0), incPtSet[i](1), incPtSet[i](2),
            offset, offset, offset, 0, nx-1, 0, ny - 1, 0, nz -1);
        VectorVec3i curGrowPtSet;
        for (int ii = xMin; ii <= xMax; ++ii){
            for(int jj = yMin; jj <= yMax; ++jj){
                for(int ij = zMin; ij <= zMax; ++ij){
                    tmp << ii,jj,ij;
                    //int test1 = signalMatrix(ii, jj, ij);
                    //int test2 = subOrigImg(ii, jj, ij);
                    if (signalMatrix(ii, jj, ij) == 0 && subOrigImg(ii, jj, ij) > threv) 
                        curGrowPtSet.push_back(tmp);
                }
            }
        }
        if (curGrowPtSet.size() > ptNumThrev){
            std::copy(curGrowPtSet.begin(), curGrowPtSet.end(), std::back_inserter(resultGrowPtSet));
            std::copy(curGrowPtSet.begin(), curGrowPtSet.end(), std::back_inserter(tmpGrowPtSet));
            for (size_t j = 0;j < curGrowPtSet.size(); ++j){
                signalMatrix(curGrowPtSet[j](0),curGrowPtSet[j](1),curGrowPtSet[j](2)) = 1;
            }
        }
    }
    resultGrowPtSet.swap(growPtSet);
    tmpGrowPtSet.swap(incPtSet);
}

void SVMTraceFilter::SVMClassfier( const MatXd& foreSamples, const MatXd& backSamples, double regularParam, VecXd& w, double &b )
{
    int foreSamplesNum = int(foreSamples.cols());
    int backSamplesNum = int(backSamples.cols());
    MatXd allSamples;
    allSamples.setZero(foreSamples.rows(),foreSamples.cols()+backSamples.cols());
    allSamples << foreSamples , -backSamples;//X^T * Y
    VecXd Y;//standard answer
    Y.setOnes(foreSamplesNum+backSamplesNum);
    Y.tail(backSamplesNum) *=-1;
    MatXd matrixS;
    MatXd tmp = MatXd::Identity(foreSamplesNum+backSamplesNum,foreSamplesNum+backSamplesNum);
    matrixS.setZero(foreSamplesNum+backSamplesNum+1,foreSamplesNum+backSamplesNum+1);
    matrixS.block(1,1,foreSamplesNum+backSamplesNum, foreSamplesNum+backSamplesNum) << allSamples.transpose() * allSamples + 1 /regularParam * tmp;
    matrixS.block(0,1,1,foreSamplesNum+backSamplesNum) = -Y.transpose();
    matrixS.block(1,0,foreSamplesNum+backSamplesNum,1) = Y;
    VecXd B;
    B.setOnes(foreSamplesNum+backSamplesNum+1);
    B.head(1) *= 0; 
    VecXd solutionMatrix = matrixS.inverse() * B;
    b = solutionMatrix(0);
    w = allSamples * solutionMatrix.segment(1,foreSamplesNum+backSamplesNum);
}


bool SVMTraceFilter::IsNearBoundary(const Vec3d& arg, int x, int y, int z)
{
    return arg(0)  <= paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(0) > double(x) - paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(1) <= paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(1) > double(y) - paramPack->boundaryDistanceThreshold_ / 2.0
        || arg(2) <= 3.0
        || arg(2) > double(z - 3);
}

//initPt and initDir are passed by value
void SVMTraceFilter::MeanShiftGrowTraceTotal( Vec3d initPt, Vec3d initDir, const SVolume& origImg, const SVolume& backImg, 
                                             int curEndNodeIndex, 
                                             VectorVec5d& resultCurve, int& collideIndex )
{
    if (IsNearBoundary(initPt, origImg.x(), origImg.y(), origImg.z())) return;
    Vec5d tmpPt; tmpPt << initPt(0),initPt(1),initPt(2),0,0;
    resultCurve.push_back(tmpPt);
    bool traceFlag = true;
    int i(0);
    //Vec3d tmpDir = initDir;
    collideIndex = 0;
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    Vec3d nextNode,nextDir,curPt, tmpDir; nextDir.setZero();
    VecXd tmpValueList,signalIntensity,backIntensity;
    MatXd curFeature,tmp;
    VectorVec3d tmpCurve;
    double maxCriterion,internalDist;
    while (traceFlag) {
        //LOG(INFO)<<"MeanShiftGrowTrace";
        MeanShiftGrowTrace(initPt,initDir,origImg,backImg,nextNode,nextDir);
        if ( (nextDir -Vec3d(0,0,0)).norm() - 0.0 > 0.001) {
            ++i;
            tmpPt.block(0,0,3,1) = nextNode;
            resultCurve.push_back(tmpPt);
            tmpDir = initDir;
            initDir = nextDir;
            initPt = nextNode;
            if ( initPt.minCoeff() < 1 || initPt(0) > nx-2 || initPt(1) > ny-2 || initPt(2) > nz-2) {//
                traceFlag = false;
            }
            if (i > 3) {
                VectorVec3d tmpResultCurve;
                NGUtility::GetPartVectorVec3d(resultCurve, i - 2, i, tmpResultCurve);
                WeighRayValue(tmpResultCurve, origImg, signalIntensity);
                WeighRayValue(tmpResultCurve, backImg, backIntensity);

                tmpValueList = signalIntensity - backIntensity - backIntensity.cwiseSqrt() *  3.0;
                if (tmpValueList.maxCoeff() < 0 ) {
                    GetSVMSampleFeatures(origImg,tmpResultCurve,curFeature);
                    //tmp = w.transpose() * curFeature;
                    //tmp.array() += b;
                    libsvmclassifier.Predict(curFeature, tmp);
                    maxCriterion = tmp.maxCoeff();
                    //tmpValueList1 = signalIntensity - backIntensity - backIntensity.cwiseSqrt() *  0.1;
                    if (maxCriterion < 0 ) {//|| tmpValueList1.maxCoeff() < 0 
                        //std::cout <<"svm 1:" << tmp.array() <<"\n feature:\n" << curFeature <<std::endl;
                        traceFlag = 0;
                        i -= 3;
                    }
                }

            }
            if (i>1) {
                internalDist = (resultCurve[i].block(0,0,3,1)-resultCurve[i-1].block(0,0,3,1)).norm();
                if (internalDist < 0.02) {
                    printf("SVM 2\n");
                    traceFlag = false;	
                }
            }

            if (i > maxSteps_ - 1) {
                traceFlag = false;
            }

            if (tmpDir.transpose() * initDir < -0.5 ) {
                printf("svm 3\n");
                traceFlag = false;
            }

            if (i > 1) {
                curPt = resultCurve[i].block(0,0,3,1);
                if (traceFlag) {
                    DetectCollisionByMeanShift(curPt,indexImg_,origImg,curEndNodeIndex,traceFlag,collideIndex);
                }
            }	
        }else{
            printf("svm 4\n");
            traceFlag = 0;
        }
    }

    {
        int n = i+1;
        VectorVec5d tmpResultCurve(n);
        std::copy(resultCurve.begin(), resultCurve.begin() + n, tmpResultCurve.begin());
        resultCurve.swap(tmpResultCurve);
    }

    {
        std::vector<double> r,v;
        VectorVec3d tmpResultCurve;
        NGUtility::GetPartVectorVec3d(resultCurve, 0, int(resultCurve.size() - 1), tmpResultCurve);
        //LOG(INFO)<<"CalcParmOfCurveNodeListSVM";
        CalcParmOfCurveNodeListSVM(origImg,backImg, tmpResultCurve, r,v);
        for (size_t j = 0; j < resultCurve.size(); ++j) {
            resultCurve[j](3) = r[j];
            resultCurve[j](4) = v[j];
        }
    }
    if (i > 1) {
        VectorVec5d interplCurve;
        Vec5d half;
        std::vector<double> fillRadius;
        for (size_t x = 0; x < resultCurve.size() - 1u; ++x){
            interplCurve.push_back(resultCurve[x]);
            half = 0.5 * (resultCurve[x] + resultCurve[x+1]);
            interplCurve.push_back(half);
        }
        interplCurve.push_back(resultCurve.back());

        for (size_t xx = 0;xx < interplCurve.size();++xx) {
            fillRadius.push_back(std::min(std::max(NGUtility::Round(interplCurve[xx](3) + 0.25), 3), 2));
        }

        //fillRadius = std::min(std::max(fillRadius,3),2);
        int id1(0), id2(0), id3(0);
        int xMin, xMax, yMin, yMax, zMin, zMax;
        for (size_t ik =1; ik <interplCurve.size();++ik) {

            id1 = NGUtility::Round(interplCurve[ik](0));
            id2 = NGUtility::Round(interplCurve[ik](1));
            id3 = NGUtility::Round(interplCurve[ik](2));

            Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, id1, id2, id3,
                std::max(0, (int)fillRadius[ik]), std::max(0, (int)fillRadius[ik]), std::max(0, (int)fillRadius[ik]),
                0, (int)nx - 1, 0, (int)ny - 1, 0, (int)nz - 1);

            for (int ii = xMin; ii <= xMax; ++ii){
                for (int jj = yMin; jj <= yMax; ++jj){
                    for (int kk = zMin; kk <= zMax; ++kk){
                        if (indexImg_(ii, jj, kk) == 0){
                            indexImg_(ii, jj, kk) = curEndNodeIndex;
                            traceLabelMatrix_(ii,jj,kk) = 1;//2016-8-1
                        }
                    }
                }
            }
        }


    }
    else resultCurve.clear();

}

void SVMTraceFilter::MeanShiftGrowTrace( const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, 
                                        Vec3d& nextNode, Vec3d& nextDir )
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    VectorMatXd samplePlaneSet;
    std::vector<VectorVec4d> extractPtSet;
    VectorVec3d centerPoints,newPositions;
    std::vector<double> centerPtValue;
    ExtractSampleConePtSet(initPt,initDir,origImg,backImg,samplePlaneSet,extractPtSet,centerPoints);
    WeighRayValue(centerPoints, origImg, centerPtValue);
    bool flag(true);
    double centerPtMin(100000000.0), centerPtXMax(0.0), centerPtYMax(0.0), centerPtZMax(0.0);
    for (size_t i = 0; i < centerPoints.size(); ++i) {
        centerPtMin = std::min(centerPtMin, centerPoints[i].minCoeff());
        centerPtXMax = std::max(centerPtXMax, centerPoints[i](0));
        centerPtYMax = std::max(centerPtYMax, centerPoints[i](1));
        centerPtZMax = std::max(centerPtZMax, centerPoints[i](2));
    }

    if (centerPtMin < 1.0 || int(centerPtXMax) > nx-2 || int(centerPtYMax) > ny-2 || int(centerPtZMax) > nz-2) 
        flag = false;

    if (std::accumulate(centerPtValue.begin(), centerPtValue.end(), 0.0) / double(centerPtValue.size()) > 20.0 && flag ) {
        CalcPositionByMeanshift(extractPtSet, centerPoints, newPositions);
        nextNode = newPositions[1];
        Vec3d nextDir1;
        Principald(newPositions,nextDir1);
        nextDir = 0.5 * (initDir + nextDir1);
        nextDir /= nextDir.norm()+0.0001;
    }
    else{
        nextDir.setZero();
        nextNode.setZero();
    }
}
void SVMTraceFilter::WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,
                                   std::vector<double> &rayNodeWet)
{
    typedef double spheredist;
    int nxss = int(rayNode.size());
    rayNodeWet.clear();

    //坐标
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

void SVMTraceFilter::ExtractSampleConePtSet( const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, 
                                            VectorMatXd& samplePlaneSet, std::vector<VectorVec4d> &extractPtSet, VectorVec3d& centerPoints )
{
    Vec3d x2, x3;
    CalcOrthoBasis(initDir, x2, x3);
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    samplePlaneSet.resize(4);
    extractPtSet.resize(4);
    centerPoints.resize(5);
    centerPoints[0] = initPt + initDir;
    Vec3d curPt;
    for (int jj = 1; jj < 5; ++jj) {
        MatXd perpendicularPlane = MatXd::Zero(3+2*jj, 3+2*jj);
        Vec3d centerPt = initPt + (1+0.5*double(jj)) * initDir;
        VectorVec4d curExtractPtSet;
        centerPoints[jj] = centerPt;
        for (int i = -jj-1; i <= jj+1; ++i) {
            for (int j = -jj - 1; j <= jj+1; ++j) {
                curPt = centerPt + double(i) * x2 + double(j) * x3;
                curPt(0) = std::max( std::min(curPt(0), double(nx - 1) ), 0.0);
                curPt(1) = std::max( std::min(curPt(1), double(ny - 1) ), 0.0);
                curPt(2) = std::max( std::min(curPt(2), double(nz - 1) ), 0.0);
                double signalIntensity = std::max(0.0, WeighRayValue(curPt, origImg) - WeighRayValue(curPt, backImg) );
                curExtractPtSet.push_back(Vec4d(curPt(0), curPt(1), curPt(2), signalIntensity));
                perpendicularPlane( i+jj+1 , j+jj+1) = signalIntensity;
            }
        }
        extractPtSet[jj-1].swap(curExtractPtSet);
        samplePlaneSet[jj-1].swap(perpendicularPlane);
    }
}

//centerPoints is passed by value
void SVMTraceFilter::CalcPositionByMeanshift( const std::vector<VectorVec4d> &extractPtSet, VectorVec3d centerPoints,//centerPoints will be modified
                                             VectorVec3d& newPositions )
{
    size_t centerPointsNum = centerPoints.size();
    Vec3d pos0, pos1, pos2, pos3;
    Vec3d pos0bak, pos1bak, pos2bak, pos3bak;
    pos0bak.setZero();
    pos1bak.setZero();
    pos2bak.setZero();
    pos3bak.setZero();
    for (size_t kkj = 0u; kkj < 100u; ++kkj) {
        pos0bak = pos0;
        pos1bak = pos1;
        pos2bak = pos2;
        pos3bak = pos3;

        CalcOnePosByMeanShift(extractPtSet[0], centerPoints[1], pos0);
        pos0 = 0.8 * pos0 + 0.1 * (centerPoints[0] + centerPoints[2]);

        CalcOnePosByMeanShift(extractPtSet[1], centerPoints[2], pos1);
        pos1 = 0.8 * pos1 + 0.1 * (centerPoints[1] + centerPoints[3]);

        CalcOnePosByMeanShift(extractPtSet[2], centerPoints[3], pos2);
        pos2 = 0.8 * pos2+ 0.1 * (centerPoints[2] + centerPoints[4]);

        CalcOnePosByMeanShift(extractPtSet[3], centerPoints[4], pos3);

        centerPoints[1] = pos0;
        centerPoints[2] = pos1;
        centerPoints[3] = pos2;
        centerPoints[4] = pos3;

        if(kkj > 10u){
            if ((pos0bak - pos0).norm() < 0.001 && (pos1bak - pos1).norm() < 0.001 && 
                (pos2bak - pos2).norm() < 0.001 && (pos3bak - pos3).norm() < 0.001 ) {
                    break;
            }
        }
    }
    newPositions.resize(centerPointsNum - 1);
    std::copy(centerPoints.begin() + 1, centerPoints.begin() + centerPointsNum, newPositions.begin());
}

void SVMTraceFilter::CalcOnePosByMeanShift( const VectorVec4d& curExtractPtSet, const Vec3d& centerPoint, Vec3d& resultPos )
{
    size_t ptSetNum = curExtractPtSet.size();
    resultPos.setZero();
    double den(0.0);
    double curDist, curPixel, gaussianDist;
    double wet;
    for (size_t i = 0; i < ptSetNum; ++i) {
        curDist = (curExtractPtSet[i].block(0,0,3,1) - centerPoint).norm();
        gaussianDist = std::exp(- std::pow(curDist, 2.0));
        curPixel = curExtractPtSet[i](3) / 100.0;
        wet = curPixel * gaussianDist;
        resultPos += wet * curExtractPtSet[i].block(0,0,3,1);
        den += wet;
    }
    if (den > 0.01) {
        resultPos = resultPos / den;
    } else {
        resultPos = centerPoint;
    }
}

void SVMTraceFilter::DetectCollisionByMeanShift( const Vec3d &curPt, const IVolume& indexImg, const SVolume& origImg, int curEndNodeIndex, 
                                                bool& traceFlag, int& collideIndex )
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    collideIndex = 0;
    traceFlag = true;
    int id1 = NGUtility::Round(curPt(0));
    int id2 = NGUtility::Round(curPt(1));
    int id3 = NGUtility::Round(curPt(2));
    int xMin = std::max(0, id1);
    int xMax = std::min(nx - 1, id1);
    int yMin = std::max(0, id2);
    int yMax = std::min(ny - 1, id2);
    int zMin = std::max(0, id3);
    int zMax = std::min(nz - 1, id3);
    VectorVec5d collideNodeSet;
    Vec5d tmp5d;
    int curIndex;
    double dist;
    for (int ii = xMin; ii <= xMax; ++ii) {
        for (int jj = yMin; jj <= yMax; ++jj) {
            for (int ij = zMin; ij <= zMax; ++ij) {
                curIndex = indexImg(ii,jj,ij);
                if (curIndex != 0 && curIndex != curEndNodeIndex) {
                    tmp5d << double(ii), double(jj), double(ij), 0.0, 0.0;
                    dist = (tmp5d.block(0,0,3,1) - curPt).norm();
                    tmp5d(3) = dist;
                    tmp5d(4) = curIndex;
                    collideNodeSet.push_back(tmp5d);
                }
            }
        }
    }
    if (collideNodeSet.size() > 0) {
        VectorVec5d::iterator it = std::min_element(collideNodeSet.begin(), collideNodeSet.end(), Vec5d_3th_less());
        if ((*it)(3) < 1.5) {
            printf("svm 5\n");
            traceFlag = false;
            collideIndex = int((*it)(4));
        }
    }
}

void SVMTraceFilter::CalcParmOfCurveNodeListSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
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

void SVMTraceFilter::CalcParmOfCurveNodeSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
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
    Vec4d procDistWet = distWet.row(1).array() / ( distWet.row(0).array() + 0.0001);
    for (int i = 1; i <= 4; ++i){
        if (procDistWet(4-i) > 0.5){
            radius = 5.0 - (double)i;
            break;
        }
    }

    if (radius > 2.1){//warning!!
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
bool SVMTraceFilter::GetSVMClassifier( const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, bool isSoma )//, VecXd& w, double &b
{

    MatXd foreSamples,backSamples;
    if (isSoma) GetSVMSampleBySoma(rawDendList, origImg, foreSamples, backSamples);
    else GetSVMSample(rawDendList,origImg,foreSamples,backSamples);
    //time_t end = time(0);
    //printf("GetSVMSample time %d\n",(end-beg));
    if (foreSamples.cols() == 0) {
        return false;
    }
    /*LOG(INFO) << "fore samples:\n";
    for (int i = 0; i < foreSamples.cols(); ++i) {
        LOG(INFO) << foreSamples.col(i).transpose();
    }
    LOG(INFO) << "back samples:\n";
    for (int i = 0; i < foreSamples.cols(); ++i) {
        LOG(INFO) << backSamples.col(i).transpose();
    }*/
    //beg = time(0);
    //SVMClassfier(foreSamples, backSamples, 2, w, b);
    //end = time(0);
    //printf("SVMClassfier time %d\n",(end-beg));
    libsvmclassifier.CreateSVMProblem(foreSamples, backSamples);
    libsvmclassifier.Train();
    //MatXd res;
    //libsvmclassifier.Predict(backSamples, res);
    //std::cout << res << std::endl;
    //std::cout << "foresample " <<foreSamples << std::endl;
    /*VecXd result; result = foreSamples.transpose() * w; 
    result = result.array() + double(b);
    std::cout << "fore " <<result.transpose() << std::endl;
    result = backSamples.transpose() * w; 
    result = result.array() + double(b);
    std::cout << "back " <<result.transpose() << std::endl;*/
    /*FILE *fp = fopen("F:/Linhuimin/back.txt", "w");
    for (int k = 0; k < backSamples.rows(); ++k) {
    for (int m = 0; m < backSamples.cols(); ++m) {
    fprintf(fp, "%lf ", backSamples(k, m));
    }
    fprintf(fp, "\n");
    }fclose(fp);*/
    return true;

}
void SVMTraceFilter::WeighRayValue( const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, VecXd &rayNodeWet )
{
    typedef double spheredist;
    int nxss = int(rayNode.size());
    rayNodeWet.setZero(nxss);

    spheredist x,y,z;
    spheredist segmentDense, segmentWet;//, w1;
    for (int i = 0; i < nxss; ++i){
        x = rayNode[i](0);
        y = rayNode[i](1);
        z = rayNode[i](2);
        segmentDense = segmentWet = 0.0;
        ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
        rayNodeWet(i) = segmentDense / (segmentWet + 0.0001);
    }
}

VecXd SVMTraceFilter::WeighRayValue( const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg )
{
    typedef double spheredist;
    int nxss = int(rayNode.size());
    VecXd rayNodeWet; rayNodeWet.setZero(nxss);

    spheredist x,y,z;
    spheredist segmentDense, segmentWet;//, w1;
    for (int i = 0; i < nxss; ++i){
        x = rayNode[i](0);
        y = rayNode[i](1);
        z = rayNode[i](2);
        segmentDense = segmentWet = 0.0;
        ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
        rayNodeWet(i) = segmentDense / (segmentWet + 0.0001);
    }
    return rayNodeWet;
}

double SVMTraceFilter::WeighRayValue( const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg )
{
    typedef double spheredist;
    spheredist x,y,z;
    spheredist segmentDense, segmentWet;//, w1;
    x = rayNode(0);
    y = rayNode(1);
    z = rayNode(2);
    segmentDense = segmentWet = 0.0;
    ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
    return segmentDense / (segmentWet + 0.0001);
}
void SVMTraceFilter::Principald( const VectorVec3d &dataL, Vec3d &x1 )
{
    int nx = int(dataL.size());
    Vec3d nn;
    nn.setZero();
    double tmp(0.0);
    Vec3d tmpVec;
    for (int i = 0; i < nx - 1; ++i){
        tmpVec = dataL[i + 1] - dataL[i];
        tmp = tmpVec.norm();
        nn += tmp * tmpVec;
    }
    x1 = nn.normalized();
}

void SVMTraceFilter::TileTreeRemoveExtraCurveEnd( const SVolume &origImg, const SVolume &backImg, const VectorMat2i& rawDendConInfo, 
                                                 std::vector<VectorVec5d>& rawDendList )
{
    size_t curvesNum = rawDendList.size();
    //LOG(INFO)<<curvesNum;
    for (size_t i = 0lu; i < curvesNum; ++i) {
        const Mat2i& curConInfo = rawDendConInfo[i];
        //head
        if (0==curConInfo(0,0) && i >0lu) {//do not remove extra start end  of root curve
            VectorVec5d& curDendList = rawDendList[i];
            size_t shortenIndex = RemoveExtraCurveEnd(curDendList, 0, origImg, backImg);
            if (shortenIndex != 0llu) {
                VectorVec5d tmp;
                std::copy(curDendList.begin() + shortenIndex, curDendList.end(), std::back_inserter(tmp));
                curDendList.swap(tmp);
            }
            /*Vec5d tmp5d;
            tmp5d << curDendList[0](0), curDendList[0](1), curDendList[0](2), 0.0, double(i + 1);
            allEndNodes.push_back(tmp5d);
            VectorVec3d tmpDirCurve;
            TraceUtil::GetPartVectorVec3d(curDendList, 0, std::max(int(curDendList.size()) - 9, 1), tmpDirCurve);
            std::reverse(tmpDirCurve.begin(), tmpDirCurve.end());
            Vec3d curDir;
            CalcPtCurveDirection(tmpDirCurve,curDir);
            allEndNodesDirection.push_back(curDir);*/
        }
        //tail
        if (0==curConInfo(0,1)) {
            VectorVec5d& curDendList = rawDendList[i];
            size_t shortenIndex = RemoveExtraCurveEnd(curDendList, 1, origImg, backImg);
            if (shortenIndex < curDendList.size() - 1llu) {
                curDendList.erase(curDendList.begin() + shortenIndex + 1, curDendList.end());
                //VectorVec5d tmp;
                //std::copy(curDendList.begin(), curDendList.begin() + shortenIndex + 1, std::back_inserter(tmp));
                //curDendList.swap(tmp);
            }
            /*Vec5d tmp5d;
            tmp5d << curDendList[shortenIndex](0), curDendList[shortenIndex](1), curDendList[shortenIndex](2), 1.0, double(i + 1);
            allEndNodes.push_back(tmp5d);
            VectorVec3d tmpDirCurve;
            TraceUtil::GetPartVectorVec3d(curDendList, std::max(int(shortenIndex) - 8, 1), int(shortenIndex), tmpDirCurve);
            Vec3d curDir;
            CalcPtCurveDirection(tmpDirCurve,curDir);
            allEndNodesDirection.push_back(curDir);*/
        }
    }
}

void SVMTraceFilter::TileTreeRemoveExtraCurveEndBySVM(const SVolume &origImg, const VectorMat2i& rawDendConInfo, std::vector<VectorVec5d>& rawDendList, VectorVec5d& allEndNodes, VectorVec3d& allEndNodesDirection)
{
    size_t curvesNum = rawDendList.size();
    /*int nxx[2] = { 0, origImg.x() - 1 };
    int nyy[2] = { 0, origImg.y() - 1 };
    int nzz[2] = { 0, origImg.z() - 1 };
    int rpx, rpy, rpz;
    int innerBoundX, outerBoundX, innerBoundY, outerBoundY, innerBoundZ, outerBoundZ;
    innerBoundX = nxx[0] - paramPack->boundaryDistanceThreshold_;
    outerBoundX = nxx[1] + paramPack->boundaryDistanceThreshold_;
    innerBoundY = nyy[0] - paramPack->boundaryDistanceThreshold_;
    outerBoundY = nyy[1] + paramPack->boundaryDistanceThreshold_;
    innerBoundZ = nzz[0] - paramPack->boundaryDistanceThreshold_;
    outerBoundZ = nzz[1] + paramPack->boundaryDistanceThreshold_;*/
    for (size_t i = 0; i < curvesNum; ++i) {
        if (rawDendList[i].empty()) continue;
        const Mat2i& curConInfo = rawDendConInfo[i];
        //head
        if (0 == curConInfo(0, 0) && i > 0lu) {//cannot calculate root node
            VectorVec5d& curDendList = rawDendList[i];
            size_t shortenIndex = RemoveExtraCurveEndBySVM(curDendList, 0, origImg);
            if (shortenIndex != 0llu) {
                /*int curId = int(i + 1);
                for (int x = 0; x < indexImg.x(); ++x) {
                for (int y = 0; y < indexImg.y(); ++y) {
                for (int z = 0; z < indexImg.z(); ++z) {
                if (indexImg(x, y, z) == curId) indexImg(x, y, z) = 0;
                }
                }
                }*/
                //FillCurveLabel(origImg, backImg, indexImg, curDendList, 0, int(i+1));
                VectorVec5d tmp;
                std::copy(curDendList.begin() + shortenIndex, curDendList.end(), std::back_inserter(tmp));
                curDendList.swap(tmp);
                //FillCurveLabel(origImg, backImg, indexImg, curDendList, int(i)+1);
            }
            Vec5d tmp5d;
            tmp5d << curDendList[0](0), curDendList[0](1), curDendList[0](2), 0.0, double(i + 1);
            allEndNodes.push_back(tmp5d);
            VectorVec3d tmpDirCurve;
            NGUtility::GetPartVectorVec3d(curDendList, 0, std::min(int(curDendList.size()) - 1, 9), tmpDirCurve);
            std::reverse(tmpDirCurve.begin(), tmpDirCurve.end());
            Vec3d curDir;
            CalcPtCurveDirection(tmpDirCurve, curDir);
            //curDir = tmpDirCurve.back() - tmpDirCurve[tmpDirCurve.size() - 2lu];
            //curDir.normalize();
            //curDir.setZero();
            allEndNodesDirection.push_back(curDir);
        }
        //tail
        if (0 == curConInfo(0, 1)) {
            VectorVec5d& curDendList = rawDendList[i];
            //Vec3d pt = curDendList.back().block(0, 0, 3, 1);
            //rpx = NGUtility::Round(pt(0));
            //rpy = NGUtility::Round(pt(1));
            //rpz = NGUtility::Round(pt(2));
            //if (AbsDiff(rpx, nxx[0]) <= paramPack->boundaryDistanceThreshold_ ||
            //    AbsDiff(rpx, nxx[1]) <= paramPack->boundaryDistanceThreshold_ ||
            //    AbsDiff(rpy, nyy[0]) <= paramPack->boundaryDistanceThreshold_ ||
            //    AbsDiff(rpy, nyy[1]) <= paramPack->boundaryDistanceThreshold_ ||
            //    AbsDiff(rpz, nzz[0]) <= paramPack->boundaryDistanceThreshold_ ||
            //    AbsDiff(rpz, nzz[1]) <= paramPack->boundaryDistanceThreshold_){//near the boundary
            //    continue;
            //}
            size_t shortenIndex = RemoveExtraCurveEndBySVM(curDendList, 1, origImg);
            if (shortenIndex < curDendList.size() - 1llu) {
                curDendList.erase(curDendList.begin() + shortenIndex + 1, curDendList.end());
            }
            //VectorVec5d tmp;
            //std::copy(curDendList.begin(), curDendList.begin() + shortenIndex + 1, std::back_inserter(tmp));
            //curDendList.swap(tmp);
            Vec5d tmp5d;
            tmp5d << curDendList[shortenIndex](0), curDendList[shortenIndex](1), curDendList[shortenIndex](2), 1.0, double(i + 1);
            allEndNodes.push_back(tmp5d);
            VectorVec3d tmpDirCurve;
            NGUtility::GetPartVectorVec3d(curDendList, std::max(int(shortenIndex) - 8, 1), int(shortenIndex), tmpDirCurve);
            Vec3d curDir;
            CalcPtCurveDirection(tmpDirCurve, curDir);
            //curDir.setZero();
            allEndNodesDirection.push_back(curDir);
        }
    }
}

size_t SVMTraceFilter::RemoveExtraCurveEnd( const VectorVec5d& dendList, int headTailFlag, const SVolume &origImg, const SVolume &backImg )
{
    int shortenIndex = 0;
    int dendListLen = int(dendList.size());
    int rmbegin, rmend;
    if (headTailFlag == 0) {
        shortenIndex = 0;
        rmbegin = 0;
        rmend = std::min(dendListLen,10);
        for (int i = rmbegin; i < rmend; ++i) {
            VectorVec3d curPt ;
            std::vector<double> v1, v2;
            curPt.push_back(dendList[i].block(0,0,3,1));
            WeighRayValue(curPt, origImg, v1);
            WeighRayValue(curPt, backImg, v2);
            if (v1[0] > v2[0] + 1.5 * std::sqrt(v2[0])) {
                shortenIndex = i;
                break;
            }
        }
    }else{
        shortenIndex = dendListLen - 1;
        rmbegin = dendListLen - 1;
        rmend = std::max(-1, dendListLen - 10);
        for (int i = rmbegin; i > rmend; --i) {
            VectorVec3d curPt ;
            std::vector<double> v1, v2;
            curPt.push_back(dendList[i].block(0,0,3,1));
            WeighRayValue(curPt, origImg, v1);
            WeighRayValue(curPt, backImg, v2);
            if (v1[0] > v2[0] + 1.5 * std::sqrt(v2[0])) {
                shortenIndex = size_t(i);
                break;
            }
        }
    }
    return shortenIndex;
}

size_t SVMTraceFilter::RemoveExtraCurveEndBySVM(const VectorVec5d& dendList, int headTailFlag, const SVolume &origImg)
{
    size_t shortenIndex = 0;
    size_t dendListLen = dendList.size();
    MatXd featureMat;
    MatXd svmRes;
    VectorVec3d curPt;
    std::vector<double> pro;

    if (headTailFlag == 0) {
        shortenIndex = 0;
        size_t end;
        if (dendListLen < 7lu) end = dendListLen;
        else end = dendListLen;
        for (size_t i = 0; i < end; ++i) {
            curPt.clear();
            curPt.push_back(dendList[i].block(0, 0, 3, 1));
            GetSVMSampleFeatures(origImg, curPt, featureMat);
            pro.clear();
            libsvmclassifier.Predict(featureMat, svmRes, pro);
            if (svmRes(0, 0) > 0 && pro[0] > 0.3) {
                shortenIndex = i;
                break;
            }
        }
    }
    else{
        shortenIndex = dendListLen - 1;
        size_t beg;
        if (dendListLen < 7lu) beg = 0llu;
        else beg = dendListLen - 6llu;
        for (size_t i = dendListLen - 1; i > beg; --i) {
            curPt.clear();
            curPt.push_back(dendList[i].block(0, 0, 3, 1));
            GetSVMSampleFeatures(origImg, curPt, featureMat);
            pro.clear();
            libsvmclassifier.Predict(featureMat, svmRes, pro);
            if (svmRes(0, 0) > 0 && pro[0] > 0.3) {
                shortenIndex = size_t(i);
                break;
            }
        }
    }
    return shortenIndex;
}

void SVMTraceFilter::SampleCylinderDataSetsExtra( const VectorVec3d curves, int radius, MatXd &matrixS )
{
    VectorVec6d dirVecs;
    VectorVec3d tmpCurve(5);
    Vec3d x1;
    Vec6d tmpdir;
    for (size_t i = 3; i < curves.size() - 3; ++i) {
        std::copy(curves.begin() + i - 2, curves.begin() + i + 3, tmpCurve.begin());
        CalcPtCurveDirection(tmpCurve, x1);
        tmpdir << curves[i], x1;
        dirVecs.push_back(tmpdir);
    }
    std::vector<double> ssStep;
    GenerateSerialVector(0, 2 * M_PI - 0.2, 0.2, ssStep);
    size_t dirNum = dirVecs.size();
    //std::vector<Vectorvec4d> sampleDataSet;
    matrixS.resize(ssStep.size(), dirNum);
    //VectorVec4d curPlanePts;
    Vec3d position, x2, x3;
    for (size_t jj = 0; jj < ssStep.size(); ++jj) {
        double ax1 = radius * std::cos(ssStep[jj]);
        double ay1 = radius * std::sin(ssStep[jj]);
        //curPlanePts.clear();
        for (size_t kk = 0; kk < dirNum; ++kk) {
            position = dirVecs[kk].block(0,0,3,1);
            CalcOrthoBasis(dirVecs[kk].block(3,0,3,1), x2, x3);
            position = position + ax1 * x2 + ay1 * x3;
            position(0) = std::min(std::max(0.0, position(0)), double(origImgPointer->x()));
            position(1) = std::min(std::max(0.0, position(1)), double(origImgPointer->y()));
            position(2) = std::min(std::max(0.0, position(2)), double(origImgPointer->z()));
            double v1 = WeighRayValue(position, *origImgPointer);
            //curPlanePts.push_back(Vec4d(position(0), position(1), position(2), v1));
            matrixS(jj,kk) = v1;
        }
    }
}

double SVMTraceFilter::CalcBifurcationThrev( const std::vector<VectorVec5d> rawDendList )
{
    VectorVec3d curves;
    //guorui turn rawDendList into curves
    size_t sum = 0;
    for (size_t i = 0; i < rawDendList.size(); ++i) sum += rawDendList[i].size();
    curves.resize(sum);
    size_t index = 0;
    for (size_t i = 0; i < rawDendList.size(); ++i) {
        for (size_t j = 0; j < rawDendList[i].size(); ++j) {
            curves[index] = rawDendList[i][j].block(0,0,3,1);
            ++index;
        }
    }
    //
    if (curves.size() < 7) {
        printf("curves size too small, will return self\n ");
        LOG(ERROR) <<"curves size error, error in CalcBifurcationThrev, will return self.";
        return paramPack->bifurcationValue_;
    }
    MatXd matrixS1, matrixS2;
    //LOG(INFO)<<"SampleCylinderDataSetsExtra";
    SampleCylinderDataSetsExtra(curves, 1, matrixS1);
    SampleCylinderDataSetsExtra(curves, 6, matrixS2);
    VecXd matrixS1MaxCol = matrixS1.colwise().maxCoeff();
    //std::cout << "matrixS1MaxCol : " << matrixS1MaxCol.transpose() << std::endl;
    std::vector<double>  origWet;
    VecXd S(matrixS1MaxCol.rows());
    WeighRayValue(curves, *origImgPointer, origWet);
    for (int i = 0; i < matrixS1MaxCol.rows(); ++i) {
        S[i] = std::max(origWet[i + 3], matrixS1MaxCol(i));
    }
    VecXd Q05, Q50;
    MatXdPrctile(matrixS2, 0.05, 0, Q05);
    MatXdPrctile(matrixS2, 0.5, 0, Q50);
    VecXd L1 = (S.array() - Q50.array()) / Q50.array().sqrt();
    VecXd L2 = (Q50.array() - Q05.array()) / Q05.array().sqrt();
    //std::cout << "L1 : " << L1.transpose() << std::endl;
    //std::cout << "L2 : " << L2.transpose() << std::endl;
    VecXd Lx1;
    {
        double tmp = Prctile(L1, 0.15);
        Lx1 = (L1.array() < tmp).select(L1, tmp);
    }
    double thre = 0.5 * (Lx1.mean() + L2.mean());
    return std::min(thre, 2.0);
}

void SVMTraceFilter::MatXdPrctile(const MatXd& mat, double prc, int axis, VecXd& vec)
{
    switch(axis){
        case 0:
            vec.resize(mat.cols());
            for (int i = 0; i < mat.cols(); ++i) {
                VecXd tmpVec = mat.col(i);
                vec(i) = Prctile(tmpVec, prc);
            }
            break;
        case 1:
            vec.resize(mat.rows());
            for (int i = 0; i < mat.rows(); ++i) {
                VecXd tmpVec = mat.row(i).transpose();
                vec(i) = Prctile(tmpVec, prc);
            }
            break;
        default:
            printf("error in matXdPrctile\n");
            LOG(ERROR) << "error in matXdPrctile.";
            break;
    }
}

double SVMTraceFilter::Prctile(VecXd vec, double prc)
{
    double val;
    std::sort(vec.data(), vec.data() + vec.size());
    double ind = vec.size() * prc;
    if ( std::abs(std::floor(ind) - ind - 0.0) < 0.1) {//inte
        int index = std::max(int(ind) - 1, 0);
        if (index == vec.size() - 1) {
            val = vec[index];
        }else{
            val = double(vec[index] + vec[index + 1]) / 2.0;
        }
    }else{//fraction
        int index = int(std::ceil(ind)) - 1;
        val = vec[index];
    }
    return val;
}

void SVMTraceFilter::CubeLabel( const VectorVec3d &curDendrite, VectorVec3i& cubeLabelPts)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    const SVolume &backImg = *backImgPointer;
    int nxx = origImg.x();
    int nyy = origImg.y();
    int nzz = origImg.z();
    cubeLabelPts.clear();

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

        /*printf("chazhi:\n");
        for (size_t i =0; i < resultCurveCopy.size(); ++i) {
            std::cout << resultCurveCopy[i].block(0,0,3,1).transpose() << std::endl;
        }*/

        std::vector<double> three;
        for (VectorVec5d::size_type i = 0; i < resultCurveCopy.size(); ++i){
            three.push_back((double)std::max(3, NGUtility::Round(resultCurveCopy[i](3) + 1.0)));//2016-5-6
        }
        int id1(0), id2(0), id3(0);
        int xMin, xMax, yMin, yMax, zMin, zMax;
        for (int ik = 1; ik < 2*int(rawSomaCurve.size())-2; ++ik){
            id1 = NGUtility::Round(resultCurveCopy[ik](0));
            id2 = NGUtility::Round(resultCurveCopy[ik](1));
            id3 = NGUtility::Round(resultCurveCopy[ik](2));

            xMin = std::max(0, id1 - std::max(0, (int)three[ik]));
            xMax = std::min(nxx - 1, id1 + std::max(0, (int)three[ik]));
            yMin = std::max(0, id2 - std::max(0, (int)three[ik]));
            yMax = std::min(nyy - 1, id2 + std::max(0, (int)three[ik]));
            zMin = std::max(0, id3 - std::max(0, (int)three[ik]));
            zMax = std::min(nzz - 1, id3 + std::max(0, (int)three[ik]));

            double threBack(0);

            for (int ii = xMin; ii <= xMax; ++ii){
                for (int jj = yMin; jj <= yMax; ++jj){
                    for (int kk = zMin; kk <= zMax; ++kk){
                        //2015-6-8
                        //threBack = double(backImg(ii, jj, kk)) + bifBinThrev_ * std::sqrt(double(backImg(ii, jj, kk)));//2016-3-18
                        threBack = double(backImg(ii, jj, kk)) + paramPack->bifurcationValue_ * std::sqrt(double(backImg(ii, jj, kk)));//2016-3-18
                        if (traceLabelMatrix_(ii, jj, kk) == 0 && double(origImg(ii, jj, kk)) > threBack){
                            indexImg_(ii, jj, kk) = globalID_;
                            traceLabelMatrix_(ii, jj, kk) = 1;
                            cubeLabelPts.push_back(Vec3i(ii,jj,kk));
                        }
                    }
                }
            }//for
        }//for
    }
}

void SVMTraceFilter::DetectBifurcation(const VectorVec3i& cubeLabelPts, 
                                       VectorVec3d& bifurcationDirections, std::vector<VectorVec3d>& bifurcationClustre )
{
    //traceLabelMatrix is traceLabelMatrixCp;
    CVolume traceLabelMatrixCp;
    traceLabelMatrixCp.QuickCopy(traceLabelMatrix_);
    std::vector<VectorVec3i> growPointSet;
    ExtractCubeBoundary(traceLabelMatrixCp, cubeLabelPts, growPointSet);
    {
        /*char fileline[256];
        int id = 0;
        for (auto &it : growPointSet) {
        sprintf(fileline, "F:/bug/DetectBifurcation%d.swc", ++id);
        NGUtility::WriteVectorVec3d(fileline, "dot", it);
        }*/

    }
    
    bool isContinue = true;
    for(size_t i = 0; i < growPointSet.size(); ++i){
        if(growPointSet[i].empty())
            isContinue = false;
    }
    if(isContinue){
        std::vector<VectorVec3i> clustrePtSet;
        std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> > clustreLabel;
        ClusterGrowPtSet(growPointSet, traceLabelMatrixCp, clustrePtSet, clustreLabel);
        size_t flagNum = 0;
        for(size_t i =0; i < clustreLabel.size(); ++i){
            flagNum = std::max(flagNum, clustreLabel[i].size());
        }
        if(flagNum > 10000){//2015-6-8
            bifurcationDirections.clear();
            bifurcationClustre.clear();
        } else{
            //LOG(INFO)<<flagNum;
            std::vector<std::vector<int> > subRegionCon;
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

            //2015-6-8
            for(int i = 0; i < pathMatrix.rows(); ++i){
                if(redundancyLabel[i]==1){
                    VectorVec3d curDendrite;
                    NGUtility::GetReversePartVectorVec3d(connectPoints[i], 0, 5, curDendrite);//2016-5-5
                    VectorVec3d tmp;
                    //tmp.push_back(soma);
                    //std::copy(curDendrite.begin() + 2, curDendrite.end(), back_inserter(tmp));
                    std::copy(curDendrite.begin(), curDendrite.end(), back_inserter(tmp));
                    bifurcationClustre.push_back(VectorVec3d());
                    bifurcationClustre.back().swap(tmp);

                    Vec3d curDir;
                    //VectorVec3d tmp;
                    //tmp.clear();
                    //std::copy(curDendrite.begin(), curDendrite.end(), back_inserter(tmp));//2015-6-8
                    //CurveDirection(tmp, curDir);
                    CurveDirection(curDendrite, curDir);
                    bifurcationDirections.push_back(curDir);
                }
            }
        }
    }else{
        bifurcationClustre.clear();
        bifurcationDirections.clear();
    }
}

void SVMTraceFilter::ExtractCubeBoundary(CVolume& traceLabelMatrix, const VectorVec3i& cubeLabelPts,
                                         std::vector<VectorVec3i>& growPointSet)
{
    //initial
    const SVolume &origImg = *origImgPointer;
    //const SVolume &backImg = *backImgPointer;
    int nxx = origImg.x() - 1;
    int nyy = origImg.y() - 1;
    int nzz = origImg.z() - 1;//warning!!
    VectorVec3i boundaryLabel;
    int pointRadius = 2;
    int xMin, xMax, yMin, yMax, zMin, zMax, ss;
    for (size_t i = 0; i < cubeLabelPts.size(); ++i) {
        const Vec3i& curPt = cubeLabelPts[i];
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPt(0), curPt(1), curPt(2), pointRadius, pointRadius, pointRadius, 0, nxx, 0, nyy, 0, nzz );
        GetAreaSum(traceLabelMatrix, xMin, xMax, yMin, yMax, zMin, zMax, ss);
        if (double(ss) < 0.5 * double((xMax - xMin + 1) * (yMax - yMin + 1)*(zMax - zMin + 1))) {
            boundaryLabel.push_back(curPt);
        }
    }
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
        RegionInflationModifyV2(boundaryLabelCp, traceLabelMatrix, 1, curPoints);//2016-3-18
        if(!curPoints.empty()){
            growPointSet[i] = curPoints;
            boundaryLabelCp.swap(curPoints);
        } else{
            break;
        }
    }
}
void SVMTraceFilter::ClusterGrowPtSet( const std::vector<VectorVec3i>& growPointSet, CVolume& growLabelMatrix, 
                                      std::vector<VectorVec3i>& clustrePtSet, std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> >& clustreLabel )
{
    //initialize
    size_t nxx = growPointSet.size();
    clustreLabel.clear();
    clustreLabel.resize(nxx);
    clustrePtSet = growPointSet;
    /*
    * cluster each level
    */
    VectorVec3i inflatPtSet;
    VectorVec3d curConnectPt;
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
        std::vector<SVMTraceFilter::CLUSTRELABEL> curPtSetFeature(subClustrePtSetNumLength);
        for(size_t i = 0; i < curPtSetFeature.size(); ++i)
            curPtSetFeature[i].index = subClustrePtSetNum[i];
        for(size_t i = 0; i < subClustrePtSetNumLength; ++i){
            curConnectPt.clear();
            /*std::copy(subClustrePtSet.begin() + subClustrePtSetCumNum[i],
            subClustrePtSet.begin() + subClustrePtSetCumNum[i + 1],
            std::back_inserter(curConnectPt));*/
            {
                Vec3d tmpVec3d;
                for (size_t k = subClustrePtSetCumNum[i]; k < subClustrePtSetCumNum[i + 1]; ++k) {
                    NGUtility::Vec3i2Vec3d(subClustrePtSet[k], tmpVec3d);
                    curConnectPt.push_back(tmpVec3d);
                }
            }
            //the same as isempty function
            if(subClustrePtSetCumNum[i + 1] - 1 > subClustrePtSetCumNum[i]){
                curCenter.setZero();
                // get center//2016-4-25
                std::vector<double> curConnectPtW;
                WeighRayValue(curConnectPt, *origImgPointer, curConnectPtW);
                Vec3d maxCurConnectPt = curConnectPt[std::distance(curConnectPtW.begin(), std::max_element(curConnectPtW.begin(), curConnectPtW.end()))];
                VectorVec4d tmpVectorVec4d(curConnectPtW.size());
                for (size_t k = 0; k < curConnectPtW.size(); ++k) {
                    tmpVectorVec4d[k] << curConnectPt[k](0), curConnectPt[k](1), curConnectPt[k](2), curConnectPtW[k];
                }
                CalcOnePosByMeanShift(tmpVectorVec4d, maxCurConnectPt, curCenter);
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

void SVMTraceFilter::RegionInflationModifyV2( const VectorVec3i &curPoints, CVolume &growLabelMatrix, double, VectorVec3i &growPoints )
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
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VectorVec3i tmpGrowPts;
    //LOG(INFO) << "RegionInflationModifyV2";
    for(size_t i = 0; i < nxx; ++i){

        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, curPoints[i](0), curPoints[i](1), curPoints[i](2),
            2,2,2, 0, xBdy, 0, yBdy, 0, zBdy);//2015-2-8
        tmpGrowPts.clear();
        tmpGrowPts.reserve(75);
        for(int ii = xMin; ii <= xMax; ++ii){
            for(int jj = yMin; jj <= yMax; ++jj){
                for(int ij = zMin; ij <= zMax; ++ij){
                    //2016-4-28
                    //flag1 = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + bifBinThrev_ * std::sqrt(double(backImg(ii,jj,ij)));
                    flag1 = double(origImg(ii,jj,ij)) > double(backImg(ii,jj,ij)) + paramPack->bifurcationValue_ * std::sqrt(double(backImg(ii,jj,ij)));
                    if(0 == growLabelMatrix(ii,jj,ij) && flag1){
                        tmpGrowPts.push_back(Vec3i(ii,jj,ij));
                        //LOG(INFO) << origImg(ii,jj,ij) << " " << backImg(ii,jj,ij) << "true";
                    }
                    //else if(0 == growLabelMatrix(ii,jj,ij)) LOG(INFO) << origImg(ii,jj,ij) << " " << backImg(ii,jj,ij) << "false"; 
                }
            }
        }//for
        
        if(tmpGrowPts.size() > 0){
            std::copy(tmpGrowPts.begin(), tmpGrowPts.end(), std::back_inserter(growPoints));
            for(size_t j = 0; j < tmpGrowPts.size(); ++j)
                growLabelMatrix(tmpGrowPts[j](0), tmpGrowPts[j](1), tmpGrowPts[j](2)) = 1;
        }
    }
}
void SVMTraceFilter::CurveDirection(const VectorVec3d &curDendrite, Vec3d &curDir)
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
void SVMTraceFilter::ConnectGrowClustreFromSoma(const std::vector<std::vector<int> > &subRegionCon,
                                                const std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> > &clustreLabel,
                                                std::vector<VectorVec4d> &connectPoints,
                                                std::vector<VectorVec4d> &shortConnectPoints)
{
    size_t nxx = clustreLabel.size();
    /*
    * clustreLabel is consist of index, three dimension coordinate of center,
    * three dimension coordinate of standard deviation, of all 7 elements.
    */
    std::vector<SVMTraceFilter::CLUSTRELABEL> terminalPtFeature = clustreLabel.back();
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
            curRecursiveIndex = subRegionCon[nxx - i - 2][curRecursiveIndex] - 1lu;
            connectPoints[j][i +1].block(0,0,3,1) = clustreLabel[nxx - i -2][curRecursiveIndex].center;
            const Vec3d &maxCoordinate = clustreLabel[nxx - i -2][curRecursiveIndex].standardDeviation;
            if(maxCoordinate.maxCoeff() < 15.1)//2015-6-8 3.1
                connectPoints[j][i +1](3) = 1.0;
            else
                connectPoints[j][i +1](3) = 0.0;
        }
    }
    shortConnectPoints.clear();//2015-6-8

}
void SVMTraceFilter::DeleteRedundancyPathV2(const MatXi &pathMatrix, std::vector<int> &redundancyLabel)
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
void SVMTraceFilter::ClustersExtractionPoints(const VectorVec3i &inflatPtSet, CVolume &growLabelMatrix,
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
void SVMTraceFilter::RegionConnections(const std::vector<VectorVec3i> &clustrePtSet,
                                       const std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> > &clustreLabel,
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
    //LOG(INFO)<<clustreNum;
    for(size_t i = 0; i < clustreNum; ++i){
        const VectorVec3i& curClustrePtSet = clustrePtSet[i];
        const std::vector<SVMTraceFilter::CLUSTRELABEL>& curClustreFeature = clustreLabel[i];
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
    //LOG(INFO)<<clustreNum;
    for(size_t i = 1; i < clustreNum; ++i){
        connectionFlag = 0;
        const VectorVec3i& curClustrePtSet = clustrePtSet[i];
        const std::vector<SVMTraceFilter::CLUSTRELABEL>& curClustreFeature = clustreLabel[i];
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
            // region id has plus one
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
void SVMTraceFilter::RegionConnectionsSub(const VectorVec3i &curConPtSet, Volume<int> &growLabelMatrix,
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
void SVMTraceFilter::ClustersExtraction(CVolume& growLabelMatrix, const VectorVec3i& inflatPtSet,
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
void SVMTraceFilter::ExtractSubRegionOP(const Vec3i &currentPoint, VectorVec3i &extractedPtSet, CVolume &growLabelMatrix)
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
void SVMTraceFilter::ExtractSubRegionAboveThrevModify(const VectorVec3i &currentPtSet, VectorVec3i &extractedPtSet,
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

        //if (id1 > 1 && id1 < nx-2 &&id2 > 1 && id2 < ny-2 && id3 > 1 && id3 < nz-2){
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, id1, id2, id3, 1,1,1, 0, nx - 1, 0, ny -1, 0, nz -1);

        //
        t_Sum = 0;
        for (i = xMin; i <= xMax; ++i)
            for (ij = zMin; ij <= zMax; ++ij)
                for (j = yMin; j <= yMax; ++j){
                    t_Sum += indexImg(i , j , ij );
                }

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
        //}
    }//for
}
void SVMTraceFilter::ConnectPathInTree(const std::vector<std::vector<int> > &subRegionCon, MatXi &pathMatrix)
{
    std::vector<std::vector<int> > subRegionConCp(subRegionCon);
    int nxx = int(subRegionConCp.size());
    int pathNum = int(subRegionConCp.back().size());
    std::reverse(subRegionConCp.begin(), subRegionConCp.end());
    pathMatrix = MatXi(pathNum, nxx +1);
    pathMatrix.setZero();
    for(int i = 0; i < pathNum; ++i){
        pathMatrix(i,0) = i+1;//warning!
        for(size_t j = 1; j < nxx + 1; ++j){
            const std::vector<int>& curData = subRegionConCp[j-1];
            pathMatrix(i,j) = curData[pathMatrix(i,j-1)-1];//C++
        }
    }
    MatXi pathMatrixCp = pathMatrix.rowwise().reverse();//resolve eigen alias
    pathMatrix = pathMatrixCp;
}

void SVMTraceFilter::CalcParmOfCurveNodeList(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
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

void SVMTraceFilter::CalcParmOfCurveNode(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
                                      const Vec3d &curNode, double &radius, double &wet)
{
    int vx = backImg.x();
    int vy = backImg.y();
    int vz = backImg.z();

    int xMin, xMax, yMin, yMax, zMin, zMax;

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

    //2015-8-13
    Vec4d level;
    level << 1.0, 2.0, 3.0, 4.0;

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
    Vec4d procDistWet = distWet.row(1).array() / ( distWet.row(0).array() + 0.0001);
    for (int i = 1; i <= 4; ++i){
        if (procDistWet(4-i) > 0.5){
            radius = 5.0 - (double)i;
            break;
        }
    }

    if (radius > 2.1){
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

void SVMTraceFilter::Train( std::vector<VectorVec5d>& arg)
{
    //VecXd  w(16); double b;//double tmpw;
    if(!GetSVMClassifier( arg, *origImgPointer)){//, paramPack->w, paramPack->b
        printf("SVM failed and skip.\n");
        LOG(ERROR) <<"SVM failed and skip .";
        return;
    }
}

ProcStatPointer SVMTraceFilter::GPSTreeUpdate()
{
    printf("begin SVM.\n");
    if (!m_Input || !m_Back || !m_Bin || !paramPack){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if (m_Input->GetProcessObject()){//|| !m_Soma
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }

    origImgPointer =std::dynamic_pointer_cast<const Volume<unsigned short>>(m_Input);
    backImgPointer =  std::dynamic_pointer_cast<const Volume<unsigned short>>(m_Back);
    binImgPointer = std::dynamic_pointer_cast<const Volume<NGCHAR>>(m_Bin);

    if (!origImgPointer || !backImgPointer || !binImgPointer) {
        printf("error occurred in %s\n", className_.c_str());
        return false;
    }
    //2016-3-22
    VectorVec5d allEndNodes;
    VectorVec3d allEndNodesDirection;

    TileTreeRemoveExtraCurveEnd(*origImgPointer, *backImgPointer, rawDendConInfo_, rawDendList_);

    GetSVMClassifier(rawDendList_, *origImgPointer, true);//soma
    //back to remove the extra curves.
    allEndNodes.clear();
    allEndNodesDirection.clear();
    TileTreeRemoveExtraCurveEndBySVM(*origImgPointer, rawDendConInfo_, rawDendList_, allEndNodes, allEndNodesDirection);
    //allEndNodes.clear();//TODO: test

    //trace
    printf("allEndNodes: %lu %lu\n", allEndNodes.size(), allEndNodesDirection.size());
    for (size_t ii = 0; ii < allEndNodes.size(); ++ii) {
        printf("                                   \r");
        printf("svm: %lu of %lu", ii, allEndNodes.size());
        int curEndNodeIndex = int(allEndNodes[ii](4));
        VectorVec5d resultCurve;
        int collideIndex = 0;
        //resultCurve.push_back(MakeVec5d(allEndNodes[ii].block(0,0,3,1) + 3.0 * allEndNodesDirection[ii]));
        MeanShiftGrowTraceTotal(allEndNodes[ii].block(0, 0, 3, 1), allEndNodesDirection[ii], *origImgPointer, *backImgPointer, curEndNodeIndex,
            resultCurve, collideIndex);
        if (!resultCurve.empty()) {
            if (int(allEndNodes[ii](3)) == 0) {
                VectorVec5d curCurve;// = rawDendList[curEndNodeIndex];
                std::copy(resultCurve.rbegin(), resultCurve.rend(), std::back_inserter(curCurve));
                std::copy(rawDendList_[curEndNodeIndex - 1].begin(), rawDendList_[curEndNodeIndex - 1].end(), std::back_inserter(curCurve));
                rawDendList_[curEndNodeIndex - 1].swap(curCurve);
                rawDendConInfo_[curEndNodeIndex - 1](0, 0) = collideIndex == std::numeric_limits<int>::max() ? collideIndex : 0;
            }
            if (int(allEndNodes[ii](3)) == 1) {
                VectorVec5d curCurve;// = rawDendList[curEndNodeIndex];

                std::copy(rawDendList_[curEndNodeIndex - 1].begin(), rawDendList_[curEndNodeIndex - 1].end(), std::back_inserter(curCurve));
                std::copy(resultCurve.begin(), resultCurve.end(), std::back_inserter(curCurve));
                rawDendList_[curEndNodeIndex - 1].swap(curCurve);
                rawDendConInfo_[curEndNodeIndex - 1](0, 1) = collideIndex == std::numeric_limits<int>::max() ? collideIndex : 0;
            }
        }
        else if (collideIndex != 0){//
            if (int(allEndNodes[ii](3)) == 0) {
                rawDendConInfo_[curEndNodeIndex - 1](0, 0) = collideIndex == std::numeric_limits<int>::max() ? collideIndex : 0;
            }
            if (int(allEndNodes[ii](3)) == 1) {
                rawDendConInfo_[curEndNodeIndex - 1](0, 1) = collideIndex == std::numeric_limits<int>::max() ? collideIndex : 0;
            }
        }
    }
    VectorVec3d newSeed;
    SelectSeedForTrace(*binImgPointer, *origImgPointer, indexImg_, newSeed);
    VectorVec5d seedEndNodes;
    VectorVec3d seedEndNodesDirection;
    int id = int(rawDendList_.size());
    for (auto &it : newSeed) {
        if (it.minCoeff() < 3 || it(0) > indexImg_.x() - 3 || it(1) > indexImg_.y() - 3 || it(2) > indexImg_.z() - 3) continue;
        int maxVal = 0;
        GetAreaMaxValue(indexImg_, NGUtility::Round(it(0)) - 1, NGUtility::Round(it(0)) + 1, NGUtility::Round(it(1)) - 1,
            NGUtility::Round(it(1)) + 1, NGUtility::Round(it(2)) - 1, NGUtility::Round(it(2)) + 1, maxVal);
        if (maxVal != 0) continue;
        ++id;//ii+1
        Vec5d tmp; tmp << it(0), it(1), it(2), 0, id;
        seedEndNodes.push_back(tmp);
        tmp << it(0), it(1), it(2), 1, id;
        seedEndNodes.push_back(tmp);
        Vec3d initDirection;
        CalcInitDirectionOnTraceSeed(*origImgPointer, *backImgPointer, it, 12.0, initDirection);
        seedEndNodesDirection.push_back(initDirection);
        seedEndNodesDirection.push_back(-initDirection);
    }
    printf(" %lu svm seed\n", seedEndNodes.size() / 2lu);
    //seedEndNodes.clear();

    rawDendConInfo_.resize(rawDendList_.size() + seedEndNodes.size() / 2llu);
    std::for_each(rawDendConInfo_.begin() + rawDendList_.size(), rawDendConInfo_.end(), [](Mat2i & arg){arg.setZero(); });
    rawDendList_.resize(rawDendList_.size() + seedEndNodes.size() / 2llu);
    for (size_t ii = 0; ii < seedEndNodes.size(); ii += 2) {
        printf("                                \r");
        printf("svm: %lu of %lu", ii / 2, seedEndNodes.size() / 2lu);
        int maxVal = 0;
        GetAreaMaxValue(indexImg_, NGUtility::Round(seedEndNodes[ii](0)) - 1, NGUtility::Round(seedEndNodes[ii](0)) + 1,
            NGUtility::Round(seedEndNodes[ii](1)) - 1, NGUtility::Round(seedEndNodes[ii](1)) + 1, 
            NGUtility::Round(seedEndNodes[ii](2)) - 1, NGUtility::Round(seedEndNodes[ii](2)) + 1, maxVal);
        if (maxVal != 0) continue;
        //if (indexImg(NGUtility::Round(seedEndNodes[ii](0)), NGUtility::Round(seedEndNodes[ii](1)), NGUtility::Round(seedEndNodes[ii](2))) != 0) continue;
        int curEndNodeIndex1 = int(seedEndNodes[ii](4));
        int curEndNodeIndex2 = int(seedEndNodes[ii + 1](4));
        VectorVec5d resultCurve1, resultCurve2;
        int collideIndex1 = 0, collideIndex2 = 0;
        MeanShiftGrowTraceTotal(seedEndNodes[ii].block(0, 0, 3, 1), seedEndNodesDirection[ii], *origImgPointer, *backImgPointer, curEndNodeIndex1,
            resultCurve1, collideIndex1);
        MeanShiftGrowTraceTotal(seedEndNodes[ii + 1].block(0, 0, 3, 1), seedEndNodesDirection[ii + 1], *origImgPointer, *backImgPointer, curEndNodeIndex2,
            resultCurve2, collideIndex2);
        VectorVec5d curCurve;
        if (!resultCurve1.empty()) {
            if (collideIndex1 != 0) {
                Vec5d head = FindNearestNode(resultCurve1.back().block(0, 0, 3, 1), rawDendList_[collideIndex1 - 1]);
                //curCurve.push_back(head);
            }
            std::copy(resultCurve1.rbegin(), resultCurve1.rend(), std::back_inserter(curCurve));

        }
        if (!resultCurve2.empty()) {//2017-5-12
            std::copy(resultCurve2.begin(), resultCurve2.end(), std::back_inserter(curCurve));
            if (collideIndex2 != 0){
                Vec5d tail = FindNearestNode(resultCurve2.back().block(0, 0, 3, 1), rawDendList_[collideIndex2 - 1]);
                //curCurve.push_back(tail);
            }
        }
        if (curCurve.size() > 10lu || (curCurve.size() > 3lu && (collideIndex1 != 0 || collideIndex2 != 0))) {

            if (rawDendList_.size() > curEndNodeIndex1 && rawDendList_[curEndNodeIndex1 - 1].size() > curCurve.size()){
                printf("nima\n"); //system("pause"); 
            }
            printf("new line: %lu\n", curCurve.size());
            rawDendList_[curEndNodeIndex1 - 1].swap(curCurve);
            rawDendConInfo_[curEndNodeIndex1 - 1](0, 0) = collideIndex1;
            rawDendConInfo_[curEndNodeIndex1 - 1](0, 1) = collideIndex2;
        }
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

Vec5d SVMTraceFilter::FindNearestNode(const Vec3d& pt, const VectorVec5d& curve)
{
    double minDist = 10000.0;
    Vec5d tmp;
    for (auto &it : curve) {
        if ((pt - it.block(0, 0, 3, 1)).norm() < minDist){
            minDist = (pt - it.block(0, 0, 3, 1)).norm();
            tmp = it;
        }
    }
    return tmp;
}