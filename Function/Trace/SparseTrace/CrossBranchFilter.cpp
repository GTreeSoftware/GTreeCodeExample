
#include <cmath>
#include <set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <Eigen/Core>
#include <iostream>
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <glog/logging.h>
#include "CrossBranchFilter.h"
#include "../../volumealgo.h"
#include "../../contourutil.h"
#include "../traceutil.h"
#include "../../NGUtility.h"

CrossBranchFilter::CrossBranchFilter()
{
    className_ = std::string("CrossFiberFilter");
}


CrossBranchFilter::~CrossBranchFilter()
{
}

ProcStatPointer CrossBranchFilter::Update()
{
    NG_DYNAMIC_CAST(const SVolume, origImgPointer, paramPack->OrigImage);
    if (!origImgPointer || !tracedCurve_ || !bifurcationClustre_ || !bifurcationDirections_) {
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    /*
    * remove cross curves
    */
    std::vector<VectorVec3d> bifurcationClustreCp;
    VectorVec3d bifurcationDirectionsCp;
    VectorVec3i boundaryLabelCp;
    std::set<size_t> removeInd;
    RemoveCrossBifurcation(*tracedCurve_, *bifurcationClustre_, removeInd);
    for (size_t i = 0; i < (*bifurcationClustre_).size(); ++i) {
        Vec3d tmpVec = (*bifurcationClustre_)[i].front().block(0, 0, 3, 1) - *soma_;
        if (tmpVec.norm() < 5.0 || (tmpVec.norm() < 10.0 &&
            (tmpVec.normalized()).dot(*initDir_) < -0.5)) continue;
        if (removeInd.count(i) != 0) continue;
        bifurcationClustreCp.push_back((*bifurcationClustre_)[i]);
        bifurcationDirectionsCp.push_back((*bifurcationDirections_)[i]);
        if(boundaryLabel_) boundaryLabelCp.push_back((*boundaryLabel_)[i]);
    }
    (*bifurcationClustre_).swap(bifurcationClustreCp);
    (*bifurcationDirections_).swap(bifurcationDirectionsCp);
    if (boundaryLabel_) boundaryLabelCp.swap((*boundaryLabel_));

    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

/*
* C++11
*/
void CrossBranchFilter::RemoveCrossBifurcation(const std::vector<VectorVec3d> &tmpDendrites,
    const std::vector<VectorVec3d>&bifurClustre, std::set<size_t> &removeCrossPairIndList)
{
    if (bifurClustre.size() < 2u) return;
    removeCrossPairIndList.clear();
    /*find the cross candidate point pairs*/
    std::vector<std::vector<VectorVec3d> >crossCandidataCurvesList;
    std::vector<VectorVec3d> crossCandidateList;
    std::vector<double> crossCandidateListThirdPointValue;
    std::vector<double> crossCandidateListFourthPointValue;
    std::vector<std::pair<size_t, size_t> > crossPairInd;
    size_t curInd1(1000000), curInd2(10000000), ptInd1(1000000), ptInd2(1000000);
    for (auto it1 = bifurClustre.begin(); it1 != bifurClustre.end(); ++it1) {
        for (auto it2 = it1 + 1; it2 != bifurClustre.end(); ++it2) {
            if ((it1->front() - it2->front()).norm() < paramPack->crossDistanceThrev_) {
                /*find the parent curves correspond to cross candidate point pairs*/
                FindNearestParentCurve(it1->front(), tmpDendrites, curInd1, ptInd1);
                FindNearestParentCurve(it2->front(), tmpDendrites, curInd2, ptInd2);
                if (curInd1 == curInd2 && tmpDendrites[curInd2].size() > 5llu && int(ptInd1 + ptInd2) / 2 - 2 > 0
                    && int(ptInd1 + ptInd2) / 2 + 2  < tmpDendrites[curInd1].size()) {
                    crossPairInd.push_back(std::make_pair(std::distance(bifurClustre.begin(), it1), std::distance(bifurClustre.begin(), it2)));
                    std::vector<double> rayValue;
                    VectorVec3d curCrossPtList;
                    std::vector<VectorVec3d> curCrossCandidataCurvesList;
                    const VectorVec3d &parentCurve = tmpDendrites[curInd2];
                    //cross curve 1
                    VectorVec3d childCurve1;
                    for (size_t i = 0; i < it1->size(); ++i) {
                        if (((*it1)[i] - it1->front()).norm() < paramPack->crossDistanceThrev_) {
                            childCurve1.push_back((*it1)[i]);
                        }
                    }
                    WeighRayValue(childCurve1, *origImgPointer, rayValue);
                    auto maxIt = std::max_element(rayValue.begin(), rayValue.end());
                    curCrossPtList.push_back(childCurve1[std::distance(rayValue.begin(), maxIt)]);
                    curCrossCandidataCurvesList.push_back(childCurve1);
                    //cross curve 2
                    VectorVec3d childCurve2;
                    for (size_t i = 0; i < it2->size(); ++i) {
                        if (((*it2)[i] - it2->front()).norm() < paramPack->crossDistanceThrev_) {
                            childCurve2.push_back((*it2)[i]);
                        }
                    }
                    WeighRayValue(childCurve2, *origImgPointer, rayValue);
                    maxIt = std::max_element(rayValue.begin(), rayValue.end());
                    curCrossPtList.push_back(childCurve2[std::distance(rayValue.begin(), maxIt)]);
                    curCrossCandidataCurvesList.push_back(childCurve2);
                    //parent curve head
                    VectorVec3d pCurvePre;
                    size_t preCurveEnd = (ptInd2 + ptInd1) / 2llu - 2llu;
                    size_t preCurveHead = 0;
                    for (size_t i = preCurveEnd; i > 0; --i) {//donot need head
                        if ((parentCurve[i] - parentCurve[preCurveEnd]).norm() > paramPack->crossDistanceThrev_) {
                            preCurveHead = i;
                            break;
                        }
                    }
                    std::copy(parentCurve.begin() + preCurveHead, parentCurve.begin() + preCurveEnd, std::back_inserter(pCurvePre));
                    std::reverse(pCurvePre.begin(), pCurvePre.end());
                    WeighRayValue(pCurvePre, *origImgPointer, rayValue);
                    maxIt = std::max_element(rayValue.begin(), rayValue.end());
                    curCrossPtList.push_back(pCurvePre[std::distance(rayValue.begin(), maxIt)]);
                    crossCandidateListThirdPointValue.push_back((*maxIt));
                    curCrossCandidataCurvesList.push_back(pCurvePre);
                    //parent curve tail
                    VectorVec3d pCurvePost;
                    size_t postCurveHead = (ptInd2 + ptInd1) / 2u + 2u;
                    size_t postCurveEnd = parentCurve.size();
                    for (size_t i = postCurveHead; i < postCurveEnd; ++i) {//donot need head
                        if ((parentCurve[i] - parentCurve[postCurveHead]).norm() > paramPack->crossDistanceThrev_) {
                            postCurveEnd = i;
                            break;
                        }
                    }
                    std::copy(parentCurve.begin() + postCurveHead, parentCurve.begin() + postCurveEnd, std::back_inserter(pCurvePost));
                    WeighRayValue(pCurvePost, *origImgPointer, rayValue);
                    maxIt = std::max_element(rayValue.begin(), rayValue.end());
                    curCrossPtList.push_back(pCurvePost[std::distance(rayValue.begin(), maxIt)]);
                    crossCandidateList.push_back(curCrossPtList);
                    crossCandidateListFourthPointValue.push_back(*maxIt);
                    curCrossCandidataCurvesList.push_back(pCurvePost);
                    crossCandidataCurvesList.push_back(std::vector<VectorVec3d>());
                    crossCandidataCurvesList.back().swap(curCrossCandidataCurvesList);
                }
            }
        }
    }
    //filter cross point
    CVolume signalMatrix;
    signalMatrix.SetSize(origImgPointer->x(), origImgPointer->y(), origImgPointer->z());
    for (size_t i = 0; i < crossCandidateList.size(); ++i) {
        VectorVec3d &curCrossPtList = crossCandidateList[i];//3rd point
        //start diffuse point;
        Vec3i startPoint;
        double startThrev;
        /*Vec3i curve to Intersect*/
        VectorVec3i iCurvePre, iCurvePost, iParentCurve;
        std::transform(crossCandidataCurvesList[i][0].begin(), crossCandidataCurvesList[i][0].end(),
            std::back_inserter(iCurvePre), [](const Vec3d& arg){return Vec3i(NGUtility::Round(arg(0)), NGUtility::Round(arg(1)), NGUtility::Round(arg(2))); });
        std::transform(crossCandidataCurvesList[i][1].begin(), crossCandidataCurvesList[i][1].end(),
            std::back_inserter(iCurvePost), [](const Vec3d& arg){return Vec3i(NGUtility::Round(arg(0)), NGUtility::Round(arg(1)), NGUtility::Round(arg(2))); });
        if (crossCandidateListThirdPointValue[i] > crossCandidateListFourthPointValue[i]) {
            startPoint = Vec3i(NGUtility::Round(curCrossPtList[2](0)), NGUtility::Round(curCrossPtList[2](1)), NGUtility::Round(curCrossPtList[2](2)));//3th
            startThrev = crossCandidateListThirdPointValue[i];
            std::transform(crossCandidataCurvesList[i][3].begin(), crossCandidataCurvesList[i][3].end(),
                std::back_inserter(iParentCurve), [](const Vec3d& arg){return Vec3i(NGUtility::Round(arg(0)), NGUtility::Round(arg(1)), NGUtility::Round(arg(2))); });
        }
        else {
            startPoint = Vec3i(NGUtility::Round(curCrossPtList[3](0)), NGUtility::Round(curCrossPtList[3](1)), NGUtility::Round(curCrossPtList[3](2)));//4th
            startThrev = crossCandidateListFourthPointValue[i];
            std::transform(crossCandidataCurvesList[i][2].begin(), crossCandidataCurvesList[i][2].end(),
                std::back_inserter(iParentCurve), [](const Vec3d& arg){return Vec3i(NGUtility::Round(arg(0)), NGUtility::Round(arg(1)), NGUtility::Round(arg(2))); });
        }
        //test
        /*WriteVectorVec3d("F:/qiao/iCurvePre1.swc", "dot", crossCandidataCurvesList[i][0]);
        WriteVectorVec3d("F:/qiao/iCurvePost1.swc", "dot", crossCandidataCurvesList[i][1]);
        WriteVectorVec3d("F:/qiao/iParentCurve1.swc", "dot", crossCandidataCurvesList[i][2]);
        WriteVectorVec3d("F:/qiao/iParentCurve2.swc", "dot", crossCandidataCurvesList[i][3]);*/
        //diffuse threshold
        std::vector<double > diffuseThrev = { 1.0, 0.975, 0.95, 0.925, 0.9, 0.875, 0.85, 0.825, 0.8, 0.775, 0.75, 0.725, 0.7 }; //
        std::transform(diffuseThrev.begin(), diffuseThrev.end(), diffuseThrev.begin(), [&](const double& arg){ return arg * startThrev; });
        /*diffuse to check*/
        int ind = -1;
        MatXi flagMat(diffuseThrev.size(), 3); flagMat.setZero();
        int flag = 0;
        for (auto &it : diffuseThrev) {
            ++ind;
            int allPtNum, validSignalSum;
            GetSVMSampleFeatureByRadiusThrev(*origImgPointer, startPoint, it, 9, allPtNum, validSignalSum);
            double ratio = double(validSignalSum) / double(allPtNum);
            if (ratio > 0.3) {
                printf("fog area\n");
                break;
            }
            /*diffuse*/
            VectorVec3i diffusePtSet;
            VectorVec3i growPtSet;
            growPtSet.push_back(startPoint);
            diffusePtSet.push_back(startPoint);
            signalMatrix.SetZero();
            signalMatrix(growPtSet[0](0), growPtSet[0](1), growPtSet[0](2)) = 1;
            for (int k = 1; k <= 100; ++k) {//for safe use int
                SignalPointsIdentifyAreaGrowing(growPtSet, *origImgPointer, signalMatrix, it, 2, 0);
                if (growPtSet.empty()) {
                    break;
                }
                std::copy(growPtSet.begin(), growPtSet.end(), std::back_inserter(diffusePtSet));
            }
            /*intersect*/
           /* char fileline[256];
            sprintf(fileline, "F:/qiao/diffusePtSet%03.3f.swc", it);
            WriteVectorVec3d(fileline, "dot", diffusePtSet);*/
            if (IsIntersect(iCurvePre, diffusePtSet)) flagMat(ind, 0) = 1;
            if (IsIntersect(iCurvePost, diffusePtSet)) flagMat(ind, 1) = 1;
            if (IsIntersect(iParentCurve, diffusePtSet)) flagMat(ind, 2) = 1;
            if (flagMat(ind, 2) == 1 && flagMat(ind, 0) == 0 && flagMat(ind, 1) == 0) flag = 1;
            if (flagMat(ind, 2) == 1 && (flagMat(ind, 0) == 1 || flagMat(ind, 1) == 1) && flag == 1) {
                flag = 2;
                break;
            }
            if (ind == int(diffuseThrev.size())-1  && flagMat(ind, 0) == 0 && flagMat(ind, 1) == 0) flag = 2;//&& flagMat(ind, 2) == 1
        }

        if (flag == 2) {
            removeCrossPairIndList.insert(crossPairInd[i].first);
            removeCrossPairIndList.insert(crossPairInd[i].second);
        }
    }
}

void CrossBranchFilter::FindNearestParentCurve(const Vec3d& head, const std::vector<VectorVec3d> &pCurves, size_t &curInd, size_t &ptInd)
{
    double dist = 10000.0;
    curInd = 100000000;
    ptInd = 100000000;
    for (size_t i = 0; i < pCurves.size(); ++i) {
        for (size_t j = 0; j < pCurves[i].size(); ++j) {
            double tmpdist = (pCurves[i][j] - head).norm();
            if (tmpdist < dist) {
                dist = tmpdist;
                curInd = i;
                ptInd = j;
            }
        }
    }
}

bool CrossBranchFilter::IsIntersect(const VectorVec3i &arg1, const VectorVec3i &arg2)
{
    for (auto &it : arg1) {//arg1 is less
        auto res = std::find(arg2.begin(), arg2.end(), it);
        if (res != arg2.end()) {
            return true;
        }
    }
    return false;
}

void CrossBranchFilter::SetData(const std::vector<VectorVec3d> &arg1, const Vec3d &arg2, const Vec3d &arg3, std::vector<VectorVec3d> &arg4, VectorVec3d &arg5, VectorVec3i &arg6)
{
    tracedCurve_ = &arg1;
    soma_ = &arg2;
    initDir_ = &arg3;
    bifurcationClustre_ = &arg4;
    bifurcationDirections_ = &arg5;
    boundaryLabel_ = &arg6;
}

void CrossBranchFilter::SetData(const std::vector<VectorVec3d> &arg1, const Vec3d &arg2, const Vec3d &arg3, std::vector<VectorVec3d> &arg4, VectorVec3d &arg5)
{
    tracedCurve_ = &arg1;
    soma_ = &arg2;
    initDir_ = &arg3;
    bifurcationClustre_ = &arg4;
    bifurcationDirections_ = &arg5;
    boundaryLabel_ = nullptr;
}

void CrossBranchFilter::WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,
    std::vector<double> &rayNodeWet)
{
    typedef double spheredist;
    int nxss = int(rayNode.size());
    rayNodeWet.clear();

    //зјБъ
    typedef double spheredist;
    spheredist x, y, z;
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

void CrossBranchFilter::GetSVMSampleFeatureByRadiusThrev(const SVolume &origImg, Vec3i seedPt, double extractThrev,
    int radius, int& allPtNum, int& validSignalSum)
{
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();

    seedPt(0) = std::min(std::max(seedPt(0), 0), nx - 1); 
    seedPt(1) = std::min(std::max(seedPt(1), 0), ny - 1);
    seedPt(2) = std::min(std::max(seedPt(2), 0), nz - 1);
    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, seedPt(0), seedPt(1), seedPt(2), radius, radius, radius, 0, nx - 1, 0, ny - 1, 0, nz - 1);

    int xLen = xMax - xMin + 1, yLen = yMax - yMin + 1, zLen = zMax - zMin + 1;

    //SVolume subOrigImg;
    ExtractArea(origImg, xMin, xMax, yMin, yMax, zMin, zMax, subOrigImg_);
    allPtNum = xLen * yLen * zLen;
    VectorVec3i growPtSet;
    growPtSet.push_back(Vec3i(seedPt(0) - xMin, seedPt(1) - yMin, seedPt(2) - zMin));

    //CVolume signalMatrix;
    signalMatrix_.SetSize(xLen, yLen, zLen);
    signalMatrix_(growPtSet[0](0), growPtSet[0](1), growPtSet[0](2)) = 1;
    validSignalSum = 1;
    for (int i = 1; i <= 8; i++) {//for safe use int
        SignalPointsIdentifyAreaGrowing(growPtSet, subOrigImg_, signalMatrix_, extractThrev);
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

void CrossBranchFilter::SignalPointsIdentifyAreaGrowing(VectorVec3i& growPtSet, const SVolume& subOrigImg, CVolume& signalMatrix, double threv, int offset, int ptNumThrev)
{
    int nx = subOrigImg.x();
    int ny = subOrigImg.y();
    int nz = subOrigImg.z();
    size_t nxx = growPtSet.size();
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VectorVec3i resultGrowPtSet;
    Vec3i tmp;
    for (size_t i = 0; i < nxx; i++) {
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, growPtSet[i](0), growPtSet[i](1), growPtSet[i](2),
            offset, offset, offset, 0, nx - 1, 0, ny - 1, 0, nz - 1);
        VectorVec3i tmpGrowPtSet;
        for (int ii = xMin; ii <= xMax; ++ii){
            for (int jj = yMin; jj <= yMax; ++jj){
                for (int ij = zMin; ij <= zMax; ++ij){
                    tmp << ii, jj, ij;
                    //int test1 = signalMatrix(ii, jj, ij);
                    //int test2 = subOrigImg(ii, jj, ij);
                    if (signalMatrix(ii, jj, ij) == 0 && subOrigImg(ii, jj, ij) > threv)
                        tmpGrowPtSet.push_back(tmp);
                }
            }
        }
        if (tmpGrowPtSet.size() > ptNumThrev){
            std::copy(tmpGrowPtSet.begin(), tmpGrowPtSet.end(), std::back_inserter(resultGrowPtSet));
            for (size_t j = 0; j < tmpGrowPtSet.size(); ++j){
                signalMatrix(tmpGrowPtSet[j](0), tmpGrowPtSet[j](1), tmpGrowPtSet[j](2)) = 1;
            }
        }
    }
    resultGrowPtSet.swap(growPtSet);
}