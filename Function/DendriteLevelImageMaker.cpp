#include "DendriteLevelImageMaker.h"
#include "NGUtility.h"
#include "volumealgo.h"

DendriteLevelImageMaker::DendriteLevelImageMaker()
{
    className_ = "DendriteLevelImageMaker";
}

DendriteLevelImageMaker::~DendriteLevelImageMaker()
{

}

ProcStatPointer DendriteLevelImageMaker::Update()
{
    if (!paramPack || !paramPack->activeTree || !paramPack->OrigImage) {
        NG_ERROR_MESSAGE("error in DendriteLevelImageMaker");
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if (layerList_.empty() || branchList_.empty() || layerList_[0] < 0 || branchList_[0] < 0) {
        MAKEPROCESSSTATUS(resSta, false, className_, "no input branch data.");
        return resSta;
    }
    //initial
    NG_CREATE_DYNAMIC_CAST(SVolume, origImg, paramPack->OrigImage);
    if (!paramPack->LayerImage) paramPack->LayerImage = std::make_shared<SVolume>();
    NG_CREATE_DYNAMIC_CAST(SVolume, layerImg, paramPack->LayerImage);
    layerImg->SetSize(origImg->x(), origImg->y(), origImg->z());
    layerImg->SetResolution(origImg->XResolution(), origImg->YResolution(), origImg->ZResolution());
    //get level curves
    auto &curTopo = paramPack->activeTree->m_Topo;
    auto &curCurveList = paramPack->activeTree->m_curveList;
    std::vector<size_t> levelCurveID;
    std::vector<decltype(curTopo.begin())> firstLevel;
    for (auto it = curTopo.begin() + 1; it != curTopo.end(); ++it) {//
        if (layerList_[0] == curTopo.depth(it)) 
            firstLevel.push_back(it);
    }
    //search level subtree 
    for (size_t i = 0; i < branchList_.size(); ++i) {
        if (branchList_[i] >= firstLevel.size()) break;
        auto subTree = curTopo.subtree(firstLevel[branchList_[i]], curTopo.next_sibling(firstLevel[branchList_[i]]));//level minus layerList [0]
        for (auto iter = subTree.begin(); iter != subTree.end(); ++iter) 
            if (NGUtility::IsInVector(layerList_, curTopo.depth(iter)+layerList_[0]))
                levelCurveID.push_back(iter->id);
    }
    std::sort(levelCurveID.begin(), levelCurveID.end());
    //get valid points
    int xMin, xMax, yMin, yMax, zMin, zMax;
    int xLen = paramPack->xMax_ - paramPack->xMin_-1;
    int yLen = paramPack->yMax_ - paramPack->yMin_-1;
    int zLen = paramPack->zMax_ - paramPack->zMin_-1;
    //VectorVec3i validPtSet;
    Vec3i tmpPt;
    for (auto &id : levelCurveID) {//all valid curve id
        for (auto &it : curCurveList[id]) {//vec5d
            if (NGUtility::IsPointInRange(it, paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_)){
                tmpPt << NGUtility::Round(it(0)), NGUtility::Round(it(1)), NGUtility::Round(it(2));
                tmpPt(0) -= paramPack->xMin_; tmpPt(1) -= paramPack->yMin_; tmpPt(2) -= paramPack->zMin_;
                Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, tmpPt(0), tmpPt(1), tmpPt(2), radius_, radius_, radius_, 0, xLen, 0, yLen, 0, zLen);
                for (int i = xMin; i <= xMax; ++i) {
                    for (int j = yMin; j <= yMax; ++j) {
                        for (int ij = zMin; ij <= zMax; ++ij) {
                            layerImg->operator()(i, j, ij) = origImg->operator()(i, j, ij);
                        }
                    }
                }
            }
        }
        //expand at the end of each curve
        VectorVec3d tmpEnd;
        int sz = int(curCurveList[id].size());
        NGUtility::GetPartVectorVec3d(curCurveList[id], std::max(0, sz - 10), sz - 1, tmpEnd);
        Vec3d tailDir;
        NGUtility::CalcPtCurveDirection(tmpEnd, tailDir);
        tailDir *= 2.0;
        tmpEnd.clear();
        Vec3d step;
        for (size_t k = 0; k < 10; ++k) {
            step = curCurveList[id].back().block(0, 0, 3, 1) + tailDir*double(k);
            tmpEnd.push_back(step);
        }
        for (auto &it : tmpEnd) {//vec5d
            if (NGUtility::IsPointInRange(it, paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_)){
                tmpPt << NGUtility::Round(it(0)), NGUtility::Round(it(1)), NGUtility::Round(it(2));
                tmpPt(0) -= paramPack->xMin_; tmpPt(1) -= paramPack->yMin_; tmpPt(2) -= paramPack->zMin_;
                Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, tmpPt(0), tmpPt(1), tmpPt(2), radius_, radius_, radius_, 0, xLen, 0, yLen, 0, zLen);
                for (int i = xMin; i <= xMax; ++i) {
                    for (int j = yMin; j <= yMax; ++j) {
                        for (int ij = zMin; ij <= zMax; ++ij) {
                            layerImg->operator()(i, j, ij) = origImg->operator()(i, j, ij);
                        }
                    }
                }
            }
        }
    }
    

    if (hasSoma_ && NGUtility::IsInVector(layerList_, 1)) {//has soma
        auto id = (curTopo.begin() + 1)->id;
        auto &pt = curCurveList[id].front();//soma
        tmpPt << NGUtility::Round(pt(0)), NGUtility::Round(pt(1)), NGUtility::Round(pt(2));
        tmpPt(0) -= paramPack->xMin_; tmpPt(1) -= paramPack->yMin_; tmpPt(2) -= paramPack->zMin_;
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, tmpPt(0), tmpPt(1), tmpPt(2), 5 * radius_, 5 * radius_, 2 * radius_, 0, xLen, 0, yLen, 0, zLen);
        for (int i = xMin; i <= xMax; ++i) {
            for (int j = yMin; j <= yMax; ++j) {
                for (int ij = zMin; ij <= zMax; ++ij) 
                    layerImg->operator()(i, j, ij) = origImg->operator()(i, j, ij);
            }
        }
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(DendriteLevelImageMaker, TreeCurve)
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(DendriteLevelImageMaker)
