
#include <QMessageBox>
#include "CrossFiberFilter.h"
#include "../traceutil.h"
#include "../../contourutil.h"
#include "../../NGUtility.h"


CrossFiberFilter::CrossFiberFilter()
{
    className_ = std::string("CrossFiberFilter");
    m_Source = std::shared_ptr<TreeCurve>(new TreeCurve(this));
}

CrossFiberFilter::~CrossFiberFilter()
{
}

ProcStatPointer CrossFiberFilter::Update()
{
    //initial
    if (!paramPack) {
        printf( "Error in CrossFiberFilter! Please Set param pack\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "Please Set param pack.");
        return resSta;
    }
    NG_DYNAMIC_CAST(const SVolume, origImgPointer, paramPack->OrigImage);
    if (!origImgPointer || !traceLine_) {
        printf("Error in CrossFiberFilter! Please Set right original img or trace line\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "Please Set right original img or trace line.");
        return resSta;
    }
    //
    //ReviseBifurcationCurves();
    if (traceLine_->empty()){
        //return false;
        MAKEPROCESSSTATUS(resSta, false, className_, "traced Lines is empty.");
        return resSta;
    }
    BuildConInfo();
    //get connect info. the connect just occur between upper level and lower level
    //collect all connect nodes
	
    std::list<Vec3d, Eigen::aligned_allocator<Vec3d>> connectInfoList;
    for (size_t i = 0; i < conInfo_.size(); ++i) {
        if (conInfo_[i](0, 0) != 0) {
            const Vec5d& head = (*traceLine_)[i].front();
            if (std::abs(soma_(0) - head(0)) < 1.0 && std::abs(soma_(1) - head(1)) < 1.0 && std::abs(soma_(2) - head(2)) < 1.0) {
                continue;//soma are not considered.
            }
            for (size_t k = 0; k < traceLine_->size(); ++k) {
                if (k == i) continue;
                auto &curLine = (*traceLine_)[k];
                for (size_t n = 0; n < curLine.size(); ++n)  {
                    auto &pt = curLine[n];
                    if (std::abs(pt(0) - head(0)) < 1.0 && std::abs(pt(1) - head(1)) < 1.0 && std::abs(pt(2) - head(2)) < 1.0) {
                        connectInfoList.push_back(pt.block(0, 0, 3, 1));
                    }
                }
            }
        }
        if (conInfo_[i](0, 1) != 0) {//in ordinary situation , just head has connect info
            const Vec5d& tail = (*traceLine_)[i].front();
            if (std::abs(soma_(0) - tail(0)) < 1.0 && std::abs(soma_(1) - tail(1)) < 1.0 && std::abs(soma_(2) - tail(2)) < 1.0) {
                continue;//soma are not considered.
            }
            for (size_t k = 0; k < traceLine_->size(); ++k) {
                if (k == i) continue;
                auto &curLine = (*traceLine_)[k];
                for (size_t n = 0; n < curLine.size(); ++n)  {
                    auto &pt = curLine[n];
                    if (std::abs(pt(0) - tail(0)) < 1.0 && std::abs(pt(1) - tail(1)) < 1.0 && std::abs(pt(2) - tail(2)) < 1.0) {
                        connectInfoList.push_back(pt.block(0, 0, 3, 1));
                    }
                }
            }
        }
    }
    if (connectInfoList.empty()) {
        MAKEPROCESSSTATUS(resSta, true, className_, "There is no lines to check.");
        return resSta;
    }
    //cluster connect nodes and corresponding line id by distance
    //auto connectInfoListCp = connectInfoList;//back up for debug
    std::vector<VectorVec3d > clustreConnectPt;
    VectorVec3d singleConnectPt;
    int cid = 0;
    while (!connectInfoList.empty()) {
        ++cid;
        VectorVec3d clustre;
        clustre.push_back(connectInfoList.front());
        const Vec3d& curPt = clustre.front();
        connectInfoList.pop_front();
        for (auto it = connectInfoList.begin(); it != connectInfoList.end();) {
            if ((curPt - *it).norm() < 5.0) {
                clustre.push_back(*it);
                it = connectInfoList.erase(it);
            }
            else ++it;
        }
        if (clustre.size() > 1lu) {
            clustreConnectPt.emplace_back();
            clustreConnectPt.back().swap(clustre);
        }
        else singleConnectPt.push_back(clustre[0]);
    }
    auto clustreConnectPtCp = clustreConnectPt;
    auto &rawTraceLine = (*traceLine_);
    //for loop for 2 connect node sets
    for (auto &curConnectPtSet : clustreConnectPt) {
        if (curConnectPtSet.size() < 2lu) continue;
        bool isFind = false;
        //find corresponding lines and break
        size_t oldsz = rawTraceLine.size();
        for (size_t i = 0; i < oldsz; ++i) {// line
            for (size_t j = 0; j < rawTraceLine[i].size(); ++j) {//line point id
                const Vec5d& pt = rawTraceLine[i][j];
                if ((pt.block(0, 0, 3, 1) - curConnectPtSet[0]).norm() < 1.0) {
                    if (j > 3lu && j < rawTraceLine[i].size() - 4lu) {//break
                        rawTraceLine.emplace_back();
                        std::copy(rawTraceLine[i].begin() + j, rawTraceLine[i].end(), std::back_inserter(rawTraceLine.back()));
                        rawTraceLine[i].erase(rawTraceLine[i].begin() + j + 1, rawTraceLine[i].end());//save connect pt
                        isFind = true;
                    }
                }
                if (isFind) break;
            }
            if (isFind) break;
        }
        if (!isFind){
            curConnectPtSet.clear();
            continue;
        }
        //build tree topology
        KPTREE::SetTreeRoot(rawTraceLine, topo_, soma_, 1.0, false, false, false);//because the trace methods decide that the curve is traced outward, so T style is impossible.
        //label curve, of which the head is connectPt
        //for (int k = 0; k < int(clustreConnectPt.size()); ++k) {
        //auto &curConnectPtSet = clustreConnectPt[k];
        //if (curConnectPtSet.size() < 2lu) continue;
        //find the top id. This curve end is connect point
        bool topFlag = false;
        tree<LineID>::pre_order_iterator topIt;
        for (size_t i = 0; i < rawTraceLine.size(); ++i) {// line
            for (size_t j = 0; j < curConnectPtSet.size(); ++j) {
                for (int num = 1; num < 4; ++num) {//the last N point is the cross node.
                    const Vec5d& pt = *(rawTraceLine[i].end() - num);
                    if ((pt.block(0, 0, 3, 1) - curConnectPtSet[j]).norm() < 1.0) {
                        tree<LineID>::iterator curIt = std::find_if(topo_.begin() + 1, topo_.end(), [&](const LineID& arg){return (arg.id % prefix_) == i; });
                        topIt = curIt; topFlag = true;
                        break;
                    }
                }
                if (topFlag) break;
            }
            if (topFlag) break;
        }
        if (!topFlag) {
            //printf("there are many connect points connect to end of curves. will be ignored by crossfiberfilter\n");
            NG_ERROR_MESSAGE("Please Debug");
            system("pause");
            //continue;
        }
        //find corresponding lines and set top id as common parent, and ensure connecting to parent curve
        for (size_t i = 0; i < rawTraceLine.size(); ++i) {// line
            for (size_t j = 0; j < curConnectPtSet.size(); ++j) {//connect point
                const Vec5d& pt = rawTraceLine[i][0];//head
                if ((pt.block(0, 0, 3, 1) - curConnectPtSet[j]).norm() < 1.0) {
                    tree<LineID>::iterator curIt = std::find_if(topo_.begin() + 1, topo_.end(), [&](const LineID& arg){return (arg.id % prefix_) == i; });
                    if (curIt.node == topIt.node) continue;
                    curIt->id += prefix_*(cid + 1);//label represent the connect point cluster num !!!!!
                    topo_.reparent(topIt, curIt, topo_.next_sibling(curIt));
                    //KPTREE::print_tree_bracketed(topo_); std::cout << std::endl;
                    //ensure to connect to the end of parent curve
                    if ((pt.block(0, 0, 3, 1) - rawTraceLine[topIt->id % prefix_].back().block(0, 0, 3, 1)).norm() > 1.0) {
                        rawTraceLine[i][0] = rawTraceLine[topIt->id % prefix_].back();
                    }
                }
            }
        }
        //KPTREE::print_tree_bracketed(topo_); std::cout << std::endl;
        //traverse label tree node
        bool traFlag = true;
        while (traFlag) {
            traFlag = false;
            if (topo_.size() < 4lu) break;//just 2 curves
            for (auto it = topo_.begin() + 2; it != topo_.end(); ++it) {//head curve can not be false curves.
                if (it->id >= prefix_) {
                    auto pIt = topo_.parent(it);
                    if (pIt->id >= prefix_) {
                        printf("debug here\n");
                        system("pause");
                    }
                    auto &pCurve = rawTraceLine[pIt->id % prefix_];
                    auto &cCurve = rawTraceLine[it->id % prefix_];
                    //connect node of parent curve in multi-connect point situation , must be the end part;
                    //calculate main direction
                    Vec3d pDir, cDir;
                    {
                        VectorVec3d tmpPCurve; NGUtility::GetPartVectorVec3d(pCurve, std::max(0, int(pCurve.size()) - 11), int(pCurve.size()) - 1, tmpPCurve);
                        CalcPtCurveDirection(tmpPCurve, pDir);
                        VectorVec3d tmpCCurve; NGUtility::GetPartVectorVec3d(cCurve, 2, std::min(int(cCurve.size()) - 2, 15), tmpCCurve);
                        CalcPtCurveDirection(tmpCCurve, cDir);
                    }
                    // whether obtuse angle
                    //double test = cDir.dot(pDir);
                    if (cDir.dot(pDir) < 0.3) {
                        bool bigObtuseFlagCheck = true;
                        //find symmetry curve
                        for (auto fit = topo_.begin() + 2; fit != topo_.end(); ++fit) {
                            if (topo_.parent(fit) == topo_.parent(it) && fit->id % prefix_ != it->id % prefix_) {//the same parent connect point clustre, just compare the head node
                                //calc main direction
                                auto &sCurve = rawTraceLine[fit->id % prefix_];
                                VectorVec3d tmpCurve; NGUtility::GetPartVectorVec3d(sCurve, 2, std::min(int(sCurve.size()) - 2, 15), tmpCurve);
                                Vec3d sDir;  CalcPtCurveDirection(tmpCurve, sDir);
                                double cosVal = sDir.dot(cDir);
                                if (cosVal < -0.9) {//is symmetry, 25 degree
                                    //clear curves corresponding to 2 branches. traverse branches
                                    //traverse fit
                                    auto subtree = topo_.subtree(fit, topo_.next_sibling(fit));
                                    for (auto sit = subtree.begin(); sit != subtree.end(); ++sit)
                                        rawTraceLine[sit->id % prefix_].clear();
                                    //traverse it
                                    subtree = topo_.subtree(it, topo_.next_sibling(it));
                                    for (auto sit = subtree.begin(); sit != subtree.end(); ++sit)
                                        rawTraceLine[sit->id % prefix_].clear();
                                    //delete 2 branches
                                    topo_.erase(fit);
                                    topo_.erase(it);//in fact just 2 curves can be symmetry
                                    //
                                    traFlag = true;
                                    bigObtuseFlagCheck = false;
                                    break;//fit
                                }
                            }
                        }
                        //if (bigObtuseFlagCheck && cDir.dot(pDir) < -0.7) {//large obtuse angle curve must be deleted.
                        //    //rawTraceLine[it->id % prefix_].clear();//delete curves 
                        //    auto subtree = topo_.subtree(it, topo_.next_sibling(it));
                        //    for (auto sit = subtree.begin(); sit != subtree.end(); ++sit)
                        //        rawTraceLine[sit->id % prefix_].clear();
                        //    topo_.erase(it); //delete this branch.
                        //    traFlag = true;
                        //}
                    }
                    else it->id %= prefix_;
                }
                //std::for_each(topo_.begin(), topo_.end(), [&](LineID& arg){ arg.id = arg.id % prefix_; });
                if (traFlag) // some curves deleted, need to restart checking
                    break;//break for it
            }//for
        }//while
        /*connect the connected lines*/
        bool searchFlag = true;
        while (searchFlag) {
            searchFlag = false;
            for (size_t i = 0; i < rawTraceLine.size(); ++i) {
                if (rawTraceLine[i].empty()) continue;
                const Vec5d &tail = rawTraceLine[i].back();
                for (int j = int(rawTraceLine.size()) - 1; j > -1; --j) {//search from back, because break curves are put into the back of curves set.
                    if (i == j || rawTraceLine[j].empty()) continue;
                    if ((rawTraceLine[j].front().block(0, 0, 3, 1) - tail.block(0, 0, 3, 1)).norm() < 1.0) {
                        searchFlag = true;
                        size_t oldsz = rawTraceLine[i].size();
                        rawTraceLine[i].resize(rawTraceLine[i].size() + rawTraceLine[j].size());
                        std::copy(rawTraceLine[j].begin(), rawTraceLine[j].end(), rawTraceLine[i].begin() + oldsz);
                        rawTraceLine[j].clear();
                        break;//reset search
                    }
                }
                if (searchFlag) break;
            }
        }
        /*remove empty lines*/
        rawTraceLine.erase(std::remove_if(rawTraceLine.begin(), rawTraceLine.end(), [](const Line5d& arg){return arg.empty(); }), rawTraceLine.end());//delete empty curves
        //}
    }
        //clustreConnectPt.erase(std::remove_if(clustreConnectPt.begin(), clustreConnectPt.end(), [](VectorVec3d &arg){return arg.empty(); }), clustreConnectPt.end());
    /*remove single obtuse angle curves*/
    //for (auto &it : singleConnectPt) {
    //    for (size_t i = 0; i < rawTraceLine.size(); ++i) {// line
    //        if (rawTraceLine[i].empty()) continue;
    //        const Vec5d& pt = rawTraceLine[i][0];//head
    //        if ((pt.block(0, 0, 3, 1) - it).norm() < 1.0) {
    //            tree<LineID>::iterator curIt = std::find(topo_.begin() + 1, topo_.end(), LineID(i));
    //            auto pIt = topo_.parent(curIt);
    //            if (pIt->id < prefix_) {
    //                auto &pCurve = rawTraceLine[pIt->id % prefix_];
    //                auto &cCurve = rawTraceLine[curIt->id % prefix_];
    //                Vec3d pDir, cDir;
    //                {
    //                    VectorVec3d tmpPCurve; NGUtility::GetPartVectorVec3d(pCurve, std::max(0, int(pCurve.size()) - 11), int(pCurve.size()) - 1, tmpPCurve);
    //                    CalcPtCurveDirection(tmpPCurve, pDir);
    //                    VectorVec3d tmpCCurve; NGUtility::GetPartVectorVec3d(cCurve, 0, std::min(int(cCurve.size()) - 1, 10), tmpCCurve);
    //                    CalcPtCurveDirection(tmpCCurve, cDir);
    //                }
    //                if (cDir.dot(pDir) < -0.3) {//large obtuse angle curve must be deleted.
    //                    //rawTraceLine[it->id % prefix_].clear();//delete curves 
    //                    auto subtree = topo_.subtree(curIt, topo_.next_sibling(curIt));
    //                    for (auto sit = subtree.begin(); sit != subtree.end(); ++sit)
    //                        rawTraceLine[sit->id % prefix_].clear();
    //                    topo_.erase(curIt); //delete this branch.
    //                }
    //            }
    //        }
    //    }
    //}
    
    KPTREE::SetTreeRoot(rawTraceLine, topo_, soma_);
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(CrossFiberFilter, TreeCurve)
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(CrossFiberFilter)

void CrossFiberFilter::CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir)
{
    if (ptCurve.empty()) {
        fir.setZero();
        return;
    }
    auto nx = ptCurve.size();
    Vec3d tmpFir;
    tmpFir.setZero();
    double tmp(0.0);
    Vec3d tmpVec;
    for (auto i = 0; i < nx - 1llu; ++i){
        tmpVec = ptCurve[i + 1] - ptCurve[i];
        tmp = tmpVec.norm();
        tmpFir += tmp * tmpVec;
    }
    fir = tmpFir.normalized();
}

void CrossFiberFilter::AddCollideConnectNode(const std::vector<VectorVec5d> &dendCurves, VectorMat2i &dendConInfo,
    const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList)
{
    std::vector<VectorVec5d>::size_type kk = dendCurves.size();
    resultDendList.clear();
    resultDendList = dendCurves;

    for (std::vector<VectorVec5d>::size_type ij = 0; ij < kk; ++ij){
        Mat2i &currentConInfo = dendConInfo[ij];
        Mat2i conCp;
        do {
            conCp = currentConInfo;
            AddCollideConnectNodeSub( dendConInfo, currentConInfo, ij, somaList, resultDendList);
        } while (conCp(0, 0) != currentConInfo(0, 0) || conCp(0, 1) != currentConInfo(0, 1));
    }
}

void CrossFiberFilter::AddCollideConnectNodeSub( VectorMat2i &dendConInfo, Mat2i& currentConInfo, size_t ij,
    const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList)
{

    if (currentConInfo(0, 0) > 0){
        ConnectNeighborCurveToHead(resultDendList, dendConInfo, ij, currentConInfo);
    }//if

    if (currentConInfo(0, 0) < 0){
        const VectorVec5d& currentDend = resultDendList[ij];
        int conInfo = -currentConInfo(0, 0) - 1;
        VectorVec5d tmp;
        Vec5d myce;
        myce << somaList[conInfo](0), somaList[conInfo](1), somaList[conInfo](2),
            1.0, 0.0;
        if ((currentDend[0].block(0, 0, 3, 1) - myce.block(0, 0, 3, 1)).norm() > 0.001) {
            tmp.reserve(1llu + currentDend.size());
            tmp.push_back(myce);
            std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
            resultDendList[ij].clear();
            resultDendList[ij].swap(tmp);
        }
    }

    if (currentConInfo(0, 1) < 0) {
        const VectorVec5d& currentDend = resultDendList[ij];
        int conInfo = -currentConInfo(0, 1) - 1;
        Vec5d myce;
        myce << somaList[conInfo](0), somaList[conInfo](1), somaList[conInfo](2),
            1.0, 0.0;
        if ((currentDend.back().block(0, 0, 3, 1) - myce.block(0, 0, 3, 1)).norm() > 0.001) {
            resultDendList[ij].push_back(myce);
        }
    }

    if (currentConInfo(0, 1) > 0){
        ConnectNeighborCurveToTail(resultDendList, dendConInfo, ij, currentConInfo);
    }
}


void CrossFiberFilter::ConnectNeighborCurveToHead(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t ij, Mat2i &currentConInfo)
{
    const VectorVec5d& currentDend = resultDendList[ij];
    Vec3d firstDendPt;
    firstDendPt << currentDend[0](0), currentDend[0](1), currentDend[0](2);

    VectorVec3d data16;
    Vec3d tmpData;
    //Change:2014-3-21-10-28
    for (int i = std::min(5, int(currentDend.size())); i > -1; --i){
        tmpData << currentDend[i](0), currentDend[i](1), currentDend[i](2);
        data16.push_back(tmpData);
    }
    Vec3d mainDir;
    CalcPtCurveDirection(data16, mainDir);

    VectorVec5d &conDendCurve = resultDendList[currentConInfo(0, 0) - 1];

    VectorVec5d::size_type nxx = conDendCurve.size();

    std::vector<std::pair<int, double> > conNodeDistList;
    Vec3d datap2;
    Vec3d datap22;
    for (VectorVec5d::size_type ii = 0; ii < nxx; ++ii){
        datap2 << conDendCurve[ii](0), conDendCurve[ii](1), conDendCurve[ii](2);
        datap22 = firstDendPt - datap2;
        conNodeDistList.push_back(std::pair<int, double>(int(ii) + 1, datap22.norm()));
    }

    std::sort(conNodeDistList.begin(), conNodeDistList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){ return lhs.second < rhs.second; });

    //Change:2014-3-21-10-32
    //int isShortestExist = std::min(1, (int)nxx);//1 or 0 // why 1 ? It just choose one closed points;
    if ((conNodeDistList[0].second - 0.0) < 0.001) {
        return; // has done connect
    }
    int isShortestExist = std::min(20, (int)nxx);//2017-5-18
    std::vector<int> shortestNodeId;
    for (int ii = 0; ii < isShortestExist; ++ii){
        if (ii > 0 && conNodeDistList[ii].second > 10.0) {
            isShortestExist = ii;
            break;
        }
        shortestNodeId.push_back(conNodeDistList[ii].first);
    }
    if (isShortestExist < 1) {
        printf("Lima\n");
    }
    if (conNodeDistList[0].second < 2.0) {
        currentConInfo(1, 0) = conNodeDistList[0].first;
    }
    else{
        //2017-5-26 sort by signal
        const SVolume &origImg = *origImgPointer;
        std::vector<std::pair<int, double> > conNodeSignalList;
        Vec3d midPt; double signal_;
        for (auto &it : shortestNodeId) {
            Vec5d &conNode = conDendCurve[it - 1];
            midPt = 2.0 / 3.0*firstDendPt + 1.0 / 3.0*conNode.block(0, 0, 3, 1);
            signal_ = NGUtility::WeighRayValue(midPt, origImg);
            //midPt = 1.0 / 3.0*firstDendPt + 2.0 / 3.0*conNode.block(0, 0, 3, 1);
            //signal_ += WeighRayValue(midPt, origImg);
            conNodeSignalList.emplace_back(it, signal_);
        }
        std::sort(conNodeSignalList.begin(), conNodeSignalList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){return lhs.second > rhs.second; });
        currentConInfo(1, 0) = conNodeSignalList[0].first;
        if (conNodeSignalList.size() > 3lu) {
            std::sort(conNodeSignalList.begin(), conNodeSignalList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){return lhs.second > rhs.second; });
            conNodeSignalList.erase(conNodeSignalList.begin() + 3lu, conNodeSignalList.end());
        }
        std::vector<double> conNodeAngleList;
        Vec3d data2p;
        for (int ii = 0; ii < conNodeSignalList.size(); ++ii){
            data2p = conDendCurve[conNodeSignalList[ii].first - 1].block(0, 0, 3, 1);
            double den = (data2p - firstDendPt).norm();
            double num = (data2p - firstDendPt).transpose() * mainDir;
            conNodeAngleList.push_back(std::abs(num / den));
        }

        std::vector<double>::const_iterator maxItem = std::max_element(conNodeAngleList.begin(), conNodeAngleList.end());
        currentConInfo(1, 0) = conNodeSignalList[maxItem - conNodeAngleList.begin()].first;
        //if parallel ? 
        {
            Vec3d conMainDir;
            int id = conNodeDistList[0].first - 1;
            VectorVec3d tmpCurve;
            NGUtility::GetPartVectorVec3d(conDendCurve, std::max(0, id - 3), std::min(int(conDendCurve.size()) - 1, id + 3), tmpCurve);
            CalcPtCurveDirection(tmpCurve, conMainDir);
            double tmpA = std::abs(conMainDir.dot(mainDir));
            if (tmpA > 0.9) {
                if (currentConInfo(1, 0) < 12) currentConInfo(1, 0) = 1;
                if (currentConInfo(1, 0) > nxx - 12) currentConInfo(1, 0) = int(nxx);
            }
        }
    }


    //2017-5-16 merge
    VectorVec5d tmp;
    if ((currentConInfo(1, 0) < 5 && dendConInfo[currentConInfo(0, 0) - 1](0, 0) == 0)
        || (currentConInfo(1, 0) > nxx - 4 && dendConInfo[currentConInfo(0, 0) - 1](0, 1) == 0)) {
        tmp.reserve(conDendCurve.size() + currentDend.size());
        int conDendId = currentConInfo(0, 0);
        if (currentConInfo(1, 0) < 5) {// connect to head
            currentConInfo(1, 0) = 1;
            std::copy(conDendCurve.rbegin(), conDendCurve.rend() - 3, std::back_inserter(tmp));
            std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
            currentConInfo(0, 0) = dendConInfo[currentConInfo(0, 0) - 1](0, 1);//tail
        }
        else {//connect to tail
            currentConInfo(1, 0) = int(nxx);
            std::copy(conDendCurve.begin(), conDendCurve.end() - 3, std::back_inserter(tmp));
            std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
            currentConInfo(0, 0) = dendConInfo[currentConInfo(0, 0) - 1](0, 0);//head
        }
        //clear dendConInfo
        conDendCurve.clear();
        dendConInfo[conDendId - 1].setZero();
        for (size_t kk = 0; kk < dendConInfo.size(); ++kk){
            if (kk == ij) {//loop
                if (dendConInfo[kk](0, 1) == conDendId) dendConInfo[kk](0, 1) = 0;
            }
            if (dendConInfo[kk](0, 0) == conDendId) dendConInfo[kk](0, 0) = int(ij) + 1;
            if (dendConInfo[kk](0, 1) == conDendId) dendConInfo[kk](0, 1) = int(ij) + 1;
        }
    }
    else{
        tmp.push_back(conDendCurve[currentConInfo(1, 0) - 1]);
        std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
        if (dendConInfo[ij](0, 1) == dendConInfo[ij](0, 0)) dendConInfo[ij](0, 1) = 0;//current curve loop
        if (dendConInfo[currentConInfo(0, 0) - 1](0, 0) == int(ij) + 1) dendConInfo[currentConInfo(0, 0) - 1](0, 0) = 0;//conj curve loop
        if (dendConInfo[currentConInfo(0, 0) - 1](0, 1) == int(ij) + 1) dendConInfo[currentConInfo(0, 0) - 1](0, 1) = 0;//conj curve loop
    }
    resultDendList[ij].swap(tmp);
}



void CrossFiberFilter::ConnectNeighborCurveToTail(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t ij, Mat2i &currentConInfo)
{
    const VectorVec5d& currentDend = resultDendList[ij];
    Vec3d firstDendPt;
    int nxx1 = int(currentDend.size());
    //ChangeVec3d(data1[nxx1 - 1], datap);
    firstDendPt << currentDend[nxx1 - 1](0), currentDend[nxx1 - 1](1), currentDend[nxx1 - 1](2);

    VectorVec3d data16;
    Vec3d tmpData;
    //Change:2014-3-21-10-42 8->6
    for (int i = std::max(nxx1 - 6, 0); i < nxx1; ++i){
        //ChangeVec3d(data1[i], tmpData);
        tmpData << currentDend[i](0), currentDend[i](1), currentDend[i](2);
        data16.push_back(tmpData);
    }
    Vec3d mainDir;
    //Vec3d directions2;
    //Vec3d directions3;
    CalcPtCurveDirection(data16, mainDir);

    VectorVec5d &conDendCurve = resultDendList[currentConInfo(0, 1) - 1];

    VectorVec5d::size_type nxx = conDendCurve.size();

    std::vector<std::pair<int, double> > conNodeDistList;
    Vec3d datap2;
    Vec3d datap22;
    for (VectorVec5d::size_type ii = 0; ii < nxx; ++ii){
        //ChangeVec3d(data2[ii], datap2);
        datap2 << conDendCurve[ii](0), conDendCurve[ii](1), conDendCurve[ii](2);
        datap22 = firstDendPt - datap2;
        conNodeDistList.push_back(std::pair<int, double>(int(ii) + 1, datap22.norm()));
    }

    std::sort(conNodeDistList.begin(), conNodeDistList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){ return lhs.second < rhs.second;  });

    if ((conNodeDistList[0].second - 0.0) < 0.001) {
        return; // has done connect
    }
    //Change:5->1
    //int isShortestExist = std::min(1, int(nxx));
    //int isShortestExist = std::min(1, (int)nxx);//1 or 0 // why 1 ? It just choose one closed points;
    size_t isShortestExist = std::min(20llu, nxx);//2017-5-18
    std::vector<int> shortestNodeId;
    for (size_t ii = 0; ii < isShortestExist; ++ii){
        if (conNodeDistList[ii].second > 10.0) break;
        shortestNodeId.push_back(conNodeDistList[ii].first);
    }

    isShortestExist = shortestNodeId.size();
    if (isShortestExist < 1llu) {
        printf("Lima\n");
    }
    if (conNodeDistList[0].second < 2.0) {
        currentConInfo(1, 1) = conNodeDistList[0].first;
    }
    else{
        //2017-5-26 sort by signal
        //const SVolume &origImg = *origImgPointer;
        std::vector<std::pair<int, double> > conNodeSignalList;
        Vec3d midPt; double signal_;
        for (auto &it : shortestNodeId) {
            Vec5d &conNode = conDendCurve[it - 1];
            midPt = 2.0 / 3.0*firstDendPt + 1.0 / 3.0*conNode.block(0, 0, 3, 1);
            signal_ = NGUtility::WeighRayValue(midPt, *origImgPointer);
            //midPt = 1.0 / 3.0*firstDendPt + 2.0 / 3.0*conNode.block(0, 0, 3, 1);
            //signal_ += WeighRayValue(midPt, origImg);
            conNodeSignalList.emplace_back(it, signal_);
        }
        std::sort(conNodeSignalList.begin(), conNodeSignalList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){return lhs.second > rhs.second; });
        currentConInfo(1, 1) = conNodeSignalList[0].first;
        if (conNodeSignalList.size() > 3lu) {
            std::sort(conNodeSignalList.begin(), conNodeSignalList.end(), [](const std::pair<int, double>& lhs, const std::pair<int, double>& rhs){return lhs.second > rhs.second; });
            conNodeSignalList.erase(conNodeSignalList.begin() + 3lu, conNodeSignalList.end());
        }

        std::vector<double> conNodeAngleList;
        Vec3d data2p;
        for (int ii = 0; ii < conNodeSignalList.size(); ++ii){
            data2p = conDendCurve[conNodeSignalList[ii].first - 1].block(0, 0, 3, 1);
            double den = (data2p - firstDendPt).norm();
            double num = (data2p - firstDendPt).transpose() * mainDir;
            conNodeAngleList.push_back(std::abs(num / den));
        }

        std::vector<double>::const_iterator maxItem = std::max_element(conNodeAngleList.begin(), conNodeAngleList.end());
        currentConInfo(1, 1) = conNodeSignalList[maxItem - conNodeAngleList.begin()].first;
        //if parallel ? 
        {
            Vec3d conMainDir;
            int id = conNodeDistList[0].first - 1;
            VectorVec3d tmpCurve;
            NGUtility::GetPartVectorVec3d(conDendCurve, std::max(0, id - 3), std::min(int(conDendCurve.size()) - 1, id + 3), tmpCurve);
            CalcPtCurveDirection(tmpCurve, conMainDir);
            double tmpA = std::abs(conMainDir.dot(mainDir));
            if (tmpA > 0.9) {
                if (currentConInfo(1, 1) < 12) currentConInfo(1, 1) = 1;
                if (currentConInfo(1, 1) > nxx - 12) currentConInfo(1, 1) = int(nxx);
            }
        }
    }

    //2017-5-16 merge
    VectorVec5d tmp;
    if ((currentConInfo(1, 1) < 5 && dendConInfo[currentConInfo(0, 1) - 1](0, 0) == 0)
        || (currentConInfo(1, 1) > nxx - 4 && dendConInfo[currentConInfo(0, 1) - 1](0, 1) == 0)) {
        int conDendId = currentConInfo(0, 1);
        tmp.reserve(conDendCurve.size() + currentDend.size());
        if (currentConInfo(1, 1) < 5) {// connect to head
            currentConInfo(1, 1) = 1;
            std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
            std::copy(conDendCurve.begin() + 3, conDendCurve.end(), std::back_inserter(tmp));
            currentConInfo(0, 1) = dendConInfo[currentConInfo(0, 1) - 1](0, 1);//tail
        }
        else{//connect to tail
            currentConInfo(1, 1) = int(nxx);
            std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
            std::copy(conDendCurve.rbegin() + 3, conDendCurve.rend(), std::back_inserter(tmp));
            currentConInfo(0, 1) = dendConInfo[currentConInfo(0, 1) - 1](0, 0);//head
        }
        //clear dendConInfo
        conDendCurve.clear();
        dendConInfo[conDendId - 1].setZero();
        for (size_t kk = 0; kk < dendConInfo.size(); ++kk){
            if (kk == ij) {//loop
                if (dendConInfo[kk](0, 0) == conDendId) dendConInfo[kk](0, 0) = 0;
            }
            if (dendConInfo[kk](0, 0) == conDendId) dendConInfo[kk](0, 0) = int(ij) + 1;
            if (dendConInfo[kk](0, 1) == conDendId) dendConInfo[kk](0, 1) = int(ij) + 1;
        }
    }
    else{
        tmp.reserve(currentDend.size() + 1llu);
        std::copy(currentDend.begin(), currentDend.end(), std::back_inserter(tmp));
        //for (VectorVec5d::size_type k = 0; k < data1.size(); ++k){
        //    tmp.push_back(data1[k]);
        //}
        tmp.push_back(conDendCurve[currentConInfo(1, 1) - 1]);
        //resultDendList[ij].clear();
        if (dendConInfo[ij](0, 0) == dendConInfo[ij](0, 1)) dendConInfo[ij](0, 0) = 0;//current curve loop
        if (dendConInfo[currentConInfo(0, 1) - 1](0, 0) == int(ij) + 1) dendConInfo[currentConInfo(0, 1) - 1](0, 0) = 0;//conj curve loop
        if (dendConInfo[currentConInfo(0, 1) - 1](0, 1) == int(ij) + 1) dendConInfo[currentConInfo(0, 1) - 1](0, 1) = 0;//conj curve loop
    }
    resultDendList[ij].swap(tmp);
}

void CrossFiberFilter::ReviseBifurcationCurves()
{
    //construct dendConInfo 
    auto &rawTraceLine = *traceLine_;
    KPTREE::SetTreeRoot(rawTraceLine, topo_, soma_, 1.0, false);
    VectorMat2i dendConInfo(rawTraceLine.size());
    std::for_each(dendConInfo.begin(), dendConInfo.end(), [](Mat2i& arg){arg.setZero(); });
    //int cFlag = 0;
    for (size_t i = 0; i < rawTraceLine.size(); ++i) {
        auto curIt = std::find(topo_.begin() + 1, topo_.end(), LineID(i));
        auto pIt = topo_.parent(curIt);
        if (pIt->id < prefix_ && pIt->id >= 0llu) {
                dendConInfo[i](0, 0) = int(pIt->id) + 1;
                //rawTraceLine[i].erase(rawTraceLine[i].begin());
        }
    }
    //revise connect part
    VectorVec3d somaList; somaList.push_back(soma_); std::vector<Line5d> tmpDendCurves;
    AddCollideConnectNode(rawTraceLine, dendConInfo, somaList, tmpDendCurves);
    std::swap(tmpDendCurves, rawTraceLine);
    rawTraceLine.erase(std::remove_if(rawTraceLine.begin(), rawTraceLine.end(), [](const Line5d& arg){return arg.empty(); }), rawTraceLine.end());//delete empty curves
}

void CrossFiberFilter::BuildConInfo()
{
    auto &rawTraceLine = *traceLine_;
    conInfo_.resize(rawTraceLine.size());
    std::for_each(conInfo_.begin(), conInfo_.end(), [](Mat2i &arg){arg.setZero(); });
    int sz = int(rawTraceLine.size());
    auto &root = rawTraceLine[0].front();
    for (int i = 0; i < sz; ++i) {// root curve cannot be
        const Vec5d &head = rawTraceLine[i].front();
        const Vec5d &tail = rawTraceLine[i].back();
        for (int j = 0; j < sz; ++j) {
            if (i == j) continue;
            for (auto &it : rawTraceLine[j]) {
                if ((head.block(0, 0, 3, 1) - root.block(0, 0, 3, 1)).norm() < 3.0 ) 
                    continue;
                if ( (head.block(0,0,3,1) -it.block(0,0,3,1)).norm() < 1.0 ) 
                    conInfo_[i](0, 0) = j + 1;
                if ((tail.block(0, 0, 3, 1) - it.block(0, 0, 3, 1)).norm() < 1.0) 
                    conInfo_[i](0, 1) = j + 1;
            }
        }
    }
}

//void CrossFiberFilter::WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,
//    std::vector<double> &rayNodeWet)
//{
//    typedef double spheredist;
//    int nxss = int(rayNode.size());
//    rayNodeWet.clear();
//    typedef double spheredist;
//    spheredist x, y, z;
//    spheredist segmentDense, segmentWet;//, w1;
//    for (int i = 0; i < nxss; ++i){
//        x = rayNode[i](0);
//        y = rayNode[i](1);
//        z = rayNode[i](2);
//        segmentDense = segmentWet = 0.0;
//        ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
//        rayNodeWet.push_back(segmentDense / (segmentWet + 0.0001));
//    }
//}
//
//double CrossFiberFilter::WeighRayValue(const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg)
//{
//    typedef double spheredist;
//    spheredist x, y, z;
//    spheredist segmentDense, segmentWet;//, w1;
//    x = rayNode(0);
//    y = rayNode(1);
//    z = rayNode(2);
//    segmentDense = segmentWet = 0.0;
//    ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
//    return segmentDense / (segmentWet + 0.0001);
//}
