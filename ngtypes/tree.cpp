#include "tree.h"
#include "../Function/NGUtility.h"


TreeConnect::TreeConnect(INeuronProcessObject *p)
{
    m_ProcessObject =  p;
    identifyName =  std::string("TreeConnect");
}

TreeConnect::~TreeConnect()
{

}

bool TreeConnect::IsEmpty() const
{
    return m_Connect.empty();
}

void TreeConnect::Swap(VectorMat2i &arg)
{
    m_Connect.swap(arg);
}

void TreeConnect::SetConnect(const VectorMat2i &arg)
{
    m_Connect = arg;
}


TreeCurve::TreeCurve(INeuronProcessObject *p)
{
    m_ProcessObject =  p;
    identifyName =  std::string("TreeCurve");
}

TreeCurve::~TreeCurve()
{

}

bool TreeCurve::IsEmpty() const
{
    return m_Curve.empty();
}

void TreeCurve::Swap(std::vector<VectorVec5d> &arg)
{
    m_Curve.swap(arg);
}

void TreeCurve::SetCurve(std::vector<VectorVec5d> &arg)
{
    m_Curve = arg;
}


namespace KPTREE{

    size_t FindNearestID(const Vec3d& pt, const Line5d& line, double threv)
    {
        size_t id = std::numeric_limits<size_t>::max();
        double mindist = 100000, curDist;
        double xDist, yDist, zDist;
        for (size_t i = 0; i < line.size(); ++i) {
            xDist = std::abs(pt(0) - line[i](0));
            yDist = std::abs(pt(1) - line[i](1));
            zDist = std::abs(pt(2) - line[i](2));
            if (xDist< threv && yDist < threv && zDist< threv)  {
                curDist = xDist + yDist + zDist;
                if (curDist < mindist) {
                    mindist = curDist;
                    id = i;
				//if (curDist < threv) return id;
                }
            }
        }
        return id;
    }


    bool SearchLineIDInTree(const NeuronTree& nTree, size_t id, tree<LineID>::iterator &resit)
    {

        if (nTree.m_curveList.empty() || nTree.m_Topo.size() <= 1){
            printf("error in SearchLineIDInTree.\n");
            return false;
        }
        if ((resit = std::find(nTree.m_Topo.begin() + 1, nTree.m_Topo.end(), LineID(id))) != nTree.m_Topo.end()) return true;
        return false;
    }

    bool SearchLineIDInPopulation(const NeuronPopulation& treeList, size_t id, tree<LineID>::iterator &resit)
    {
        if (treeList.IsEmpty()) {
            printf("error in SearchLineIDInTreeSearchLineIDInPopulation.\n");
            return false;
        }
        bool flag = false;
        for (auto &it : treeList.m_pop) {
            if (SearchLineIDInTree(*it, id, resit)) {
                flag = true;
                break;
            }
        }
        return flag;
    }

    void UpdateDeleteTreeIDList(NeuronTree& nTree, size_t deletedID)
    {
        for (auto iter = nTree.m_Topo.begin(); iter != nTree.m_Topo.end(); ++iter) {
            if ((*iter).id >= deletedID){//normally it cannot be equal
                --(*iter).id;
            }
        }
    }

    void UpdateDeleteTreeIDList(NeuronTree& nTree, std::vector<size_t>& deletedIDList)
    {
        if (deletedIDList.empty()) {
            printf("UpdateDeleteTreeIDList: why the deletedIDList is empty?\n");
            return;
        }
        std::sort(deletedIDList.begin(), deletedIDList.end());
        UpdateDeleteSortedTreeIDList(nTree, deletedIDList);
    }

    void UpdateDeleteSortedTreeIDList(NeuronTree& nTree, std::vector<size_t>& deletedIDList)
    {
        if (deletedIDList.empty()) {
            printf("UpdateDeleteSortedTreeIDList: why the deletedIDList is empty?\n");
            return;
        }
        if (deletedIDList.size() == 1lu) {
            UpdateDeleteTreeIDList(nTree, deletedIDList.front());
            return;
        }
        for (auto iter = nTree.m_Topo.begin(); iter != nTree.m_Topo.end(); ++iter) {
            if ((*iter).id < deletedIDList.front()) continue;
            size_t tmpID = (*iter).id;
            for (size_t j = 0; j < deletedIDList.size(); ++j) {
                if (tmpID >= deletedIDList[j]){//normally it cannot be equal
                    if (tmpID == 0llu)
                        NG_ERROR_MESSAGE("tree id is 0, cannot minus");
                    --(*iter).id;
                }
            }
        }
    }

    bool SetTreeRoot(NeuronTree&nTree, const Vec3d& root, double threv, bool checkFlag, bool castrationFlag, bool refreshFlag)
    {
        //return SetTreeRoot(nTree.m_curveList, nTree.m_Topo, root, threv);
    //}

    //bool SetTreeRoot(std::vector<Line5d>& curve, tree<LineID>& topo, const Vec3d& root, double threv /*= 0.5*/, bool checkFlag)
    //{
        auto &curve = nTree.m_curveList;
        auto &topo = nTree.m_Topo;
        tree<LineID> resTree;
        //search corresponding tree
        bool flag = false;
        for (auto &it : curve) {
            if ((it.front().block(0, 0, 3, 1) - root).norm() < 3.0 || (it.back().block(0, 0, 3, 1) - root).norm() < 3.0) {
                flag = true;
                break;
            }
        }
        //check
        if (!flag) {
            //printf("error occurred, cannot find corresponding soma for tree\n");
            return false;
        }
        //collect all id
        std::vector<size_t > curCurveIDList(curve.size());
        std::vector<size_t > curCurveIDListCp(curve.size());
        {
            size_t i = 0;
            std::for_each(curCurveIDList.begin(), curCurveIDList.end(), [&](size_t& arg){arg = i; ++i; });
        }
        //iterative search tree level. curCurveIDList save current tree curve list. currentTreeId is corresponding to current tree id
        resTree.clear();
        resTree.set_head(LineID(std::numeric_limits<size_t>::max()));//root null
        //1-level
        double headDist, tailDist;
        curCurveIDListCp = curCurveIDList;
        NEURITETYPE curNeuriteType = NEURITETYPE::DENDRITE;
        for (auto &it : curCurveIDList) {
            flag = false;
            headDist = tailDist = 10000000000.0;
            headDist = (curve[it].front().block(0, 0, 3, 1) - root).norm();
            tailDist = (curve[it].back().block(0, 0, 3, 1) - root).norm();
            if (topo.size() > 1) {
                auto oldTreeNode = std::find(topo.begin() + 1, topo.end(), LineID(it));
                if (oldTreeNode != topo.end())
                    curNeuriteType = oldTreeNode->type;
            }
            if (headDist < threv) {
                resTree.append_child(resTree.begin(), LineID(it, curNeuriteType));
                flag = true;
            }
            else if (tailDist < threv) {
                std::reverse(curve[it].begin(), curve[it].end());
                resTree.append_child(resTree.begin(), LineID(it, curNeuriteType));
                flag = true;
            }else{
                auto nearestID = KPTREE::FindNearestID(root, curve[it],2.0);
                if (nearestID == 0lu) {
                    resTree.append_child(resTree.begin(), LineID(it, curNeuriteType));
                    flag = true;
                }else if (nearestID == curve[it].size() - 1lu) {
                    std::reverse(curve[it].begin(), curve[it].end());
                    resTree.append_child(resTree.begin(), LineID(it, curNeuriteType));
                    flag = true;
                }
                else if (nearestID < curve[it].size() - 1lu) {
                    size_t newID = curve.size();
                    curve.push_back(VectorVec5d());
                    std::move(curve[it].begin() + nearestID, curve[it].end(), std::back_inserter(curve.back()));
                    curve[it].erase(curve[it].begin() + nearestID, curve[it].end());
                    std::reverse(curve[it].begin(), curve[it].end());
                    resTree.append_child(resTree.begin(), LineID(it, curNeuriteType));
                    resTree.append_child(resTree.begin(), LineID(newID, curNeuriteType));
                    flag = true;
                }
            }
            if(flag)
                curCurveIDListCp.erase(std::remove(curCurveIDListCp.begin(), curCurveIDListCp.end(), it), curCurveIDListCp.end());
        }
        //KPTREE::print_tree_bracketed(resTree);
        //printf("\n");
        curCurveIDListCp.swap(curCurveIDList);
        size_t nid;
        size_t childID;
        std::vector<tree<LineID>::leaf_iterator> iterList;
        curNeuriteType = NEURITETYPE::DENDRITE;
        for (auto itLeaf = resTree.begin_leaf(); itLeaf != resTree.end_leaf(); ++itLeaf)
            iterList.push_back(itLeaf);//valid leaf iterator
        while (!curCurveIDList.empty()) {
            size_t oldsz = curCurveIDList.size();
            std::vector<tree<LineID>::leaf_iterator> newIterList;//next leaf
            for (auto it : iterList){//each leaf
                if (curCurveIDListCp.size()!=curCurveIDList.size()) curCurveIDListCp = curCurveIDList;
                for (size_t itID = 0; itID < curCurveIDList.size(); ++itID) {//each left node
                    childID = curCurveIDList[itID];
                    if (topo.size() > 1) {
                        auto oldTreeNode = std::find(topo.begin() + 1, topo.end(), LineID(childID));
                        if (oldTreeNode != topo.end())
                            curNeuriteType = oldTreeNode->type;
                    }
                    if ((nid = FindNearestID(curve[childID].front().block(0,0,3,1), curve[(*it).id])) != std::numeric_limits<size_t>::max()){
                        newIterList.push_back(resTree.append_child(it, LineID(childID, curNeuriteType)));
                    }
                    else if ((nid = FindNearestID(curve[childID].back().block(0,0,3,1), curve[(*it).id])) != std::numeric_limits<size_t>::max()){
                        std::reverse(curve[childID].begin(), curve[childID].end());
                        newIterList.push_back(resTree.append_child(it, LineID(childID, curNeuriteType)));
                    }
                    else if ((nid = FindNearestID(curve[(*it).id].back().block(0,0,3,1), curve[childID])) != std::numeric_limits<size_t>::max()){
                        if (nid == 0 || nid >= curve[childID].size() - 1) {
                            printf("please check the T style curve.\n");
                            return false;
                        }
                        curve.push_back(Line5d());
                        //curCurveIDList.push_back(lineList.size() - 1lu);// though push back, it still be deleted later
                        std::copy(curve[childID].begin() + nid, curve[childID].end(), std::back_inserter(curve.back()));
                        curve[childID].erase(curve[childID].begin() + nid + 1, curve[childID].end());
                        std::reverse(curve[childID].begin(), curve[childID].end());
                        newIterList.push_back(resTree.append_child(it, LineID(childID, curNeuriteType)));
                        newIterList.push_back(resTree.append_child(it, LineID(curve.size() - 1lu, curNeuriteType)));
                    }
                    else if ((nid = FindNearestID(curve[(*it).id].front().block(0, 0, 3, 1), curve[childID])) != std::numeric_limits<size_t>::max()){
                        //very rare situation, third curve connect to the cross point of two curves.
                        //break the third curve and connect to the highest node.
                        if(nid == 0 || nid >= curve[childID].size() - 1) {
                            printf("please check the T style curve.\n");
                            return false;
                        }
                        curve.push_back(Line5d());
                        //curCurveIDList.push_back(lineList.size() - 1lu);// though push back, it still be deleted later
                        std::copy(curve[childID].begin() + nid, curve[childID].end(), std::back_inserter(curve.back()));
                        curve[childID].erase(curve[childID].begin() + nid + 1, curve[childID].end());
                        std::reverse(curve[childID].begin(), curve[childID].end());
                        newIterList.push_back(resTree.append_child(resTree.parent(it), LineID(childID, curNeuriteType)));
                        newIterList.push_back(resTree.append_child(resTree.parent(it), LineID(curve.size() - 1lu, curNeuriteType)));
                    }
                    else  continue;
                    curCurveIDListCp.erase(std::remove(curCurveIDListCp.begin(), curCurveIDListCp.end(), childID), curCurveIDListCp.end());
                }
                if (curCurveIDListCp.size() != curCurveIDList.size()) curCurveIDListCp.swap(curCurveIDList);
            }
            newIterList.swap(iterList);
            //KPTREE::print_tree_bracketed(resTree);
            //printf("\n");
            if (curCurveIDList.size() - oldsz == 0lu && curCurveIDList.size() > 0lu){
                /*for (size_t k = 0; k < curve.size(); ++k) {
                    if (k == curCurveIDList.front()) continue;
                    nid = FindNearestID(curve[curCurveIDList.front()].front().block(0,0,3,1), curve[k]);
                    if (nid != std::numeric_limits<size_t>::max()) {
                    NG_ERROR_MESSAGE("");
                    printf("nima\n");
                    }
                    nid = FindNearestID(curve[curCurveIDList.front()].back().block(0, 0, 3, 1), curve[k]);
                    if (nid != std::numeric_limits<size_t>::max()) {
                    NG_ERROR_MESSAGE("");
                    printf("nima\n");
                    }
                    }*/

				if (checkFlag) {
                    NG_ERROR_MESSAGE("error occurred in sorting dendrites: cannot iterative rebuild tree.");
                    //debug
       //             {
       //                 FILE* fp = fopen("F:/chenyijun/nima.swc", "w");
       //                 int id = 0;
       //                 auto &it = curve[curCurveIDList.front()];
       //                 for (auto &it : curve) {
							//auto &nodee = it[0];
							//fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, nodee(0), nodee(1), nodee(2));
       //                     for (auto &node : it) {
       //                         fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", ++id, node(0), node(1), node(2),id-1);
       //                     }
       //                 }
       //                 fclose(fp);
       //                 //system("pause");
       //             }
                    for (size_t k = 0; k < curve.size(); ++k) {
                        nid = FindNearestID(curve[curCurveIDList.front()].front().block(0, 0, 3, 1), curve[k], 5.0);
                        if (nid != std::numeric_limits<size_t>::max()) {
                            double dist = (curve[curCurveIDList.front()].front().block(0, 0, 3, 1) - curve[k][nid].block(0, 0, 3, 1)).norm();
                            printf("may connect to %d curve %d point, dist : %lf\n", k, nid, dist);
                        }
                        nid = FindNearestID(curve[curCurveIDList.front()].back().block(0, 0, 3, 1), curve[k], 5.0);
                        if (nid != std::numeric_limits<size_t>::max()) {
                            double dist = (curve[curCurveIDList.front()].front().block(0, 0, 3, 1) - curve[k][nid].block(0, 0, 3, 1)).norm();
                            printf("may connect to %d curve %d point, dist : %lf\n", k, nid, dist);
                        }
                    }
                    KPTREE::print_tree_bracketed(resTree);
                    system("pause");
                    //break;
                    return false;
                }
                else{
                    NG_ERROR_MESSAGE("just continue to read ignore the break lines.");
                    break;
                }
            }
        }
        topo = resTree;
        //print_tree_bracketed(topo);
        //conjunct continuous curves.
        if (refreshFlag) {
            if (castrationFlag)
                KPTREE::RefreshTree(nTree);
            else
                KPTREE::RefreshTreeCastration(nTree);
        }
        //
        return true;
    }

    bool SetTreeRoot(std::vector<Line5d>& curve, tree<LineID>& topo, const Vec3d& root, double threv /*= 0.5*/, bool flag /*= true*/, bool castrationFlag, bool refreshFlag)
    {
        NeuronTree tmp;
        tmp.m_curveList.swap(curve);
        tmp.m_Topo=topo;
        bool res = SetTreeRoot(tmp, root, threv, flag, castrationFlag, refreshFlag);
        tmp.m_curveList.swap(curve);
        topo = tmp.m_Topo;
        return res;
    }

    bool MergeTwoLine(NeuronTree& nTree, size_t mergeID1, size_t mergeID2)
    {
        tree<LineID> &topo = nTree.m_Topo;
        std::vector<Line5d> &allCurves = nTree.m_curveList;
        //assert whether the 1st line is in subtree of the 2nd
        auto it1 = std::find(topo.begin() + 1, topo.end(), LineID(mergeID1));
        auto it2 = std::find(topo.begin() + 1, topo.end(), LineID(mergeID2));
        if (nTree.m_Topo.is_in_subtree(it2, it1, nTree.m_Topo.next_sibling(it1))) {
            NG_ERROR_MESSAGE("line2 is in the subtree of merge line1");
            return false;
        }
        /*merge or connect   line1 to line2*/
        size_t deleteid = mergeID1;
        double dist = (allCurves[deleteid].front().block(0, 0, 3, 1) - allCurves[mergeID2].back().block(0, 0, 3, 1)).norm();
        if (dist < 2.0) {//merge 2 lines as one line
            if (dist < 1.0) std::copy(allCurves[deleteid].begin() + 1, allCurves[deleteid].end(), std::back_inserter(allCurves[mergeID2]));
            else std::copy(allCurves[deleteid].begin(), allCurves[deleteid].end(), std::back_inserter(allCurves[mergeID2]));
            allCurves.erase(allCurves.begin() + deleteid);
            //move all children of iter1 to iter2 and delete iter1
            //auto it1 = std::find(topo.begin() + 1, topo.end(), LineID(mergeID1));
            //auto it2 = std::find(topo.begin() + 1, topo.end(), LineID(mergeID2));
            topo.reparent(it2, it1);//jump root node
            auto delIt = std::find(topo.begin() + 1, topo.end(), LineID(mergeID1));// a guarentee
            if (delIt == topo.end()) {
                printf("error in MergeTwoLine: %d\n", __LINE__); return false;
            }
            topo.erase(delIt);
            //modify number of tree node
            KPTREE::UpdateDeleteTreeIDList(nTree, deleteid);
        }
        else{ // connect two lines
            size_t nearID = KPTREE::FindNearestID(allCurves[deleteid].front().block(0, 0, 3, 1), allCurves[mergeID2], 5.0);//connect far away lines
            if ((allCurves[deleteid].front().block(0, 0, 3, 1) - allCurves[mergeID2][nearID].block(0, 0, 3, 1)).norm() > 0.1) {
                allCurves[deleteid].insert(allCurves[deleteid].begin(), allCurves[mergeID2][nearID]);//add head node
            }
            auto piter1 = topo.parent(std::find(topo.begin() + 1, topo.end(), LineID(mergeID1)));
            auto piter2 = std::find(topo.begin() + 1, topo.end(), LineID(mergeID2));
            topo.reparent(piter2, piter1);
        }
        return true;
    }

    size_t FastSearchNearestLine(const Vec3d& pt, const std::vector<Line5d>& lines, const size_t interval, const double threv)
    {
        size_t vertID;
        return FastSearchNearestLine(pt, lines, vertID, interval, threv);
    }

    size_t FastSearchNearestLine(const Vec3d& pt, const std::vector<Line5d>& lines, size_t &vertID, const size_t interval/*= 3lu*/, const double threv/*= 10.0*/)
    {
        size_t nearestID = std::numeric_limits<size_t>::max();
        vertID = std::numeric_limits<size_t>::max();
        double minDist = 100000.0, tmpDist;
        //bool flag = false;
        for (size_t i = 0; i < lines.size(); ++i) {//set a rough search range
            if (lines[i].empty()) {
                printf("big bug !!!!! How could traced lines be empty?\n");
                return std::numeric_limits<size_t>::max();
            }
            //flag = false;
            const Line5d& curLine = lines[i];
            //if ((curLine[0].block(0, 0, 3, 1) - pt).array().abs().sum() < threv) {
            //    nearestID = i;
            //    vertID = 0;
            //    //return nearestID;
            //    //flag = true;
            //}
            //else if ((curLine.back().block(0, 0, 3, 1) - pt).array().abs().sum() < threv){
            //    nearestID = i;
            //    vertID = curLine.size()-1lu;
            //    //return nearestID;
            //}
            //else{
                for (size_t j = 0; j < curLine.size(); ++j) {
                    if (std::abs(curLine[j](0) - pt(0)) > threv || std::abs(curLine[j](1) - pt(1)) > threv || std::abs(curLine[j](2) - pt(2)) > threv)
                        continue;
                    tmpDist = (curLine[j].block(0, 0, 3, 1) - pt).array().abs().sum();
                    if (tmpDist < 1.0) {
                        nearestID = i;
                        vertID = j;
                        return nearestID;
                    }
                    if (tmpDist < threv && tmpDist < minDist) {
                        minDist = tmpDist;
                        nearestID = i;
                        vertID = j;
                        //return nearestID;
                    }
                }
            }
        //}
        return nearestID;
    }

    void FindSymmericDiff(const NeuronTree &tree1, const NeuronTree &tree2 , int rad, int threv, Line3d& resPts, std::vector<int> &resState)
    {
        resPts.clear();
        resState.clear();
        int disconNum = 0;
        bool disconFlag = false, tmpflag;
        Vec3d disconPt;
        //int conNum = 0, allNum = 0;
        auto topo2 = tree2.m_Topo;//must copy
        auto topo1 = tree1.m_Topo;
        auto &curveList1 = tree1.m_curveList;
        auto &curveList2 = tree2.m_curveList;
        auto &hash1 = tree1.m_hashList;
        auto &hash2 = tree2.m_hashList;
        auto &preSusPoints = tree1.m_suspiciousPositions_;
        auto &preSusState = tree1.m_suspiciousCheckState_;
        //along the topology
        for (auto iter = topo1.begin() + 1; iter != topo1.end(); ++iter) {
            if (iter->id >= 1000000) continue;
            disconNum = 0;
            disconFlag = false;
            auto &curCurve = curveList1[iter->id];
            size_t curID = 0lu;
            for (size_t i = 0; i < curCurve.size(); ++i) {
                auto &v = curCurve[i];
                tmpflag = !hash2.isExistInRange(Vec3i(int(v(0)), int(v(1)), int(v(2))), rad);
                if (tmpflag){//do not registration
                    if (!disconFlag) disconPt = v.block(0, 0, 3, 1);
                    ++disconNum;
                    disconFlag = true;
                    curID = i;
                }
                else if (!tmpflag) {
                    if (!disconFlag) {
                        continue;
                    }
                    if (disconNum > threv){
                        resPts.push_back(disconPt);
                    }
                    else{
                        //conNum += disconNum;
                    }
                    disconFlag = false;
                    disconNum = 0;
                }
                if (disconNum > threv)
                {
                    //search previous suspicious points and check states. If tolerable, then continue compare. If new, checked or no then stop
                    //preSusPoints preSusState
                    resPts.push_back(disconPt);
                    auto stateIt = std::find_if(preSusPoints.begin(), preSusPoints.end(), [&](const Vec3d& arg){return (arg - disconPt).norm() < rad; });
                    int stateFlag = 0;
                    if ( stateIt != preSusPoints.end()) {
                        size_t dist = std::distance(preSusPoints.begin(), stateIt);
                        stateFlag = preSusState[dist];
                    }
                    if (stateFlag == 2) {//tolerable
                        resState.push_back(stateFlag);
                        continue; //continue to compare branch
                    }
                    resState.push_back(stateFlag);
                    disconFlag = false;
                    disconNum = 0;
                    //remove downstream curve
                    if (iter.number_of_children() > 0) {
                        auto subtree = topo1.subtree(iter, topo1.next_sibling(iter));
                        //calculate connect node id
                        for (auto subIt = subtree.begin(); subIt != subtree.end(); ++subIt){
                            auto &conCurve = curveList1[subIt->id];
                            auto &headPt = conCurve.front();
                            size_t conId = 0;
                            for (size_t k = 0; k < curCurve.size(); ++k) {
                                if ( (headPt.block(0,0,3,1) - curCurve[k].block(0,0,3,1)).norm() < 1.0) {
                                    conId = k;
                                    break;
                                }
                            }//for
                            if (conId > curID) {
                                auto findIt = std::find(topo1.begin(), topo1.end(), LineID(subIt->id));
                                auto subsubtree = topo1.subtree(findIt, topo1.next_sibling(findIt));
                                for (auto subiter = subsubtree.begin(); subiter != subsubtree.end(); ++subiter) {
                                    auto modifyIter = std::find(topo1.begin() + 1, topo1.end(), LineID(subiter->id));
                                    modifyIter->id += 1000000;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
        //std::for_each(topo1.begin() + 1, topo1.end(), [](tree<LineID>::pre_order_iterator it){ if (it->id > 1000000) it->id -= 1000000;});
        //printf("conNum:%d\n  allNum:%d\n", conNum, allNum);
        for (auto iter = topo2.begin() + 1; iter != topo2.end(); ++iter) {
            if (iter->id >= 1000000) continue;
            disconNum = 0;
            disconFlag = false;
            auto &curCurve = curveList2[iter->id];
            size_t curID = 0lu;
            for (size_t i = 0; i < curCurve.size(); ++i) {
                auto &v = curCurve[i];
                //int a = v(0); int b = v(1);
                //++allNum;
                tmpflag = !hash1.isExistInRange(Vec3i(int(v(0)), int(v(1)), int(v(2))), rad);
                if (tmpflag){//do not registration
                    if (!disconFlag) disconPt = v.block(0, 0, 3, 1);
                    ++disconNum;
                    disconFlag = true;
                    curID = i;
                }
                else if (!tmpflag) {
                    if (!disconFlag) {
                        //++conNum;
                        continue;
                    }
                    if (disconNum > threv){
                        resPts.push_back(disconPt);
                    }
                    else{
                        //conNum += disconNum;
                    }
                    disconFlag = false;
                    disconNum = 0;
                }
                if (disconNum > threv)
                {
                    //search previous suspicious points and check states. If tolerable, then continue compare. If new, checked or no then stop
                    //preSusPoints preSusState
                    resPts.push_back(disconPt);
                    auto stateIt = std::find_if(preSusPoints.begin(), preSusPoints.end(), [&](const Vec3d& arg){return (arg - disconPt).norm() < rad; });
                    int stateFlag = 0;
                    if (stateIt != preSusPoints.end()) {
                        size_t dist = std::distance(preSusPoints.begin(), stateIt);
                        stateFlag = preSusState[dist];
                    }
                    if (stateFlag == 2) {//tolerable
                        resState.push_back(stateFlag);
                        continue; //continue to compare branch
                    }
                    resState.push_back(stateFlag);
                    disconFlag = false;
                    disconNum = 0;
                    //remove downstream curve
                    if (iter.number_of_children() > 0) {
                        auto subtree = topo2.subtree(iter, topo2.next_sibling(iter));
                        //calculate connect node id
                        for (auto subIt = subtree.begin()+1; subIt != subtree.end(); ++subIt){//the root is parent node
                            auto &conCurve = curveList2[subIt->id];
                            auto &headPt = conCurve.front();
                            size_t conId = 0;
                            for (size_t k = 0; k < curCurve.size(); ++k) {
                                if ((headPt.block(0, 0, 3, 1) - curCurve[k].block(0, 0, 3, 1)).norm() < 1.0) {
                                    conId = k;
                                    break;
                                }
                            }//for
                            if (conId > curID) {
                                auto findIt = std::find(topo2.begin(), topo2.end(), LineID(subIt->id));
                                auto subsubtree = topo2.subtree(findIt, topo2.next_sibling(findIt));
                                for (auto subiter = subsubtree.begin(); subiter != subsubtree.end(); ++subiter) {
                                    auto modifyIter = std::find(topo2.begin() + 1, topo2.end(), LineID(subiter->id));
                                    modifyIter->id += 1000000;
                                }
                            }
                        }
                    }
                    break;
                }
            }
        }
        //std::for_each(topo2.begin() + 1, topo2.end(), [](tree<LineID>::pre_order_iterator it){ if (it->id > 1000000) it->id -= 1000000; });
        //remove repeat point
        for (size_t i = 0; i < resPts.size()-1; ++i) {
            if (resPts[i](0) < 1.0) continue;
            for (size_t j = i + 1; j < resPts.size(); ++j) {
                if (resPts[j](0) < 1.0) continue;
                if ( (resPts[i]-resPts[j]).norm() < 5.0) {
                    resPts[j] << 0, 0, 0;
                    resState[j] = -1;
                }
            }
        }
        resPts.erase(std::remove_if(resPts.begin(), resPts.end(), [](Vec3d &arg){return arg(0) < 1.0; }), resPts.end());
        resState.erase(std::remove(resState.begin(), resState.end(), -1), resState.end());
    }
    void FindCurrentImagePoints(const std::vector<Line5d> &lines, std::vector<Line5d> &Currentlines, std::vector<VectorVec2i> &CurrentIDs,
        size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max)
    {
        Currentlines.clear();
        Line5d tmp;
        VectorVec2i points;
        for (size_t i = 0; i < lines.size(); ++i){
            tmp.clear();
            points.clear();
            for (size_t j = 0; j < lines[i].size(); ++j){

                if (lines[i][j][0] >= x_min&&lines[i][j][0] <= x_max&&lines[i][j][1] >= y_min&&lines[i][j][1] <= y_max&&lines[i][j][2] >= z_min&&lines[i][j][2] <= z_max){
                    tmp.push_back(lines[i][j]);
                    Vec2i point;
                    point[0] = int(i);
                    point[1] = int(j);
                    points.push_back(point);
                }
            }
            if (!points.empty()){
                CurrentIDs.push_back(points);
                Currentlines.push_back(tmp);
            }
        }
    }
    void FindCurrentImagePoints(const std::vector<Line5d> &lines, const std::vector<size_t> curBranchID, std::vector<Line5d> &Currentlines, std::vector<VectorVec2i> &CurrentIDs,
        size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max)
    {
        Currentlines.clear();
        Line5d tmp;
        VectorVec2i points;
        for each(size_t i in curBranchID){
            tmp.clear();
            points.clear();
            for (size_t j = 0; j < lines[i].size(); ++j){

                if (lines[i][j][0] >= x_min&&lines[i][j][0] <= x_max&&lines[i][j][1] >= y_min&&lines[i][j][1] <= y_max&&lines[i][j][2] >= z_min&&lines[i][j][2] <= z_max){
                    tmp.push_back(lines[i][j]);
                    Vec2i point;
                    point[0] = int(i);
                    point[1] = int(j);
                    points.push_back(point);
                }
            }
            if (points.size() > 3){
                CurrentIDs.push_back(points);
                Currentlines.push_back(tmp);
            }

        }
    }

    void InsertCurrentLine(std::vector<Line5d> &lines, const CellVec3d &InsetLines, std::vector<VectorVec2i> &CurrentIDs){
        assert(InsetLines.size() == CurrentIDs.size());
        for (size_t i = 0; i < InsetLines.size(); ++i){
            //assert(InsetLines[i].size()>= CurrentIDs.size());
            Line5d newline;
            size_t begin = CurrentIDs[i][0][1];
            size_t nx = CurrentIDs[i][0][0];
            size_t index = 0;
            for (index = 0; index < begin; ++index)
            {
                newline.push_back(lines[nx][index]);
            }
            for (size_t j = 0; j < InsetLines[i].size(); ++j){
                Vec5d tmp;
                tmp[0] = InsetLines[i][j][0];
                tmp[1] = InsetLines[i][j][1];
                tmp[2] = InsetLines[i][j][2];
                tmp[3] = 1;
                tmp[4] = 0;
                newline.push_back(tmp);
            }
            for (index += CurrentIDs[i].size(); index < lines[nx].size(); ++index)
            {
                newline.push_back(lines[nx][index]);
            }
            lines[nx].swap(newline);
        }

    }

    size_t FindNearestID2(const Vec3d& pt, const Line5d& line, double &mindist){
        size_t id = std::numeric_limits<size_t>::max();
        mindist = 100000;
        double curDist;
        for (size_t i = 0; i < line.size(); ++i) {
            curDist = std::abs(pt(0) - line[i](0)) + std::abs(pt(1) - line[i](1)) + std::abs(pt(2) - line[i](2));
            if (curDist < mindist) {
                mindist = curDist;
                id = i;
                if (curDist < 1.0) return id;
            }
        }
        return id;
    }

    bool ConjunctChildCurve(NeuronTree &tree, size_t id)
    {
        auto &topo = tree.m_Topo;
        auto topoIt = std::find(topo.begin() + 1, topo.end(), LineID(id));
        if (topo.number_of_children(topoIt) != 0lu) {
            bool conjFlag = false;
            auto childIt = topo.child(topoIt, 0);
            if (childIt.node == NULL) return false;
            std::vector<Line5d> &oldLines = tree.m_curveList;
            Line5d& curLine = oldLines[topoIt->id];
            Vec3d curLineEnd = curLine.back().block(0, 0, 3, 1);
            for (; childIt != topo.end() && childIt.node != NULL; childIt = topo.next_sibling(childIt)) {
                Line5d& childLine = oldLines[childIt->id];
                if ((childLine[0].block(0, 0, 3, 1) - curLineEnd).norm() < 0.01) {
                    size_t id = childIt->id;
                    std::vector<size_t> deleteID; deleteID.push_back(id);
                    //modify lines
                    curLine.reserve(curLine.size() + childLine.size());
                    std::copy(childLine.begin(), childLine.end(), std::back_inserter(curLine));
                    oldLines.erase(oldLines.begin() + id);
                    //modify topo
                    auto subtree = topo.subtree(childIt, topo.next_sibling(childIt));
                    topo.reparent(topoIt, subtree.begin());
                    topo.erase(childIt);
                    KPTREE::UpdateDeleteTreeIDList(tree, deleteID);
                    conjFlag = true;
                    break;
                }
            }
            return conjFlag;
        }
        else return false;
    }

    int RefreshTreeCastration(NeuronTree& tree, int xMax, int yMax, int zMax, READMODE readmode)
    {
        auto &topo = tree.m_Topo;
        //conjunct all continuous curves
        int changeStatus = 0;
        bool conjFlag = true, statusFlag = false;
        while (conjFlag){
            auto treeNode = topo.begin() + 1;
            statusFlag = false;
            for (; treeNode != topo.end() && treeNode.node != NULL; ++treeNode) {
                if (KPTREE::ConjunctChildCurve(tree, treeNode->id)){
                    statusFlag = true;
                    ++changeStatus;
                    break;
                }
            }
            if (statusFlag) continue;
            else break;
        }
        //remove redundant point, not include connect points.
        {
            auto &allCurve = tree.m_curveList;
            //bool removeFlag = false;
            double xM = double(xMax);
            double yM = double(yMax);
            double zM = double(zMax);
            for (auto &it : allCurve) {
                size_t len;
                /*removeFlag = false;
                len = it.size() - 1llu;
                for (size_t i = 0; i < len; ++i){
                if (std::abs(it[i](0) - it[i + 1](0)) < 0.1 && std::abs(it[i](1) - it[i + 1](1)) < 0.1 && std::abs(it[i](2) - it[i + 1](2)) < 0.1) {
                it[i + 1](0) = 0.0;
                removeFlag = true;
                }
                }
                if (removeFlag)
                it.erase(std::remove_if(it.begin(), it.end(), [](const Vec5d &arg){ return std::abs(arg(0) - 0.0) < 0.1; }), it.end());*/
                //boundary limit
                if (readmode == READMODE::TIFF3D) {
                    len = it.size();
                    if (xM > 0 && yM > 0 && zM > 0) {
                        for (size_t i = 0; i < len; ++i){
                            if (it[i](0) < 0) it[i](0) = 0.0;
                            if (it[i](0) > xM) it[i](0) = xM;
                            if (it[i](1) < 0) it[i](1) = 0.0;
                            if (it[i](1) > yM) it[i](1) = yM;
                            if (it[i](2) < 0) it[i](2) = 0.0;
                            if (it[i](2) > zM) it[i](2) = zM;
                        }
                    }
                    else{
                        for (size_t i = 0; i < len; ++i){
                            if (it[i](0) < 0) it[i](0) = 0.0;
                            if (it[i](1) < 0) it[i](1) = 0.0;
                            if (it[i](2) < 0) it[i](2) = 0.0;
                        }
                    }
                }

            }
        }
        return changeStatus;
    }

    int RefreshTree(NeuronTree& tree, int xMax , int yMax , int zMax, READMODE readmode )
    {
        auto &topo = tree.m_Topo;
        //conjunct all continuous curves
        int changeStatus = 0;
        bool conjFlag = true, statusFlag = false;
        while (conjFlag){
            auto treeNode = topo.begin() + 1;
            statusFlag = false;
            for (; treeNode != topo.end() && treeNode.node != NULL; ++treeNode) {
                if (KPTREE::ConjunctChildCurve(tree, treeNode->id)){
                    statusFlag = true;
                    ++changeStatus;
                    break;
                }
            }
            if (statusFlag) continue;
            else break;
        }
        //remove redundant point, not include connect points.
        {
            auto &allCurve = tree.m_curveList;
            bool removeFlag = false;
            double xM = double(xMax);
            double yM = double(yMax);
            double zM = double(zMax);
            for (auto &it : allCurve) {
                removeFlag = false;
                size_t len = it.size() - 1llu;
                for (size_t i = 0; i < len; ++i){
                    if (std::abs(it[i](0) - it[i + 1](0)) < 0.1 && std::abs(it[i](1) - it[i + 1](1)) < 0.1 && std::abs(it[i](2) - it[i + 1](2)) < 0.1) {
                        it[i + 1](0) = 0.0;
                        removeFlag = true;
                    }
                }
                if (removeFlag)
                    it.erase(std::remove_if(it.begin(), it.end(), [](const Vec5d &arg){ return std::abs(arg(0) - 0.0) < 0.1; }), it.end());
                //boundary limit
                if (readmode == READMODE::TIFF3D) {
                    len = it.size();
                    if (xM > 0 && yM > 0 && zM > 0) {
                        for (size_t i = 0; i < len; ++i){
                            if (it[i](0) < 0) it[i](0) = 0.0;
                            if (it[i](0) > xM) it[i](0) = xM;
                            if (it[i](1) < 0) it[i](1) = 0.0;
                            if (it[i](1) > yM) it[i](1) = yM;
                            if (it[i](2) < 0) it[i](2) = 0.0;
                            if (it[i](2) > zM) it[i](2) = zM;
                        }
                    }
                    else{
                        for (size_t i = 0; i < len; ++i){
                            if (it[i](0) < 0) it[i](0) = 0.0;
                            if (it[i](1) < 0) it[i](1) = 0.0;
                            if (it[i](2) < 0) it[i](2) = 0.0;
                        }
                    }
                }
                
            }
        }
        return changeStatus;
    }

    NEURITETYPE Num2NeuriteType(int num)
    {
        switch (num)
        {
        case 1:
            return NEURITETYPE::SOMA;
        case 2:
            return NEURITETYPE::AXON;
        case 3:
            return NEURITETYPE::DENDRITE;
        case 4:
            return NEURITETYPE::APICAL;
        default:
            return NEURITETYPE::UNKNOWN;
            break;
        }
    }

    int NeuriteType2Num(NEURITETYPE arg)
    {
        switch (arg)
        {
        case NEURITETYPE::SOMA:
            return 1;
            break;
        case NEURITETYPE::DENDRITE:
            return 3;
            break;
        case NEURITETYPE::APICAL:
            return 4;
            break;
        case NEURITETYPE::AXON:
            return 2;
            break;
        case NEURITETYPE::UNKNOWN:
            return 0;
            break;
        default:
            return 0;
            break;
        }
    }


}


void NeuronTree::BuildHash()
{
    m_hashList.Reset();
    m_hashList.insert(this->m_curveList);
}

void NGHashList::insert(const Vec3i& arg)
{
    _insert(arg);
    std::sort(v[mod_hash(arg(0))][mod_hash(arg(1))]->begin(), v[mod_hash(arg(0))][mod_hash(arg(1))]->end(), sortFunc);
}


void NGHashList::insert(const std::vector<Line5d> &arg)
{
    //int n;
    std::set<std::pair<int, int> > modifyList;
    Vec3i tmp;
    int sz = int(arg.size());
    for (int i = 0; i < sz; ++i){
        auto & it1 = arg[i];
        for (auto &it2 : it1){
            tmp << NGUtility::Round((std::max)(it2(0), 0.0)), NGUtility::Round((std::max)(it2(1), 0.0)), NGUtility::Round((std::max)(it2(2), 0.0));
            _insert(tmp);
            modifyList.insert(std::make_pair(mod_hash(tmp(0)), mod_hash(tmp(1))));
        }
    }
    for (auto &it : modifyList)  std::sort(v[it.first][it.second]->begin(), v[it.first][it.second]->end(), sortFunc);
}

bool NGHashList::isExistInRange(const Vec3i &val, int thre /*= 2*/) const
{
    int xMin = std::max(val(0) - thre, 0); int xMax = val(0) + thre;
    int yMin = std::max(val(1) - thre, 0); int yMax = val(1) + thre;
    Vec3i minVec = -thre + val.array(); minVec = (minVec.array() > -1).select(minVec, 0); Vec3i maxVec = thre + val.array();
    for (int i = xMin; i <= xMax; ++i){
        VectorVec3i** cv = v[mod_hash(i)];
        if (cv == NULL) continue;
        for (int j = yMin; j <= yMax; ++j) {
            VectorVec3i *ccv = cv[mod_hash(j)];
            if (ccv == NULL) continue;
            auto lb3 = std::lower_bound(ccv->begin(), ccv->end(), minVec, comp3);
            auto hb3 = std::upper_bound(ccv->begin(), ccv->end(), maxVec, comp3);
            if (std::distance(lb3, hb3) > 0){
                return true;
            }
        }
    }
    return false;
}


void NGHashList::Clear()
{
    if (NULL == v) return;
    for (int i = 0; i < modn; ++i) {
        if (v[i] == NULL) continue;
        for (int j = 0; j < modn; ++j) {
            if (v[i][j] == NULL) continue;
            v[i][j]->clear();
            delete v[i][j];
            v[i][j] = NULL;
        }
        delete[] v[i];
        v[i] = NULL;
    }
    delete[] v;
    v = NULL;
}

void NGHashList::Reset()
{
    for (int i = 0; i < modn; ++i) {
        if (v[i] == NULL) continue;
        for (int j = 0; j < modn; ++j){
            if (v[i][j] == NULL) continue;
            v[i][j]->clear();
        }
    }
}

void NGHashList::_insert(const Vec3i& arg)
{
    int i = mod_hash(arg(0));
    int j = mod_hash(arg(1));
    if (v[i] == NULL){
        v[i] = new VectorVec3i *[modn];
        for (int k = 0; k < modn; ++k)
            v[i][k] = NULL;
    }
    if (v[i][j] == NULL) v[i][j] = new VectorVec3i;
    v[i][j]->push_back(arg);
}
