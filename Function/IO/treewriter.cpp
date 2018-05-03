#include <stdio.h>
#include <stdlib.h>
#include "treewriter.h"
#include "../../ngtypes/tree.h"
#include "../../Function/NGUtility.h"
#pragma warning(disable:4996)

TreeWriter::TreeWriter()
{
    className_ = std::string("TreeWriter");
}

TreeWriter::~TreeWriter()
{

}


ProcStatPointer TreeWriter::Update()
{
    if (!m_Input || m_Input->GetIdentifyName() != std::string("NeuronPopulation")){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if (m_Input->GetProcessObject()){
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    if (!param) {
        printf("please set param.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "please set param.");
        return resSta;
    }
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Input);
    if (tmpTree->m_pop.size() != fileName.size()) {
        printf("error in treewriter: the save name list cannot correspond to the trees.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "the save name list cannot correspond to the trees.");
        return resSta;
    }
    for (size_t i = 0; i < tmpTree->m_pop.size(); ++i) {//traverse tree list
        //traverse tree curves
        int changNum=0;
        if (param->OrigImage) {
            NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, param->OrigImage);
            changNum = KPTREE::RefreshTreeCastration(*(tmpTree->m_pop[i]), tmpOrig->x() - 1, tmpOrig->y() - 1, tmpOrig->z() - 1, param->dataMode_);
        }
        else changNum = KPTREE::RefreshTreeCastration(*(tmpTree->m_pop[i]));
        printf("tree %d has been refresh for %d places before saving.\n", i, changNum);
        if (!SaveATree(*(tmpTree->m_pop[i]), fileName[i])){
            printf("error in treewriter, loop : %llu\n", i);
            char line[256];
            sprintf(line, "error in treewriter, loop : %llu", i);
            MAKEPROCESSSTATUS(resSta, false, className_, line);
            return resSta;
        }
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

bool TreeWriter::SaveATree(NeuronTree& nTree, std::string& path)
{
    if (nTree.m_curveList.empty()) {
        printf("error in SaveATree. empty tree list.\n");
        return false;
    }
    
    FILE* fp = fopen(path.c_str(), "w");
    if (!fp) return false;
    int id = 0;
    size_t curLineID, parentLineID, nearestID, pid;
    std::vector<std::vector<int>> idList(nTree.m_curveList.size());//save 1 column id
    for (size_t i = 0; i < idList.size(); ++i) 
        idList[i].resize(nTree.m_curveList[i].size());//warning !!! : VS2013 set default value as 0
    //1st line
    {
        curLineID = (nTree.m_Topo.begin() + 1)->id;
        const Line5d &line = nTree.m_curveList[curLineID];
        if (line.empty()){
            printf("error in SaveATree. empty 1st tree.\n");
            return false;
        }
        idList[curLineID][0] = ++id;
        fprintf(fp, "%d %d %lf %lf %lf 1.0 -1\n", id, 1, 
            line[0](0) *param->xRes_, line[0](1) *param->yRes_, line[0](2) *param->zRes_);//root
        int curType = KPTREE::NeuriteType2Num((nTree.m_Topo.begin() + 1)->type);
        for (size_t i = 1; i < line.size(); ++i){
            idList[curLineID][i] = ++id;
            fprintf(fp, "%d %d %lf %lf %lf 1.0 %d\n", id, curType, 
                line[i](0) *param->xRes_, line[i](1) *param->yRes_, line[i](2) *param->zRes_, id - 1);
        }
    }
    //left lines
    for (auto it = nTree.m_Topo.begin() + 2; it != nTree.m_Topo.end(); ++it) {
        curLineID = it->id;
        const Line5d &line = nTree.m_curveList[curLineID];
        if (nTree.m_Topo.parent(it) == nTree.m_Topo.begin()){
            //set the pid as 1st node
            idList[curLineID][0] = ++id;
            int curType = KPTREE::NeuriteType2Num(it->type);
            fprintf(fp, "%d %d %lf %lf %lf 1.0 1\n", id, curType,
                line[0](0) *param->xRes_, line[0](1) *param->yRes_, line[0](2) *param->zRes_);
            for (size_t i = 1; i < line.size(); ++i){
                idList[curLineID][i] = ++id;
                fprintf(fp, "%d %d %lf %lf %lf 1.0 %d\n", id, curType,
                    line[i](0) *param->xRes_, line[i](1) *param->yRes_, line[i](2) *param->zRes_, id - 1);
            }
        }
        else{//find parent line-vertex id
            parentLineID = nTree.m_Topo.parent(it)->id;
            const Vec3d &headNode = line.front().block(0, 0, 3, 1);
			nearestID = KPTREE::FindNearestID(headNode, nTree.m_curveList[parentLineID], 1);
            //if (nearestID > 1000000) {
				//nearestID = KPTREE::FindNearestID(headNode, nTree.m_curveList[parentLineID], 10);
			if (nearestID > 1000000) {
				NG_ERROR_MESSAGE("how can pid so large? please debug.");
				fclose(fp);
				return false;
			}
            pid = idList[parentLineID][nearestID];
            //TODO: debug
            if (pid > 1000000) {
                NG_ERROR_MESSAGE("how can pid so large? please debug.");
            }
            if (pid == 0) {
                printf("error in SaveATree. the parent vertex id is 0, please debug.\n");
                return false;
            }
            idList[curLineID][0] = ++id;
            int curType = KPTREE::NeuriteType2Num(it->type);
            fprintf(fp, "%d %d %lf %lf %lf 1.0 %d\n", id, curType,
                line[0](0) *param->xRes_, line[0](1) *param->yRes_, line[0](2) *param->zRes_, pid);
            for (size_t i = 1; i < line.size(); ++i){
                idList[curLineID][i] = ++id;
                fprintf(fp, "%d %d %lf %lf %lf 1.0 %d\n", id, curType,
                    line[i](0) *param->xRes_, line[i](1) *param->yRes_, line[i](2) *param->zRes_, id - 1);
            }
        }
    }
    fclose(fp);
    return true;
}



bool TreeWriter::SetOutputFileName(const std::vector<std::string> & str)
{
    if (str.empty()) return false;
    fileName = str;
    return true;
}

ProcStatPointer TreeWriter::OldUpdate()
{
    if(!m_Input || m_Input->GetIdentifyName() != std::string("Tree")){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if(m_Input->GetProcessObject()){
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    //output dd1 into many text and dd2 in one text
#ifdef _WIN32
    std::tr1::shared_ptr<const TreeCurve> tmptree = std::tr1::dynamic_pointer_cast<const TreeCurve>(m_Input);
#else
    std::shared_ptr<const TreeCurve> tmptree = std::dynamic_pointer_cast<const TreeCurve>(m_Input);
#endif
    if(tmptree && dd2ptr){
        std::string testPath("/home/zhouhang/lijing.txt");
        FILE *fp = fopen(testPath.c_str(), "w");
        if(fp){
            for(int i = 0; i < tmptree->size(); ++i){
                const VectorVec5d& cur = tmptree->GetCurve()[i];
                fprintf(fp,"%lf %lf %lf %lf %lf -1\n", cur[0](0), cur[0](1), cur[0](2),
                        cur[0](3), cur[0](4) );
                for(VectorVec5d::size_type j = 1; j < cur.size(); ++j){
                    fprintf(fp,"%lf %lf %lf %lf %lf 0\n", cur[j](0), cur[j](1), cur[j](2),
                            cur[j](3), cur[j](4) );
                }
            }
            fclose(fp);
        }
        
        testPath = "/home/zhouhang/gcd2.txt";
        FILE *fp1 = fopen(testPath.c_str(), "w");
        if(fp1){
            for(VectorMat2i::size_type i = 0; i < dd2ptr->size(); ++i){
                fprintf(fp1, "%d %d %d %d\n", dd2ptr->operator [](i)(0,0), dd2ptr->operator [](i)(0,1),
                        dd2ptr->operator [](i)(1,0), dd2ptr->operator [](i)(1,1) );
            }
            fclose(fp1);
        }
        MAKEPROCESSSTATUS(resSta, true, className_, "");
        return resSta;
    }
    printf("error occurred in %s\n", className_.c_str());
    MAKEPROCESSSTATUS(resSta, false, className_, "data is wrong.");
    return resSta;
}

void TreeWriter::SetDD2(const VectorMat2i &dd2)
{
    dd2ptr = &dd2;
}

ProcStatPointer TreeWriter::UpdateLayer(const std::vector<int>& layerDisplay, const std::vector<int>& branchDisplay)
{
    if (!m_Input || m_Input->GetIdentifyName() != std::string("NeuronPopulation")){
        printf("error occurred in %s\n", className_.c_str());
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
    if (m_Input->GetProcessObject()){
        ProcStatPointer res = m_Input->GetProcessObject()->Update();
        if (!res->success()){
            printf("error occurred in %s\n", className_.c_str());
            //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
            return res;
        }
    }
    if (!param) {
        printf("please set param.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "Please set parameter.");
        return resSta;
    }
    if (fileName.size() > 1lu) {
        NG_ERROR_MESSAGE("just save one tree.");
        MAKEPROCESSSTATUS(resSta, false, className_, "just save one tree.");
        return resSta;
    }
    if (!param->activeTree) {
        NG_ERROR_MESSAGE("there is no active tree.");
        MAKEPROCESSSTATUS(resSta, false, className_, "there is no active tree.");
        return resSta;
    }
    auto &activeTree = *(param->activeTree);
    auto &curTopo = activeTree.m_Topo;
    auto &curCurveList = activeTree.m_curveList;
    FILE* fp = fopen(fileName[0].c_str(), "w");
    if (!fp) return false;
    //find correspond curves
    std::vector<size_t> levelCurveID;
    std::vector<decltype(curTopo.begin())> firstLevel;
    for (auto it = curTopo.begin() + 1; it != curTopo.end(); ++it) {//
        if (layerDisplay[0] == curTopo.depth(it)){
            firstLevel.push_back(it);
        }
    }
    for (size_t i = 0; i < branchDisplay.size(); ++i) {
        if (branchDisplay[i] >= firstLevel.size()) break;
        auto subTree = curTopo.subtree(firstLevel[branchDisplay[i]], curTopo.next_sibling(firstLevel[branchDisplay[i]]));//level minus layerList [0]
        for (auto iter = subTree.begin(); iter != subTree.end(); ++iter)
            if (NGUtility::IsInVector(layerDisplay, curTopo.depth(iter) + layerDisplay[0]))
                levelCurveID.push_back(iter->id);
    }
    //write file
    int id = 0;// , pid = -1;
    for (auto &curID : levelCurveID) {
        auto &curve = curCurveList[curID];
        fprintf(fp, "%d %d %lf %lf %lf %lf %d\n", ++id, 0, curve[0](0), curve[0](1), curve[0](2), 1.0, -1);
        for (size_t k = 1; k < curve.size(); ++k) {
            ++id;
            fprintf(fp, "%d %d %lf %lf %lf %lf %d\n", id, 0, curve[k](0), curve[k](1), curve[k](2), 1.0, id - 1);
        }
    }
    fclose(fp);
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

