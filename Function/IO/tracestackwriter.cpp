/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include <cstdlib>
#include "tracestackwriter.h"
#include "../../ngtypes/ParamPack.h"
#include "../../ngtypes/tree.h"
#include "inifile.h"
#pragma warning(disable:4996)

TraceStackWriter::TraceStackWriter()
{
    className_ = std::string("TraceStackWriter");
}

TraceStackWriter::~TraceStackWriter()
{

}

bool TraceStackWriter::SetOutputFileName(const std::string &path)
{
    fileName = path;
    return true;
}

ProcStatPointer TraceStackWriter::Update()
{
    //FILE* fp = fopen(fileName.c_str(), "w");
    //if(!fp) {
   //     printf("error occurred in %s\n", identifyName.c_str());
   //     return false;
    //}
    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Input);
    if (!tmpTree) {
        NG_ERROR_MESSAGE("seperateTree is null.");
        MAKEPROCESSSTATUS(resSta, false, className_, "no input tree data.");
        return resSta;
    }
    inifile::IniFile ini;
    char line[256];
    sprintf_s(line, "%d", tmpTree->m_pop.size());
    ini.setValue("common", "tree_num", line);
    sprintf_s(line, "%d", paramPack->xMin_);
    ini.setValue("common", "xMin", line);
    sprintf_s(line, "%d", paramPack->xMax_);
    ini.setValue("common", "xMax", line);
    sprintf_s(line, "%d", paramPack->yMin_);
    ini.setValue("common", "yMin", line);
    sprintf_s(line, "%d", paramPack->yMax_);
    ini.setValue("common", "yMax", line);
    sprintf_s(line, "%d", paramPack->zMin_);
    ini.setValue("common", "zMin", line);
    sprintf_s(line, "%d", paramPack->zMax_);
    ini.setValue("common", "zMax", line);
    sprintf_s(line, "%d", paramPack->xExtract_);
    ini.setValue("common", "x_Extract", line);
    sprintf_s(line, "%d", paramPack->yExtract_);
    ini.setValue("common", "y_Extract", line);
    sprintf_s(line, "%d", paramPack->zExtract_);
    ini.setValue("common", "z_Extract", line);
    sprintf_s(line, "%lf", paramPack->binThreshold_);
    ini.setValue("common", "Tree_Binary_Value", line);
    sprintf_s(line, "%lf", paramPack->axonBinaryThreshold_);
    ini.setValue("common", "Axon_Binary_Value", line);
    sprintf_s(line, "%lf", paramPack->diffuseValue_);
    ini.setValue("common", "Tree_Diffuse_Value", line);
    sprintf_s(line, "%lf", paramPack->traceValue_);
    ini.setValue("common", "Tree_Trace_Value", line);
    sprintf_s(line, "%lf", paramPack->axonDiffuseValue_);
    ini.setValue("common", "Axon_Diffuse_Value", line);
    sprintf_s(line, "%lf", paramPack->axonTraceValue_);
    ini.setValue("common", "Axon_Trace_Value", line);
    sprintf_s(line, "%d", paramPack->maxBoundNum_);
    ini.setValue("common", "Max_Bound_Num", line);
    sprintf_s(line, "%d", paramPack->lowOpac_);
    ini.setValue("common", "Opacity_Low", line);
    sprintf_s(line, "%d", paramPack->highOpac_);
    ini.setValue("common", "Opacity_High", line);
    sprintf_s(line, "%d", paramPack->xScale_);
    ini.setValue("common", "XScale", line);
    sprintf_s(line, "%d", paramPack->yScale_);
    ini.setValue("common", "YScale", line);
    sprintf_s(line, "%d", paramPack->lowAdjustOpac_);
    ini.setValue("common", "Opacity_Adjust_Low", line);
    sprintf_s(line, "%d", paramPack->highAdjustOpac_);
    ini.setValue("common", "Opacity_Adjust_High", line);
    sprintf_s(line, "%d", paramPack->lowDestOpac_);
    ini.setValue("common", "Opacity_Dest_Low", line);
    sprintf_s(line, "%d", paramPack->highDestOpac_);
    ini.setValue("common", "Opacity_Dest_High", line);
    sprintf_s(line, "%llu", paramPack->runningMinutes_);
    ini.setValue("common", "Running_Minutes", line);
    char section[128];
    char arrow[128];
    char boundString[256];
    char susStr[128];
    for (size_t i = 0; i < tmpTree->m_pop.size(); ++i) {
        const NeuronTree &nTree = *(tmpTree->m_pop[i]);
        sprintf_s(section, "tree%d", i + 1);
        sprintf_s(line, "%s", swcName[i].c_str());
        ini.setValue(section, "path", line);
        sprintf_s(line, "%lf %lf %lf", nTree.m_curInitPt(0), nTree.m_curInitPt(1), nTree.m_curInitPt(2));
        ini.setValue(section, "curInit", line);
        sprintf_s(line, "%d", nTree.m_traceInitInfo.size());
        ini.setValue(section, "arrow_num", line);
        for (size_t j = 0; j < nTree.m_traceInitInfo.size(); ++j){
            const Vec3d& pos = nTree.m_traceInitInfo[j].first;
            const Vec3d& dir = nTree.m_traceInitInfo[j].second;
            sprintf_s(arrow, "arrow_%03d", j + 1);
            sprintf_s(boundString, "%lf %lf %lf %lf %lf %lf", pos(0), pos(1), pos(2), dir(0), dir(1), dir(2));
            ini.setValue(section, arrow, boundString);
        }
        //suspicious
        if (!nTree.m_suspiciousPositions_.empty()) {
            sprintf_s(line, "%lu", nTree.m_suspiciousPositions_.size());
            ini.setValue(section, "suspicious_num", line);
            for (size_t k = 0; k < nTree.m_suspiciousPositions_.size(); ++k) {
                auto &pos = nTree.m_suspiciousPositions_[k];
                sprintf_s(line, "sus_%03d", k+1);
                sprintf_s(susStr, "%lf %lf %lf %lf", pos(0), pos(1), pos(2), double(nTree.m_suspiciousCheckState_[k]));
                ini.setValue(section, line, susStr);
            }
        }
        //traverse
        bool traNum = false;
        for (auto &it1 : nTree.m_curveList) {
            for (auto &it2 : it1) {
                if (it2(4) > 0.5) {
                    traNum = true;
                    break;
                }
            }
            if (traNum == true) break;
        }
        if (traNum) {
            std::vector<std::string> dest;
            std::string curTraverseFile;
            Split(fileName, "/", dest);
            for (size_t k = 0; k < dest.size() - 1; ++k) {
                curTraverseFile += dest[k];
                curTraverseFile += "/";
            }
            std::vector<std::string> fileNamePre;
            Split(dest.back(), ".", fileNamePre);
            if (fileNamePre.empty()) 
                NG_ERROR_MESSAGE("");
            sprintf_s(line, "_traverse_%d.txt", i);
            curTraverseFile += fileNamePre[0];
            curTraverseFile += line;
            ini.setValue(section, "traverse", curTraverseFile);
            if (!WriteTraverse(curTraverseFile, nTree)){
                NG_ERROR_MESSAGE("WriteTraverse failed.");
                MAKEPROCESSSTATUS(resSta, false, className_, "WriteTraverse failed.");
                return resSta;
            }
        }
    }

    ini.saveas(fileName);
    /*fprintf(fp, "tree_num=%lu\n", tmpTree->m_pop.size());
    fprintf(fp, "image range : %d %d %d %d %d %d\n", paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
    fprintf(fp, "image extract range : %d %d %d\n", paramPack->xExtract_, paramPack->yExtract_, paramPack->zExtract_);
    fprintf(fp, "Binary Value of Tree and Axon : %lf %lf\n", paramPack->binThreshold_, paramPack->axonBinaryThreshold_);
    fprintf(fp, "Diffuse and Trace Value of Tree : %lf %lf\n", paramPack->diffuseValue_, paramPack->traceValue_);
    fprintf(fp, "Diffuse and Trace Value of Axon : %lf %lf\n", paramPack->axonDiffuseValue_, paramPack->axonTraceValue_);
    fprintf(fp, "Max Bound Num : %d\n", paramPack->maxBoundNum_);
    fprintf(fp, "Opacity : %d %d\n", paramPack->lowOpac_, paramPack->highOpac_);
    for (size_t i = 0; i < tmpTree->m_pop.size(); ++i) {
        const NeuronTree &nTree = *(tmpTree->m_pop[i]);
        fprintf(fp, "[tree]\n");
        fprintf(fp, "path:%s\n", swcName[i].c_str());
        fprintf(fp, "curInit:%lf %lf %lf\n", nTree.m_curInitPt(0), nTree.m_curInitPt(1), nTree.m_curInitPt(2));
        for (size_t j = 0; j < nTree.m_traceInitInfo.size(); ++j){
            const Vec3d& pos = nTree.m_traceInitInfo[j].first;
            const Vec3d& dir = nTree.m_traceInitInfo[j].second;
            fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", pos(0), pos(1), pos(2),
                dir(0), dir(1), dir(2));
        }
    }
    fclose(fp);*/
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

bool TraceStackWriter::WriteTraverse(const std::string &fileName, const NeuronTree &nTree)
{
    if (nTree.m_curveList.empty()) {
        NG_ERROR_MESSAGE("");
        return false;
    }
    FILE *fp = fopen(fileName.c_str(), "w");
    if (!fp) return false;
    //1st line
    size_t curLineID, parentLineID;// , nearestID, pid;
    int id = 0;
    {
        curLineID = (nTree.m_Topo.begin() + 1)->id;
        const Line5d &line = nTree.m_curveList[curLineID];
        if (line.empty()){
            printf("error in SaveATree. empty 1st tree.\n");
            return false;
        }
        ++id;
        fprintf(fp, "%d %d -1\n", id, line[0](4) > 0.5 ? 1 : 0);//root
        for (size_t i = 1; i < line.size(); ++i){
            ++id;
            fprintf(fp, "%d %d 0\n", id, line[i](4) > 0.5 ? 1 : 0);//root
        }
    }
    //left lines
    for (auto it = nTree.m_Topo.begin() + 2; it != nTree.m_Topo.end(); ++it) {
        curLineID = it->id;
        const Line5d &line = nTree.m_curveList[curLineID];
        if (nTree.m_Topo.parent(it) == nTree.m_Topo.begin()){
            //set the pid as 1st node
            ++id;
            fprintf(fp, "%d %d -1\n", id, line[0](4) > 0.5 ? 1 : 0);//root
            for (size_t i = 1; i < line.size(); ++i){
                ++id;
                fprintf(fp, "%d %d 0\n", id, line[i](4) > 0.5 ? 1 : 0);//root
            }
        }
        else{//find parent line-vertex id
            parentLineID = nTree.m_Topo.parent(it)->id;
            //const Vec3d &headNode = line.front().block(0, 0, 3, 1);
            //nearestID = KPTREE::FindNearestID(headNode, nTree.m_curveList[parentLineID]);
            //TODO: debug
            //if (nearestID > 1000000) {
            //    NG_ERROR_MESSAGE("how can pid so large? please debug.\n");
            //}
            ++id;
            
            fprintf(fp, "%d %d -1\n", id, line[0](4) > 0.5 ? 1 : 0);//root
            for (size_t i = 1; i < line.size(); ++i){
                ++id;
                fprintf(fp, "%d %d 0\n", id, line[i](4) > 0.5 ? 1 : 0);//root
            }
        }
    }
    /*int id = 0;
    for (auto &it : allCurve) {
    fprintf(fp, "%d %d %d\n", ++id, it[0](4) > 0.5 ? 1 : 0, -1);
    for (size_t k = 1; k < it.size(); ++k) {
    fprintf(fp, "%d %d %d\n", ++id, it[k](4) > 0.5 ? 1 : 0, 0);
    }
    }*/
    fclose(fp);
    return true;
}

void TraceStackWriter::Split(const std::string& src, const std::string& separator, std::vector<std::string>& dest)
{
    string str = src;
    string substring;
    string::size_type start = 0, index;
    do
    {
        index = str.find_first_of(separator, start);
        if (index != string::npos)
        {
            substring = str.substr(start, index - start);
            dest.push_back(substring);
            start = str.find_first_not_of(separator, index);
            if (start == string::npos) return;
        }
    } while (index != string::npos);

    //the last token  
    substring = str.substr(start);
    dest.push_back(substring);
}
