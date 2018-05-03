/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/

#include <cstdlib>
#include <cmath>
#include "tracestackreader.h"
#include "../../ngtypes/soma.h"
#include "../../ngtypes/tree.h"
#include "Function/IO/treeReader.h"
#include "inifile.h"
#ifdef _WIN32
#pragma warning(disable : 4996)
#endif

TraceStackReader::TraceStackReader()
{
    className_ = std::string("TraceStackReader");
    isFileValid = false;
}

TraceStackReader::~TraceStackReader()
{

}

bool TraceStackReader::SetInputFileName(const std::string &path)
{
    if(!path.empty()){
        filename = path;
        isFileValid = true;
        return true;
    }
    return false;
}

ProcStatPointer TraceStackReader::Update()
{
    Cell tempCell;
    //char line[256];

    if (isFileValid){
        if (!m_Input) {
            NG_ERROR_MESSAGE("no input seperateTree"); 
            MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
            return resSta;
        }
        
        FILE* fp = fopen(filename.c_str(), "r");
        if (!fp){
            printf("error occurred in %s, the configure file is invalid\n", className_.c_str());
            MAKEPROCESSSTATUS(resSta, false, className_, "File is invalid.");
            return resSta;
        }
        fclose(fp);

        if (!m_Input) {
            NG_ERROR_MESSAGE("no input seperateTree"); 
            MAKEPROCESSSTATUS(resSta, false, className_, "there is no tree data.");
            return resSta;
        }
        if (paramPack->activeTree) paramPack->activeTree.reset();//no active

        inifile::IniFile ini;
        ini.load(filename);
        int ret = 0;
        int tree_num = 0;
        tree_num = ini.getIntValue("common", "tree_num", ret);
        if (ret < 0){ 
            NG_ERROR_MESSAGE(""); 
            MAKEPROCESSSTATUS(resSta, false, className_, "There are no common items.");
            return resSta;
        }
        paramPack->xMin_ = ini.getIntValue("common", "xMin", ret);
        paramPack->xMax_ = ini.getIntValue("common", "xMax", ret);
        paramPack->yMin_ = ini.getIntValue("common", "yMin", ret);
        paramPack->yMax_ = ini.getIntValue("common", "yMax", ret);
        paramPack->zMin_ = ini.getIntValue("common", "zMin", ret);
        paramPack->zMax_ = ini.getIntValue("common", "zMax", ret);
        paramPack->xExtract_ = ini.getIntValue("common", "x_Extract", ret);
        if (ret < 0) paramPack->xExtract_ = 60;
        paramPack->yExtract_ = ini.getIntValue("common", "y_Extract", ret);
        if (ret < 0) paramPack->yExtract_ = 60;
        paramPack->zExtract_ = ini.getIntValue("common", "z_Extract", ret);
        if (ret < 0) paramPack->zExtract_ = 54;
        paramPack->binThreshold_ = ini.getDoubleValue("common", "Tree_Binary_Value", ret);
        if (ret < 0) paramPack->binThreshold_ = 20;
        paramPack->axonBinaryThreshold_ = ini.getDoubleValue("common", "Axon_Binary_Value", ret);
        if (ret < 0) paramPack->axonBinaryThreshold_ = 3;
        paramPack->diffuseValue_ = ini.getDoubleValue("common", "Tree_Diffuse_Value", ret);
        if (ret < 0) paramPack->diffuseValue_ = 1;
        paramPack->traceValue_ = ini.getDoubleValue("common", "Tree_Trace_Value", ret);
        if (ret < 0) paramPack->traceValue_ = 4096;
        paramPack->axonDiffuseValue_ = ini.getDoubleValue("common", "Axon_Diffuse_Value", ret);
        if (ret < 0) paramPack->axonDiffuseValue_ = 1;
        paramPack->axonTraceValue_ = ini.getDoubleValue("common", "Axon_Trace_Value", ret);
        if (ret < 0) paramPack->axonTraceValue_ = 4096;
        paramPack->maxBoundNum_ = ini.getIntValue("common", "Max_Bound_Num", ret);
        if (ret < 0) paramPack->maxBoundNum_ = 5;
        paramPack->lowOpac_ = ini.getIntValue("common", "Opacity_Low", ret);
        if (ret < 0) paramPack->lowOpac_ = 0;
        paramPack->highOpac_ = ini.getIntValue("common", "Opacity_High", ret);
        if (ret < 0) paramPack->highOpac_ = 1000;
        paramPack->xScale_ = ini.getIntValue("common", "XScale", ret);
        if (ret < 0) paramPack->xScale_ = 1;
        paramPack->yScale_ = ini.getIntValue("common", "YScale", ret);
        if (ret < 0) paramPack->yScale_ = 1;
        paramPack->lowAdjustOpac_ = ini.getIntValue("common", "Opacity_Adjust_Low", ret);
        if (ret < 0) paramPack->lowAdjustOpac_ = 0;
        paramPack->highAdjustOpac_ = ini.getIntValue("common", "Opacity_Adjust_High", ret);
        if (ret < 0) paramPack->highAdjustOpac_ = 0;
        paramPack->lowAdjustOpac_ = ini.getIntValue("common", "Opacity_Dest_Low", ret);
        if (ret < 0) paramPack->lowDestOpac_ = 0;
        paramPack->highAdjustOpac_ = ini.getIntValue("common", "Opacity_Dest_High", ret);
        if (ret < 0) paramPack->highDestOpac_ = 0;
        paramPack->runningMinutes_ = ini.getIntValue("common", "Running_Minutes", ret);
        if (ret < 0) paramPack->runningMinutes_ = 0;
        char line[10];
        std::string tmpFile;
        swcNames.clear();
        stringstream Oss;
        VectorVec3d curInitList; curInitList.reserve(tree_num);
        Vec3d tmpPt, tmpInitPt, tmpInitDir, tmpSus;
        std::vector<BOUNDINFOLIST> tmpBoundList; tmpBoundList.resize(tree_num);
        std::vector<std::string> traverseFileList; traverseFileList.resize(tree_num);
        std::vector<VectorVec3d> allSusList; allSusList.resize(tree_num);
        std::vector<std::vector<double> > allSusFlagList; allSusFlagList.resize(tree_num);
        for (int k = 0; k < tree_num; ++k) {
            sprintf_s(line, "%d", k + 1);
            tmpFile = "tree" + std::string(line);
            std::string tmpSwc = ini.getStringValue(tmpFile, "path", ret);
            if (ret < 0) {
                tree_num = k;
                break;
            }
            swcNames.push_back(tmpSwc);
            std::string curInitStr = ini.getStringValue(tmpFile, "curInit", ret);
            if (ret!=-1) {
                std::vector<std::string> resList;
                Split(curInitStr, " ", resList);
                if (resList.size() != 3lu) {
                    tree_num = k;
                    break;
                }
                for (size_t i = 0; i < resList.size(); ++i) {
                    Oss << resList[i];
                    Oss >> tmpPt(i);
                    Oss.clear();
                    Oss.str("");
                }
            }
            else
                tmpPt << 0, 0, 0;
            
            curInitList.push_back(tmpPt);
            std::string tmpArrowName;
            int arrowNum = ini.getIntValue(tmpFile, "arrow_num", ret);
            auto &boundInfo = tmpBoundList[k];
            for (int i = 0; i < arrowNum; ++i) {
                sprintf_s(line, "%03d", i + 1);
                tmpArrowName = "arrow_" + std::string(line);
                std::string curArrow = ini.getStringValue(tmpFile, tmpArrowName, ret);
                std::vector<std::string> arrowList;
                Split(curArrow, " ", arrowList);
                if (arrowList.size() != 6lu) {
                    MAKEPROCESSSTATUS(resSta, false, className_, "Arrow list info is valid.");
                    return resSta;
                }
                for (size_t j = 0; j < 3lu; ++j) {
                    Oss << arrowList[j];
                    Oss >> tmpInitPt(j);
                    Oss.clear();
                    Oss.str("");
                }
                for (size_t j = 3lu; j < 6lu; ++j) {
                    Oss << arrowList[j];
                    Oss >> tmpInitDir(j - 3);
                    Oss.clear();
                    Oss.str("");
                }
                boundInfo.push_back(std::make_pair(tmpInitPt, tmpInitDir));
            }
            //read suspicious point
            int spNum = ini.getIntValue(tmpFile, "suspicious_num", ret);
            std::string tmpSpName;
            if (ret != -1) {
                for (int i = 0; i < spNum; ++i) {
                    sprintf_s(line, "%03d", i + 1);
                    tmpSpName = "sus_" + std::string(line);
                    std::string curSus = ini.getStringValue(tmpFile, tmpSpName, ret);
                    std::vector<std::string> susList;
                    Split(curSus, " ", susList);
                    if (susList.size() != 4lu) {
                        NG_ERROR_MESSAGE("");
                        MAKEPROCESSSTATUS(resSta, false, className_, "Suspicious info is invalid.");
                        return resSta;
                    }
                    for (size_t j = 0; j < 3lu; ++j) {
                        Oss << susList[j];
                        Oss >> tmpSus(j);
                        Oss.clear();
                        Oss.str("");
                    }
                    double tmpState;
                    Oss << susList[3];
                    Oss >> tmpState;
                    Oss.clear();
                    Oss.str("");
                    allSusList[k].push_back(tmpSus);
                    allSusFlagList[k].push_back(tmpState);
                }
            }
            //read traverse file
            std::string traverseFile = ini.getStringValue(tmpFile, "traverse", ret);
            if (ret != -1) {
                traverseFileList[k] = traverseFile;
            }
            
        }
        NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Input);
        //read old SWC
        if (!swcNames.empty()){
            NGTreeReader treeReader = TreeReader::New();
            treeReader->SetInputFileNameList(swcNames);
            treeReader->SetParam(paramPack);
            treeReader->SetTraverseFile(traverseFileList);
            auto tres = treeReader->Update();
            if (!tres->success()){
                NG_ERROR_MESSAGE(""); 
                return tres;
            }
            auto tmp1 = treeReader->ReleaseData();//seperateTree = is wrong, it is value transfer.
            NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmp2, tmp1);
            tmpTree->m_pop.swap(tmp2->m_pop);
            paramPack->activeTree = tmpTree->m_pop[0];
        }

        //set trace init info
        for (size_t kk = 0; kk < tmpTree->m_pop.size(); ++kk) {
            auto &curTree = tmpTree->m_pop[kk];
            curTree->m_curInitPt.swap(curInitList[kk]);
            curTree->m_traceInitInfo.swap(tmpBoundList[kk]);
            curTree->m_suspiciousPositions_.swap(allSusList[kk]);
            curTree->m_suspiciousCheckState_.resize(curTree->m_suspiciousPositions_.size(), 0);
            auto &curSusState = curTree->m_suspiciousCheckState_;
            auto &curReadSusState = allSusFlagList[kk];
            for (size_t n = 0; n < curReadSusState.size(); ++n) {
                if (curReadSusState[n] > 0.5 &&curReadSusState[n] < 1.5) curSusState[n]=1;
                else if (curReadSusState[n] < 0.5) curSusState[n] = 0;
                else curSusState[n] = 2;
            }
        }
        MAKEPROCESSSTATUS(resSta, true, className_, "");
        return resSta;
    }
    MAKEPROCESSSTATUS(resSta, false, className_, "File is invalid.");
    return resSta;
}

void TraceStackReader::Split(const std::string& src, const std::string& separator, std::vector<std::string>& dest)
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

