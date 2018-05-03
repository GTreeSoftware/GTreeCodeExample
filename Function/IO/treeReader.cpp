#include "TreeReader.h"
#include "../ngtypes/tree.h"

TreeReader::TreeReader(void)
{
    className_ = std::string("TreeReader");
    m_Source = std::shared_ptr<NeuronPopulation>(new NeuronPopulation(this));
    isFileValid = false;
}

TreeReader::~TreeReader(void)
{
}

void TreeReader::SetInputFileNameList( const std::vector<std::string>& path)
{
    if(!path.empty()){
        fileNameList = path;
        isFileValid = true;
    }
}


ProcStatPointer TreeReader::Update()
{
    if(!isFileValid){
        printf("no input tree files.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "no input file.");
        return resSta;
    }

    NG_CREATE_DYNAMIC_CAST(NeuronPopulation, tmpTree, m_Source);
    std::vector<NeuronTreePointer>& neuronTreeList = tmpTree->m_pop;
    //std::vector<std::vector<VectorVec5d> >& readTree = tmpSeperateTree->GetTree();
    //std::vector<int> &typeList = tmpSeperateTree->GetTypeList();
    int offset = 0;
    if (traverseNameList.size() != fileNameList.size()) {
        printf("The traverse file path is invalid.\n");
        traverseNameList.clear();
    }
    for(size_t i = 0; i < fileNameList.size();++i){
        FILE* fp = fopen((const char*)(fileNameList[i].c_str()), "r");
        if(!fp) {
            printf("error write tree!\n");
            MAKEPROCESSSTATUS(resSta, false, className_, "File is invalid.");
            return resSta;
        }
        //One file One tree
        neuronTreeList.push_back(std::make_shared<NeuronTree>());
        std::vector<Line5d> &localTree = neuronTreeList[i]->m_curveList;
        auto& localTopo = neuronTreeList[i]->m_Topo;
        localTopo.clear();
        localTopo.set_head(LineID(std::numeric_limits<size_t>::max()));//root null
        Vec5d tmpVec5d;
        tmpVec5d.setZero();
        VectorVec5d tmpVectorVec5d;
        VectorVec5d tmpVectorVec5d1;
        int id(0), type(0), pid(0);
        char line[256];
        neuronTreeList[i]->m_type = 2;
        size_t lineID = 0;
        //VectorVec5d pidNodeList;//the parent id for add head node of each curve
        while (!feof(fp)) {
            if (fgets(line, 256, fp) == 0) continue;
            if (line[0] == '#' || line[0] == '\n' || line[0] == '/' || line[0] == '\\') continue;
            offset = 0;
            if (line[0] == ' ') offset = 1;
            sscanf(line + offset, "%d %d %lf %lf %lf %lf %d", &id, &type, &tmpVec5d(0), &tmpVec5d(1), &tmpVec5d(2), &tmpVec5d(3), &pid);
            tmpVec5d(0) /= param->xRes_;
            tmpVec5d(1) /= param->yRes_;
            tmpVec5d(2) /= param->zRes_;
            //tmpVec5d(3) = id;
            tmpVectorVec5d1.push_back(tmpVec5d);
            if (pid == -1 || id != pid + 1) {
                if (tmpVectorVec5d.size() > 1) {
                    localTree.push_back(tmpVectorVec5d);
                    localTopo.append_child(localTopo.begin(), LineID(lineID, KPTREE::Num2NeuriteType(type)));
                    ++lineID;
                }
                tmpVectorVec5d.clear();
                if(pid != -1){//force add head node
                    Vec5d tmpParentVec5d = tmpVectorVec5d1[pid - 1];
                    tmpVectorVec5d.push_back(tmpParentVec5d); 
                }
                else{
                    tmpVectorVec5d.push_back(tmpVec5d);
                }
                tmpVectorVec5d.push_back(tmpVec5d);
            }else{
                //if (preFlag) {
                    //preFlag = false;
                    //if ((tmpVec5d.block(0, 0, 3, 1) - tmpVectorVec5d.front().block(0, 0, 3, 1)).norm() < 1.0 )
                    //    continue;
                //}
                tmpVectorVec5d.push_back(tmpVec5d);
            }
        }
        if (!tmpVectorVec5d.empty()) {
            localTree.push_back(tmpVectorVec5d);
            tmpVectorVec5d.clear();
        }
        fclose(fp);
        if (traverseNameList.size() > i) {
            ReadTraverseFile(traverseNameList[i], neuronTreeList[i]->m_curveList);
        }
        
    }

    size_t maxLen;
    for (size_t i = 0; i < neuronTreeList.size(); ++i) {
        //remove repeat
        auto &curTreeCurve = neuronTreeList[i]->m_curveList;
        for (auto &curve: curTreeCurve) {
            maxLen = std::min(curve.size(), 10llu);
            curve.erase(std::unique(curve.begin(), curve.begin() + maxLen, [](const Vec5d& lhs, const Vec5d &rhs){return (lhs.block(0, 0, 3, 1) - rhs.block(0, 0, 3, 1)).cwiseAbs().sum() < 1.0; }), curve.begin() + maxLen);
        }
        curTreeCurve.erase(std::remove_if(curTreeCurve.begin(), curTreeCurve.end(), [](const Line5d &arg){return arg.size() < 2lu; }), curTreeCurve.end());
    }
    //rebuild topology
    for (size_t i = 0; i < neuronTreeList.size(); ++i) {
        Vec3d root = neuronTreeList[i]->m_curveList.front().front().block(0, 0, 3, 1);
        KPTREE::SetTreeRoot(*(neuronTreeList[i]), root);
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(TreeReader)
//IDataPointer TreeReader::ReleaseData()
//{
//    m_Source->ReleaseProcessObject();
//    IDataPointer tData = m_Source;
//    m_Source.reset();
//    return tData;
//}

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(TreeReader, NeuronPopulation)
//ConstIDataPointer TreeReader::GetOutput()
//{
//    if(!m_Source)
//#ifdef _WIN32
//        m_Source = std::tr1::shared_ptr<NeuronPopulation>(new NeuronPopulation(this));
//#else
//        m_Source = std::shared_ptr<NeuronPopulation>(new NeuronPopulation(this));
//#endif
//    return m_Source;
//}

bool TreeReader::ReadTraverseFile(const std::string &file, std::vector<Line5d> &allCurve)
{
    if (file.empty()) return false;
    FILE* fp = fopen(file.c_str(), "r");
    if (!fp) return false;
    char line[64];
    int curLineID = -1, curPos = -1;
    int id, state, pid;//debug
    while (!feof(fp)){// -1 denote next line
        if (fgets(line, 64, fp) == 0) continue;
        if (line[0] == '#' || line[0] == '\n' || line[0] == '/' || line[0] == '\\') continue;
        sscanf(line, "%d %d %d", &id, &state, &pid);
        if (-1 == pid) {
            //attention for the split curves, of which should be one curves
            if (curLineID == -1 || curPos == int(allCurve[curLineID].size()) - 1) {
                ++curLineID;
                curPos = 0;
                allCurve[curLineID][curPos](4) = double(state);
            }
        }
        ++curPos;
        allCurve[curLineID][curPos](4) = double(state);
    }
    fclose(fp);
    return true;
}

void TreeReader::SetTraverseFile(const std::vector<std::string>& arg)
{
    if (!arg.empty()) {
        traverseNameList = arg;
    }
}
