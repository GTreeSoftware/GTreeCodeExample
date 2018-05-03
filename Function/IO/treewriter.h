#ifndef TREEWRITER_H
#define TREEWRITER_H
#include <string>
#include <vector>
#include "../ineuronio.h"
#include "../../ngtypes/basetypes.h"
#include "../../ngtypes/tree.h"
#include "../../ngtypes/ParamPack.h"
class TreeWriter;
typedef std::shared_ptr<TreeWriter> NGTreeWriter;
class TreeWriter : public INeuronIO
{
public:
    static NGTreeWriter New(){return NGTreeWriter(new TreeWriter());}
    TreeWriter();
    ~TreeWriter();
    bool SetOutputFileName(const std::vector<std::string> &);
    virtual ProcStatPointer Update();
    ProcStatPointer UpdateLayer(const std::vector<int>& layerDisplay, const std::vector<int>& branchDisplay);
    void SetParam(NGParamPack arg){ param = arg; }
    ProcStatPointer OldUpdate();
    void SetDD2(const VectorMat2i&);
protected:
    bool SaveATree(NeuronTree&, std::string&);
    IDataPointer ReleaseData(){return IDataPointer(new TreeCurve());}
    bool SetInputFileName(const std::string &){ return false; }
    ConstIDataPointer GetOutput(){return m_Source;}
    bool SetOutputFileName(const std::string &){ return true; }
    std::vector<std::string> fileName;
private:
    NGParamPack param;
    const VectorMat2i* dd2ptr;
};//TODO

#endif // TREEWRITER_H
