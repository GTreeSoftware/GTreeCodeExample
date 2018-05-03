/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef TRACESTACKWRITER_H
#define TRACESTACKWRITER_H
#include "../ineuronio.h"
#include "../ngtypes/basetypes.h"
#include "../ngtypes/tree.h"
#include "ngtypes/soma.h"
class TraceStackWriter;
struct ParamPack;
typedef std::shared_ptr<TraceStackWriter> NGTraceStackWriter;
typedef std::shared_ptr<ParamPack> NGParamPack;
class TraceStackWriter : public INeuronIO
{
public:

    static NGTraceStackWriter New(){return NGTraceStackWriter(new TraceStackWriter());}
    TraceStackWriter();
    ~TraceStackWriter();
    bool SetOutputFileName(const std::string &);
    virtual ProcStatPointer Update();
	void SetParam(NGParamPack arg){ paramPack = arg; }
    void SetSwcFileNames(const std::vector<std::string> &f){ swcName = f; }
protected:
    bool WriteTraverse(const std::string &fileName, const NeuronTree &traList);
    IDataPointer ReleaseData(){return m_Source;}
    bool SetInputFileName(const std::string &){ return false; }
    ConstIDataPointer GetOutput(){return m_Source;}
    void Split(const std::string& src, const std::string& separator, std::vector<std::string>& dest);
    std::string fileName;
    std::vector<std::string> swcName;
//    std::string identifyName;
//    ConstDataPointer m_Input;
//    DataPointer m_Source;//output data
	NGParamPack paramPack;
};

#endif // SOMAWRITER_H
