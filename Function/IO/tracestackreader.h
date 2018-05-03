/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef TRACESTACKREADER_H
#define TRACESTACKREADER_H
#include <memory>
#include <string>
#include "../ineuronio.h"
#include "ngtypes/soma.h"
#include "ngtypes/ParamPack.h"
//#include "../../ngtypes/ineurondataobject.h"
class INeuronDataObject;
//class Soma;
class TraceStackReader;
typedef std::shared_ptr<TraceStackReader> NGTraceStackReader;

class TraceStackReader : public INeuronIO
{
public:
    static NGTraceStackReader New(){return NGTraceStackReader(new TraceStackReader());}
    TraceStackReader();
    ~TraceStackReader();
    bool SetInputFileName(const std::string&);
    virtual ProcStatPointer Update();
    void SetParam(NGParamPack arg){paramPack = arg;}
    const std::vector<std::string>& GetSwcFileName(){ return swcNames; }

protected:
    bool SetOutputFileName(const std::string&){ return false; }
    IDataPointer ReleaseData(){ return m_Source; }
    ConstIDataPointer GetOutput(){ return m_Source; }
    void Split(const std::string& src, const std::string& separator, std::vector<std::string>& dest);
   

private:
    bool isFileValid;
    std::string filename;
    std::vector<std::string> swcNames;
    NGParamPack paramPack;
    
};

#endif // SOMAREADER_H
