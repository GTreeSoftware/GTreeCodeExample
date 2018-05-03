/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef SOMAWRITER_H
#define SOMAWRITER_H
#include "../ineuronio.h"
#include "../ngtypes/soma.h"
class SomaWriter;
#ifdef _WIN32
typedef std::tr1::shared_ptr<SomaWriter> NGSomaWriter;
#else
typedef std::shared_ptr<SomaWriter> NGSomaWriter;
#endif
class SomaWriter : public INeuronIO
{
public:

    static NGSomaWriter New(){return NGSomaWriter(new SomaWriter());}
    SomaWriter();
    ~SomaWriter();
    bool SetOutputFileName(const std::string &);
    ProcStatPointer Update();
    //void SetInput(ConstDataPointer input)

protected:
    IDataPointer ReleaseData(){return IDataPointer(new Soma());}
    bool SetInputFileName(const std::string &){ return false; }
    ConstIDataPointer GetOutput(){return m_Source;}
    std::string fileName;
//    std::string identifyName;
//    ConstDataPointer m_Input;
//    DataPointer m_Source;//output data

};

#endif // SOMAWRITER_H
