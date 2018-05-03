/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef SOMAREADER_H
#define SOMAREADER_H
#include <memory>
#include <string>
#include "../ineuronio.h"
#include "ngtypes/ParamPack.h"
//#include "../../ngtypes/ineurondataobject.h"
class INeuronDataObject;
class Soma;
class SomaReader;
#ifdef _WIN32
typedef std::tr1::shared_ptr<SomaReader> NGSomaReader;
#else
typedef std::shared_ptr<SomaReader> NGSomaReader;
#endif

class SomaReader : public INeuronIO
{
public:
#ifdef _WIN32
    typedef std::tr1::shared_ptr<Soma> SomaPointer;
#else
    typedef std::shared_ptr<Soma> SomaPointer;
#endif
    static NGSomaReader New(){return NGSomaReader(new SomaReader());}
    SomaReader();
    ~SomaReader();
    bool SetInputFileName(const std::string&);
    INEURONPROCESSOBJECT_DEFINE
    void SetParam(NGParamPack arg){paramPack = arg;}
    //void SetImageResolution(double x, double y, double z);
private:
    void SetInput(ConstIDataPointer){}
    bool SetOutputFileName(const std::string&){ return false; }
    bool isFileValid;
    std::string filename;
    //double xres_, yres_, zres_;
    double scale_;
    NGParamPack paramPack;
};

#endif // SOMAREADER_H
