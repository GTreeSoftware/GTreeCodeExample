#ifndef INEURONIO_H
#define INEURONIO_H

#include "ineuronprocessobject.h"
struct ParamPack;
NG_SMART_POINTER_TYPEDEF(ParamPack, NGParamPack);
class INeuronIO : public INeuronProcessObject
{
public:
    virtual ~INeuronIO(){}
    virtual bool SetInputFileName(const std::string&)=0;
    virtual bool SetOutputFileName(const std::string&)=0;
};

class INeuronBigReader : public INeuronIO
{
public:
    virtual ~INeuronBigReader(){}
    virtual void SetParam(NGParamPack arg)=0;//TODO:
    virtual void GetMaxRange(int& x, int& y, int& z)=0;
    virtual DATATYPE GetImageType()const=0;
};

#endif // INEURONIO_H
