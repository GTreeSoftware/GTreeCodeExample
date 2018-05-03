#ifndef INEURONDATAOBJECT_H
#define INEURONDATAOBJECT_H
#include <memory>
#include <string>
#include "basetypes.h"
class INeuronProcessObject;
class INeuronDataObject;

//#ifdef _WIN32
//#define NG_SMART_POINTER_TYPEDEF(type, name) typedef std::tr1::shared_ptr<type> name
//#else
#define NG_SMART_POINTER_TYPEDEF(type, name) typedef std::shared_ptr<type> name
//#endif 
NG_SMART_POINTER_TYPEDEF(INeuronDataObject, IDataPointer);
NG_SMART_POINTER_TYPEDEF(INeuronDataObject, ConstIDataPointer);

class INeuronDataObject
{
public:
    virtual bool IsEmpty()const=0;
    virtual void SetProcessObject(INeuronProcessObject* p){m_ProcessObject = p;}
    virtual INeuronProcessObject *GetProcessObject()const{return m_ProcessObject;}
    virtual void ReleaseProcessObject(){m_ProcessObject = NULL;}
    const std::string& GetIdentifyName()const{return identifyName;}
    DATATYPE DataType()const{ return dataType; }
    void SetDataType(DATATYPE arg){ dataType = arg; }
    virtual ~INeuronDataObject(){}
protected:
    std::string identifyName;
    DATATYPE dataType;
    INeuronProcessObject* m_ProcessObject;//father process

};

#endif // NEURONDATAOBJECT_H
