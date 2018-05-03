#ifndef INEURONPROCESSOBJECT_H
#define INEURONPROCESSOBJECT_H
#include "../ngtypes/ineurondataobject.h"
class ProcessStatus
{
public:
    ProcessStatus(){ isSuccess = false;  }
    ProcessStatus(bool arg, const std::string& str, const std::string &infostr):isSuccess(arg), className(str), errorStr(infostr){}
    ~ProcessStatus(){}
    bool success() const{ return isSuccess; }
    std::string ErrorClassName() const { return className; }
    std::string ErrorInfo() const { return errorStr; }
    void SetSuccess(bool arg) {isSuccess = arg; }
    void SetClassName(const std::string &name) { className = name; }
    void SetErrorStr(const std::string &str) { errorStr = str; }
private:
    bool isSuccess;
    std::string className;
    std::string errorStr;
};
typedef std::shared_ptr<ProcessStatus> ProcStatPointer;

class INeuronProcessObject
{
public:
    //INeuronProcessObject();
    virtual ProcStatPointer Update() = 0;
    virtual void SetInput(ConstIDataPointer input){m_Input = input;}
    virtual ConstIDataPointer GetOutput()=0;//{return m_Source;}
    virtual IDataPointer ReleaseData()=0;
    virtual ~INeuronProcessObject(){}
    std::string ClassName()const { return className_; }
protected:
    std::string className_;
    ConstIDataPointer m_Input;
    IDataPointer m_Source;//output data
};

typedef std::shared_ptr<INeuronProcessObject> INeuronProcessPointer;

#define MAKEPROCESSSTATUS(statusName, flag, className, infoStr) std::shared_ptr<ProcessStatus> statusName = \
    std::shared_ptr<ProcessStatus>(new ProcessStatus(flag, className, infoStr)); 

#define INEURONPROCESSOBJECT_DEFINE virtual ProcStatPointer Update(); \
    virtual ConstIDataPointer GetOutput(); \
    virtual IDataPointer ReleaseData();

#define INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(className) IDataPointer className::ReleaseData() \
{ \
    m_Source->ReleaseProcessObject(); \
    IDataPointer tData(m_Source); \
    m_Source.reset(); \
    return tData; } 

#define INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(className, typeName) ConstIDataPointer className::GetOutput() \
{    if (!m_Source) \
        m_Source = std::shared_ptr<typeName>(new typeName(this)); \
        return m_Source; }
        //NG_SMART_POINTER_NEW(typeName, m_Source, this); \
    

#endif // INEURONPROCESSOBJECT_H
