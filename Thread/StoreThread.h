

#ifndef STORETHREAD_H
#define STORETHREAD_H
#include <QApplication> 
#include <QThread>
//#include <QMutex>
//#include <QMutexLocker>
#include <memory>
#include "Function/ineuronprocessobject.h"
#include "ngtypes/ParamPack.h"
#include "Function/ineuronio.h"
//#include "Function/IO/MOSTDReader.h"

NG_SMART_POINTER_TYPEDEF(INeuronBigReader, NGNeuronBigReader);
class StoreThread : public QThread
{
	Q_OBJECT
public:
	StoreThread(QObject *parent = nullptr);
	~StoreThread();
    //void SetMutex(QMutex *arg){ mutex = arg; }
	void setParam(NGParamPack &arg, NGNeuronBigReader arg2);
	void stop();
	bool RetVal()const{ return retVal; }
	IDataPointer Imagedata;

protected:
	void run();

signals:

private:
    //
    //QMutex *mutex;
    //
	NGNeuronBigReader mostdReader;
	NGParamPack paramPack;
	volatile bool stopped;//no use
	///--------------data------------------------------
	bool			isInit;//
	INeuronProcessPointer filter;
	//run status----//
	bool retVal;
};

#endif