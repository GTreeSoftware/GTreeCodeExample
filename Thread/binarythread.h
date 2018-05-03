
/*
*图像二值化操作的线程类
*在AutoNeuronGPS中没有
*目的是在二值化的同时，不会使得窗口假死
*
*作者：周航
*时间：2013-8-26
*武汉国家光电实验室
*/

#ifndef BINARYTHREAD_H
#define BINARYTHREAD_H
#include <QThread>
#include <memory>
#include "../Function/ineuronprocessobject.h"

class BinaryThread : public QThread
{
	Q_OBJECT
public:
	BinaryThread(QObject *parent);
	~BinaryThread();
    void setParam(INeuronProcessPointer);
	void stop();
    ProcStatPointer RetVal()const{ return retVal; }

protected:
	void run();

signals:

private:
	volatile bool stopped;//no use
	///--------------data------------------------------
	bool			isInit;//
    INeuronProcessPointer filter;
	//run status----//
    ProcStatPointer retVal;
};

#endif