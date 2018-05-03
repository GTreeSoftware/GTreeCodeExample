
//#include ¡°stdafx.h¡±
#include "binarythread.h"

BinaryThread::BinaryThread(QObject *parent)
: QThread(parent)
{
    stopped = false;
    isInit		= false;
}

BinaryThread::~BinaryThread()
{
    wait();
}

void BinaryThread::run()
{
    MAKEPROCESSSTATUS(resSta, true, filter->ClassName(), "");
    if(isInit){
        retVal = filter->Update();
    }
}

void BinaryThread::stop()
{
    /*stopped = true;*/
}

void BinaryThread::setParam(INeuronProcessPointer arg)
{
    filter = arg;
    isInit = true;
}
