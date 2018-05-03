#include "StoreThread.h"

StoreThread::StoreThread(QObject *parent)
: QThread(parent)
{
	stopped = false;
	isInit = false;
}

StoreThread::~StoreThread()
{
	wait();
}

void StoreThread::run()
{
	if (isInit){
        retVal = true;
		//printf("nima\n");
		//mostdReader->SetParam(paramPack);
		printf("INThread %d %d %d %d %d %d\n", paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
        /*if (!mostdReader->UpdateROI()){
            retVal = false;
            printf("nima1\n");
            return;
            }*/
		if (!mostdReader->Update()->success()){
            retVal = false;
			printf("nima2\n");
            return;
		}
		//Imagedata = mostdReader->ReleaseData();
	}
}

void StoreThread::stop()
{
	/*stopped = true;*/
}

void StoreThread::setParam(NGParamPack &arg, NGNeuronBigReader arg2)
{
	paramPack = arg;
	mostdReader = arg2;
	isInit = true;
}
