#ifndef PRESTORE_H
#define PRESTORE_H
#include <QObject>
#include <algorithm>
#include "../../ngtypes/basetypes.h"
#include "../../ngtypes/tree.h"
#include "Function/ineuronprocessobject.h"
#include "../../ngtypes/ParamPack.h"
#include "Thread/binarythread.h"
#include "Function/ineuronio.h"
#include "ngtypes/volume.h"
#include "../../Thread/StoreThread.h"
#include <tuple>

NG_SMART_POINTER_TYPEDEF(INeuronBigReader, NGNeuronBigReader);
class Prestore;
NG_SMART_POINTER_TYPEDEF(Prestore, NGPrestore);

class Prestore : public QObject, public INeuronProcessObject
{
	Q_OBJECT
public:
	static NGPrestore New(){ return NGPrestore(new Prestore()); };
	Prestore();
	~Prestore();
    virtual ProcStatPointer Update(){
        MAKEPROCESSSTATUS(resSta, true, className_, "");
        return resSta;
    }//get image and start read thread ??? useless?, size_t arg

	
	bool GetNextImage(IDataPointer &arg);//get check start thread
    bool GetPreviousImage(IDataPointer &arg);//get check start thread

	void AppendNewROI();//check start thread

	size_t GetWaitingROINum(){ return waitingROIDeque_.size(); }
    void SetParam(NGNeuronBigReader arg2, NGParamPack arg3);

	IDataPointer previousImgptr;
	//Vec6i previousCoor;
	//VectorVec6i wrongCoor;
	//void showWrongFrame();
	void BackTrace();

    bool Initial();
    bool Finish();
    void Reset();

    void SetTraversePosition(int, int);

protected slots:
	void EndThread_Slot();
	
signals:
    void CacheComplete_Signal();

protected:
	ConstIDataPointer GetOutput(){ return m_Source; }
	IDataPointer ReleaseData(){ return IDataPointer(new TreeCurve()); }
	bool CheckThreadStart();
	void PopROI();//
	bool IsImageDequeValid();
    void StartThread();
	void applyread();
	void buildFrame();
	void FowardTrace();
	//void AddCurLineID();
	void calculateRange();
private:
    struct RangeIterator{
        RangeIterator(const tree<LineID>* t, const std::vector<VectorVec5d> *c)
            :treePointer(t), curvePointer(c){}
        ~RangeIterator(){}
        int curBeg = std::numeric_limits<int>::max(), curEnd = 0, curLineID = 0;
        const tree<LineID>* treePointer = NULL;
        const std::vector<VectorVec5d> *curvePointer = NULL;
        tree<LineID>::pre_order_iterator treeIter;
        void Next(int arg);
        void Back(int arg);
    };
    RangeIterator  *traverseIter = NULL;
    RangeIterator  *traverseIterCp= NULL;
    bool forceThreadStop_ = false;
	Vec6i FrameCoo;
	//size_t curLineID, preLineID;
	int numPoints_ = 40;
    int preBeg = std::numeric_limits<int>::max(), preEnd = 0, preLineID = 0;
	bool firstWork_ = true;
	bool firstWork_2 = true;
    bool directionFlag_ = true;
    int direction_;
	IDataPointer ImageData;
	size_t numBox, readingBox = 0;
	int storageImageMaxNum_ = 3;
	NGParamPack paramPack;
	NGParamPack paramPackCp;
	//VectorVec6i FrameCoordinate;
    NGNeuronBigReader mostdReader;
	std::deque<size_t> waitingROIDeque_;
	std::deque<std::tuple<std::shared_ptr<SVolume>, Vec6i, Vec3i> > imageDeque_;
	//std::deque<std::pair<std::shared_ptr<SVolume>, size_t>> imageDeque_;
	std::shared_ptr<SVolume> recentdata;
	StoreThread* worker;
};

#endif // NEURONPROCESS
