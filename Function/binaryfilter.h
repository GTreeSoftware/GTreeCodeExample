#ifndef BINARYFILTER_H
#define BINARYFILTER_H
#include "../ngtypes/ParamPack.h"
#include "ineuronprocessobject.h"
#include "../ngtypes/basetypes.h"
#include "../ngtypes/volume.h"
class BinaryFilter;
NG_SMART_POINTER_TYPEDEF(BinaryFilter, NGBinaryFilter);
NG_SMART_POINTER_TYPEDEF(ParamPack, NGParamPack);

class BinaryFilter : public INeuronProcessObject
{
public:
    NG_SMART_POINTER_TYPEDEF(VectorVec3i, BinPtSetPointer);
    static NGBinaryFilter New(){return NGBinaryFilter(new BinaryFilter());}
    BinaryFilter();
    ~BinaryFilter();
    INEURONPROCESSOBJECT_DEFINE
    /*bool Update();
    ConstIDataPointer GetOutput();
    IDataPointer ReleaseData();*/
    ConstIDataPointer GetBackNoiseImage();
    VectorVec3i &GetBinPtSet();
    IDataPointer ReleaseBackNoiseImage();
    //BinPtSetPointer ReleaseBinPtSet();
    void SetParam(NGParamPack arg){paramPack = arg;}
    //void SetThreadNum(int);
    void SetThreshold(double);
    //void SetthreValue(double);
    //void SetRatio(double);
	void CalcBackImage(const SVolume&, SVolume&, int);
private:
    bool Binary();
    void CalcRandomMean(const SVolume&, double&);
    //int threadNum ;
    double binThreshold_ ;
    NGParamPack paramPack;
    //double threValue;
    //double ratio;
    //m_Source binary image
    IDataPointer m_Back;//backnoise
    VectorVec3i m_BinPtSet;
};

#endif // BINARYFILTER_H
