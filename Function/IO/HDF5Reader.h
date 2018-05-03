#ifndef HDF5READER_H
#define HDF5READER_H

#include "../ineuronio.h"
#include "ngtypes/volume.h"
#include "ngtypes/basetypes.h"
class HDF5Reader;
typedef std::shared_ptr<HDF5Reader> NGHDF5Reader;

class HDF5Reader : public INeuronBigReader
{
public:
    static NGHDF5Reader New(){ return NGHDF5Reader(new HDF5Reader()); }
    HDF5Reader();
    virtual ~HDF5Reader();

    INEURONPROCESSOBJECT_DEFINE

    bool SetInputFileName(const std::string& arg);
    virtual void SetParam(NGParamPack arg){ paramPack_ = arg; }
    virtual void GetMaxRange(int& x, int& y, int& z);/*the size of mostd dataq*/
    virtual DATATYPE GetImageType()const{ return dataType; }
protected:
    virtual bool SetOutputFileName(const std::string&){ return false; }
    void Transfer1DPointerTo2DPointer_USHORT(unsigned short *data, int x, int y, unsigned short*** dst);
    DATATYPE dataType = DATATYPE::IMAGE16;
    std::string filename;
    NGParamPack paramPack_;
};

#endif