#ifndef HDF5WRITER_H
#define HDF5WRITER_H
#include "../ineuronio.h"
#include "ngtypes/ParamPack.h"
#include "ngtypes/volume.h"

#include <QProgressDialog>
#include <QApplication>
#include <QObject>

class HDF5Writer;
typedef std::shared_ptr<HDF5Writer> NGHDF5Writer;

class HDF5Writer : public QObject, public INeuronIO
{
    Q_OBJECT
public:
    static NGHDF5Writer New(){ return NGHDF5Writer(new HDF5Writer()); }
    HDF5Writer();
    virtual ~HDF5Writer();

    virtual ProcStatPointer Update();
    bool SetOutputFileName(const std::string& arg){ filename = arg; return true; }
    void SetTIFFFileList(const std::vector<std::string>& arg){ tiffList = arg; }
    void SetParam(NGParamPack arg){ paramPack_ = arg; }
    volatile bool stopFlag_ = false;

signals:
    void Progress_Signal(int);

protected:
    //private function
    IDataPointer ReleaseData();
    ConstIDataPointer GetOutput();
    bool SetInputFileName(const std::string&){ return false; }
    //void Transfer1DPointerTo2DPointer_UCHAR(unsigned char *data, int x, int y, unsigned char** dst);
    void Transfer1DPointerTo2DPointer_USHORT(unsigned short *data, int x, int y, unsigned short*** dst);
    bool MakePyramid(const SVolume& orig, int level, SVolume &dst);

    int subblockSz_ = 256;
    NGParamPack paramPack_;
    std::string filename;
    std::vector<std::string> tiffList;
    
};

#endif // !NGHDF5WRITER_H
