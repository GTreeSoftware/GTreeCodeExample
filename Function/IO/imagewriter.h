/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef IMAGEWRITER_H
#define IMAGEWRITER_H
#include <memory>
#include "../ineuronio.h"
#include <string>
template<typename T> class Volume;
typedef Volume<unsigned char> CVolume;
typedef Volume<unsigned short> SVolume;
class INeuronDataObject;
class ImageWriter;
typedef std::shared_ptr<ImageWriter> NGImageWriter;

class ImageWriter : public INeuronIO
{
public:
    typedef std::shared_ptr<const SVolume> CSVolumePointer;
    typedef std::shared_ptr<SVolume> SVolumePointer;
    typedef std::shared_ptr<CVolume> CVolumePointer;
    
    static NGImageWriter New(){return NGImageWriter(new ImageWriter());}
    ImageWriter(void);
    ~ImageWriter(void);

    //void SetInput(ConstDataPointer&);derived function
    bool SetOutputFileName(const std::string&);
    ProcStatPointer Update();

private:
    //private function
    IDataPointer ReleaseData();
    ConstIDataPointer GetOutput();
    bool SetInputFileName(const std::string&){ return false; }
    bool WriteImage(const char *path);
    //**//
    bool WriteTiff(const char* path, SVolumePointer &pData);
    bool WriteTiff(const char* path, CVolumePointer &data);
    //output jpg
    //bool WriteJPEG(const char* path, CVolumePointer &pData);
    //
    std::string filename;
    bool isFileValid;
};
#endif
