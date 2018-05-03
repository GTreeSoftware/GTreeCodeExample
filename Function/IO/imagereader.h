/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef IMAGEREADER_H
#define IMAGEREADER_H
#include <memory>
#include <string>
#include "../ineuronio.h"
#include "../../ngtypes/ParamPack.h"
template<typename T> class Volume;
typedef Volume<unsigned char> CVolume;
typedef Volume<unsigned short> SVolume;
class INeuronDateObject;
class ImageReader;

typedef std::shared_ptr<ImageReader> NGImageReader;

class ImageReader : public INeuronIO
{
public:
    typedef std::shared_ptr<SVolume> VolumePointer;
    static NGImageReader New(){return NGImageReader(new ImageReader());}

    ImageReader();
    ~ImageReader();
    INEURONPROCESSOBJECT_DEFINE
    bool SetInputFileName(const std::string&);
    void SetParam(NGParamPack arg){paramPack = arg;}
    bool GetImageInfo(const std::string&);//save information into paramPack
    DATATYPE GetImageType()const{ return dataType; }
    
protected:
    //--- do not use--//
    void SetInput(ConstIDataPointer&){}
    bool SetOutputFileName(const std::string&){ return false; }
    //---//
    bool ReadImage(const char*);
    
    //2015-6-8
    //void GetHist(int z);
    //void FillBlackArea(int z);
private:
    bool isFileValid;
    std::string filename; //read path
    DATATYPE dataType;
    //int histogram[1021]; //2015-6-8
    //double xres_, yres_, zres_;
    NGParamPack paramPack;

};

#endif // IMAGEREADER_H
