/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include "imagereader.h"
#include <tiffio.h>
#include <cstdlib>
#include "../../ngtypes/volume.h"
#include "Function/volumealgo.h"
#include "Function/NGUtility.h"

ImageReader::ImageReader()
{
   className_ = std::string("ImageReader");
   m_Source = std::shared_ptr<SVolume>(new SVolume(this));//2015-6-8
   isFileValid = false;
   /*xres_= 0.5;
   yres_ = 0.5;
   zres_ = 2.0;*/
}

ImageReader::~ImageReader()
{
    //printf("hello");
}

bool ImageReader::SetInputFileName(const std::string &path)
{
    TIFF* fp = TIFFOpen(path.c_str(), "r");
    if(NULL == fp){
        printf("can not open file : %s.\n error from %s", path.c_str(), className_.c_str());
        return false;
    }
    else {
        filename = path;
        isFileValid = true;
    }
    TIFFClose(fp);
    return true;
}

ProcStatPointer ImageReader::Update()
{
    if(isFileValid){
        if(!ReadImage(filename.c_str())){
            printf("error occurred in %s\n", className_.c_str());
            MAKEPROCESSSTATUS(resSta, false, className_, "File is invalid.");
            return resSta;
        }
            //printf("can not read image file.this is %s", identifyName.c_str());
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(ImageReader)
INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(ImageReader, SVolume)

bool ImageReader::ReadImage(const char *path)
{
    TIFFSetWarningHandler(0);//do not show warning dialog
    TIFF* in = TIFFOpen(path, "r");
    if(!in){
        return false;
    }
    else{
        int x, y, z;
        //int rx, ry, rz;
        //int xScale, yScale;
        uint16 bitspersample = 8;
        uint16 photo=1;//1 : 0 is dark, 255 is white
        uint16 samplesperpixel = 1;//8bit denote a pixel
        TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
        TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photo); // PHOTOMETRIC_MINISWHITE=0, PHOTOMETRIC_MINISBLACK=1
        TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
        if(bitspersample != 8){
            dataType = DATATYPE::IMAGE16;
            printf("Image is above 8 bpp.\n");
            //return false;
        }
        else{
            dataType = DATATYPE::IMAGE8;
            paramPack->filterNum_ = 20;
        }
        //-----get image size
        TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &x);
        TIFFGetField(in, TIFFTAG_IMAGELENGTH,&y);
        z = TIFFNumberOfDirectories(in);
       
        //double xRes = paramPack->xRes_;
        //double yRes = paramPack->yRes_;
        //double zRes = paramPack->zRes_;
        //zRes *= double(paramPack->zScale_);
        if (!m_Source) m_Source = std::shared_ptr<SVolume>(new SVolume(this));
        m_Source->SetProcessObject(this);//Attach m_Source/OrigImage to this module.
        VolumePointer image = std::dynamic_pointer_cast<SVolume>( m_Source );//2015-6-8 type conversion
        image->SetDataType(dataType);
        image->SetSize(x,y,z);
        image->SetResolution(paramPack->xRes_, paramPack->yRes_, paramPack->zRes_);      
        /*-----------read data from file--------*/
        TIFFSetDirectory(in,0);
        //int temp(0);
        //int maxValue= 0;
        if (bitspersample == 8) {
            //double ratio = 4.0 / double(paramPack->xScale_ * paramPack->yScale_);//2015-6-8 //the max value of read image must be 1020 = 255 * 4.
            //---------pre read buffer------------//
            uint32* raster = new uint32[x * y];//here must be x y z
            for(int nCur = 0; nCur < z; ++nCur){
                    //-----advanced read function is able to read image of any type into a 24bit image, maybe cost many memory----//
                    TIFFReadRGBAImageOriented(in, x, y, raster, ORIENTATION_TOPLEFT);
                    //if (paramPack->xScale_ * paramPack->yScale_ != 1 ){//re sample
                    //    for (int i = 0; i < x; ++i)//
                    //    for (int j = 0; j < y; ++j){
                    //        temp = 0;
                    //        //sum neighborhood value, do not divide
                    //        for (int ii = 0; ii < paramPack->xScale_; ++ii){
                    //            for (int jj = 0; jj < paramPack->yScale_; ++jj){
                    //                temp += TIFFGetG(raster[i*paramPack->xScale_ + ii + (j * paramPack->yScale_ + jj) * x]);
                    //            }
                    //        }
                    //        image->operator ()(i, j, nCur) = (unsigned short)(double(temp) * ratio);//2015-6-8//int(double(temp) * ratio);
                    //    }
                    //}
                    //else{//no re sample
                        for (int i = 0; i < x; ++i)
                        for (int j = 0; j < y; ++j){
                            image->operator ()(i, j, nCur) = TIFFGetG(raster[i + j * x]);//2015-6-8
                            
                            //image->operator ()(i, j, nCur) = 4*TIFFGetG(raster[i + j * x]);//2015-6-8
                            //temp = TIFFGetG(raster[i + j * x]);
                            //image->operator ()(i,j,nCur) = 4 * temp;
                        }
                    //}
                    TIFFReadDirectory(in);//read next page
             }
            delete[] raster;//clear stack memory
        }  else{
            //double ratio = 1.0 / double(paramPack->xScale_ * paramPack->yScale_);//2015-6-8 //the max value of read image must be 1020 = 255 * 4.
            uint16 *data = new uint16[x*y];
            int width(0), height(0);
            TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &width);
            TIFFGetField(in, TIFFTAG_IMAGELENGTH, &height);
            memset(data, 0, sizeof(short)* width * height);
            TIFFSetDirectory(in, 0);
            int StripSize = TIFFStripSize(in);
            //int Vstripsize = TIFFVStripSize(in, height);
            //int rowsperpix = TIFFComputeStrip(in, height, 0);
            int NumberofStrips = TIFFNumberOfStrips(in);
            for (int nCur = 0; nCur < z; ++nCur) {
                for (int i = 0; i < NumberofStrips; ++i) {
                    TIFFReadEncodedStrip(in, i, (data + i * StripSize / sizeof(short)), StripSize);//
                }
                //if (paramPack->xScale_ * paramPack->yScale_ != 1){//re sample
                //    for (int i = 0; i < x; ++i)//
                //    for (int j = 0; j < y; ++j){
                //        temp = 0;
                //        //sum neighborhood value, do not divide
                //        for (int ii = 0; ii < paramPack->xScale_; ++ii){
                //            for (int jj = 0; jj < paramPack->yScale_; ++jj){
                //                temp += data[i*paramPack->xScale_ + ii + (j * paramPack->yScale_ + jj) * x];
                //            }
                //        }
                //        image->operator ()(i, j, nCur) = (unsigned short)(double(temp) * ratio);//2015-6-8//int(double(temp) * ratio);
                //    }
                //}
                //else{//no re sample
                    for (int i = 0; i < x; ++i)
                    for (int j = 0; j < y; ++j){
                        image->operator ()(i, j, nCur) = data[i + j * x];//2015-6-8
                    }
                //}
                TIFFReadDirectory(in);//read next page
            }
            delete[] data;
        }
        
        TIFFClose(in);

        //----- fill black area
        //2015-6-8
        //        for(int k = 0; k < z; ++k){
        //            GetHist(k);//calc background value
        //            FillBlackArea(k);//fill k-th page
        //        }
        return true;
    }
}

bool ImageReader::GetImageInfo(const std::string& string)
{
    if (!paramPack) return false;
    TIFFSetWarningHandler(0);//do not show warning dialog
    TIFF* in = TIFFOpen(string.c_str(), "r");
    if (!in)
        return false;
    uint16 bitspersample;
    int x, y, z;
    TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &x);
    TIFFGetField(in, TIFFTAG_IMAGELENGTH, &y);
    z = TIFFNumberOfDirectories(in);
    TIFFClose(in);
    paramPack->xRangeMax_ = x;
    paramPack->yRangeMax_ = y;
    paramPack->zRangeMax_ = z;
    dataType = bitspersample == 8 ? DATATYPE::IMAGE8 : DATATYPE::IMAGE16;
    return true;
}

//the following function is no use

//--calc background value--//
//void ImageReader::GetHist(int z)
//{
//    //initial
//#ifdef _WIN32
//	VolumePointer image = std::tr1::dynamic_pointer_cast<SVolume>( m_Source );
//#else
//	VolumePointer image = std::dynamic_pointer_cast<SVolume>( m_Source );
//#endif
    
//    int x = image->x();
//    int y = image->y();
//    memset(histogram, 0, sizeof(int) * 1021);//the num of pixel.
//    //make histgram
//    for(int i = 0; i < x; ++i){
//        for (int j = 0; j < y; ++j){
//            ++histogram[image->GetPixel(i,j,z)];
//        }
//    }
//}

//-- fill z - th pages black area--//
//void ImageReader::FillBlackArea(int z)
//{
//    //initial
//#ifdef _WIN32
//	VolumePointer image = std::tr1::dynamic_pointer_cast<SVolume>( m_Source );
//#else
//	VolumePointer image = std::dynamic_pointer_cast<SVolume>( m_Source );
//#endif
    
//    int x = image->x();
//    int y = image->y();
//    //int z = image->z();
//    //get fill color
//    int sumIndex = 0;
//    int maxPixel = 0;
//    int maxPixelSum = 0;
//    int newMaxPixel = 255 * 4;
//    for(int i = 1; i < newMaxPixel; ++i){//donot care about 0
//        if(sumIndex > 10) break;
//        if(histogram[i] > maxPixel){
//            ++sumIndex;
//            maxPixelSum = histogram[i];
//            maxPixel = i;
//        }
//    }
//    //fill black area
//    for(int i = 0; i < x; ++i){
//        for (int j = 0; j < y; ++j){
//            //for (int ij = 0; ij < z; ++ij){
//                if(image->GetPixel(i,j,z) == 0)
//                    image->GetPixel(i,j,z) = maxPixel;
//            //}
//        }
//    }
//    printf("black area of z: %d is painted with color: %d\n", z, maxPixel);
//}
