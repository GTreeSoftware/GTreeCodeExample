/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#include <tiffio.h>
#include <cstdlib>
#include "imagewriter.h"
#include "../ngtypes/basetypes.h"
#include "../ngtypes/volume.h"

ImageWriter::ImageWriter(void):isFileValid(false)
{
	className_ = std::string("ImageWriter");
	m_Source = std::shared_ptr<SVolume>(new SVolume(this));
}

ImageWriter::~ImageWriter(void)
{
}

bool ImageWriter::SetOutputFileName(const std::string& path)
{
	filename = path;
	isFileValid = true;
    return true;
}

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(ImageWriter, SVolume)
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(ImageWriter)

ProcStatPointer ImageWriter::Update()
{
    if(isFileValid && m_Input ){
        if (m_Input->GetProcessObject()){
            ProcStatPointer res = m_Input->GetProcessObject()->Update();
            if (!res->success()){
                printf("error occurred in %s\n", className_.c_str());
                //LOG(ERROR) << "error occurred in "<<className_<<" when check m_Input.";
                return res;
            }
        }
        if (WriteImage(filename.c_str())){
            MAKEPROCESSSTATUS(resSta, true, className_, "");
            return resSta;
        }
        //printf("can not read image file.this is %s", identifyName.c_str());
    }
	NG_ERROR_MESSAGE("error occurred in ImageWriter.");
    MAKEPROCESSSTATUS(resSta, false, className_, "Cannot write image.");
    return resSta;
}

bool ImageWriter::WriteImage(const char *path)
{
    CVolumePointer imageC = std::const_pointer_cast<CVolume>(std::dynamic_pointer_cast<const CVolume>(m_Input));
    if (imageC) return WriteTiff(path, imageC);
    else {
        SVolumePointer imageS = std::const_pointer_cast<SVolume>(std::dynamic_pointer_cast<const SVolume>(m_Input));
        return WriteTiff(path, imageS);
    }
}

bool ImageWriter::WriteTiff(const char* path, CVolumePointer &data)
{
    TIFFSetWarningHandler(0);
    TIFF* out = TIFFOpen(path, "w");
    if (out)
    {
        size_t nWidth = data->x();
        size_t nHeight = data->y();
        size_t nTotalFrame = data->z();
        //size_t wh = nWidth * nHeight;
        unsigned char *temp = new unsigned char[nWidth * nHeight];
        size_t nCur = 0;
        do
        {
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(out, TIFFTAG_PAGENUMBER, nTotalFrame);
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32)nWidth);
            TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32)nHeight);
            //TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, 2);
            /*TIFFSetField(out, TIFFTAG_YRESOLUTION, 196.0f);
            TIFFSetField(out, TIFFTAG_XRESOLUTION, 204.0f);*/
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            //
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, nHeight);

            for (size_t i = 0; i < nWidth; ++i){
                for (size_t j = 0; j < nHeight; ++j){
                    temp[i + j * nWidth] = data->operator ()(i, j, nCur);// / den;//不支持sint
                }
            }

            TIFFWriteEncodedStrip(out, 0, temp, tsize_t(nWidth * nHeight));
            ++nCur;
        } while (TIFFWriteDirectory(out) && nCur < nTotalFrame);
        delete[] temp;
        TIFFClose(out);
        return true;
    }
    else
    {
        printf("cant create tiff!\n");
        return false;
    }
}

bool ImageWriter::WriteTiff(const char* path, SVolumePointer &data)
{
    TIFFSetWarningHandler(0);
    TIFF* out = TIFFOpen(path, "w");
    if (out)
    {
        size_t nWidth = data->x();
        size_t nHeight = data->y();
        size_t nTotalFrame = data->z();
        //size_t wh = nWidth * nHeight;
        unsigned short *temp = new unsigned short[nWidth * nHeight];
        size_t nCur = 0;
        do
        {
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(out, TIFFTAG_PAGENUMBER, nTotalFrame);
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, (uint32)nWidth);
            TIFFSetField(out, TIFFTAG_IMAGELENGTH, (uint32)nHeight);
            //TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, 2);
            /*TIFFSetField(out, TIFFTAG_YRESOLUTION, 196.0f);
            TIFFSetField(out, TIFFTAG_XRESOLUTION, 204.0f);*/
            TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            //
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 16);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);
            TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(out, TIFFTAG_ROWSPERSTRIP,nHeight);

            for (size_t i = 0; i < nWidth;++i){
                for(size_t j = 0; j < nHeight; ++j){
                    temp[i + j * nWidth] = data->operator ()(i,j,nCur);// / den;//不支持sint
                }
            }

            TIFFWriteEncodedStrip(out, 0, temp, tsize_t(2 * nWidth * nHeight));
            ++nCur;
        }
        while(TIFFWriteDirectory(out) && nCur < nTotalFrame);
        delete[] temp;
        TIFFClose(out);
        return true;
    }
    else
    {
        printf("cant create tiff!\n");
        return false;
    }
}

