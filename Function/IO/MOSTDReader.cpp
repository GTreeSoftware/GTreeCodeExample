/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang, liyuxin
*	2015-10-28
*/
#include "MOSTDReader.h"
#include "mostd/most_data.h"
#include "mostd/most_data16.h"
#include "ngtypes/volume.h"
#pragma warning (disable: 4996)

MOSTDReader::MOSTDReader()
{
    reader = new MostImage;
    reader16 = new MostImage16;
    //dataType_ = IMAGE8;
    //level_ = 1;
}

MOSTDReader::~MOSTDReader()
{
    if(reader->IORdata)
        delete[] reader->IORdata;
    if(reader16->IORdata)
        delete[] reader16->IORdata;
    delete reader;
    delete reader16;
}

bool MOSTDReader::SetInputFileName( const std::string& fileName)
{
    if (!reader || !reader16 || !paramPack) {
        printf("error mostd reader.\n");
        return false;
    }
    FILE *fp = fopen(fileName.c_str(), "r");
    if(!fp) return false;
    int id, datatype;
    if (!feof(fp)){
        fscanf(fp, "%d\n",&id);
        fscanf(fp, "%d\n", &datatype);//second data is data type
    }
    fclose(fp);
    dataType = datatype == 1 ? DATATYPE::IMAGE8 : DATATYPE::IMAGE16;
    if (dataType == DATATYPE::IMAGE8)
    {
        reader->loadDataInfo(fileName.c_str());
        SetResolution(reader->rez_x, reader->rez_y, reader->rez_z );
        paramPack->xRangeMax_ = reader->sz0;
        paramPack->yRangeMax_ = reader->sz1;
        paramPack->zRangeMax_ = reader->sz2;
        paramPack->maxValue_ = 255;
    }
    else if (dataType == DATATYPE::IMAGE16){
        reader16->loadDataInfo(fileName.c_str());
        paramPack->xRangeMax_ = reader16->sz0;
        paramPack->yRangeMax_ = reader16->sz1;
        paramPack->zRangeMax_ = reader16->sz2;
        SetResolution(reader16->rez_x, reader16->rez_y, reader16->rez_z );
        paramPack->maxValue_ = 4096;
    }
    return true;
}

void MOSTDReader::SetResolution( double x, double y, double z)
{
    //xScale_ = 1;// SetScale(x);
    //yScale_ = 1;// SetScale(y);
    //zScale_ = 1;// SetScale(z);
    paramPack->xRes_ = x;
    paramPack->yRes_ = y ;
    paramPack->zRes_ = z ;
}

ProcStatPointer MOSTDReader::Update()
{
    if (!UpdateROI()){
        printf("error mostd reader in UpdateROI.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "error mostd reader.");
        return resSta;
    }
    //no most image or cannot read
    // choose most reader or mostreader16 to read image
    if ((dataType == DATATYPE::IMAGE8 && (!reader || !reader->IORdata)) && (dataType != DATATYPE::IMAGE8 && (!reader16 || !reader16->IORdata))) {
        printf("error mostd reader.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "error mostd reader.");
        return resSta;
    }
    //
    if(!TransferMOSTD2Volume()){//must be transfer into volume structure
        printf("cannot transfer mostd data into volume.\n");
        MAKEPROCESSSTATUS(resSta, false, className_, "cannot transfer mostd data into volume.");
        return resSta;
    }
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

IDataPointer MOSTDReader::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData(m_Source);
    m_Source.reset();
    return tData;
}

ConstIDataPointer MOSTDReader::GetOutput()
{
    if(!m_Source)
        m_Source = std::shared_ptr<SVolume>(new SVolume(this));
    return m_Source;
}

bool MOSTDReader::UpdateROI()
{
    if (!reader && !reader16) {
        printf("error mostd reader.\n");
        return false;
    }
    // check if out of boundary
    if(paramPack->xMax_ <= paramPack->xMin_ || paramPack->yMax_ <= paramPack->yMin_ 
        || paramPack->zMax_ <= paramPack->zMin_ || paramPack->xMax_ <0 || paramPack->xMin_ < 0 || paramPack->yMax_ < 0 || 
        paramPack->yMin_ < 0 || paramPack->zMax_ < 0 || paramPack->zMin_ < 0){
            printf("error range.\n");
            return false;
    }
    //reader ROI image
    if (dataType == DATATYPE::IMAGE8){
        reader->select_IOR_file(paramPack->xMin_, paramPack->xMax_, 
            paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_, paramPack->mostdLevel_);//always the original image level
        SetResolution(reader->rez_x, reader->rez_y, reader->rez_z);
    }else{
        reader16->select_IOR_file(paramPack->xMin_, paramPack->xMax_, 
            paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_, paramPack->mostdLevel_);//always the original image level
        SetResolution(reader16->rez_x, reader16->rez_y, reader16->rez_z);
    }
    return true;
}

bool MOSTDReader::TransferMOSTD2Volume()
{
    // check original data
    if(!reader->IORdata && !reader16->IORdata){
        printf("error mostd reader.\n");
        return false;
    }
    int xLen, yLen, zLen;
    //check output
    if (dataType == DATATYPE::IMAGE8) {
        if(!m_Source) NG_SMART_POINTER_NEW(SVolume, m_Source, this);
        xLen = reader->IOR_sz0;
        yLen = reader->IOR_sz1;
        zLen = reader->IOR_sz2;
    }else{
        if(!m_Source) NG_SMART_POINTER_NEW(SVolume, m_Source, this);
        xLen = reader16->IOR_sz0;
        yLen = reader16->IOR_sz1;
        zLen = reader16->IOR_sz2;
    }

    //2015-8-13
    /* int rx, ry, rz;
     rx = xLen / paramPack->xScale_;
     ry = yLen / paramPack->yScale_;
     rz = zLen;*/

    //double ratio = 1.0 / double(paramPack->xScale_ * paramPack->yScale_);//2015-6-8 //the max value of output image must be 255.

    //int temp(0);
    //int xyLen = xLen * yLen;
    double levelRes = std::pow(2.0, paramPack->mostdLevel_-1);
    if (dataType == DATATYPE::IMAGE8) {//read 8bit
        NG_CREATE_DYNAMIC_CAST(SVolume, image, m_Source);
        image->SetResolution(paramPack->xRes_ * levelRes, paramPack->yRes_ * levelRes, paramPack->zRes_* levelRes);
        image->Swap(&(reader->IORdata), xLen, yLen, zLen);
        image->SetDataType(DATATYPE::IMAGE8);
        if(reader->IORdata) delete[] reader->IORdata;//clear mostd data for memory
        reader->IORdata = NULL;
    }else{//read 16bit, must set max value, because the max value maybe 4096 instead of 65535
        //ratio *=  255.0 / maxValue_;
        NG_CREATE_DYNAMIC_CAST(SVolume, image, m_Source);
        image->SetResolution(paramPack->xRes_ * levelRes, paramPack->yRes_ * levelRes, paramPack->zRes_* levelRes);
        image->Swap(&(reader16->IORdata), xLen, yLen, zLen);
        image->SetDataType(DATATYPE::IMAGE16);
        if(reader16->IORdata) delete[] reader16->IORdata;//clear mostd data for memory
        reader16->IORdata = NULL;
    }
    return true;
}

//the size of mostd data
void MOSTDReader::GetMaxRange( int& x, int& y, int& z )
{
    if (dataType == DATATYPE::IMAGE8){
        x = reader->sz0;
        y = reader->sz1;
        z = reader->sz2;
    }else{
        x = reader16->sz0;
        y = reader16->sz1;
        z = reader16->sz2;
    }
}

//void MOSTDReader::SetRealResolution( double res, int origSz, int& scale, int& sz )
//{
//	//2015-8-13
//	if(res > 0.25 && res <= 0.5 ) {
//		scale = 2;
//		sz = origSz / scale;
//	}else if(res <= 0.25 ) {
//		scale = 3;
//		sz = origSz / scale;
//	} else{
//		scale = 1;
//		sz = origSz;
//	}
//}

//int MOSTDReader::SetScale( double res )
//{
//    //2015-8-13
//    int scale;
//    if(res > 0.25 && res <= 0.5 ) {
//        scale = 2;
//    }else if(res <= 0.25 ) {
//        scale = 3;
//    } else{
//        scale = 1;
//    }
//    return scale;
//}

void MOSTDReader::ClearCache()
{
    if (!reader && !reader16) {
        printf("error mostd reader.\n");
        return;
    }

    if (dataType == DATATYPE::IMAGE8){
        //reader->
        reader->imageCaches_.ClearCache();
    }
    else{
        reader16->imageCaches_.ClearCache();
    }
}


