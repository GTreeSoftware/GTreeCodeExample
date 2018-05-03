
#include "H5Cpp.h"
#include "HDF5Reader.h"
#include "imagewriter.h"
#include "ngtypes/ParamPack.h"
using namespace H5;

HDF5Reader::HDF5Reader()
{
    className_ = std::string("HDF5Reader");
}


HDF5Reader::~HDF5Reader()
{
}

ProcStatPointer HDF5Reader::Update()
{
    H5File *file = new H5File(filename.c_str(), H5F_ACC_RDONLY);
    if (!file) {
        MAKEPROCESSSTATUS(resSta, false, className_, "Cannot create H5");
        return resSta;
    }
	if (paramPack_->xMin_<0 || paramPack_->yMin_<0 || paramPack_->zMin_<0 )
	{
		MAKEPROCESSSTATUS(resSta, false, className_, "Image coordinate is less than zero");
		return resSta;
	}
	if (paramPack_->xMax_ - paramPack_->xMin_<10 || paramPack_->yMax_ - paramPack_->yMin_<10 || paramPack_->zMax_ - paramPack_->zMin_<10  )
	{
		MAKEPROCESSSTATUS(resSta, false, className_, "Image is too small");
		return resSta;
	}
	if (paramPack_->xMax_ >= paramPack_->xRangeMax_ || paramPack_->yMax_ >= paramPack_->yRangeMax_ || paramPack_->zMax_ >= paramPack_->zRangeMax_)
	{
		MAKEPROCESSSTATUS(resSta, false, className_, "Image Size is out of range");
		return resSta;
	}
    //DSetCreatPropList cparms;
    //FileAccPropList fapl = FileAccPropList();
    //fapl.setCache(100, 521, 256 * 256, 0.5);
    
    H5Pset_cache(file->getAccessPlist().getId(), 0, 521, 256 * 256, 0.5);
    char line[256];
    ////read size
    //{
    //    DataSet dataset = file->openDataSet("size");
    //    DataSpace filespace = dataset.getSpace();
    //    int rank = filespace.getSimpleExtentNdims();
    //    hsize_t *dims = new hsize_t[rank];
    //    rank = filespace.getSimpleExtentDims(dims);
    //    DataSpace mspace1(rank, dims);
    //    int sz[3];
    //    dataset.read(sz, PredType::NATIVE_INT, mspace1, filespace);
    //    paramPack_->xRangeMax_ = sz[0];
    //    paramPack_->yRangeMax_ = sz[1];
    //    paramPack_->zRangeMax_ = sz[2];
    //}
    //{//read resolution
    //    DataSet dataset = file->openDataSet("resolution");
    //    DataSpace filespace = dataset.getSpace();
    //    int rank = filespace.getSimpleExtentNdims();
    //    hsize_t *dims = new hsize_t[rank];
    //    rank = filespace.getSimpleExtentDims(dims);
    //    DataSpace mspace1(rank, dims);
    //    double res[3];
    //    dataset.read(res, PredType::NATIVE_DOUBLE, mspace1, filespace);
    //    paramPack_->xRes_ = res[0];
    //    paramPack_->yRes_ = res[1];
    //    paramPack_->zRes_ = res[2];
    //}
    
    {//read data
        if (!m_Source) NG_SMART_POINTER_NEW_DEFAULT(SVolume, m_Source);
        NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, m_Source);
        int den = int(std::pow(2.0, paramPack_->mostdLevel_ - 1));
        //tmpOrig->SetSize((paramPack_->xMax_ - paramPack_->xMin_ + 1)/den, (paramPack_->yMax_ - paramPack_->yMin_ + 1)/den,
        //    (paramPack_->zMax_ - paramPack_->zMin_ + 1)/den);
		tmpOrig->SetSize((paramPack_->xMax_ - paramPack_->xMin_ + 1) / den, (paramPack_->yMax_ - paramPack_->yMin_ + 1) / den,
			(paramPack_->zMax_ - paramPack_->zMin_ + 1) );
        tmpOrig->SetDataType(dataType);
        //tmpOrig->SetResolution(paramPack_->xRes_*double(den), paramPack_->yRes_*double(den), paramPack_->zRes_*double(den));
		tmpOrig->SetResolution(paramPack_->xRes_*double(den), paramPack_->yRes_*double(den), paramPack_->zRes_);
        //int zBeg = paramPack_->zMin_ / den, zEnd = zBeg + tmpOrig->z();
		int zBeg = std::min(paramPack_->zMin_, paramPack_->zRangeMax_ - 1), zEnd = std::max(zBeg,std::min(zBeg + tmpOrig->z(), paramPack_->zRangeMax_ - 1));
        for (int z = zBeg; z < zEnd; ++z) {
            sprintf_s(line, "/level%d", paramPack_->mostdLevel_);
            if (!file->exists(line)){
                tmpOrig->Clear();
                MAKEPROCESSSTATUS(resSta, false, className_, "There is no this level data in H5");
                return resSta;
            }
            sprintf_s(line, "/level%d/%d", paramPack_->mostdLevel_, z);
            DataSet dataset = file->openDataSet(line);
            /*DataSpace filespace = dataset.getSpace();
            int rank = filespace.getSimpleExtentNdims();
            hsize_t *dims = new hsize_t[rank];
            rank = filespace.getSimpleExtentDims(dims);
            DataSpace mspace1(rank, dims);*/
            hsize_t offset[2] = { paramPack_->yMin_ / den, paramPack_->xMin_ / den };
            hsize_t count[2] = { tmpOrig->y(), tmpOrig->x() };
            //hsize_t stride[2] = { 1, 1 };
            //hsize_t block[2] = { 1, 1 };
            DataSpace memspace(2, count, NULL);
            DataSpace dataspace = dataset.getSpace();
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);//, stride, block);
            unsigned short **data = NULL;
            unsigned short *origPointer = tmpOrig->GetPointer() + (z - paramPack_->zMin_) *tmpOrig->x() * tmpOrig->y();
            Transfer1DPointerTo2DPointer_USHORT(origPointer, tmpOrig->x(), tmpOrig->y(), &data);
            try{
                dataset.read(data[0], PredType::NATIVE_USHORT, memspace, dataspace);
            }
            catch (DataSetIException &error){
                error.printError();
            }
        }
    }
    //
    file->close();
    delete file;
    //test
    /*NGImageWriter writer = ImageWriter::New();
    writer->SetInput(paramPack_->OrigImage);
    writer->SetOutputFileName("E:/hehe.tif");
    writer->Update();*/
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

void HDF5Reader::Transfer1DPointerTo2DPointer_USHORT(unsigned short *data, int x, int y, unsigned short*** dst)
{
    *dst = new unsigned short*[y];
    for (int k = 0; k < y; ++k) {
        (*dst)[k] = data + k * x;
    }
}

INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(HDF5Reader)

INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(HDF5Reader, SVolume)

bool HDF5Reader::SetInputFileName(const std::string& arg)
{
    filename = arg; 
    H5File *file = new H5File(filename.c_str(), H5F_ACC_RDONLY);
    if (!file) {
        return false;
    }
    //char line[256];
    //data type
    {
        DataSet dataset = file->openDataSet("dataType");
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();
        hsize_t *dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);
        DataSpace mspace1(rank, dims);
        int datatype[1];
        dataset.read(datatype, PredType::NATIVE_INT, mspace1, filespace);
        dataType = datatype[0] == 1 ? DATATYPE::IMAGE8 : DATATYPE::IMAGE16;
    }
    //read size
    {
        DataSet dataset = file->openDataSet("size");
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();
        hsize_t *dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);
        DataSpace mspace1(rank, dims);
        int sz[3];
        dataset.read(sz, PredType::NATIVE_INT, mspace1, filespace);
        paramPack_->xRangeMax_ = sz[0];
        paramPack_->yRangeMax_ = sz[1];
        paramPack_->zRangeMax_ = sz[2];
    }
    {//read resolution
        DataSet dataset = file->openDataSet("resolution");
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();
        hsize_t *dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);
        DataSpace mspace1(rank, dims);
        double res[3];
        dataset.read(res, PredType::NATIVE_DOUBLE, mspace1, filespace);
        paramPack_->xRes_ = res[0];
        paramPack_->yRes_ = res[1];
        paramPack_->zRes_ = res[2];
    }
    file->close();
    return true;
}

void HDF5Reader::GetMaxRange(int& x, int& y, int& z)
{
    x = paramPack_->xRangeMax_;
    y = paramPack_->yRangeMax_;
    z = paramPack_->zRangeMax_;
}
