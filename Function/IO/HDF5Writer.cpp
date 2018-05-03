
#include "H5Cpp.h"
#include "HDF5Writer.h"
#include "ngtypes/volume.h"
#include "imagereader.h"
using namespace H5;


HDF5Writer::HDF5Writer()
{
    className_ = std::string("HDF5Writer");
}


HDF5Writer::~HDF5Writer()
{
}

ProcStatPointer HDF5Writer::Update()
{
    
    char line[256];
    H5File *file = new H5File(filename.c_str(), H5F_ACC_TRUNC);//std::string has bug
    if (!file) {
        MAKEPROCESSSTATUS(resSta, false, className_, "Cannot create H5");
        return resSta;
    }
    //
    NGImageReader reader = ImageReader::New();
    reader->SetParam(paramPack_);
    reader->GetImageInfo(tiffList[0]);
    //DATATYPE dataType = reader->GetImageType();
    int xR = paramPack_->xRangeMax_;
    int yR = paramPack_->yRangeMax_;
    int zR = int(tiffList.size());
    
    //calculate pyramid level num
    int levelNum = 1;
	int maxLength = std::max(xR, yR);//xR > yR ? (xR > zR ? xR: zR) : (yR > zR ? yR : zR);
	int minLength = std::min(xR, yR);//xR < yR ? (xR < zR ? xR : zR) : (yR < zR ? yR : zR);
    //int maxLengthCp = maxLength;
    //int minLengthCp = minLength;
    while (true){
        if (minLength/2 <subblockSz_ || (minLength < 50 && maxLength < 2*subblockSz_) || minLength < 10 ) 
            break;
        maxLength /= 2;
        minLength /= 2;
        ++levelNum;
    }
    //datatype
    {
        int datatype[1] = { reader->GetImageType() == DATATYPE::IMAGE8 ? 1 :2 };
        hsize_t     dimsf[1];              // dataset dimensions
        dimsf[0] = 1;
        DataSpace dataspace(1, dimsf);
        DataSet dataset = file->createDataSet("dataType", PredType::NATIVE_INT, dataspace);
        dataset.write(datatype, PredType::NATIVE_INT);
    }
    //xyz
    {
        int xyz[3] = { xR, yR, zR };
        hsize_t     dimsf[1];              // dataset dimensions
        dimsf[0] = 3;
        DataSpace dataspace(1, dimsf);
        DataSet dataset = file->createDataSet("size", PredType::NATIVE_INT, dataspace);
        dataset.write(xyz, PredType::NATIVE_INT);
    }
    //xyz resolution
    {
        double xyz[3] = { paramPack_->xRes_, paramPack_->yRes_, paramPack_->zRes_ };
        hsize_t     dimsf[1];              // dataset dimensions
        dimsf[0] = 3;
        DataSpace dataspace(1, dimsf);
        DataSet dataset = file->createDataSet("resolution", PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(xyz, PredType::NATIVE_DOUBLE);
    }
    //level
     {
         int level[1] = { levelNum };
         hsize_t     dimsf[1];              // dataset dimensions
         dimsf[0] = 1;
         DataSpace dataspace(1, dimsf);
         DataSet dataset = file->createDataSet("level", PredType::NATIVE_INT, dataspace);
         dataset.write(level, PredType::NATIVE_INT);
     }
    // Create groups
     Group *groupList = new Group[levelNum];
     for (int k = 0; k < levelNum; ++k) {
         sprintf_s(line, "/level%d", k+1);
         groupList[k] = file->createGroup(line);
     }
     for (int k = 0; k < levelNum; ++k)
         groupList[k].close();
     //delete[] groupList;
     file->close();
     //file = new H5File(filename.c_str(), H5F_ACC_RDWR);
    //chunk
    {
        DSetCreatPropList cparms;
        hsize_t      chunk_dims[2] = { subblockSz_, subblockSz_ };
        cparms.setChunk(2, chunk_dims);
        int fill_val = 0;
        cparms.setFillValue(PredType::NATIVE_INT, &fill_val);
        cparms.setDeflate(6);
        int id = 0;

		/*QProgressDialog progressDialog;
		progressDialog.setWindowTitle("Creat HDF5");
		progressDialog.setLabelText("Creating...");
		progressDialog.setRange(0, zR);
		progressDialog.setModal(true);
		progressDialog.setCancelButtonText("Cancel");
		QApplication::processEvents();
		progressDialog.show();*/

        for (auto &tiffName : tiffList) {
            if (stopFlag_) {
                MAKEPROCESSSTATUS(resSta, false, className_, "Cancelled");
                stopFlag_ = false;
                return resSta;
            }
            emit Progress_Signal(id);
			/////////////////////////////////////////////////////////////////////////////////////////////////
			//progressDialog.setValue(id + 1);
			//QApplication::processEvents();
			//progressDialog.update();
            /*if (progressDialog.wasCanceled())
                break;*/

            printf("HDF5 process %d%% : %s\r", int(double(id+1) / double(zR) * 100), tiffName.c_str());
            reader->SetInputFileName(tiffName);
            reader->SetParam(paramPack_);
            reader->Update();
            IDataPointer img = reader->ReleaseData();
            std::shared_ptr<SVolume> tmpImg = std::dynamic_pointer_cast<SVolume>(img);
            if (!tmpImg) {
                MAKEPROCESSSTATUS(resSta, false, className_, "Image is invalid");
                return resSta;
            }
            file = new H5File(filename.c_str(), H5F_ACC_RDWR);
            
            for (int k = 0; k < levelNum; ++k) {
                sprintf_s(line, "/level%d", k + 1);
                groupList[k] = file->openGroup(line);
            }
            sprintf_s(line, "%d", id);
            unsigned short **data = NULL;
            for (int k = 1; k <= levelNum; ++k) {
                data = NULL;
                if (k == 1){
                    hsize_t      dims[2] = { tmpImg->y(),  tmpImg->x() };  // dataset dimensions at creation
                    DataSpace mspace(2, dims);
                    DataSet dataset = groupList[k - 1].createDataSet(line, PredType::NATIVE_USHORT, mspace, cparms);
                    Transfer1DPointerTo2DPointer_USHORT(tmpImg->GetPointer(), tmpImg->x(), tmpImg->y(), &data);
                    dataset.write(data[0], PredType::NATIVE_USHORT);
                    groupList[k - 1].close();
                }
                else{//pyramid
                    SVolume levelImg;
                    MakePyramid(*tmpImg, k, levelImg);
                    hsize_t      dims[2] = { levelImg.y(), levelImg.x() };  // dataset dimensions at creation
                    DataSpace mspace(2, dims);
                    DataSet dataset = groupList[k - 1].createDataSet(line, PredType::NATIVE_USHORT, mspace, cparms);
                    Transfer1DPointerTo2DPointer_USHORT(levelImg.GetPointer(), levelImg.x(), levelImg.y(), &data);
                    dataset.write(data[0], PredType::NATIVE_USHORT);
                    groupList[k - 1].close();
                }
                //delete[] data;
            }
            ++id;
            file->close();
        }
        printf("\n");
    }
    //file->close();
    delete[] groupList;
    delete file;

    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

void HDF5Writer::Transfer1DPointerTo2DPointer_USHORT(unsigned short *data, int x, int y, unsigned short*** dst)
{
    *dst = new unsigned short*[y];
    for (int k = 0; k < y; ++k) {
        (*dst)[k] = data + k * x;
    }
}

bool HDF5Writer::MakePyramid(const SVolume& orig, int level, SVolume &dst)
{
	if (level == 1)
		return false;
	int den = std::pow(2, level - 1);
	int ox = orig.x() / den;
	int oy = orig.y() / den;
	dst.SetSize(ox, oy, 1);
	int i, j;
	for (i = 0; i < dst.x(); ++i) {
		for (j = 0; j < dst.y(); ++j) {
			dst(i, j, 0) = orig(i*den, j*den, 0);
		}
	}
	return true;
}


INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(HDF5Writer, SVolume)
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(HDF5Writer)
