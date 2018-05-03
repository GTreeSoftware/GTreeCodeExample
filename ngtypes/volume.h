//time : 2013-11-20
//this file provide a structure to save image data and temp image data
//author : zhouhang, WLNO

#ifndef VOLUME_H
#define VOLUME_H
#include "ineurondataobject.h"
#include "volumecreator.h"
typedef unsigned char NGCHAR;
template<typename T> class Volume;
typedef Volume<NGCHAR> CVolume;
typedef Volume<unsigned short> SVolume;
typedef Volume<int> IVolume;
typedef Volume<double> DVolume;
template<typename T>
class Volume : public INeuronDataObject
{
public:
    Volume(INeuronProcessObject* p = NULL);
    //Volume& operator=(const VolumeBlock&);
    ~Volume();

    bool IsEmpty()const{return dataptr_ == NULL;}
    /*T& operator()(int x, int y, int z){return dataptr_[x + y * x_ + z * xy_];}
    T  operator()(int x, int y, int z)const{return dataptr_[x + y * x_ + z * xy_];}
    T& GetPixel(int x, int y, int z){return dataptr_[x + y * x_ + z * xy_];}
    T GetPixel(int x, int y, int z)const{return dataptr_[x + y * x_ + z * xy_];}*/
    T& operator()(size_t x, size_t y, size_t z){ return dataptr_[x + y * x_ + z * xy_]; }
    T  operator()(size_t x, size_t y, size_t z)const{ return dataptr_[x + y * x_ + z * xy_]; }
    T& GetPixel(size_t x, size_t y, size_t z){ return dataptr_[x + y * x_ + z * xy_]; }
    T GetPixel(size_t x, size_t y, size_t z)const{ return dataptr_[x + y * x_ + z * xy_]; }
    void SetSize(int,int, int);//original data will lost
    void SetResolution(double, double, double);
    void Clear();
    void QuickCopy(const Volume<T>&);
    void SetZero();
    T* GetPointer(){return dataptr_;}
	const T* GetPointer()const{ return dataptr_; }
    void Swap(T**, int, int, int);
    void ComputeVolumeMinMax();
    void Adjust(int, int, int, int);
//    void SetProcessObject(INeuronProcessObject*);
//    void ReleaseProcessObject();

    int x()const{return x_;}
    int y()const{return y_;}
    int z()const{return z_;}
    T MinVal()const{ return minVal_; }
    T MaxVal()const{ return maxVal_; }
    T& MinVal(){ return minVal_; }
    T& MaxVal(){ return maxVal_; }
    double XResolution()const{return xRes_;}
    double YResolution()const{return yRes_;}
    double ZResolution()const{return zRes_;}

private:
    //data
    int x_ = 0 , y_ = 0 , z_ = 0 , xy_ = 0 ;//size
    T minVal_ = 0, maxVal_ = 0;
    double xRes_ = 1.0, yRes_ = 1.0, zRes_ = 1.0;
    T *dataptr_ = NULL;

    //private function
    Volume(const Volume&);//no copy assign
    Volume& operator=(const Volume&);//no copy
};

template<typename T>
void Volume<T>::Adjust(int origLow, int origHigh, int destLow, int destHigh)
{
    int index = 0;
    int origLen = origHigh - origLow;
    int destLen = destHigh - destLow;
    double ratio = destLen / origLen;
    for (int i = 0; i < x_; ++i) {
        for (int j = 0; j < y_; ++j) {
            for (int ij = 0; ij < z_; ++ij) {
                index = i + j * x_ + ij * xy_;
                if (dataptr_[index] <= origLow) dataptr_[index] = destLow;
                else if (dataptr_[index] >= origHigh) dataptr_[index] = destHigh;
                else{
                    dataptr_[index] = unsigned short(double(dataptr_[index] - origLow) * ratio) + destLow;
                }
            }
        }
    }
}

template<class T>
Volume<T>::Volume(INeuronProcessObject* p)
{
    m_ProcessObject =  p;
    identifyName =  std::string("Volume");
    x_ = 0;
    y_ = 0;
    z_ = 0;
    xy_ = 0;//size
    xRes_=1.0;
    yRes_=1.0;
    zRes_=1.0;
    if (sizeof(T) == 2) dataType = DATATYPE::IMAGE16;
    else if (sizeof(T) == 1) dataType = DATATYPE::IMAGE1;
    dataptr_ = NULL;
}

template<class T>
Volume<T>::~Volume()
{
    Destroy2DArray<T>(&dataptr_, z_, xy_);
}

template<class T>
void Volume<T>::Swap(T** p, int x, int y, int z)
{
    std::swap(*p, dataptr_);
    x_ = x;
    y_ = y;
    z_ = z;
    xy_ = x_ * y_;
}


template<class T>
void Volume<T>::SetSize(int x, int y, int z)
{
    if (x < 0 || y < 0 || z < 0) return;
    if (x_ == x && y_ == y && z_ == z) {
        memset(dataptr_, 0, sizeof(T) * xy_*z_);//2016-8-24
        return;
    }
    if (0 == (x * y * z)){
        if (dataptr_)
            Destroy2DArray<T>(&dataptr_, z_, xy_);
        x_ = x; y_ = y; z_ = z; xy_ = x_ * y_;
    }
    else{
        if (dataptr_)
            Destroy2DArray<T>(&dataptr_, z_, xy_);
        x_ = x; y_ = y; z_ = z; xy_ = x_ * y_;
        Create2DArray<T>(&dataptr_, z_, xy_);
    }
}

template<class T>
void Volume<T>::SetResolution(double x, double y, double z)
{
    xRes_ = x;
    yRes_ = y;
    zRes_ = z;
}

template<class T>
void Volume<T>::Clear()
{
    SetSize(0,0,0);
}

template<class T>
void Volume<T>::QuickCopy(const Volume<T> &src)
{
    SetSize(src.x(), src.y(), src.z());
    dataType = src.dataType;
    T* sptr = src.dataptr_;
    T* dptr = dataptr_;
    memcpy(dptr, sptr, sizeof(T) * src.x() * src.y() * src.z());
}

template<class T>
void Volume<T>::SetZero()
{
    Set2DArray<T>(dataptr_, z_, xy_, 0);
}

template<class T>
void Volume<T>::ComputeVolumeMinMax()
{
    int index = 0;
    for (int i = 0; i < x_; ++i) {
        for (int j = 0; j < y_; ++j) {
            for (int ij = 0; ij < z_; ++ij) {
                index = i + j * x_ + ij * xy_;
                minVal_ = (std::min)(minVal_, dataptr_[index]);
                maxVal_ = (std::max)(maxVal_, dataptr_[index]);
            }
        }
    }
}

//2015-6-18
struct VolumeBoundingBox{
    int xMin, xMax, yMin, yMax, zMin, zMax;
    VolumeBoundingBox(int x1, int x2, int y1, int y2, int z1, int z2):
    xMin(x1), xMax(x2), yMin(y1), yMax(y2), zMin(z1), zMax(z2){}
};

#endif // VOLUME_H
