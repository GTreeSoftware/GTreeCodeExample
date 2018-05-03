#include "ImageScaleFiter.h"


ImageScaleFiter::ImageScaleFiter()
{
    className_ = std::string("ImageScaleFiter");
    NG_SMART_POINTER_NEW(SVolume, m_Source, this);
}


ImageScaleFiter::~ImageScaleFiter()
{
    
}

ProcStatPointer ImageScaleFiter::Update()
{
    if (!m_Input) {
        NG_ERROR_MESSAGE("no input image data.");
        MAKEPROCESSSTATUS(resSta, false, className_, "There is no image data.");
        return resSta;
    }
    if (paramPack_->xScale_ < 1 && paramPack_->yScale_ < 1){
        NG_ERROR_MESSAGE("Up sampling function is not developed.");
        MAKEPROCESSSTATUS(resSta, false, className_, "Up sampling function is not developed.");
        return resSta;
    }
    //
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpSource, m_Source);
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpImg, m_Input);
    if (paramPack_->xScale_ == 1 && paramPack_->yScale_ == 1) {
        if (m_Input->DataType() == DATATYPE::IMAGE16) m_Source = m_Input;//do nothing
        else{//8bit image will multiply 4
            tmpSource->QuickCopy(*tmpImg);
            int rx, ry, rz;
            rx = tmpImg->x();
            ry = tmpImg->y();
            rz = tmpImg->z();
            for (int ij = 0; ij < rz; ++ij) {
                for (int i = 0; i < rx; ++i) {
                    for (int j = 0; j < ry; ++j) {
                        tmpSource->operator()(i, j, ij) *= 4;
                    }
                }
            }
        }
        tmpSource->SetResolution(tmpImg->XResolution(), tmpImg->YResolution(), tmpImg->ZResolution());
    }
    else{
        int rx, ry, rz;
        rx = tmpImg->x() / paramPack_->xScale_;
        ry = tmpImg->y() / paramPack_->yScale_;
        rz = tmpImg->z();
        int temp(0);
        int xLen = tmpImg->x();
        int xyLen = tmpImg->x()  *  tmpImg->y();
        tmpSource->SetSize(rx, ry, rz);
        tmpSource->SetResolution(tmpImg->XResolution() * paramPack_->xScale_, tmpImg->YResolution() * paramPack_->yScale_, tmpImg->ZResolution());
        if (tmpImg->DataType() == DATATYPE::IMAGE16) {
            double ratio = 1.0 / double(paramPack_->xScale_ * paramPack_->yScale_);
            for (int ij = 0; ij < rz; ++ij) {
                for (int i = 0; i < rx; ++i) {
                    for (int j = 0; j < ry; ++j) {
                        temp = 0;
                        for (int ii = 0; ii < paramPack_->xScale_; ++ii){
                            for (int jj = 0; jj < paramPack_->yScale_; ++jj){
                                temp += tmpImg->GetPointer()[i*paramPack_->xScale_ + ii + (j * paramPack_->yScale_ + jj) * xLen + ij * xyLen];
                            }
                        }
                        tmpSource->operator ()(i, j, ij) = unsigned short(temp * ratio);
                    }
                }
            }
        }
        else{//8 bit image
            double ratio = 4.0 / double(paramPack_->xScale_ * paramPack_->yScale_);
            for (int ij = 0; ij < rz; ++ij) {
                for (int i = 0; i < rx; ++i) {
                    for (int j = 0; j < ry; ++j) {
                        temp = 0;
                        for (int ii = 0; ii < paramPack_->xScale_; ++ii){
                            for (int jj = 0; jj < paramPack_->yScale_; ++jj){
                                temp += tmpImg->GetPointer()[i*paramPack_->xScale_ + ii + (j * paramPack_->yScale_ + jj) * xLen + ij * xyLen];
                            }
                        }
                        tmpSource->operator ()(i, j, ij) = unsigned short(temp * ratio);
                    }
                }
            }
        }
    }

    //adjust illumination
    if (paramPack_->highAdjustOpac_!=0    && paramPack_->highDestOpac_ != 0) {//adjust
        tmpSource->Adjust(paramPack_->lowAdjustOpac_, paramPack_->highAdjustOpac_, paramPack_->lowDestOpac_, paramPack_->highDestOpac_);
    }

    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(ImageScaleFiter)
INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(ImageScaleFiter, SVolume)
