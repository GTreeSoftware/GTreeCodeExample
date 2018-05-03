#pragma once
#include "ineuronprocessobject.h"
#include "../ngtypes/volume.h"
#include "../ngtypes/ParamPack.h"

class ImageScaleFiter;
typedef std::shared_ptr<ImageScaleFiter> NGImageScaleFiter;

class ImageScaleFiter :
    public INeuronProcessObject
{
public:
    static NGImageScaleFiter New(){ return NGImageScaleFiter(new ImageScaleFiter()); }
    ImageScaleFiter();
    virtual ~ImageScaleFiter();
    void SetParam(NGParamPack arg){ paramPack_ = arg; }
    void SetInputScale(NGParamPack arg){ paramPack_ = arg; }

    INEURONPROCESSOBJECT_DEFINE

private:
    NGParamPack paramPack_;

};

