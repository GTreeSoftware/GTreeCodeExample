#ifndef DENDRITELEVELIMAGEMAKER_H
#define DENDRITELEVELIMAGEMAKER_H
#include "../ngtypes/basetypes.h"
#include "../ngtypes/volume.h"
#include "../ngtypes/tree.h"
#include "../ngtypes/ParamPack.h"
#include "ineuronprocessobject.h"

class DendriteLevelImageMaker :public INeuronProcessObject{
public:
    DendriteLevelImageMaker();
    ~DendriteLevelImageMaker();
   INEURONPROCESSOBJECT_DEFINE
    void SetParam(NGParamPack arg){ paramPack = arg; }
    void SetLayer(const std::vector<int>& arg){ layerList_ = arg; }
    void SetBranch(const std::vector<int>& arg){ branchList_ = arg; }
    void SetLayerRadius(int arg){ radius_ = arg; }
    void SetHasSoma(bool arg){ hasSoma_ = arg; }
protected:
private:
    bool hasSoma_;
    NGParamPack paramPack;
    int radius_;
    std::vector<int> layerList_;
    std::vector<int> branchList_;
};

#endif