#ifndef CROSSBRANCHFILTER_H
#define CROSSBRANCHFILTER_H
#include "SVMTraceFilter.h"
#include "../../../ngtypes/basetypes.h"
#include "../../../ngtypes/volume.h"
#include "../../ngtypes/tree.h"
#include "ngtypes/ParamPack.h"

class Tree;
class CrossBranchFilter;
NG_SMART_POINTER_TYPEDEF(CrossBranchFilter, NGCrossBranchFilter);

class CrossBranchFilter :public INeuronProcessObject
{
public:
    static NGCrossBranchFilter New(){ return NGCrossBranchFilter(new CrossBranchFilter()); }
    CrossBranchFilter();
    ~CrossBranchFilter();

    virtual ProcStatPointer Update();

    void SetParam(NGParamPack val) { paramPack = val; }
    void SetData(const std::vector<VectorVec3d> &, const Vec3d&, const Vec3d&, std::vector<VectorVec3d> &, VectorVec3d &, VectorVec3i &);
    void SetData(const std::vector<VectorVec3d> &, const Vec3d&, const Vec3d&, std::vector<VectorVec3d> &, VectorVec3d &);
protected:
    ConstIDataPointer GetOutput(){ return m_Source; }
    IDataPointer ReleaseData(){ return IDataPointer(new TreeCurve()); }
    bool SetInputFileName(const std::string &){ return false; }
    void RemoveCrossBifurcation(const std::vector<VectorVec3d> &tmpDendrites, const std::vector<VectorVec3d>&bifurClustre, std::set<size_t> &removeCrossPairIndList);
    void FindNearestParentCurve(const Vec3d& head, const std::vector<VectorVec3d> &pCurves, size_t &curInd, size_t &ptInd);
    bool IsIntersect(const VectorVec3i &arg1, const VectorVec3i &arg2);
    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, std::vector<double> &rayNodeWet);
    void GetSVMSampleFeatureByRadiusThrev(const SVolume &origImg, Vec3i seedPt, double extractThrev, int radius, int& allPtNum, int& validSignalSum);
    void SignalPointsIdentifyAreaGrowing(VectorVec3i& growPtSet, const SVolume& subOrigImg, CVolume& signalMatrix, double threv, int = 2, int = 2);
private:
    NGParamPack paramPack;
    const std::vector<VectorVec3d> *tracedCurve_;
    const Vec3d* soma_;
    const Vec3d* initDir_;
    std::vector<VectorVec3d> *bifurcationClustre_;
    VectorVec3d *bifurcationDirections_;
    VectorVec3i* boundaryLabel_ = nullptr;
    NG_SMART_POINTER_DEFINE(const Volume<unsigned short>, origImgPointer);
    NG_SMART_POINTER_DEFINE(const Volume<unsigned short>, backImgPointer);
    NG_SMART_POINTER_DEFINE(const Volume<unsigned char>, binImgPointer);

    SVolume subOrigImg_;
    CVolume signalMatrix_;
};

#endif