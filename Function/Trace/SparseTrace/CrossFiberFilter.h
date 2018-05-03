#ifndef  CROSSFIBERFILTER_H
#define CROSSFIBERFILTER_H
#include "../../ineuronprocessobject.h"
#include "../../../ngtypes/basetypes.h"
#include "../../../ngtypes/volume.h"
#include "ngtypes/ParamPack.h"
#include "Function\ineuronprocessobject.h"
#include "../../../ngtypes/basetypes.h"
#include "../../../ngtypes/tree.h"

class CrossFiberFilter;
NG_SMART_POINTER_TYPEDEF(CrossFiberFilter, NGCrossFiberFilter);
class CrossFiberFilter :
    public INeuronProcessObject
{
public:
    static NGCrossFiberFilter New(){return NGCrossFiberFilter(new CrossFiberFilter());}
    CrossFiberFilter();
    ~CrossFiberFilter();

    INEURONPROCESSOBJECT_DEFINE

    void SetParam(NGParamPack val) { paramPack = val; }
    void SetTraceLine(std::vector<Line5d> &arg){ traceLine_ = &arg; }
    void SetSoma(const Vec3d& arg){ soma_ = arg; }

protected:
    struct ConnectInfo{
        int lineID1;
        int lineID2;
        int connectPtID;//lineID1 connect to lineID2
        Vec3d connectPt;
    };

    void ReviseBifurcationCurves();
    void BuildConInfo();

    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir);
    void ConnectNeighborCurveToTail(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t ij, Mat2i &currentConInfo);
    void ConnectNeighborCurveToHead(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t ij, Mat2i &currentConInfo);
    void AddCollideConnectNodeSub(VectorMat2i &dendConInfo, Mat2i& currentConInfo, size_t ij, const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList);
    void AddCollideConnectNode(const std::vector<VectorVec5d> &dendCurves, VectorMat2i &dendConInfo, const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList);
    //double WeighRayValue(const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg);
    //void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, std::vector<double> &rayNodeWet);
private:
    int prefix_ = 100000;
    tree<LineID> topo_;
    VectorMat2i conInfo_;
    std::vector<Line5d> *traceLine_ = nullptr;
    //VectorMat2i *traceLineConInfo = nullptr;
    NGParamPack paramPack;
    NG_SMART_POINTER_DEFINE(const Volume<unsigned short>, origImgPointer);
    Vec3d soma_;
};

#endif // ! CROSSFIBERFILTER_H
