/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef SPARSETRACEAXONFILTER_H
#define SPARSETRACEAXONFILTER_H
#include <set>
#include "sparsetracefilter.h"
#include "../../../ngtypes/basetypes.h"
#include "../../../ngtypes/volume.h"
#include "SVMTraceFilter.h"
class Tree;
class SparseTraceAxonFilter;
NG_SMART_POINTER_TYPEDEF(SparseTraceAxonFilter, NGSparseTraceAxonFilter);

class SparseTraceAxonFilter : public SparseTraceFilter
{
public:
    static NGSparseTraceAxonFilter New(){return NGSparseTraceAxonFilter(new SparseTraceAxonFilter());}
    SparseTraceAxonFilter();
    ~SparseTraceAxonFilter();

    INEURONPROCESSOBJECT_DEFINE

    void Train(std::vector<VectorVec5d>&);
   
protected:
    void TraceDendritesConnectedAxon(const Vec3d&, std::vector<VectorVec3d>&);
    virtual void RecursiveTrace(std::vector<std::vector<VectorVec3d> > &treeDataSet);
    virtual void TraceNextCurveNodeV2(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode, const double &threv, const Vec3d &initDir, bool flag, Vec3d &nextCurveNode, Vec3d &nextDenDir, int &isNextNodeValid);
    //2015-6-8 here the traceLabelMatrix is the copy of global traceLabelMatrix
    virtual void DetectBifurcation(CVolume &traceLabelMatrix, const std::vector<VectorVec3d> &curveCluster,
                                          VectorVec3d &bifurcationDirections, std::vector<VectorVec3d> &bifurcationClustre,
                                          VectorVec3i& boundaryLabel);
    virtual void ExtractCubeBoundary(const std::vector<VectorVec3d>& curveCluster,
        CVolume& traceLabelMatrix,
        std::vector<VectorVec3i>& growPointSet,VectorVec3i& boundaryLabel);

    virtual void RegionInflationModifyV2(const VectorVec3i &curPoints,
                                              CVolume &growLabelMatrix, 
                                              VectorVec3i &growPoints);

    virtual void CubeLabel(const VectorVec3d &curDendrite);
    
    //
    //void RemoveCrossBifurcation(const std::vector<VectorVec3d> &, const std::vector<VectorVec3d>&, std::set<size_t> &);
    //virtual void FindNearestParentCurve(const Vec3d& head, const std::vector<VectorVec3d> &pCurves, size_t &curInd, size_t &ptInd);
    //void GetSVMSampleFeatureByRadiusThrev(const SVolume &, Vec3i , double , int , int& , int& );
    //void SignalPointsIdentifyAreaGrowing(VectorVec3i& , const SVolume& , CVolume& , double, int = 2, int = 2 );
    //bool IsIntersect(const VectorVec3i&, const VectorVec3i&);
    //virtual void GetCrossCandidateCurves(size_t ptInd1, size_t ptInd2, std::vector<VectorVec3d>&bifurClustre)
    //below is save in base class
    //m_Input is origImg
    //m_Output is TreeCurve;
    /*ConstDataPointer m_Bin;
    ConstDataPointer m_Back;
    ConstDataPointer m_Soma;*/
    //output boundary points set
    //int boundaryDistanceThreshold_;//the least distance from boundary
    //VectorVec3d boundaryPtSet_;
    //VectorVec3d boundaryDirectionSet_;
    //std::vector<VectorVec3d> boundaryCurves_;
//#ifdef _WIN32
//    std::tr1::shared_ptr<const Volume<unsigned short> > origImgPointer;
//    std::tr1::shared_ptr<const Volume<unsigned short> > backImgPointer;
//    std::tr1::shared_ptr<const Volume<NGCHAR> > binImgPointer;
//#else
//    std::shared_ptr<const Volume<unsigned short> > origImgPointer;
//    std::shared_ptr<const Volume<unsigned short> > backImgPointer;
//    std::shared_ptr<const Volume<NGCHAR> > binImgPointer;
//#endif
//    Volume<unsigned short> traceLabelMatrix_;// it is needed to destroy
//    Volume<int> indexImg_;//
    int maxSteps_;
    NGSVMTraceFilter svmFilter;
    NGCrossBranchFilter cff;
};

#endif // SPARSETRACEFILTER_H
