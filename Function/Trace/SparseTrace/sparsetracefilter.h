/* SparseTracer
*	function: trace single neuron
*	Author : zhouhang
*	2015-10-28
*/
#ifndef SPARSETRACEFILTER_H
#define SPARSETRACEFILTER_H
#include "../../ineuronprocessobject.h"
#include "../../../ngtypes/basetypes.h"
//template<typename T> class Volume;
#include "../../../ngtypes/volume.h"
#include "ngtypes/ParamPack.h"
class Tree;
class SparseTraceFilter;
NG_SMART_POINTER_TYPEDEF(SparseTraceFilter, NGSparseTraceFilter);

class SparseTraceFilter : public INeuronProcessObject
{
public:
    static NGSparseTraceFilter New(){return NGSparseTraceFilter(new SparseTraceFilter());}
    SparseTraceFilter();
    ~SparseTraceFilter();

    INEURONPROCESSOBJECT_DEFINE

    void SetInputBack(ConstIDataPointer);
    void SetInputBin(ConstIDataPointer);
    void SetSoma(ConstIDataPointer);
    void SetParam(NGParamPack arg){paramPack = arg;}
    //void SetDiffuseValue(double arg){diffuseValue_ = arg;}
    //void SetTraceValue(double arg){traceValue_ = arg+5.0;}
    //void SetBinThreshold(double arg){binThrev_ = arg;}
    void SetInitialDirection(const Vec3d&);

    //const std::vector<std::vector<size_t> >& GetTreePIDLevel()const{return curTreePIDLevel_;}
    //const std::vector<std::vector<size_t> > & GetTreeLevel()const {return treeLevel_;}

    //const VectorVec3d& GetBoundaryPtSet()const{return boundaryPtSet_;}//return boundary points set
    //const VectorVec3d& GetBoundaryDirectionSet()const{return boundaryDirectionSet_;}//
    //const std::vector<VectorVec3d>& GetBoundaryCurves()const{return boundaryCurves_;}

protected:
    //structure for clustre label
    struct CLUSTRELABEL{
        int index;
        Vec3d center;
        Vec3d standardDeviation;
    };

    void SignalIdentityMeanshift(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode,
        const double &threv, const Vec3d &initDir, bool flag, int &isNextNodeValid);

    void CubeLabelL1Opt(const VectorVec3d &curDendrite);

    void BuildCurTreeLevel(const std::vector<std::vector<VectorVec3d> >&);
    //trace the curve from soma
    void TraceDendritesConnectedSoma(const Vec3d& soma, //const Volume<unsigned short>& origImg,//2015-6-8
                                     //const Volume<unsigned short>& backImg, Volume<unsigned short> indexImg,
                                     std::vector<VectorVec3d>& initialDendrites);
                                     //Volume<unsigned short>& traceLabelMatrix);
    //continue tracing the curve.
    virtual void RecursiveTrace(std::vector<std::vector<VectorVec3d> > &treeDataSet);

    void ReconstructSomaShapeQuanReviV2(const Vec3d& initialPoint,
                                      //2015-6-8
                                      VectorVec3i& innerSomaPts);

    void RegionInflationModify(const VectorVec3i& curPoints, //2015-6-8
                               CVolume& traceLabelMatrix,
                               double threv,
                               VectorVec3i& growPoints);

    //2015-6-8
    virtual void RegionInflationModifyV2(const VectorVec3i &curPoints,
                                                  CVolume &growLabelMatrix, double threv,
                                                  VectorVec3i &growPoints);

    void ClusterGrowPtSet(const std::vector<VectorVec3i>& growPointSet, CVolume& growLabelMatrix,
                          std::vector<VectorVec3i>& clustrePtSet,
                          std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> >& clustreLabel);

    void RegionConnections(const std::vector<VectorVec3i>& clustrePtSet,
                           const std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> >& clustreLabel, std::vector<std::vector<int> >& clustreConInfo);

    void ConnectGrowClustreFromSoma(const std::vector<std::vector<int> >& subRegionCon,
                                    const std::vector<std::vector<SparseTraceFilter::CLUSTRELABEL> >& clustreLabel,
                                    std::vector<VectorVec4d>& connectPoints,
                                    std::vector<VectorVec4d>& shortConnectPoints);

    void AddSomaToInitDendrites(const std::vector<VectorVec4d>& connectPoints,
                                const std::vector<VectorVec4d>& shortConnectPoints,
                                const Vec3d& soma, std::vector<VectorVec3d>& initDendrites);

    void IdentifyingRayburst(const Volume<unsigned short>& origImg,//2015-6-8
                             std::vector<VectorVec3d>& initialDendrites);

    void CurveDirection(const VectorVec3d& curDendrite, Vec3d& curDir);

    void ContrainedPCATraceSomaBifur(//const Volume<unsigned short>& origImg,//2015-6-8
                                 //const Volume<unsigned short>& backImg,//2015-6-8
                                 //Volume<unsigned short>& traceLabelMatrix,
                                 //Volume<int>& indexImg,
                                 //int &globalID,
                                 const std::vector<VectorVec3d>& initialDendritesCp,
                                 const VectorVec3d& directionSet,
                                 std::vector<VectorVec3d> &resultDendrites);

    void ContrainedPCATraceBifur(const std::vector<VectorVec3d> &initialDendritesCp,
                                 const VectorVec3d &directionSet,
                                 std::vector<std::vector<VectorVec3d> > &treeDataSet);

    void CalVectorAngle(const Vec3d& curNormPt, double& altitudeAngle, double& azimuthAngle);

    void ClustersExtractionPoints(const VectorVec3i &inflatPtSet, CVolume& growLabelMatrix,
                                  VectorVec3i& subClustrePtSet, std::vector<int>& subClustrePtSetNum);

    void RegionConnectionsSub(const VectorVec3i &curConPtSet, Volume<int>& growLabelMatrix, int threvValue,
                              int &connectionFlag);

    void RayburstShapeTrackV2(const Vec3d& end2ndPt, const Vec3d& endPt, const Volume<unsigned short>& origImg,//2015-6-8
                              int curveLen, VectorVec3d& forwardArea);

    void RayBurstShapeV2(const Vec3d& initPt, const Volume<unsigned short>& origImg,//2015-6-8
                         VectorVec3i& rayNode);


    void TraceCurveForwardFromSomaV2(//const Volume<unsigned short>& origImg,//2015-6-8
                                     //const Volume<unsigned short>& backImg,
                                     const Vec3d& seedPoint, const Vec3d& curDirection,
                                     //Volume<int>& traceLabelMatrix,
                                     VectorVec5d &traceDendrite);

    void TraceCurveForwardFromSomaV3(const Vec3d &seedPoint, const Vec3d &curDirection,
                                                        VectorVec5d &resultSeedCurve, int &collideID, int maxStep = 3000);//2015-6-8

    //2015-6-8
    virtual void CubeLabel(//SVolume &traceLabelMatrix,
                   //const Volume<unsigned short>& origImg, //2015-6-8
                   //const Volume<unsigned short>& backImg,
                   const VectorVec3d &curDendrite);//, int globalID, Volume<int> &indexImg

    //2016-13
    void CubeLabelSoma(const VectorVec3d &curDendrite);//, int globalID, Volume<int> &indexImg

    virtual void DetectBifurcation(//const Volume<unsigned short>& origImg,//2015-6-8
                           //const Volume<unsigned short>& backImg,
                           CVolume& traceLabelMatrix,
                           const std::vector<VectorVec3d>& curveCluster,
                           VectorVec3d& bifurcationDirections,
                           std::vector<VectorVec3d>& bifurcationClustre,
                           VectorVec3i& boundaryLabel);

    virtual void ExtractCubeBoundary(const std::vector<VectorVec3d>& curveCluster,
                             CVolume& traceLabelMatrix,
                             //const Volume<unsigned short>& origImg,//2015-6-8
                             //const Volume<unsigned short>& backImg,
                             std::vector<VectorVec3i>& growPointSet,VectorVec3i& boundaryLabel);

    void ThreeDimdCurveResample(const VectorVec3d& curvePoints, double incValue,
                                VectorVec3d& resampleCurve);

    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,//2015-6-8
                       std::vector<double> &rayNodeWet);

    void CalcParmOfCurveNodeList(const Volume<unsigned short> &origImg,//2015-6-8
                                 const Volume<unsigned short> &backImg,
                                 const VectorVec3d &curNode,
                                 std::vector<double> &radius, std::vector<double> &rav);

    void CalcParmOfCurveNode(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,//2015-6-8
                             const Vec3d &curNode, double &radius, double &wet);

    virtual void TraceNextCurveNodeV2(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode, //2015-6-8
                            const double &threv, const Vec3d &initDir, bool flag,
                            Vec3d &nextCurveNode, Vec3d &nextDenDir, int &isNextNodeValid);


    void CalcNeighborSignalV2(const Volume<unsigned short> &origImg, const Vec3d &curSeedNode, //2015-6-8
                            const Vec3d &initVec, const double threv, bool flag, VectorVec3d &neighborPtSet,
                            std::vector<double> &neighborWet, Vec3d &firDir, Vec3d &secDir);


    void CalcConstraintPCA(const VectorVec3d &neighborPtSet, const Vec3d &curSeedNode,
                           const Mat3d &convMat, const std::vector<double> &neighborWet,
                           double &P, Mat3d &sigmaH, Vec3d &mdata3);


    void CalcPCADirections(const Mat3d &sigmaH, const Vec3d &initVec, const Vec3d &T2, const double threv, Vec3d &vec1);


    void CalcPCAMainDirection(const Vec3d &x0, const Mat3d &sigmaH, const Vec3d &T1, const Vec3d &T2, const double threv, Vec3d &x1);

    void IsCurvesCollideV2(const Vec5d &seed, int &isNextNodeValid);//2015-6-8 const Volume<int> &traceLabelMatrix,

    void IsCurvesCollideV3(const Vec5d &seed, int& isNextNodeValid, int &collideID);//2015-6-8

    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir);

    void CalcOrthoBasis(const Vec3d &vec1, Vec3d &vec2, Vec3d &vec3);

    void ClustersExtraction(CVolume &growLabelMatrix, const VectorVec3i &inflatPtSet,
                            VectorVec3i &subClustrePtSet, std::vector<int> &subClustrePtSetNum);

    void ExtractSubRegionOP(const Vec3i &currentPoint, VectorVec3i &extractedPtSet, CVolume &growLabelMatrix);

    void ExtractSubRegionAboveThrevModify(const VectorVec3i &currentPtSet, VectorVec3i &extractedPtSet,
                                          CVolume &indexImg);

    void GetRayLimitV2(const Volume<double> &sphereRayWet, const Volume<double> &sphereBackWet,
                       const double constrictionThrev, std::vector<std::vector<double> > &rayLimit);//2015-6-8

    void ConnectPathInTree(const std::vector<std::vector<int> > &subRegionCon, MatXi& pathMatrix);//2015-6-8

    void DeleteRedundancyPath(const MatXi& pathMatrix, std::vector<int>& redundancyLabel);

    void DeleteRedundancyPathV2(const MatXi &pathMatrix, std::vector<int> &redundancyLabel);//2015-6-8

    void Ind2CellInd(std::vector<std::vector<VectorVec3d> >&, int ind, int& cellLevel, int& cellInd);
    void ExtractLocalDomainV2(const Vec3d &initPoint, const SVolume &origImg, SVolume &locOrigImg, Vec3d &locPoint);


    //Mean shift algorithm
    bool IsNearBoundary(const Vec3d& arg, int x, int y, int z);
    void MeanShiftGrowTraceTotal(Vec3d initPt, Vec3d initDir, const SVolume& origImg, const SVolume& backImg, VectorVec5d& resultCurve, int& collideIndex);
    void MeanShiftGrowTrace(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, Vec3d& nextNode, Vec3d& nextDir);

    void ExtractSampleConePtSet(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, VectorMatXd& samplePlaneSet, std::vector<VectorVec4d> &extractPtSet, VectorVec3d& centerPoints);
    void CalcPositionByMeanshift(const std::vector<VectorVec4d> &extractPtSet, VectorVec3d centerPoints,/*centerPoints will be modified */ VectorVec3d& newPositions);
    void CalcOnePosByMeanShift(const VectorVec4d& curExtractPtSet, const Vec3d& centerPoint, Vec3d& resultPos);
    void DetectCollisionByMeanShift(const Vec3d &curPt, const IVolume& indexImg, const SVolume& origImg, bool& traceFlag, int& collideIndex);
    void CalcParmOfCurveNodeListSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const VectorVec3d &curNode, std::vector<double> &radius, std::vector<double> &rav);
    void CalcParmOfCurveNodeSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const Vec3d &curNode, double &radius, double &wet);
    

    //m_Input is origImg
    //m_Output is TreeCurve;
    ConstIDataPointer m_Bin;
    ConstIDataPointer m_Back;
    ConstIDataPointer m_Soma;
    NGParamPack paramPack;
    //output boundary points set
    //int boundaryDistanceThreshold_;//the least distance from boundary
    //VectorVec3d boundaryPtSet_;
    //VectorVec3d boundaryDirectionSet_;
    //std::vector<VectorVec3d> boundaryCurves_;
    Vec3d soma_;

    //2015-6-8
    int maxSteps_ = 1000;
    int globalID_; //for convenient to use
    //double diffuseValue_;//it is a new constant value for trace
    //double traceValue_;
    //double binThrev_;
    Vec3d initialDirection_;
    //tree level
    std::vector<std::vector<size_t> > treeLevel_;
    //std::vector<std::vector<size_t> > curTreePIDLevel_;
    //std::vector<std::vector< std::vector<size_t> > > curTreeSIDLevel_;
    
#ifdef _WIN32
    std::tr1::shared_ptr<const Volume<unsigned short> > origImgPointer;
    std::tr1::shared_ptr<const Volume<unsigned short> > backImgPointer;
    std::tr1::shared_ptr<const Volume<NGCHAR> > binImgPointer;
#else
    std::shared_ptr<const Volume<unsigned short> > origImgPointer;
    std::shared_ptr<const Volume<unsigned short> > backImgPointer;
    std::shared_ptr<const Volume<NGCHAR> > binImgPointer;
#endif
    Volume<unsigned char> traceLabelMatrix_;// it is needed to destroy
    Volume<int> indexImg_;//

};

#endif // SPARSETRACEFILTER_H
