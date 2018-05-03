#ifndef SVMTRACEFILTER_H
#define SVMTRACEFILTER_H
#include <set>
#include "../tracefilter.h"
#include "../ngtypes/tree.h"
#include "../../ngtypes/basetypes.h"
#include "../../ngtypes/volume.h"
#include "../../ngtypes/ParamPack.h"
#include "LIBSVMClassifier.h"
class SVMTraceFilter;
NG_SMART_POINTER_TYPEDEF(SVMTraceFilter, NGSVMTraceFilter);
class CrossBranchFilter;
NG_SMART_POINTER_TYPEDEF(CrossBranchFilter, NGCrossBranchFilter);

class SVMTraceFilter :
    public TraceFilter
{
public:
    
    static NGSVMTraceFilter New(Volume<int>& p, Volume<unsigned char>& t){return NGSVMTraceFilter(new SVMTraceFilter(p, t));}
    SVMTraceFilter(Volume<int>&, Volume<unsigned char>&);
    ~SVMTraceFilter(void);
    virtual ProcStatPointer Update();
    void SetInputRawCurves(const std::vector<VectorVec5d>& p);
    void SetInputRawConInfo(const VectorMat2i&);
    void SetSoma(const Vec3d& arg){soma_ = arg;}
    void SetInitialDirection(const Vec3d& arg){initDir_ = arg;}
    const std::vector<VectorVec5d>& GetRawCurves()const{return rawDendList_;}
    const VectorMat2i& GetRawConInfo()const{return rawDendConInfo_;}
    void SwapRawCurves(std::vector<VectorVec5d>& arg){arg.swap(rawDendList_);}
    void SwapRawConInfo(VectorMat2i& arg){arg.swap(rawDendConInfo_);}
    void SetTraceInitPoint(bool arg){isTraceInitialPoint_ = arg;}
    void Train(std::vector<VectorVec5d>& arg);
    void SetParam(NGParamPack val) { paramPack = val; }
    void SetCrossFiberFilter(NGCrossBranchFilter arg){ cff = arg; }
    void SetSeed(const VectorVec3d &arg){ seed = arg; }
    ProcStatPointer GPSTreeUpdate();

protected:
    struct CLUSTRELABEL{
        int index;
        Vec3d center;
        Vec3d standardDeviation;
    };

    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, std::vector<double> &rayNodeWet);
    void Principald(const VectorVec3d &dataL, Vec3d &x1);
    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, VecXd &rayNodeWet);
    VecXd WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg);
    double WeighRayValue(const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg);
    bool GetSVMClassifier( const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, bool = false);//, VecXd& w, double &b 
    //function [foreSamples,backSamples]=GetSVMSample(rawDendList,origImg)
    void GetSVMSample(const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, MatXd &foreSamples, MatXd& backSamples);
    void GetSVMSampleBySoma(const std::vector<VectorVec5d>& rawDendList, const SVolume &origImg, MatXd &foreSamples, MatXd& backSamples);
    //function featureVector=GetSVMSampleFeatures(origImg,samplePtSet)
    void GetSVMSampleFeatures(const SVolume &origImg, VectorVec3d &samplePtSet, MatXd &featureVector);
    //function [allPtNum,validSignalSum]=GetSVMSampleFeatureByRadiusThrev(origImg,seedPt,extractThrev,radius)
    void GetSVMSampleFeatureByRadiusThrev(const SVolume &origImg, Vec3i seedPt, double extractThrev, int radius, int& allPtNum, int& validSignalSum);
    //function  [resultGrowPtSet,signalMatrix]=SignalPointsIdentifyAreaGrowing(growPtSet,origImg,signalMatrix,thre)
    void SignalPointsIdentifyAreaGrowing(VectorVec3i& growPtSet, VectorVec3i& incPtSet, const SVolume& subOrigImg, CVolume& signalMatrix, double threv, int = 2, int = 2);
    //function [w,b]=SVMClassfier(foreSamples,backSamples,regularParameter)
    void SVMClassfier(const MatXd& foreSamples, const MatXd& backSamples, double regularParam, VecXd& w, double &b);
    //
    void MeanShiftGrowTraceTotal( Vec3d initPt, Vec3d initDir, const SVolume& origImg, const SVolume& backImg, 
        int curEndNodeIndex, VectorVec5d& resultCurve, int& collideIndex );//const VecXd& w, double b, 
    //function [nextNode,nextDir]=MeanShiftGrowTrace(initPt,initDir,origImg,backImg)
    void MeanShiftGrowTrace(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, Vec3d& nextNode, Vec3d& nextDir);
    //function [samplePlaneSet,extractPtSet,centerPoints]=ExtractSampleConePtSet(initPt,initDir,origImg,backImg)
    void ExtractSampleConePtSet(const Vec3d& initPt, const Vec3d& initDir, const SVolume& origImg, const SVolume& backImg, VectorMatXd& samplePlaneSet,
        std::vector<VectorVec4d> &extractPtSet, VectorVec3d& centerPoints);
    //function newPositions=CalcPositionByMeanshift(extractPtSet,centerPoints)
    void CalcPositionByMeanshift(const std::vector<VectorVec4d> &extractPtSet, VectorVec3d centerPoints, VectorVec3d& newPositions);
    //function resultPos=CalcOnePosByMeanShift(curExtractPtSet,centerPoint)
    void CalcOnePosByMeanShift(const VectorVec4d& curExtractPtSet, const Vec3d& centerPoint, Vec3d& resultPos);
    //function [traceFlag,collideIndex]=DetectCollisionByMeanShift(curPt,indexImg,origImg,curEndNodeIndex)
    void DetectCollisionByMeanShift(const Vec3d &curPt, const IVolume& indexImg, const SVolume& origImg, int curEndNodeIndex, 
        bool& traceFlag, int& collideIndex);
    //
    void CalcParmOfCurveNodeSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const Vec3d &curNode, 
        double &radius, double &wet);
    void CalcParmOfCurveNodeListSVM(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const VectorVec3d &curNode,
        std::vector<double> &radius, std::vector<double> &rav);
    //
    /* void TileTreeRemoveExtraCurveEnd( const SVolume &origImg, const SVolume &backImg, const VectorMat2i& rawDendConInfo,
         std::vector<VectorVec5d>& rawDendList, VectorVec5d& allEndNodes, VectorVec3d& allEndNodesDirection );*/
    void TileTreeRemoveExtraCurveEnd(const SVolume &origImg, const SVolume &backImg, const VectorMat2i& rawDendConInfo,
        std::vector<VectorVec5d>& rawDendList);
    //
    size_t RemoveExtraCurveEnd( const VectorVec5d& dendList, int headTailFlag, const SVolume &origImg, const SVolume &backImg );
 
    //forbidden
    IDataPointer ReleaseData(){return IDataPointer(new TreeCurve());}
    ConstIDataPointer GetOutput(){return m_Source;}
    virtual void ExtractCubeBoundary( CVolume& traceLabelMatrix, const VectorVec3i& cubeLabelPts, std::vector<VectorVec3i>& growPointSet);
    virtual void ClusterGrowPtSet( const std::vector<VectorVec3i>& growPointSet, CVolume& growLabelMatrix, std::vector<VectorVec3i>& clustrePtSet, std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> >& clustreLabel );
    virtual void RegionInflationModifyV2( const VectorVec3i &curPoints, CVolume &growLabelMatrix, double threv, VectorVec3i &growPoints );
    virtual void CurveDirection(const VectorVec3d &curDendrite, Vec3d &curDir);
    virtual void ConnectGrowClustreFromSoma(const std::vector<std::vector<int> > &subRegionCon, const std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> > &clustreLabel, std::vector<VectorVec4d> &connectPoints, std::vector<VectorVec4d> &shortConnectPoints);
    virtual void RegionConnections(const std::vector<VectorVec3i> &clustrePtSet, const std::vector<std::vector<SVMTraceFilter::CLUSTRELABEL> > &clustreLabel, std::vector<std::vector<int> > &clustreConInfo);
    virtual void RegionConnectionsSub(const VectorVec3i &curConPtSet, Volume<int> &growLabelMatrix, int threvValue, int &connectionFlag);
    virtual void ClustersExtraction(CVolume& growLabelMatrix, const VectorVec3i& inflatPtSet, VectorVec3i& subClustrePtSet, std::vector<int>& subClustrePtSetNum);
    virtual void ExtractSubRegionOP(const Vec3i &currentPoint, VectorVec3i &extractedPtSet, CVolume &growLabelMatrix);
    virtual void ExtractSubRegionAboveThrevModify(const VectorVec3i &currentPtSet, VectorVec3i &extractedPtSet, CVolume &indexImg);
    virtual void ConnectPathInTree(const std::vector<std::vector<int> > &subRegionCon, MatXi &pathMatrix);
    virtual void CubeLabel( const VectorVec3d &curDendrite, VectorVec3i& cubeLabelPts);
    virtual double CalcBifurcationThrev( const std::vector<VectorVec5d> rawDendList );
    virtual void MatXdPrctile(const MatXd& mat, double prc, int axis, VecXd& vec);
    virtual double Prctile(VecXd vec, double prc);
    virtual void DetectBifurcation(const VectorVec3i& cubeLabelPts, VectorVec3d& bifurcationDirections, std::vector<VectorVec3d>& bifurcationClustre );
    virtual void DeleteRedundancyPathV2(const MatXi &pathMatrix, std::vector<int> &redundancyLabel);
    virtual void SampleCylinderDataSetsExtra( const VectorVec3d curves, int radius, MatXd &matrixS );
    virtual void ClustersExtractionPoints(const VectorVec3i &inflatPtSet, CVolume &growLabelMatrix, VectorVec3i &subClustrePtSet, std::vector<int> &subClustrePtSetNum);
    virtual void CalcParmOfCurveNode(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const Vec3d &curNode, double &radius, double &wet);
    virtual void CalcParmOfCurveNodeList(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const VectorVec3d &curNode, std::vector<double> &radius, std::vector<double> &rav);
    bool IsNearBoundary(const Vec3d& arg, int x, int y, int z);
    size_t RemoveExtraCurveEndBySVM(const VectorVec5d& dendList, int headTailFlag, const SVolume &origImg);
    void TileTreeRemoveExtraCurveEndBySVM(const SVolume &origImg, const VectorMat2i& rawDendConInfo, std::vector<VectorVec5d>& rawDendList, VectorVec5d& allEndNodes, VectorVec3d& allEndNodesDirection);
    Vec5d FindNearestNode(const Vec3d& pt, const VectorVec5d& curve);


    //void RemoveCrossBifurcation(const std::vector<VectorVec3d> &tmpDendrites, const std::vector<VectorVec3d>&bifurClustre, std::set<size_t> &removeCrossPairIndList);
    //void FindNearestParentCurve(const Vec3d& head, const std::vector<VectorVec3d> &pCurves, size_t &curInd, size_t &ptInd);
    //bool IsIntersect(const VectorVec3i &arg1, const VectorVec3i &arg2);
    //new data
private:
    bool isTraceInitialPoint_;
    int maxSteps_ = 150;
    std::vector<VectorVec5d> rawDendList_;
    VectorMat2i rawDendConInfo_;
    Volume<int> &indexImg_;
    Volume<unsigned char>& traceLabelMatrix_;
    int globalID_; //for convenient to use
    double diffuseValue_;//it is a new constant value for trace
    Vec3d soma_;
    Vec3d initDir_;
    //restore old sample
    MatXd foreFeatureVector_old_;
    MatXd backFeatureVector_old_;
    NGParamPack paramPack;
    //boost getsvmsample
    SVolume subOrigImg_;
    CVolume signalMatrix_;
    //libsvm
    LIBSVMClassifier libsvmclassifier;
    NGCrossBranchFilter cff;
    //
    VectorVec3d seed;
};

#endif


