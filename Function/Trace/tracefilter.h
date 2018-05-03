/*
 * Copyright (c)2013-2015  Zhou Hang, Shaoqun Zeng, Tingwei Quan
 * Britton Chance Center for Biomedical Photonics, Huazhong University of Science and Technology
 * All rights reserved.
 */
#ifndef TRACEFILTER_H
#define TRACEFILTER_H
#include "../ineuronprocessobject.h"
#include "../../ngtypes/basetypes.h"
#include "../../ngtypes/volume.h"
#include "../../ngtypes/ParamPack.h"
//template<typename T> class Volume;
class Tree;
class TraceFilter;
#ifdef _WIN32
typedef std::tr1::shared_ptr<TraceFilter> NGTraceFilter;
#else
typedef std::shared_ptr<TraceFilter> NGTraceFilter;
#endif
class TraceFilter : public INeuronProcessObject
{
public:

    static NGTraceFilter New(){return NGTraceFilter(new TraceFilter());}
    TraceFilter();
    ~TraceFilter();
    INEURONPROCESSOBJECT_DEFINE

    IDataPointer GetConnect();
    void SetInputBack(ConstIDataPointer);
    void SetInputBin(ConstIDataPointer);
    void SetSoma(ConstIDataPointer);
    IDataPointer &GetDendConInfo();
    void SetThreValue(double);//2015-8-13
    void SetParam(NGParamPack arg){ paramPack = arg; }

protected:
    //m_Input origimage
    ConstIDataPointer m_Bin;
    ConstIDataPointer m_Back;
    ConstIDataPointer m_Soma;
    IDataPointer rawDendConInfoPointer;
    //VectorMat2i rawDendConInfo;
    //m_Source tree
    double fillThrev, endThrev;
    //SVM
    NGParamPack paramPack;
    //VecXd w;double b;
    NG_SMART_POINTER_DEFINE(const Volume<unsigned short> , origImgPointer);
    NG_SMART_POINTER_DEFINE(const Volume<unsigned short> , backImgPointer);
    NG_SMART_POINTER_DEFINE(const Volume<unsigned char> , binImgPointer);
//#ifdef _WIN32
//	std::tr1::shared_ptr<const Volume<unsigned short> > origImgPointer;
//	std::tr1::shared_ptr<const Volume<unsigned short> > backImgPointer;
//	std::tr1::shared_ptr<const Volume<unsigned char> > binImgPointer;
//#else
//	std::shared_ptr<const SVolume > origImgPointer;
//	std::shared_ptr<const SVolume > backImgPointer;
//	std::shared_ptr<const CVolume > binImgPointer;
//#endif

    void TraceCurvesAndConInfo(Volume<int> &resultIndexImg, const Volume<unsigned short> &origImg, const Volume<NGCHAR> &binImg, const Volume<unsigned short> &backImg,
                 const VectorVec3d &somaList,
       std::vector<VectorVec5d> &rawDendList, //std::vector<mat2d > &dd2,
       VectorMat2i& rawDendConInfo);

    void SelectSeedForTrace(const Volume<NGCHAR> &binImg, const Volume<unsigned short> &origImg, const Volume<int> &indexImg, VectorVec3d &traceSeed);

    void CalcNeighborSignal(
            const Volume<unsigned short> &origImg,//XXv
            const Vec3d& curSeedNode, //data1
            const Vec3d& initVec, //x1
            const double threv,
            VectorVec3d& neighborPtSet,	//dataS0
            std::vector<double>& neighborWet,	//W
            Vec3d& firDir,//x10
            Vec3d& secDir // x11
            );
    //function aa_data=weigthvalue(dataL1,L_XX3)
    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, std::vector<double> &rayNodeWet);
    //void WeighRayValue(const VectorVec3d &dataL1, const Volume<unsigned char> &L_XX3, std::vector<double> &aa_data);

    void TraceCurvesFromSeed(const Volume<unsigned short> &orig, const Volume<unsigned short> &back,
                                const Vec3d &initialP, Volume<int> &YY,
                                const VectorVec3d &dataC,
                                const int kkms, VectorVec5d &dataLLss, Mat2i &MM);

    void CalcInitDirectionOnTraceSeed(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
                                      const Vec3d &initPoint, const double windowSize,
                                      Vec3d &vec1);

    void TraceCurvesForwardFromSeed(const Vec3d &seedInitDir,
                                  const Vec3d &initPt, const Volume<int> &indexImg, const VectorVec3d &somaList,
                                  const Volume<unsigned short> &backImg, const Volume<unsigned short> &origImg,
                                  VectorVec5d &resultSeedCurve, Vec2i &curConInfo);

    void IsCurvesCollide(const Volume<int> &indexImg, const int xMin, const int xMax, const int yMin, const int yMax,
                     const int zMin, const int zMax, int &conInfo, double &maxCollideSetWet, double &somaCollideAngle);

    void TraceCurvesFromSoma(const Volume<unsigned short> &XX3, const Volume<unsigned short> &XXb, Volume<int> &YY,
        const VectorVec5d &initialP, //const Volume<unsigned short> &XXb, const Volume<unsigned short> &XXv,
        const int kkms, const Mat3d &somaInitDir, const VectorVec3d &somaList, /*const VectorVec5d &existDendCurves, */
                             VectorVec5d &dataLLss, Mat2i &MM);

    void TraceCurveForwardFromSoma(const Vec3d &dendSomaInitDir, /*const double thre,*/
                                  const Vec3d &initPt, const Volume<int> &indexImg, /*const Volume<int> &YY,*/ const Volume<unsigned short> &backImg,
                                   const Volume<unsigned short> &origImg, const VectorVec3d &somaList,
                                  VectorVec5d &resultSomaCurve);

    void ReconstructShapeForTrace(
        const Vec3d &intialPoint, const Volume<unsigned short> &origImg, //const Volume<unsigned short> &backImg,
        //Vec3d &denSomaInitDirection, /*Vec3d &xx2, Vec3d &xx3,*/
            //Vec3d &denInitCenter,
            MatXd& resultRayLength,
            VectorVec3d &innerSomaPts);

    void InflateBoundary(const VectorVec3d &innerSomaPts, const Volume<unsigned short> &locOrigImg, const double inflateThrev,
                          VectorVec3d &inflatedArea);

    void InflateMarginalAreaAboveThrev(const VectorVec3d &initPoints, Volume<unsigned short> &locIndexImg, const double thre,
                         const int nx, const int ny, const int nz, VectorVec3d &resultInflatedArea);

//    void DetectCellSignal(const Volume<NGCHAR> &locIndexImg, const VectorVec3d &inflatedPoints,
//                                         VectorVec3d &dendInitPts, std::vector<int> &dendInitPtSetSum);

    void DetectCellSignalModify(const Volume<NGCHAR> &locIndexImg, const VectorVec3d &inflatedPoints, const int threv,
                                         VectorVec3d &dendInitPts, std::vector<int> &dendInitPtSetSum);

    void CalcOrthoBasis(const Vec3d &vec1, Vec3d &vec2, Vec3d &vec3);

    void CalcConstraintPCA(const VectorVec3d &neighborPtSet, //data1
                           const Vec3d &curSeedNode,	//data3
                           const Mat3d &convMat, //data2
                           const std::vector<double> &neighborWet, double &P,
                           Mat3d &sigmaH, //sigmaH
                           Vec3d &mdata3//mdata3
                           );

    void CalcPCADirections(const Mat3d &sigmaH,//sigmaH
                   const Vec3d &initVec,//T1
                   const Vec3d &T2,//T2
                   const double threv,
                   Vec3d &vec1
                   );

    void CalcPCAMainDirection(const Vec3d &x0,
                              const Mat3d &sigmaH,//sigmaH
                              const Vec3d &T1, //T1
                              const Vec3d &T2, //T2
                              const double threv,
                              Vec3d &x1//x1
                              );

    void CalcParmOfCurveNode(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const Vec3d &curNode,
                             double &radius, double &wet);

    void CalcParmOfCurveNodeList(const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg,
                                 const VectorVec3d &curNode,
                              std::vector<double> &radius, std::vector<double> &rav);

    //function [x1,x2,x3]=principald(dataL)
    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir /*, Vec3d &x2, Vec3d &x3*/);

    void TraceNextCurveNode(
            const Volume<unsigned short> &locOrigImg,//XXv
            //const Volume<unsigned short> &XXb,//XXb
            const Vec3d &curSeedNode,//dataL1
            const double &threv,
            const Vec3d &initDir,
            /*const Vec3d &T2,
            const Vec3d &T3,*/
            Vec3d &nextCurveNode,
            //Vec3d &vmmk,
            Vec3d &nextDenDir,//vmmk
            /*Vec3d &x2,
            Vec3d &x3,*/
            int &isNextNodeValid//idexx
            );

    void ClearShortCurvesAndInvalidConnection(const std::vector<VectorVec5d > &denPt, const VectorMat2i &denCx,
                      std::vector<VectorVec5d > &denPtNew, VectorMat2i &denCxNew/*, std::vector<int> &denNums*/);

    /*void AddCollideConnectNode(const std::vector<VectorVec5d > &denPt, const VectorMat2i &denCx, const VectorVec3d &myCell,
                      std::vector<VectorVec5d > &f_pts);*/

    void ReconstructSomaShapeForTrace(const VectorVec3d &somaList, const Volume<unsigned short>& origImg,
                                      Volume<int>& indexImg, VectorMatXd& allRayLen);

    void FindThickDendDirFromSoma(const MatXd& rayLength, const Volume<unsigned short>& origImg,
                                  const Volume<unsigned short>& backImg,
                                  const Vec3d &currentSoma, const Volume<int>& indexImg, const int seriID,
                                  Vec3d& mainDir1, Vec3d& mainDir2, Vec3d& mainDir3, Vec3d &thickDendInitPt);

    void ExtractLocalDomainAlias(const Vec3d& currentSoma, const Volume<unsigned short>& origImg,
                                 const Volume<unsigned short>& backImg,
                                 const Volume<int>& indexImg, const int seriID,
                                 Volume<unsigned short> &locOrigImg, Vec3d& locPt);

    void ForceModifyDirection(VectorMat3d& allThickDendInitDir, std::vector<int>& largeAngleFlagList);

    void DetectCollision(const VectorVec5d& somaCurve, const VectorVec3d& somaList, const Vec2i& collideInfo,
                         const Volume<int>& indexImg, const Volume<unsigned short>& origImg,
                         const Volume<unsigned short>& backImg, const int origID, const int isSoma, const int isContinue,
                         VectorVec5d& modifiedSomaCurve, int &isContinueTrace, Vec2i& resultCollideInfo,
                         int &curID, int &resultConInfo);

    void TraceForwardInKnot(const Vec3d &initMainDir,const Vec3d& initP,const Volume<unsigned short>& backImg,
                            const Volume<unsigned short>& origImg,
                            const Volume<int>& indexImg,const int collideConInfo,
                            VectorVec3d& crossCurve, int &isContinueTrace);

    void PushNodeUseRayBurst(const Vec3d& curveNode, const Volume<unsigned short>& origImg, const double threv, const Vec3d& x1,
            Vec3d& pushedNode);

    void RayBurstSampling(const Volume<double> &Sphere_XX, const double three_vs, std::vector<std::vector<double> > &Uzz);
    void RayBurstSampling1(const std::vector<double> &LLs, const double three_vs, int &arra);
    void Principald(const VectorVec3d &dataL, Vec3d &x1, Vec3d &x2, Vec3d &x3);
    bool Corrode(const Volume<unsigned char> &binImage, const VectorVec3d &binPtSet,
                 const double eroIntensity, Volume<unsigned char> &eroBinImg, VectorVec3d &eroPtSet);
    void RayBurstShape(const Vec3d &initSoma, const Volume<unsigned short> &v, VectorVec3d &rayNode, Volume<double> &smoothRay);
    void RayburstShapeTrack(const Vec3d &initPointt, const Volume<unsigned short> &origImg, const Vec3d &initDir, int len, VectorVec5d &datatt);
    void ReconstructSomaShapeQuanRevi(const Vec3d &initialPoint, MatXd& resultRayLength, VectorVec3d &innerSomaPts);
    void GetRayLimitV2(const Volume<double> &sphereRayWet, const Volume<double> &sphereBackWet, const double constrictionThrev, std::vector<std::vector<double> > &rayLimit);
    void ExtractLocalDomainV2(const Vec3d &initPoint, const SVolume &origImg, SVolume &locOrigImg, Vec3d &locPoint);
    void CalVectorAngle(const Vec3d &curNormPt, double &altitudeAngle, double &azimuthAngle);
    void CalcNeighborSignalV2(const Volume<unsigned short> &origImg, const Vec3d &curSeedNode, const Vec3d &initVec, const double threv, VectorVec3d &neighborPtSet, std::vector<double> &neighborWet, Vec3d &firDir, Vec3d &secDir);
    void CalcConstraintPCAv2(const VectorVec3d &neighborPtSet, const Vec3d &curSeedNode, const Mat3d &convMat, const std::vector<double> &neighborWet, double &P, Mat3d &sigmaH, Vec3d &mdata3);
    void TraceNextCurveNodeV2(const Volume<unsigned short> &locOrigImg, const Vec3d &curSeedNode, const double &threv, const Vec3d &initDir, Vec3d &nextCurveNode, Vec3d &nextDenDir, int &isNextNodeValid);
    void DetectCollisionV2(const VectorVec5d &somaCurve, const VectorVec3d &somaList, const Vec2i &initCollideInfo, const Volume<int> &indexImg, const Volume<unsigned short> &origImg, const Volume<unsigned short> &backImg, const int origID, const int isSoma, const int isContinue, VectorVec5d &modifiedSomaCurve, int &isContinueTrace, Vec2i &resultCollideInfo, int &curID, int& resultConInfo);
    void TraceForwardInKnotV2(const Vec3d &initMainDir, const Vec3d &initP, const Volume<unsigned short> &backImg, const Volume<unsigned short> &origImg, const Volume<int> &indexImg, const int collideConInfo, VectorVec3d &crossCurve, int &isContinueTrace);
    void IsCurvesCollideV2(const Volume<int> &indexImg, const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, int &conInfo, double &curveCollideFlag, double &somaCollideFlag);
    //function [ddk1,ddk2]=MyCrossDetection(dendCurves,dendConInfo)
    void MyCrossDetection(const std::vector<VectorVec5d>& dendCurves, const VectorMat2i &dendConInfo, 
        std::vector<VectorVec5d>& resultDendCurves, VectorMat2i &resultDendConInfo);
    //function clustreConnectDotSet = ClustreFuseConnectDotSet(connectDotSet, distThrev)
    void ClustreFuseConnectDotSet(const VectorVec5d& connectDotSet, double distThrev, std::vector<std::vector<int> >& clustreConnectDotSet);
    //[fuseFlag,paraDist,corrDir]=MyIdentifyCurvesMerge(curCurve1,curCurve2,headTailFlags);
    void MyIdentifyCurvesMerge(const VectorVec5d&, const VectorVec5d&, const Vec2i&, int&, double&, double&);
    //function [paraDist,corrDir]=MyCalcCurvesSimilarity(curve1,curve2)
    void MyCalcCurvesSimilarity(const VectorVec3d&, const VectorVec3d&, double&, double&);
    //function [paraDist, corrDir]=MyCalcSegmentSimilarity(curve1Dot1,curve1Dot2,curve2Dot1,curve2Dot2)
    void MyCalcSegmentSimilarity(const Vec3d&, const Vec3d&, const Vec3d&, const Vec3d&, double&, double&);
    //function [resultDendCurves,resultDendConInfo]= MyMergeCurves(resultDendCurves,resultDendConInfo,fuseCurveIndex1,fuseCurveIndex2,headTailFlags)
    void MyMergeCurves(std::vector<VectorVec5d>&, VectorMat2i&, int, int, const Vec2i&);

    void ConnectNeighborCurveToHead(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t k, Mat2i &currentConInfo);
    void ConnectNeighborCurveToTail(std::vector<VectorVec5d > &resultDendList, VectorMat2i &dendConInfo, size_t k, Mat2i &currentConInfo);

    void AddCollideConnectNode(const std::vector<VectorVec5d> &dendCurves, const VectorMat2i &dendConInfo, const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList);
    void AddCollideConnectNodeSub(VectorMat2i &denCx, Mat2i& currentConInfo, size_t ij, const VectorVec3d &myCell,
        std::vector<VectorVec5d > &f_pts);
    void AddCollideConnectNodeV2(const std::vector<VectorVec5d> &dendCurves, VectorMat2i &dendConInfo, const VectorVec3d &somaList, std::vector<VectorVec5d> &resultDendList);
};

#endif // TRACEFILTER_H
