/*
 * Copyright (c)2013-2015  Zhou Hang, Shaoqun Zeng, Tingwei Quan
 * Britton Chance Center for Biomedical Photonics, Huazhong University of Science and Technology
 * All rights reserved.
 */
#ifndef PROBEFILTER_H
#define PROBEFILTER_H
#include <stack>
#include "ineuronprocessobject.h"
#include "../ngtypes/basetypes.h"

typedef std::stack<Vec4d,  std::deque<Vec4d,Eigen::aligned_allocator<Vec4d>> > StackVec4d;
class Soma;
struct Cell;
template<typename T> class Volume;
typedef Volume<NGCHAR> CVolume;
typedef Volume<unsigned short> SVolume;
typedef Volume<double> DVolume;
class ProbeFilter;
typedef std::shared_ptr<VectorVec3i> BinPtSetPointer;
typedef std::shared_ptr<ProbeFilter> NGProbeFilter;

class ProbeFilter : public INeuronProcessObject
{
public:
    static NGProbeFilter New(){return NGProbeFilter(new ProbeFilter());}
    ProbeFilter();
    ~ProbeFilter();
    ProcStatPointer Update();
    ConstIDataPointer GetOutput();
    IDataPointer ReleaseData();

    void SetThreadNum(int);
    void SetInputBin(ConstIDataPointer);
    void SetInputBinPtSet(BinPtSetPointer);
    void SetConnectDomainRange(int, int);
    void SetMinRadius(double);
    double GetMinRadius()const;

protected:
    void SomaShape(const Vec3d&, const SVolume&, const CVolume&, VectorVec3i&, MatXd&, std::vector<std::vector<double> >&,
        std::vector<VectorVec3i> &, double& );

    void NeighborhoodPtSet(const VectorVec3i&, std::vector<VectorVec3i>&);

    double GetOverlapRate(const VectorVec3i&, const VectorVec3i&);

    bool WriteSwcTest(const char*, const VectorVec3i&);

private:
    int threadNum;
    int k_theta ;
    int k_phi;
    int minConnectNum_ ;
    int maxConnectNum_ ;
    int maxBoundingBoxSize_;
    double slice_;
    double minRadius_ ;
    double xVoxel_;
    double yVoxel_;
    double zVoxel_;
    ConstIDataPointer m_BinImg;
    BinPtSetPointer m_BinPtSet;
    //m_Source soma

    //function
    /*locate soma cell*/
    bool IsInVolume(const Vec3d&, int, int, int, int, int, int);

    void ExtractLocalDomain(const Vec3d &initSoma, const SVolume &v, SVolume &locOrigImg, Vec3d &locSoma);//extractLocalArea

    void GetSubDomain(const Vec3i &initPoint, VectorVec3i &connectPtSet, CVolume &eroImage);

    void SubDomainExpand(const VectorVec3i &initPoints, VectorVec3i &connectPtSet, CVolume &eroImage);

    bool IsCxDomainBoundingBoxTooLarge(const VectorVec3i& locPtSet);
    void CalculateCxDomainBoundingBox(const VectorVec3i& locPtSet, int& xMin, int &xMax, int &yMin,
                                      int &yMax, int &zMin, int &zMax);


    /*erode*/
    bool Corrode(const CVolume &binImage, const VectorVec3i &binPtSet,
                 const double eroIntensity,
                 CVolume &eroBinImg, VectorVec3i &eroPtSet);

    /*extract connected domain*/
    void GetConnectedDomain(const CVolume &eroImage, const VectorVec3i &eroPtSet,
                            CVolume &remainEroImg, VectorVec3i &remainBinPtSet,
                            VectorVec3i &connectPtSet, std::vector<int> &connectNum);

    /*find soma*/
    void SomaLocalize(const SVolume &v, const VectorVec3i &connectPtSet, const std::vector<int> &connectNum,
                       Soma &resultSoma);

    /*pre-select seed point */
    void CreateSeeds(const VectorVec3i &curConPtSet, const SVolume &origImage,
                     CVolume &dstBinImg, SVolume &dstOrigImg, VectorVec4d &locSeed, VectorVec3d &globalSeed);

    //
    void RayBurstShape(const Vec3d &initSoma, const SVolume &v, MatXd &resultRayLength, DVolume &smoothRay);

    /*private function*/
    /*total function of localize*/
    void _SomaRecognition(const CVolume &locBinImg, const VectorVec3d &locPtSet, const SVolume &locOrigImg,
                          VectorVec3d &resultLocPtSet, std::vector<double> &radius);

    /*correct the position of seed*/
    void _PostionUpdate(const SVolume &locOrigImg, const VectorVec3d &locPtSet, const std::vector<double> &radius,
                        VectorVec3d &resultLocPtSet);

    /*correct the position of seed, sub function*/
    void _CalculateNewPosition(const SVolume &locOrigImg, const Vec3d &locPtSet, const double radius, Vec3d &resultLocPos);

    void _CalculateL1Minimum(const CVolume &locBinImg, const VectorVec3d &locPtSet,
                             const std::vector<double> &rWet, const std::vector<double> &radius, std::vector<double> &rGrad);

    /*cal differential*/
    void _CalculateGrad(const CVolume &locBinImg, const VectorVec3d &locPtSet,
                        const std::vector<double> &rWet, const std::vector<double> &radius, std::vector<double> &rGrad);

    void _CalculateL1Value(const CVolume &locBinImg, const VectorVec3d &locPtSet, const std::vector<double> &L1Wet,
                           const std::vector<double> &radius, double &L1Value);

    /*cal sphere function*/
    void _CalculateSphereWet(const Vec3d &seed, const int vx, const int vy, const int vz,
                             const double radius, std::vector<double> &sphereWet);

    //long short axis, surface, volume
    void CalSomaInfo(const Vec3d &initSoma, const DVolume &smoothRay,
                     double &majAxis, double &minAxis, double &somaSurface, double &somaVolume);

    //
    void GetContour(const Vec3d &initSoma, const DVolume &smoothRay,
                    VectorVec3d &contourPtSet, MatXd &resultRayLength);

    void GetContour(const Vec3d &initSoma, const MatXd &resultRayLength,const int Theta, const int Phi,
                    VectorVec3d &contourPtSet);

    void GetAxis(const MatXd &resultRayLength, double &majAxis, double &minAxis);

    int CalcConv3dRadius(const double resolution);

    void GetInnerPtSet(const MatXd &resultRayLength, int Theta, int Phi, double slice,
                       double xVox, double yVox, double zVox,
                       VectorVec3i& innerPtSet);


    void CalcInnerRegionPCA(const VectorVec3i& innerPtSet, Vec3d& firDir,
                            double &lambda1, double &lambda2, double &lambda3);

    int JudgeRegionShapeUsingEigenValue(double lambda1, double lambda2, double lambda3);

    void CorrectRayMatrix(const MatXd &rayLen, int Theta, int Phi, MatXd &newRayLen);

    void MergeSomasByShape(const SVolume &v, StackVec4d& gloSomaStack, 
        const Vec3i& offSet, MatXd& newUx, CVolume &locBinImg );

    void FillShape(const Vec3d& , MatXd&, CVolume&);

};

#endif // PROBEFILTER_H
