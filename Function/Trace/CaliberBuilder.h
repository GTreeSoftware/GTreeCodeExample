#ifndef CALIBERBUILDER_H
#define CALIBERBUILDER_H
#include <iostream>
#include <math.h>
#include<string.h>
#include <ctime>
#include <algorithm>
#include"../Function/ineuronprocessobject.h"
#include"../../ngtypes/basetypes.h"
#include"../../ngtypes/volume.h"
#include "../../ngtypes/ParamPack.h"
#include"../binaryfilter.h"
#include"../../ngtypes/tree.h"
//using namespace Eigen;
//using namespace Eigen::internal;
//using namespace Eigen::Architecture;
class CaliberBuilder;
typedef std::shared_ptr<CaliberBuilder> NGCaliberBuilder;
class CaliberBuilder :public INeuronProcessObject
{
public:
    static NGCaliberBuilder New(){ return NGCaliberBuilder(new CaliberBuilder()); }
    CaliberBuilder();
    ~CaliberBuilder();
    //bool RunSmooth();
    INEURONPROCESSOBJECT_DEFINE
    void SetLinePoints();//
    void SetLinePoints(const std::vector<Line5d>& arg, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax);
    void SetLinePoints(const std::vector<size_t> &curBranchID);
    void SetLinePoints(size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max);
    void SetLinePoints(const std::vector<Line3d>& arg);
    void SetParam(NGParamPack arg){ paramPack = arg; }
	void NeuritesShapeReconModifyTotal(CellVec3d &CurveSet, VectorVec2d& ConnectSet, SVolume &OrigImag, CellVec3d &ShapePointSet);
    std::shared_ptr<CellVec3d> GetCaliberPointSet() { 
		if (!(paramPack->caliberBulidData))paramPack->caliberBulidData = std::make_shared<CellVec3d>();
		paramPack->caliberBulidData = caliberPointSet;
		return paramPack->caliberBulidData;
	}

protected:
    bool IsForRange = true;
    std::vector<std::pair<size_t, std::vector<double>>> radiusSet_;
    NGParamPack paramPack;
    //SVolume src;
    CellVec3d myCellCurves;//y,x,z
    VectorVec2d myConnectSet;
    std::shared_ptr<CellVec3d> caliberPointSet;//x,y,z
    std::vector<Line5d> Currentlines;
    std::vector<VectorVec2i> CurrentIDs;
    SVolume commonLabel_;
    MatXd NSIOPImag0, NSIOPImag1, NSIOPGradientImage2, NSIOPGradientImage3, NSIOPGradientImage4, NSIOPGradientImage, NSIOPtmp,
        NSIOPJJk;
    MatXd FDTImag00, FDTImag11;
    MatXd NSOFSS, NSOFImag0, NSOFImag1;
    
    bool FastSWCFormatChange(double x_min, double y_min, double z_min);
    bool SWCFormatChangeConnect(double Thres, double x_min, double y_min, double z_min);
    void SWCFormatChangeSub11(const Vec3d &data1, const Vec3d &data2, const VectorVec3d &Currdata0, double &ds1, double &ds2);
    void SWCFormatChangeSub1(const CellVec3d &Curvesets, MatXd &Distances1, MatXd &Distances2);
    void directionc(const Vec3d &x1, Vec3d &x2, Vec3d& x3);
    //int Sign(double);

    void SampleCylinderDataSetsExtraModify(VectorVec3d &Curves, SVolume & OrigIma, int ScaleSubImage, float Pixelsize, bool IsNew, VectorMatXd &SampleDataset, MatXd &DirecVecs);
    void SBGCSTotalTT(const MatXd& ImageF, MatXd &ImagPrio, double Lamda0, double Threv, double LamdaPrio, MatXd & USeg, MatXd &U, Vec2d & c1_c2);
    void FusedImageD(MatXd& RecImage, MatXd& Imag0, MatXd& Imag1);
    void FusedImageDTranspose(MatXd &Imag0, MatXd &Imag1, MatXd &GradientImage2);
    void SBGCSsub1(const MatXd &ImageF, MatXd& U, double Thres, MatXd & R, Vec2d & c1_c2);
    void NueriteShapeIterativeOP(MatXd &SegImage, MatXd &SegImageInf, MatXd &SegImageprio, double Lamda, double Mu, double Lamda0, MatXd &Ua, MatXd & Ub, MatXd & Va, MatXd &Vb, double &Objectvalue);
    double  NueriteShapeObjectFun(MatXd &SegImage, MatXd & SegImageInf, MatXd &SegImagePrio, MatXd &Va, MatXd &Vb, MatXd &Ua, MatXd & Ub, double Lamda, double Mu, double Lamda0);
    void SoftThresholdingOperation(MatXd &X, MatXd & LamdaVector, MatXd &Y);
    //void NeuritesShapeReconModifyTotal(CellVec3d &CurveSet, VectorVec2d& ConnectSet, SVolume &OrigImag, CellVec3d &ShapePointSet);
    void NeuritesShapeReconModify(CellVec3d & CurveSet, SVolume & OrigImage, SVolume & ImagLabel, int SampleImgsize, VectorVecXd &LabelSet, bool IsNew, CellVec3d &ShapePointSet, VectorVecXd &LabelVecSet);
    void principald(VectorVec3d &dataL, Vec3d&x1);
    void OrthogonalsystemMatch(Vec3d &x1, Vec3d &x2, Vec3d &y1, Vec3d &y2, Vec3d &yy1, Vec3d &yy2);
    double weigthvalue(Vec3d &dataL1, SVolume &L_XX3);
    void SegmentpointModify(VectorMatXd &Segmentdata, MatXd &DircVecs, SVolume &OrigImag, VectorVec2d &ThresMatrix, CellVec7d &dataPointSet);
    void LabelImage(MatXd&LabelImag);
    void LabelImageExtra(MatXd& BinaryImage, VectorVec2d&Seedp, VectorVec2d&Datapoints);
    void LabelImageExtra_other(MatXd& BinaryImage, VectorVec2d&Seedp, VectorVec2d&Datapoints);
    void SegmentPointSceening(CellVec7d &dataPointSet, SVolume &ImagLabel, SVolume &OrigImage, int CurveLabel, VectorVec3d & dataSet, VecXd& Labelvec);
    int  SegmentPointScreeningSub1(VectorVec3d &dataPointS, SVolume &ImagLabel, int CurveLabel);
    void FilleddataPointModify(VectorVec7d &dataPoint, SVolume & OrigImag, VectorVec3d &dataPointS);
    void CurvesDircs(VectorVec3d &Curves, MatXd&DirecVecs);
};

#endif