#ifndef CORRECTTRACE_H
#define CORRECTTRACE_H
#include"../Function/ineuronprocessobject.h"
#include"../../ngtypes/basetypes.h"
#include"../../ngtypes/volume.h"
#include "../../ngtypes/ParamPack.h"
#include"../binaryfilter.h"
#include"../../ngtypes/tree.h"
#include <iostream>
#include <math.h>
#include<string.h>
#include <ctime>
#include <algorithm>
using namespace Eigen;
using namespace Eigen::internal;
using namespace Eigen::Architecture;
using namespace std;
class CorrectTrace;
typedef std::shared_ptr<CorrectTrace> NGCorrectTrace;
class CorrectTrace :public INeuronProcessObject
{
public:
	static NGCorrectTrace New(){ return NGCorrectTrace(new CorrectTrace()); }
	CorrectTrace();
	~CorrectTrace();
	//bool RunSmooth();
	INEURONPROCESSOBJECT_DEFINE
		//bool Update();
		//virtual ConstIDataPointer GetOutput();
		//virtual IDataPointer ReleaseData();
		void SetLinePoints();//
	void SetLinePoints(const std::vector<size_t> curBranchID);
	void SetLinePoints(size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max);
	void SetLinePoints(const VectorVec5d SemiLine,int flag);//20180409
	void SetParam(NGParamPack arg){ paramPack = arg; }
	ProcStatPointer SaveImageAndSwcForTest(string fliename);
	VectorVec5d GetSmoSemiLine();//20180409
	int SemiFlag = 0;//20180409

	//setinput curve paramPack->origimg backImg
protected:
	NGParamPack paramPack;
	SVolume src;
	CellVec3d myCellCurves;//y,x,z
	VectorVec2d myConnectSet;
	CellVec3d FinalCurves;//x,y,z
	std::vector<Line5d> Currentlines;
	std::vector<VectorVec2i> CurrentIDs;
	int Sign(double);
	Vec3d Meanshiftmodifyposition(const VectorVec4d &currPositions, const Vec3d &centerPoint, const Mat3d &SigmaM);

	Vec3d MeanShiftSinglePointModify(const Vec3d &singleDataPoint, const Mat3d &SigmaM, const SVolume &origImg);


	void Softthresholdingoperation(const VectorVec3d &, const VectorVec3d &, VectorVec3d &);

	void CurvesCorrSparse(const VectorVec3d &, const VectorVec3d &, const VectorVec3d&, VectorVec3d &, VectorVec3d &);
	//The last two argument is the result

	void InitialDR_matrix(const VectorVec3d &, VectorVec3d &, VectorVec3d &);
	//The last two argument is the result

	//void Lishiwei(const VectorVec5d &, const SVolume &,VectorVec5d &);
	//The last argument is the result 

	void Get3DRegion(int &xMin, int &xMax, int &yMin, int &yMax, int &zMin, int &zMax,
		const int xCenter, const int yCenter, const int zCenter,
		const int xOffset, const int yOffset, const int zOffset,
		const int xlower, const int xupper, const int ylower, const int yupper, const int zlower, const int zupper);
	size_t findId(VectorVec3d &Curves1, Vec3d &Curves2);
	void BifurLocalCurvesExtract(const CellVec3d& CellCurves, const VectorVec2d& ConnectSet, XCellVec3d &LocalCurve, VectorVec3d &LableVec);
	//void dataprocesstotalF(SVolume& XX, double threv, SVolume &YY, SVolume &XXb);
	void RevisedCurvesBifur0Modify(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new);
	//---2018/01/18&2018/01/24---//
	void RevisedCurvesBifur0Modify_3line(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new);
	void RevisedCurvesBifur0Modify_3line(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS,
		VectorVec3d &Curve2_new, VectorVec3d &Curve1_new, Vec3d &ModifiedPoint1,
		VectorVec3d &CurveBase_bak1, VectorVec3d &CurveBir_bak1, VectorVec3d &CurveBase_bak2, VectorVec3d &CurveBir_bak2);
	void CurvesCorrforBifurp(const VectorVec3d &Curve1, const VectorVec3d &Curve2, VectorVec3d& Curve3, VectorVec3d& Curve11, VectorVec3d& Curve10);
	void distance_check(const VectorVec3d &Up_curve, const VectorVec3d &down_curve1, const VectorVec3d &down_curve2,
		VectorVec3d &CurveBase_bak1, VectorVec3d &CurveBir_bak1, VectorVec3d &CurveBase_bak2, VectorVec3d &CurveBir_bak2);
	void RevisedCurvesBifur0Modify(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new, Vec3d &ModifiedPoint);
	void weigthvalue(const VectorVec3d &dataL1, SVolume &L_XX3, std::vector<double> &aa_data);
	void SWCFormatChangeConnectModify(double Thres, double x_min, double y_min, double z_min);
	void ThreeDimdCurveResampleModify(const VectorVec3d &DataPoints, double Increase_value, VectorVec3d &ResampleData);
	void Linear_Interpolation(VectorVec3d &MM, const VectorVec3d &M);
	void FastSWCFormatChangeModify(double x_min, double y_min, double z_min);
	void CurvesCorrSubMeanShiftModify0001(const SVolume &ImageS, const VectorVec3d &dataPoints,
		const std::vector<size_t> &LabelMatrix, double Lamda, VectorVec3d &d, VectorVec3d &r, VectorVec3d& dataPoints1);
	Vec3d  MeanshiftmodifypositionS(VectorVec4d &CurrPositions, const Vec3d &CenterPoint, const Mat3d &SigmaM);
	///-----------------------------///
	void RevisedCurvesFixed(VectorVec3d & SingCurveSet, SVolume &ImageS, std::vector<size_t> &LabelMatrix, double ThreShift, double ThreL1);
	void CurvesCorrSubMeanShiftModify000(const SVolume &ImageS, const VectorVec3d &dataPoints, const std::vector<size_t> &LabelMatrix, double Lamda, VectorVec3d &d, VectorVec3d &r, VectorVec3d& dataPoints1);
	void directionc(const Vec3d &x1, Vec3d &x2, Vec3d& x3);
	void conv3d(const SVolume &src, DVolume &fliter, SVolume &output);
	void FastConv3d(const DVolume &fliter, const VectorVec3d &Curve1, const VectorVec3d &Curve2, size_t Numx, size_t Numy, size_t Numz, SVolume &src);
	void AddFliter(DVolume & src, const DVolume &fliter, size_t x, size_t y, size_t z, int radius);
	double Getsum(SVolume & src, DVolume &fliter, size_t x, size_t y, size_t z, int radius);
	void CurvesCorrforBifurp(const VectorVec3d &Curve1, const VectorVec3d &Curve2, VectorVec3d& Curve3);
	void FindSpecialPointInCurve2(const VectorVec3d &Curve2, const VectorVec3d &Curve1, double Thre, Vec3d& Point, size_t &Index);
	Vec3d RevisedCurveBifurTerminalPoint(Vec3d SinglePoint, SVolume &ReImag, Vec3d &PriorDirc, Vec3d &ReferenceP);
	void weak_signal_analysis();

	void CellCurvesaddedPoints(const CellVec3d& CellCurves, VectorVec3d &LabelVec, CellVec4d&);
	void AddingPointsTreeModify(VectorVec4d &curCurves, Vec3d &AddPoint);
	void RevisedCurvesSetFixed(CellVec4d &CurvesSet, SVolume &ImageS, double Threshift, double ThreL1);
	void MakeVec4D(const Vec3d &arg, double a4, Vec4d &aim);
	void Make4DTo3D(const Vec4d &arg, Vec3d &aim);
	void Make4DTo3Dadv(CellVec4d arg, CellVec3d &aim);
	//void SubVector3d(const Vec3d &arg1, const Vec3d &arg2, Vec3d &aim);
	bool IsSame(Vec3d &frist, Vec5d&second);
	///
	void FastSWCFormatChange(double x_min, double y_min, double z_min);
	void SWCFormatChangeConnect(double Thres, double x_min, double y_min, double z_min);
	void SWCFormatChangeSub11(const Vec3d &data1, const Vec3d &data2, const VectorVec3d &Currdata0, double &ds1, double &ds2);
	void SWCFormatChangeSub1(const CellVec3d &Curvesets, MatXd &Distances1, MatXd &Distances2);
	void ThreeDimdCurveResample(const VectorVec3d &DataPoints, double Increase_value, VectorVec3d &ResampleData);
private:

};



#endif
