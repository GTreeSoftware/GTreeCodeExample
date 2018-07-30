#include"CorrectTrace.h"
#include"../IO/imagewriter.h"
#include "../Function/Trace/CaliberBuilder.h"
#include<Eigen\Dense>
#include<Eigen\Eigenvalues>
#include <Eigen/Core>
#pragma warning(disable:4996)
/*core work*/

CorrectTrace::CorrectTrace()
{
}

CorrectTrace::~CorrectTrace()
{
}
/*
save current image and swc
*/
/*don't use*/
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(CorrectTrace)
/*don't use*/
INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(CorrectTrace, NeuronPopulation)
ProcStatPointer CorrectTrace::SaveImageAndSwcForTest(string fliename){
	if (Currentlines.empty()){
		SetLinePoints();
	}
	for (int i = 1; i < CurrentIDs.size(); ++i){
		if (CurrentIDs[i][0][1] != 0){
			//iscorrect = false;
			CurrentIDs.erase(CurrentIDs.begin() + i);
			Currentlines.erase(Currentlines.begin() + i);
		}
	}
	if (paramPack->OrigImage) {
		//QString saveName = QFileDialog::getSaveFileName(this, "save current image", paramPack->defaultDir, "TIF (*.tif)");
		//paramPack->defaultDir = saveName.section('/', 0, -2);
		NGImageWriter writer = ImageWriter::New();
		writer->SetOutputFileName(fliename + ".tif");
		writer->SetInput(paramPack->OrigImage);//OrigImage
		auto wres = writer->Update();
		if (!wres->success())
			return wres;
		std::cout << "save image: " << fliename + ".tif" << std::endl;
	}
	FILE* fp = fopen((fliename + ".swc").c_str(), "w");
	int id = 0;
	size_t curLineID, parentLineID, nearestID, pid;
	std::vector<std::vector<int>> idList(Currentlines.size());//save 1 column id
	for (size_t i = 0; i < idList.size(); ++i)
		idList[i].resize(Currentlines[i].size());//warning !!! : VS2013 set default value as 0
	//1st line
	{
		curLineID = 0;
		const Line5d &line = Currentlines[curLineID];
		if (line.empty()){
			printf("error in SaveATree. empty 1st tree.\n");
			MAKEPROCESSSTATUS(resSta, false, className_, "empty 1st tree.");
			return resSta;
		}
		idList[curLineID][0] = ++id;
		fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", id, (line[0](0) - paramPack->xMin_), (line[0](1) - paramPack->yMin_), (line[0](2) - paramPack->zMin_));//root
		for (size_t i = 1; i < line.size(); ++i){
			idList[curLineID][i] = ++id;
			fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", id, (line[i](0) - paramPack->xMin_), (line[i](1) - paramPack->yMin_), (line[i](2) - paramPack->zMin_), id - 1);
		}
	}
	//left lines
	double mindist, curDist;
	for (curLineID = 1; curLineID < Currentlines.size(); ++curLineID)
	{
		mindist = 100000;
		curDist = 99999;
		const Line5d &line = Currentlines[curLineID];
		const Vec3d &headNode = line.front().block(0, 0, 3, 1);
		tree<LineID>::iterator CurrentLineID;
		KPTREE::SearchLineIDInTree(*paramPack->activeTree, CurrentIDs[curLineID][0][0], CurrentLineID);
		auto parentID = paramPack->activeTree->m_Topo.parent(CurrentLineID)->id;
		for (size_t j = 0; j < CurrentIDs.size(); ++j){
			if (parentID == CurrentIDs[j][0][0] && curLineID != j){
				nearestID = KPTREE::FindNearestID2(headNode, Currentlines[j], curDist);
				parentLineID = j;
			}
		}
		const Vec3d &beginNode = Currentlines[parentLineID][nearestID].block(0, 0, 3, 1);
		//std::cout << "parentLineID= " << parentLineID << " nearestID= " << nearestID << std::endl;
		//std::cout << "headNode=" << headNode << std::endl << " nearestpoints=" << beginNode << std::endl;
		double DistantOfConnection = (headNode - beginNode).transpose()*(headNode - beginNode);
		std::cout << "DistantOfConnection=" << DistantOfConnection << std::endl;
		//if (DistantOfConnection > 99999) {
		//	NG_ERROR_MESSAGE("how can different so large? please debug.\n");
		//}
		if (DistantOfConnection > 30) {
			NG_ERROR_MESSAGE("two points are different,maybe not interconnected!");
			//1st line
			{
				const Line5d &line = Currentlines[curLineID];
				if (line.empty()){
					printf("error in SaveATree. empty 1st tree.\n");
					MAKEPROCESSSTATUS(resSta, false, className_, "empty 1st tree.");
					return resSta;
				}
				idList[curLineID][0] = ++id;
				fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", id, (line[0](0) - paramPack->xMin_), (line[0](1) - paramPack->yMin_), (line[0](2) - paramPack->zMin_));//root
				for (size_t i = 1; i < line.size(); ++i){
					idList[curLineID][i] = ++id;
					fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", id, (line[i](0) - paramPack->xMin_), (line[i](1) - paramPack->yMin_), (line[i](2) - paramPack->zMin_), id - 1);
				}
			}
			continue;
		}
		pid = idList[parentLineID][nearestID];
		//TODO: debug
		if (pid > 1000000) {
			NG_ERROR_MESSAGE("how can pid so large? please debug.");
		}
		if (pid == 0) {
			printf("error in SaveATree. the parent vertex id is 0, please debug.");
			MAKEPROCESSSTATUS(resSta, false, className_, "the parent vertex id is 0, please debug.");
			return resSta;
		}
		idList[curLineID][0] = ++id;
		fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", id, (line[0](0) - paramPack->xMin_), (line[0](1) - paramPack->yMin_), (line[0](2) - paramPack->zMin_), pid);
		for (size_t i = 1; i < line.size(); ++i){
			idList[curLineID][i] = ++id;
			fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", id, (line[i](0) - paramPack->xMin_), (line[i](1) - paramPack->yMin_), (line[i](2) - paramPack->zMin_), id - 1);
		}
	}
	fclose(fp);
	std::cout << "save swc: " << fliename + ".swc" << std::endl;
	MAKEPROCESSSTATUS(resSta, true, className_, "");
	return resSta;

}
/*
set the line to smooth and correct
default function choose all lines in current image
*/
void CorrectTrace::SetLinePoints(){
	KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, Currentlines, CurrentIDs,
		paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
}
/*
set the line to smooth and correct
choose  part of lines in current image
*/
void CorrectTrace::SetLinePoints(const std::vector<size_t> curBranchID){
	KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, curBranchID, Currentlines, CurrentIDs,
		paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
}
/*
set the line to smooth and correct
local range lines
*/
void CorrectTrace::SetLinePoints(size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max){
	x_min = x_min > paramPack->xMin_ ? x_min : paramPack->xMin_;
	y_min = y_min > paramPack->yMin_ ? y_min : paramPack->yMin_;
	z_min = z_min > paramPack->zMin_ ? z_min : paramPack->zMin_;
	x_max = x_max < paramPack->xMax_ ? x_max : paramPack->xMax_;
	y_max = y_max < paramPack->yMax_ ? y_max : paramPack->yMax_;
	z_max = z_max < paramPack->zMax_ ? z_max : paramPack->zMax_;
	KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, Currentlines, CurrentIDs,
		x_min, x_max, y_min, y_max, z_min, z_max);
}
/*
set the line to smooth and correct
choose  part of lines in current image
*/
void CorrectTrace::SetLinePoints(const VectorVec5d SemiLine,int SemiFlag1){
	Currentlines .push_back(SemiLine);
	SemiFlag = SemiFlag1;
}
VectorVec5d CorrectTrace::GetSmoSemiLine(){
	VectorVec5d lline;
	for (size_t i = 0; i < FinalCurves[0].size();++i){
		Vec5d po;
		po[0] = FinalCurves[0][i][0];
		po[1] = FinalCurves[0][i][1];
		po[2] = FinalCurves[0][i][2];
		po[3] = 1;
		po[4] = 0;
		lline.push_back(po);
	}
	return lline;
}
/* run smooth lines*/
ProcStatPointer CorrectTrace::Update(){
	//find inner part of target curve and connective curves. get connect information from m_Topo
	//bool iscorrect = true;
	if (Currentlines.empty()){
		MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
		return resSta;
	}
	for (int i = 1; i < CurrentIDs.size(); ++i){
		if (CurrentIDs[i][0][1] != 0){
			//iscorrect = false;
			CurrentIDs.erase(CurrentIDs.begin() + i);
			Currentlines.erase(Currentlines.begin() + i);
		}
	}
	//if (!iscorrect)return iscorrect;
	SWCFormatChangeConnect(1, paramPack->xMin_, paramPack->yMin_, paramPack->zMin_);
	//FastSWCFormatChange(paramPack->xMin_, paramPack->yMin_, paramPack->zMin_);

	NG_CREATE_DYNAMIC_CAST(SVolume, mysrc, paramPack->OrigImage);
	if (mysrc->IsEmpty())return false;
	src.QuickCopy(*mysrc);
	CellVec3d caliberPointSet;
	std::vector<Line3d> Curves;
	for (int i = 0; i < myCellCurves.size(); ++i){
		Curves.push_back(myCellCurves[i]);
	}
	NGCaliberBuilder cb = CaliberBuilder::New();
	cb->SetLinePoints(Curves);
	cb->SetParam(paramPack);
	cb->Update();
	caliberPointSet = *(cb->GetCaliberPointSet());
	size_t nx = src.x();
	size_t ny = src.y();
	size_t nz = src.z();
	SVolume ImageS_tmp;
	ImageS_tmp.SetSize(nx, ny, nz);
	ImageS_tmp.SetZero();
	for (size_t i = 0; i < caliberPointSet.size();++i)
	{
		for (size_t j = 0; j < caliberPointSet[i].size();++j)
		{
			int nxx = caliberPointSet[i][j](0);
			int nyy = caliberPointSet[i][j](1);
			int nzz = caliberPointSet[i][j](2);
			ImageS_tmp(nxx, nyy, nzz) = std::max(src.operator()(nxx, nyy, nzz), (unsigned short)8);
		}
	}
	src.QuickCopy(ImageS_tmp);
	caliberPointSet.clear();
/*	size_t nx = src.x();
	size_t ny = src.y();
	size_t nz = src.z();
	SVolume ImageS_tmp, backimage;
	ImageS_tmp.SetSize(nx, ny, nz);
	for (int nxx = 0; nxx < nx; nxx++){
		for (int nyy = 0; nyy < ny; nyy++){
			for (int nzz = 0; nzz < nz; nzz++){
				ImageS_tmp(nxx, nyy, nzz) = std::max(src.operator()(nxx, nyy, nzz), (unsigned short)8);
			}
		}
	}
	NGBinaryFilter filter = BinaryFilter::New();
	filter->CalcBackImage(ImageS_tmp, backimage, 5);
	for (size_t nxx = 0; nxx < nx; nxx++){
		for (size_t nyy = 0; nyy < ny; nyy++){
			for (size_t nzz = 0; nzz < nz; nzz++){
				src.operator()(nxx, nyy, nzz) = std::max(src.operator()(nxx, nyy, nzz) - backimage.operator()(nxx, nyy, nzz), 0);
			}
		}
	}
	backimage.Clear();*/
	ImageS_tmp.Clear();
	//Currentlines.clear();

	//FastSWCFormatChange(1, paramPack->xMin_, paramPack->yMin_, paramPack->zMin_);

	weak_signal_analysis();
	//check

	//for (int i = 0; i < FinalCurves.size(); ++i){
	//	if (!IsSame(FinalCurves[i][0],Currentlines[i][0])){
	//		std::cout << "Line[" << i << "]: change head!" << std::endl;
	//	}
	//	if (!IsSame(FinalCurves[i][FinalCurves[i].size() - 1], Currentlines[i][Currentlines[i].size()-1])){
	//		std::cout << "Line[" << i << "]: change end!" << std::endl;
	//		iscorrect = false;
	//	}
	//}
	//if (!iscorrect)return iscorrect;
	if (SemiFlag==0){
		KPTREE::InsertCurrentLine(paramPack->activeTree->m_curveList, FinalCurves, CurrentIDs);
 }
else{
		//SemiFlag = 0;
	}
	
	MAKEPROCESSSTATUS(resSta, true, className_, "");
	return resSta;
}
bool CorrectTrace::IsSame(Vec3d &frist, Vec5d&second){
	double tmp = std::abs(frist[0] - second[0]) + std::abs(frist[1] - second[1]) + std::abs(frist[2] - second[2]);
	if (tmp < 0.0001)return true;
	else
	{
		return false;
	}
}
//ConstIDataPointer CorrectTrace::GetOutput(){
//
//
//	ConstIDataPointer kk;
//	return kk;
//
//}
//IDataPointer CorrectTrace::ReleaseData(){
//
//	//paramPack.reset();
//	//src.ReleaseProcessObject();
//	//CellCurves.clear();
//	//ConnectSet.clear();
//	IDataPointer kk;
//	return kk;
//
//}

//has
int CorrectTrace::Sign(double n)
{
	if (0 == n)
		return 0;
	else if (n > 0)
		return 1;
	else
		return -1;
}



Vec3d CorrectTrace::Meanshiftmodifyposition(const VectorVec4d &currPositions, const Vec3d &centerPoint, const Mat3d &SigmaM){
	size_t Num0 = currPositions.size();
	Vec3d returnPos; returnPos.setZero();// = VecXd::Zero(3, 1);
	Vec3d temCurrPt; temCurrPt.setZero();// = VecXd::Zero(3, 1);
	double ss = 0;

	for (size_t i = 0; i < Num0; ++i) {
		for (size_t j = 0; j < 3; ++j)	temCurrPt(j) = currPositions[i](j) - centerPoint(j);
		//double sss = 0;
		double sss = temCurrPt.transpose()*SigmaM*temCurrPt;
		//currv = temCurrPt.norm();
		double weights0 = std::exp(-sss / 9.0);
		double weights1 = currPositions[i](3) / 100.0;

		for (size_t z = 0; z < 3; ++z)
			returnPos(z) += weights1*weights0*currPositions[i](z);

		ss += weights0 * weights1;
	}
	returnPos /= ss;

	return returnPos;
}

Vec3d CorrectTrace::MeanShiftSinglePointModify(const Vec3d &singleDataPoint, const Mat3d &SigmaM, const SVolume &origImg){
	int Nxx, Nyy, Nzz;
	Nxx = origImg.x();
	Nyy = origImg.y();
	Nzz = origImg.z();

	int ptX = std::round(singleDataPoint(0));
	int ptY = std::round(singleDataPoint(1));
	int ptZ = std::round(singleDataPoint(2));
	/*
	int maxxx, maxyy, maxzz, minxx, minyy, minzz;
	int kernel_x = 2 + paramPack->kernel_size * 2;
	int kernel_y = 2 + paramPack->kernel_size * 2;
	int kernel_z = 1 + paramPack->kernel_size * 2;
	Get3DRegion(minxx, maxxx, minyy, maxyy, minzz, maxzz, ptX, ptY, ptZ, kernel_x, kernel_y, kernel_z, 1, Nxx, 1, Nyy, 1, Nzz);
	*/
	//VectorVec4d currPositions;
	//currPositions.reserve((maxxx - minxx + 1)*(maxyy - minyy + 1)*(maxzz - minzz + 1));//preserve memory,but did not put any elements
	//---
	vector<int> Indexx, Indeyy, Indezz;//kernel size
	for (size_t i = std::max(ptX - 8, 0); i <= std::min(ptX + 8, Nxx-1); ++i)		Indexx.push_back(i);
	for (size_t i = std::max(ptY - 8, 0); i <= std::min(ptY + 8, Nyy-1); ++i)		Indeyy.push_back(i);
	for (size_t i = std::max(ptZ - 6, 0); i <= std::min(ptZ + 6, Nzz-1); ++i)		Indezz.push_back(i);
	//int Indexx_len = std::min(ptX + 4, Nxx) - std::max(ptX - 4, 1);
	//int Indeyy_len = std::min(ptY + 4, Nyy) - std::max(ptY - 4, 1);
	//int Indezz_len = std::min(ptZ + 3, Nzz) - std::max(ptZ - 3, 1);
	VectorVec4d CurrPositions(Indexx.size()*Indeyy.size()*Indezz.size());
	std::for_each(CurrPositions.begin(), CurrPositions.end(), [&](Vec4d &arg){arg.setZero(); });
	//---
	/*change*/
	//Vec3d returnPos; returnPos.setZero();// = VecXd::Zero(3, 1);
	//Vec3d temCurrPt; temCurrPt.setZero();// = VecXd::Zero(3, 1);
	//double ss = 0;
	/*....*/
	//clock_t beg = clock();
	size_t kk = 0;

	for (int i = 0; i < Indexx.size(); ++i) {
		for (int j = 0; j < Indeyy.size(); ++j) {
			for (int ij = 0; ij < Indezz.size(); ++ij) {
				//Vec4d tmp; tmp.setZero();
				CurrPositions[kk](0) = Indexx[i];
				CurrPositions[kk](1) = Indeyy[j];
				CurrPositions[kk](2) = Indezz[ij];
				CurrPositions[kk](3) = origImg(Indexx[i], Indeyy[j], Indezz[ij]);
				//CurrPositions.push_back(tmp);
				kk = kk + 1;
			}
		}
	}
	Vec3d Corrdatap;
	if (!CurrPositions.empty())
		Corrdatap = MeanshiftmodifypositionS(CurrPositions, singleDataPoint, SigmaM);
	else
	{
		//Corrdatap = [];
	}
	return Corrdatap;
}

Vec3d  CorrectTrace::MeanshiftmodifypositionS(VectorVec4d &CurrPositions, const Vec3d &CenterPoint, const Mat3d &SigmaM){
	size_t Num0 = CurrPositions.size();
	Vec3d Position0;
	Position0.setZero();
	double ss = 0;
	for (size_t i = 0; i < Num0; i++)
	{
		Vec3d Currv;
		Currv(0) = CurrPositions[i](0) - CenterPoint(0);
		Currv(1) = CurrPositions[i](1) - CenterPoint(1);
		Currv(2) = CurrPositions[i](2) - CenterPoint(2);
		double sss = Currv.transpose()*SigmaM*Currv;
		double weights0 = std::exp(-sss / 9.0);//exp
		double weights1 = CurrPositions[i](3) / 100.0;//p(x)
		Position0(0) += weights1*weights0*CurrPositions[i](0);
		Position0(1) += weights1*weights0*CurrPositions[i](1);
		Position0(2) += weights1*weights0*CurrPositions[i](2);
		ss = ss + weights0 * weights1;
	}
	if (ss>0.1)
		Position0 /= ss;
	else
	{
		Position0 = CenterPoint;
		//std::cout << "Error! Please Debug!" << std::endl;
	}

	return Position0;
}

void CorrectTrace::Softthresholdingoperation(const VectorVec3d &X, const VectorVec3d &lamdaVector, VectorVec3d& returnData){
	size_t Nxx = X.size(), Nyy = 3;
	returnData.resize(Nxx);
	std::for_each(returnData.begin(), returnData.end(), [&](Vec3d &arg){arg.setZero(); });

	for (size_t j = 0; j < Nyy; ++j){
		for (size_t i = 0; i < Nxx; ++i){
			double currData = X[i](j);
			double currThreshold = lamdaVector[i](j);
			if (std::abs(currData) > currThreshold)
				returnData[i](j) = Sign(currData)*(std::abs(currData) - currThreshold);//Change Sign to TraceUtil::sign
		}
	}
}

void CorrectTrace::CurvesCorrSparse(const VectorVec3d &dataPoints, const VectorVec3d &r, const VectorVec3d& LamdaMatrix, VectorVec3d &d_new, VectorVec3d &r_new){
	if (dataPoints.size() < 2){
		printf("error in Curvescorrsparse \n");
		return;
	}
	size_t numP = dataPoints.size();
	VectorVec3d currMatrix(numP - 2);
	std::for_each(currMatrix.begin(), currMatrix.end(), [&](Vec3d& arg){ arg.setZero(); });

	for (size_t i = 1; i < numP - 1; ++i)
		currMatrix[i - 1] = 2 * dataPoints[i] - dataPoints[i - 1] - dataPoints[i + 1] + r[i - 1];

	/*VectorVec3d temMatr(numP - 2);
	std::for_each(temMatr.begin(), temMatr.end(), [&](Vec3d &arg){ arg = Lamda*Vec3d::Ones(3, 1); });
	*/
	Softthresholdingoperation(currMatrix, LamdaMatrix, d_new);

	size_t countNum = 0;
	r_new.resize(currMatrix.size());
	//std::for_each(r_new.begin(), r_new.end(), [&](Vec3d){ r_new[countNum] = currMatrix[countNum] - d_new[countNum]; ++countNum; });
	std::for_each(r_new.begin(), r_new.end(), [&](Vec3d &arg){ arg = currMatrix[countNum] - d_new[countNum]; ++countNum; });

}

void CorrectTrace::InitialDR_matrix(const VectorVec3d &dataPoints, VectorVec3d &d, VectorVec3d &r)
{
	//std::cout << "InitialDR_matrix start!" << std::endl;;/////////////////////////////////////////////////////

	if (dataPoints.size() < 2llu) {
		printf("error in InitialDR_matrix\n");
		return;
	}
	size_t numPt = dataPoints.size();
	d.resize(numPt - 2llu);
	std::for_each(d.begin(), d.end(), [&](Vec3d &arg){arg.setZero(); });
	r.resize(numPt - 2llu);
	r = d;
	for (size_t i = 1; i < numPt - 1llu; ++i)
		d[i - 1] = 2 * dataPoints[i] - dataPoints[i - 1] - dataPoints[i + 1];

	//std::cout << "InitialDR_matrix over!" << std::endl;;/////////////////////////////////////////////////////
}

//has
void CorrectTrace::Get3DRegion(int &xMin, int &xMax, int &yMin, int &yMax, int &zMin, int &zMax,
	const int xCenter, const int yCenter, const int zCenter,
	const int xOffset, const int yOffset, const int zOffset,
	const int xlower, const int xupper, const int ylower, const int yupper, const int zlower, const int zupper)
{
	xMin = max(xlower, xCenter - xOffset);
	xMax = min(xupper, xCenter + xOffset);
	yMin = max(ylower, yCenter - yOffset);
	yMax = min(yupper, yCenter + yOffset);
	zMin = max(zlower, zCenter - zOffset);
	zMax = min(zupper, zCenter + zOffset);
}
void CorrectTrace::FastSWCFormatChange(double x_min, double y_min, double z_min){
	if (Currentlines.empty()){
		///error;
	}
	myCellCurves.clear();
	myConnectSet.clear();
	VectorVec3d DataPoints, ResampleData;
	for (size_t i = 0; i < Currentlines.size(); ++i){
		DataPoints.clear();
		//  DataPoints=SWCdata(CurveIndex(i):CurveIndex(i+1)-1,[4,3,5])';
		for (int j = 0; j < Currentlines[i].size(); ++j){
			/*
			if (j == 0 && i != 0 && Currentlines[i].size()>2){
			continue;
			}*/
			Vec3d tmp;
			tmp[0] = Currentlines[i][j][0] - x_min;
			tmp[1] = Currentlines[i][j][1] - y_min;
			tmp[2] = Currentlines[i][j][2] - z_min;
			DataPoints.push_back(tmp);
		}

		ThreeDimdCurveResample(DataPoints, 1, ResampleData);
		myCellCurves.push_back(ResampleData);
		//myCellCurves.push_back(DataPoints);
	}
	Vec2d ConVec;
	for (size_t i = 0; i < CurrentIDs.size(); ++i){
		if (i == 0){
			ConVec[0] = 0;
			ConVec[1] = 0;
		}
		else
		{
			bool flag = false;
			/*	auto &topo = paramPack->activeTree->m_Topo;
			auto CurrentLineID = std::find(topo.begin() + 1, topo.end(), LineID(CurrentIDs[i][0][0]));*/
			tree<LineID>::iterator CurrentLineID;
			KPTREE::SearchLineIDInTree(*paramPack->activeTree, CurrentIDs[i][0][0], CurrentLineID);
			auto parentID = paramPack->activeTree->m_Topo.parent(CurrentLineID)->id;
			for (size_t j = 0; j < CurrentIDs.size(); ++j){
				if (i != j){
					if (parentID == CurrentIDs[j][0][0])
					{
						flag = true;
						ConVec[0] = j + 1;
						ConVec[1] = 0;
					}
				}
			}
			if (!flag){
				ConVec[0] = 0;
				ConVec[1] = 0;
			}
		}
		//Vec2d ConVec;
		//ConVec[0] = i;
		//Vec3d Currentpoint = Currentlines[i][0].block(0, 0, 3, 1);
		//ConVec[1]=KPTREE::FastSearchNearestLine(Currentpoint, Currentlines);
		myConnectSet.push_back(ConVec);
		/*paramPack->activeTree->m_Topo->*/
	}

}

void CorrectTrace::FastSWCFormatChangeModify(double x_min, double y_min, double z_min){
	if (Currentlines.empty()){
		///error;
	}
	myCellCurves.clear();
	myConnectSet.clear();
	VectorVec3d DataPoints, ResampleData;
	for (size_t i = 0; i < Currentlines.size(); ++i){
		DataPoints.clear();
		//  DataPoints=SWCdata(CurveIndex(i):CurveIndex(i+1)-1,[4,3,5])';
		for (int j = 0; j < Currentlines[i].size(); ++j){
			if (j == 0 && i != 0 && Currentlines[i].size()>2){
				continue;
			}
			Vec3d tmp;
			tmp[0] = Currentlines[i][j][0] - x_min;
			tmp[1] = Currentlines[i][j][1] - y_min;
			tmp[2] = Currentlines[i][j][2] - z_min;
			DataPoints.push_back(tmp);
		}

		ThreeDimdCurveResample(DataPoints, 1, ResampleData);
		myCellCurves.push_back(ResampleData);
	}
	Vec2d ConVec;
	for (size_t i = 0; i < CurrentIDs.size(); ++i){
		if (i == 0){
			ConVec[0] = 0;
			ConVec[1] = 0;
		}
		else
		{
			bool flag = false;
			/*	auto &topo = paramPack->activeTree->m_Topo;
			auto CurrentLineID = std::find(topo.begin() + 1, topo.end(), LineID(CurrentIDs[i][0][0]));*/
			tree<LineID>::iterator CurrentLineID;
			KPTREE::SearchLineIDInTree(*paramPack->activeTree, CurrentIDs[i][0][0], CurrentLineID);
			auto parentID = paramPack->activeTree->m_Topo.parent(CurrentLineID)->id;
			for (size_t j = 0; j < CurrentIDs.size(); ++j){
				if (i != j){
					if (parentID == CurrentIDs[j][0][0])
					{
						flag = true;
						ConVec[0] = j + 1;
						ConVec[1] = 0;
					}
				}
			}
			if (!flag){
				ConVec[0] = 0;
				ConVec[1] = 0;
			}
		}
		//Vec2d ConVec;
		//ConVec[0] = i;
		//Vec3d Currentpoint = Currentlines[i][0].block(0, 0, 3, 1);
		//ConVec[1]=KPTREE::FastSearchNearestLine(Currentpoint, Currentlines);
		myConnectSet.push_back(ConVec);
		/*paramPack->activeTree->m_Topo->*/
	}
}
/*
change swc format and resample lines
*/
void CorrectTrace::SWCFormatChangeConnect(double Thres, double x_min, double y_min, double z_min)
{
	if (Currentlines.empty()){
		///error;
	}
	myCellCurves.clear();
	myConnectSet.clear();
	VectorVec3d DataPoints, ResampleData;
	for (size_t i = 0; i < Currentlines.size(); ++i){
		DataPoints.clear();
		//  DataPoints=SWCdata(CurveIndex(i):CurveIndex(i+1)-1,[4,3,5])';
		for (int j = 0; j < Currentlines[i].size(); ++j){
			//if (j == 0 && i != 0 && Currentlines[i].size()>2){
			//	continue;}
			Vec3d tmp;
			tmp[0] = Currentlines[i][j][0] - x_min;
			tmp[1] = Currentlines[i][j][1] - y_min;
			tmp[2] = Currentlines[i][j][2] - z_min;

			DataPoints.push_back(tmp);
		}
		if (DataPoints.size()<2)
		{
			ResampleData = DataPoints;
		}
		else
		{
			ThreeDimdCurveResample(DataPoints, 1, ResampleData);
		}
		myCellCurves.push_back(ResampleData);
	}
	//ConnectSet=cell(1,NumCurves);
	MatXd DistMatrixStarP, DistMatrixEndP;
	SWCFormatChangeSub1(myCellCurves, DistMatrixStarP, DistMatrixEndP);
	//cout << "DistMatrixStarP=" << DistMatrixStarP << endl;
	//cout << "DistMatrixEndP=" << DistMatrixEndP << endl;
	for (int i = 0; i < Currentlines.size(); i++){
		MatXd	CurrVec1 = DistMatrixStarP.col(i);
		MatXd	CurrVec2 = DistMatrixEndP.col(i);
		//std::cout << CurrVec1 << std::endl;
		//std::cout << CurrVec2 << std::endl;
		MatrixXf::Index minRow1, minCol1;
		double min1 = CurrVec1.minCoeff(&minRow1, &minCol1);
		MatrixXf::Index minRow2, minCol2;
		double min2 = CurrVec2.minCoeff(&minRow2, &minCol2);
		Vec2d ConVec;
		ConVec.setZero();
		if (min1<Thres)
			ConVec[0] = minRow1 + 1;
		if (min2<Thres)
			ConVec[1] = minRow2 + 1;
		myConnectSet.push_back(ConVec);
	}
}
//*****************//
//**82018/1/24**//
//change swc format and resample lines//
void CorrectTrace::SWCFormatChangeConnectModify(double Thres, double x_min, double y_min, double z_min)
{
	if (Currentlines.empty()){
		///error;
	}
	myCellCurves.clear();
	myConnectSet.clear();
	/////
	VectorVec3d Col7;
	for (int i = 0; i < CurrentIDs.size(); i++)
	{
		if (i > 0)
		{
			tree<LineID>::iterator CurrentLineID;
			KPTREE::SearchLineIDInTree(*paramPack->activeTree, CurrentIDs[i][0][0], CurrentLineID);
			auto parentID = paramPack->activeTree->m_Topo.parent(CurrentLineID)->id;
			size_t parentIDj = 0;
			for (size_t j = 0; j < CurrentIDs.size(); ++j){
				if (i != j){
					if (parentID == CurrentIDs[j][0][0])   parentIDj = j + 1;
				}
			}
			size_t id = 100000;
			double mindist = 100000, curDist;
			double xDist, yDist, zDist;
			for (size_t j = 0; j < Currentlines[parentIDj - 1].size(); j++)
			{
				xDist = std::abs(Currentlines[i][0](0) - Currentlines[parentIDj - 1][j](0));
				yDist = std::abs(Currentlines[i][0](1) - Currentlines[parentIDj - 1][j](1));
				zDist = std::abs(Currentlines[i][0](2) - Currentlines[parentIDj - 1][j](2));
				if (xDist < 1 && yDist < 1 && zDist < 1)  {
					curDist = xDist + yDist + zDist;
					if (curDist < mindist) {
						mindist = curDist;
						id = j;
					}
				}

			}
			Vec3d tmp;
			tmp[0] = parentIDj;//curve id
			tmp[1] = id;//point id
			tmp[2] = 0;//is flag 
			Col7.push_back(tmp);
		}
	}
	for (int i = 0; i < Col7.size(); i++)
	{
		for (int j = i; j < Col7.size(); j++)
		{
			if (i != j&&Col7[j][2] == 0)
			{
				if (Col7[i](0) == Col7[j](0) && Col7[i](1) == Col7[j](1))
				{
					Col7[i](2) == 1;
					Col7[j](2) == 1;
					//cishu++;
					//AnaMat.push_back([i,cishu]);
				}
			}
		}
	}



	/////
	/*
	tree<LineID>::iterator CurrentLineID;
	KPTREE::SearchLineIDInTree(*paramPack->activeTree, CurrentIDs[i][0][0], CurrentLineID);
	auto parentID = paramPack->activeTree->m_Topo.parent(CurrentLineID)->id;
	for (size_t j = 0; j < CurrentIDs.size(); ++j){
	if (i != j){
	if (parentID == CurrentIDs[j][0][0])
	{
	flag = true;
	ConVec[0] = j+1;
	ConVec[1] = 0;
	}
	}
	}*/
	VectorVec3d DataPoints, ResampleData;
	for (size_t i = 0; i < Currentlines.size(); ++i){
		DataPoints.clear();
		//  DataPoints=SWCdata(CurveIndex(i):CurveIndex(i+1)-1,[4,3,5])';
		for (int j = 0; j < Currentlines[i].size(); ++j){
			if (j == 0 && i != 0 && Currentlines[i].size()>2){
				continue;
			}
			Vec3d tmp;
			tmp[0] = Currentlines[i][j][0] - x_min;
			tmp[1] = Currentlines[i][j][1] - y_min;
			tmp[2] = Currentlines[i][j][2] - z_min;

			DataPoints.push_back(tmp);
		}

		ThreeDimdCurveResample(DataPoints, 1, ResampleData);
		myCellCurves.push_back(ResampleData);
	}
	//ConnectSet=cell(1,NumCurves);
	MatXd DistMatrixStarP, DistMatrixEndP;
	SWCFormatChangeSub1(myCellCurves, DistMatrixStarP, DistMatrixEndP);
	//cout << "DistMatrixStarP=" << DistMatrixStarP << endl;
	//cout << "DistMatrixEndP=" << DistMatrixEndP << endl;
	for (int i = 0; i < Currentlines.size(); i++){
		//CurrVec1=DistMatrixStarP(:,i);
		MatXd	CurrVec1 = DistMatrixStarP.col(i);
		MatXd	CurrVec2 = DistMatrixEndP.col(i);
		//std::cout << CurrVec1 << std::endl;
		//std::cout << CurrVec2 << std::endl;
		MatrixXf::Index minRow1, minCol1;
		double min1 = CurrVec1.minCoeff(&minRow1, &minCol1);
		MatrixXf::Index minRow2, minCol2;
		double min2 = CurrVec2.minCoeff(&minRow2, &minCol2);
		Vec2d ConVec;
		ConVec.setZero();
		if (min1<Thres)
			ConVec[0] = minRow1 + 1;
		if (min2<Thres)
			ConVec[1] = minRow2 + 1;
		myConnectSet.push_back(ConVec);
	}
}


//big bug
void CorrectTrace::ThreeDimdCurveResample(const VectorVec3d &DataPoints, double Increase_value, VectorVec3d &ResampleData)
{

	//VectorVec3d DataPoints_addlast = DataPoints;
	//DataPoints.push_back(DataPoints[DataPoints.size() - 1]);
	size_t Num_DataPoints = DataPoints.size();
	VectorVec3d AdjacentPoints;


	//Index=zeros(1,5*Num_DataPoints);
	std::vector<double>  DistanceVec, Flag_Point;
	//size_t allpoints = std::round((DataPoints[Num_DataPoints - 1] - DataPoints[0]).norm() / Increase_value)*2;
	//if (allpoints < Num_DataPoints)
	size_t allpoints = 6* Num_DataPoints;//20180403  , 5* Num_DataPoints;
	std::vector<size_t> Index(allpoints);
	ResampleData.resize(allpoints);
	std::for_each(ResampleData.begin(), ResampleData.end(), [&](Vec3d& arg){ arg.setZero(); });
	ResampleData[0] = DataPoints[0];
	std::for_each(Index.begin(), Index.end(), [&](size_t& arg){ arg = -1; });
	Index[0] = 0;
	size_t i;
	for (i = 0; i < allpoints - 1; i++){
		if (Index[i] < Num_DataPoints - 1){
			AdjacentPoints.clear();
			for (size_t j = Index[i] + 1; j < std::min(Index[i] + 4, Num_DataPoints); j++){
				AdjacentPoints.push_back(DataPoints[j]);
			}
			if (AdjacentPoints.size() > 0){//2018 04 03, >1
				DistanceVec.clear();
				for (size_t j = 0; j < AdjacentPoints.size(); j++){
					DistanceVec.push_back((AdjacentPoints[j] - ResampleData[i]).norm());
				}
				// Flag_Point=find(DistanceVec>Increase_value);
				Flag_Point.clear();
				for (size_t k = 0; k < DistanceVec.size(); k++){
					if (DistanceVec[k] > Increase_value)Flag_Point.push_back(k);
				}
				//
				if (Flag_Point.size() == 0){
					//Index(i+1)=min(Index(i)+3,Num_DataPoints);
					Index[i + 1] = std::min(Index[i] + 3, Num_DataPoints - 1);
					ResampleData[i + 1] = DataPoints[std::min(Index[i] + 3, Num_DataPoints - 1)];
				}
				else
				{
					Vec3d	CurrDirec;
					if (Flag_Point[0] == 0){
						Index[i + 1] = Index[i];
						CurrDirec = (DataPoints[Index[i] + 1] - ResampleData[i]) / DistanceVec[0];
						ResampleData[i + 1] = ResampleData[i] + 1.5*CurrDirec;
					}
					else
					{
						CurrDirec = AdjacentPoints[Flag_Point[0]] - AdjacentPoints[Flag_Point[0] - 1];
						CurrDirec = CurrDirec / CurrDirec.norm();
						double	IncreaseStep = Increase_value*Increase_value - std::pow((AdjacentPoints[Flag_Point[0] - 1] - ResampleData[i]).norm(), 2);
						ResampleData[i + 1] = AdjacentPoints[Flag_Point[0] - 1] + std::sqrt(IncreaseStep)*CurrDirec;
						Index[i + 1] = Flag_Point[0] + Index[i];
					}
				}
			}
			if (AdjacentPoints.size() == 0){//2018 04 03 , ==1 cannot enter forever
				Index[i + 1] = Num_DataPoints - 1;
				ResampleData[i + 1] = AdjacentPoints[0];
			}
		}
		else
		{
			//ResampleData = ResampleData(:, 1 : i - 1);
			//ResampleData.erase(ResampleData.begin() + i+1, ResampleData.end());
			break;
		}
	}
	//ResampleData.pop_back();
	ResampleData.erase(ResampleData.begin() + i + 1, ResampleData.end());
	//ResampleData.push_back(DataPoints[Num_DataPoints - 1]);
	//ResampleData.pop_back();//2018 04 03
	//ResampleData.erase(ResampleData.end());
}

//2018/1/24//
void CorrectTrace::ThreeDimdCurveResampleModify(const VectorVec3d &DataPoints, double Increase_value, VectorVec3d &ResampleData)
{

	//VectorVec3d DataPoints_addlast = DataPoints;
	//DataPoints.push_back(DataPoints[DataPoints.size() - 1]);
	size_t Num_DataPoints = DataPoints.size();
	VectorVec3d AdjacentPoints;


	//Index=zeros(1,5*Num_DataPoints);
	std::vector<double>  DistanceVec, Flag_Point;
	size_t allpoints = std::round((DataPoints[Num_DataPoints - 1] - DataPoints[0]).norm() / Increase_value) * 2;
	if (allpoints < Num_DataPoints)allpoints = 7 * Num_DataPoints;
	std::vector<size_t> Index(allpoints);
	ResampleData.resize(allpoints);
	std::for_each(ResampleData.begin(), ResampleData.end(), [&](Vec3d& arg){ arg.setZero(); });
	ResampleData[0] = DataPoints[0];
	std::for_each(Index.begin(), Index.end(), [&](size_t& arg){ arg = -1; });
	Index[0] = 0;
	for (size_t i = 0; i < allpoints - 1; i++){
		if (Index[i] < Num_DataPoints - 1){
			AdjacentPoints.clear();
			for (size_t j = Index[i] + 1; j < std::min(Index[i] + 4, Num_DataPoints); j++){
				AdjacentPoints.push_back(DataPoints[j]);
			}
			if (AdjacentPoints.size() > 1){
				DistanceVec.clear();
				for (size_t j = 0; j < AdjacentPoints.size(); j++){
					DistanceVec.push_back((AdjacentPoints[j] - ResampleData[i]).norm());
				}
				// Flag_Point=find(DistanceVec>Increase_value);
				Flag_Point.clear();
				for (size_t k = 0; k < DistanceVec.size(); k++){
					if (DistanceVec[k] > Increase_value)Flag_Point.push_back(k);
				}
				//
				if (Flag_Point.size() == 0){
					//Index(i+1)=min(Index(i)+3,Num_DataPoints);
					Index[i + 1] = std::min(Index[i] + 3, Num_DataPoints - 1);
					ResampleData[i + 1] = DataPoints[std::min(Index[i] + 3, Num_DataPoints - 1)];
				}
				else
				{
					Vec3d	CurrDirec;
					if (Flag_Point[0] == 0){
						Index[i + 1] = Index[i];
						CurrDirec = (DataPoints[Index[i] + 1] - ResampleData[i]) / DistanceVec[0];
						ResampleData[i + 1] = ResampleData[i] + 1.5*CurrDirec;
					}
					else
					{
						CurrDirec = AdjacentPoints[Flag_Point[0]] - AdjacentPoints[Flag_Point[0] - 1];
						CurrDirec = CurrDirec / CurrDirec.norm();
						double	IncreaseStep = Increase_value*Increase_value - std::pow((AdjacentPoints[Flag_Point[0] - 1] - ResampleData[i]).norm(), 2);
						ResampleData[i + 1] = AdjacentPoints[Flag_Point[0] - 1] + std::sqrt(IncreaseStep)*CurrDirec;
						Index[i + 1] = Flag_Point[0] + Index[i];
					}
				}
			}
			if (AdjacentPoints.size() == 1){
				Index[i + 1] = Num_DataPoints - 1;
				ResampleData[i + 1] = AdjacentPoints[0];
			}
		}
		else
		{
			//ResampleData = ResampleData(:, 1 : i - 1);
			ResampleData.erase(ResampleData.begin() + i + 1, ResampleData.end());
			break;
		}
	}
	double dis = (ResampleData[ResampleData.size() - 1] - DataPoints[DataPoints.size() - 1]).norm();
	if (dis > 1.5)
	{
		VectorVec3d Mm(2);
		Mm[0] = ResampleData[ResampleData.size() - 1];
		Mm[1] = DataPoints[DataPoints.size() - 1];
		VectorVec3d MM;
		Linear_Interpolation(MM, Mm);
		ResampleData.pop_back();
		for (int i = 0; i < MM.size(); ++i)
		{
			ResampleData.push_back(MM[i]);
		}
	}
	//ResampleData.push_back(DataPoints[Num_DataPoints - 1]);
	//ResampleData.pop_back();
	//ResampleData.erase(ResampleData.end());
}

//2018/1/24//
void CorrectTrace::Linear_Interpolation(VectorVec3d &MM, const VectorVec3d &M){
	int kkk = 0;
	VectorVec3d tmpM = M;
	for (int j = 0; j < M.size() - 2; ++j)
	{
		double long_ = (M[j + 1] - M[j]).norm();
		if (long_<10)
		{
			Vec3d dire = M[j + 1] - M[j];
			dire = dire / (dire.norm() + 0.00001);
			for (double k = 0; k <= long_; k = k + 0.2)
			{
				Vec3d P = M[j] + k*dire;
				kkk = kkk + 1;
				tmpM[kkk] = P;
			}
		}
	}
	MM = tmpM;
}

//2018/1/24
void reJoineCurves(const VectorVec2d &tagMat, const CellVec3d& CellCurves, CellVec3d& CellCurvesNew){
	int Nums = CellCurves.size();
	CellCurvesNew = CellCurves;
	for (int i = 0; i < Nums; i++)
	{
		int CurveId = tagMat[i](1);
		if (CurveId > 0)
		{
			VectorVec3d Currdata = CellCurves[i];
			int FatherId = tagMat[i](0);
			std::vector<int> idx;
			for (int j = 0; j < tagMat.size(); ++j)
			{

				if (tagMat[j][0] == FatherId)
					idx.push_back(j);
			}
			int Numidx = idx.size();
			if (FatherId == 0 && Numidx == 2)//if (FatherId==1)&&(Numidx==2)
			{
				if (CurveId != idx[1])
				{
					VectorVec3d Currdata1 = CellCurves[idx[1]];
					int j;
					//CellCurvesNew[i].erase(CellCurvesNew[i].begin(), CellCurvesNew[i].end());
					for (j = 1; j <Currdata.size(); j++)
						CellCurvesNew[i][j] = Currdata[Currdata.size() - 1 - j];
					CellCurvesNew.pop_back();
					for (int jj = 0; jj < Currdata1.size(); jj++)
						CellCurvesNew[i].insert(CellCurvesNew[i].end(), Currdata1[jj]);
					//CellCurvesNew{ idx(2) } = [];
					CellCurvesNew[idx[1]].clear();
				}
				else
					continue;
			}
			if (FatherId != 0 && Numidx>2)
			{
				if (CurveId != idx[1])
				{
					VectorVec3d Currdata1 = CellCurves[idx[1]];
					int j;
					//CellCurvesNew[i].erase(CellCurvesNew[i].begin(), CellCurvesNew[i].end());
					for (j = 1; j <Currdata.size(); j++)
						CellCurvesNew[i][j] = Currdata[Currdata.size() - 1 - j];
					CellCurvesNew.pop_back();
					for (int jj = 0; jj < Currdata1.size(); jj++)
						CellCurvesNew[i].insert(CellCurvesNew[i].end(), Currdata1[jj]);
					//CellCurvesNew{ idx(2) } = [];
					CellCurvesNew[idx[1]].clear();
				}
				else
				{
					continue;
				}
			}
		}
	}
	int k = 1;
	VectorVec3d CellCurves0;
	for (int j = 0; j < Nums; j++)
	{
		VectorVec3d currD = CellCurvesNew[j];
		if (!currD.empty()){
			if (k == 1){
				CellCurves0 = currD;
				++k;
			}
			else
			{
				//CellCurves0.insert(CellCurves0.end(),currD.begin(),currD.end());
				for (size_t p = 0; p < currD.size(); p++)
				{
					CellCurves0.push_back(currD[p]);
				}
			}
		}
	}
	CellCurvesNew.erase(CellCurvesNew.begin(), CellCurvesNew.end());
	CellCurvesNew.push_back(CellCurves0);
}

void CorrectTrace::SWCFormatChangeSub1(const CellVec3d &Curvesets, MatXd &Distances1, MatXd &Distances2){
	size_t Nums = Curvesets.size();
	double ds1, ds2;
	Distances1.setOnes(Nums, Nums);
	Distances1 = Distances1 * 100;
	Distances2.setOnes(Nums, Nums);
	Distances2 = Distances2 * 100;
	for (int i = 0; i < Nums; i++){
		VectorVec3d Currdata0 = Curvesets[i];
		for (int j = 0; j < Nums; j++){
			if (i != j){
				VectorVec3d Currdata = Curvesets[j];
				Vec3d data1 = Currdata[0];
				Vec3d data2 = Currdata[Currdata.size() - 1];
				SWCFormatChangeSub11(data1, data2, Currdata0, ds1, ds2);

				Distances1.operator()(i, j) = ds1;
				Distances2.operator()(i, j) = ds2;
			}
		}
	}
	//std::cout << Distances1 << std::endl;
	//std::cout << Distances2 << std::endl;
}
void CorrectTrace::SWCFormatChangeSub11(const Vec3d &data1, const Vec3d &data2, const VectorVec3d &Currdata0, double &ds1, double &ds2)
{
	size_t Nums = Currdata0.size();
	//std::vector<double> Dis1, Dis2;//=DisMatrix
	double ds1_tmp, ds2_tmp;
	//MatXd DisMatrix = MatXd::Zero(2, Nums);
	for (int i = 0; i < Nums; i++){
		ds1_tmp = (data1 - Currdata0[i]).norm();
		ds2_tmp = (data2 - Currdata0[i]).norm();
		if (i == 0){
			ds1 = ds1_tmp;
			ds2 = ds2_tmp;
		}
		else
		{
			if (ds1_tmp < ds1){
				ds1 = ds1_tmp;
			}
			if (ds2_tmp < ds2){
				ds2 = ds2_tmp;
			}
		}
	}
	/*
	[ds1,Index1]=min(DisMatrix(1,:));
	[ds2,Index2]=min(DisMatrix(2,:));
	//*/
	//std::vector<double>::iterator result = std::min_element(std::begin(Dis1), std::end(Dis1));
	//ds1 = Dis1[std::distance(std::begin(Dis1), result)];
	//result = std::min_element(std::begin(Dis2), std::end(Dis2));
	//ds2 = Dis2[std::distance(std::begin(Dis2), result)];
}
/*find the points which are around the  connection points  */
void CorrectTrace::BifurLocalCurvesExtract(const CellVec3d& CellCurves, const VectorVec2d& ConnectSet, XCellVec3d &LocalCurve, VectorVec3d &LableVec)
{
	size_t Numss = CellCurves.size();
	Vec3d point;
	size_t Index2;
	//int kk = 0;
	std::vector<double> Index;
	for (int ii = 0; ii < Numss; ii++){
		Vec2d CurrV;
		VectorVec3d Data0;
		VectorVec3d Data1;
		CellVec3d tmp;
		Vec3d ttt;
		VectorVec3d Data11, Data00;
		CurrV = ConnectSet[ii];
		if (CurrV[0] > 0){
			Data0 = CellCurves[ii];
			Data1 = CellCurves[CurrV[0] - 1];
			//FindSpecialPointInCurve2(Data1, Data0, 3, point, Index2);
			//if (Index2>0)Data0.erase(Data0.begin(), Data0.begin() + Index2);
			size_t b = findId(Data1, Data0[0]);
			b = std::min(std::max(b, (size_t)3), (size_t)(std::max((int)Data1.size() - 5, 0)));
			Index.clear();
			Data11.clear();
			Data00.clear();
			tmp.clear();
			//bug  
			for (size_t j = std::max((int)b - 10, 0); j < std::min((size_t)b + 11, Data1.size()); j++){
				Index.push_back(j);
				Data11.push_back(Data1[j]);
			}
			Data00.push_back(Data1[b]);
			for (size_t jj = 0; jj < std::min(15, (int)Data0.size()); jj++){
				Data00.push_back(Data0[jj]);
			}
			tmp.push_back(Data00);
			tmp.push_back(Data11);
			LocalCurve.push_back(tmp);
			ttt[0] = ii + 1;
			ttt[1] = CurrV[0];
			ttt[2] = 0;
			LableVec.push_back(ttt);
		}
		if (CurrV[1] > 0){
			Data0 = CellCurves[ii];
			// Data0=Data0(:,[size(Data0,2):-1:1]);
			std::reverse(Data0.begin(), Data0.end());
			Data1 = CellCurves[CurrV[1] - 1];
			//FindSpecialPointInCurve2(Data1, Data0, 3, point, Index2);
			//if (Index2>0)Data0.erase(Data0.begin(), Data0.begin() + Index2);
			size_t b = findId(Data1, Data0[0]);
			b = std::min(std::max(b, (size_t)3), (size_t)(std::max((int)Data1.size() - 5, 0)));
			Index.clear();
			Data11.clear();
			Data00.clear();
			tmp.clear();

			for (size_t j = std::max((int)b - 10, 0); j < std::min((size_t)(b + 11), Data1.size()); ++j){
				Index.push_back(j);
				Data11.push_back(Data1[j]);
			}
			Data00.push_back(Data1[b]);
			for (size_t jj = 0; jj < std::min(15, (int)Data0.size()); ++jj){
				Data00.push_back(Data0[jj]);
			}
			tmp.push_back(Data00);
			tmp.push_back(Data11);
			LocalCurve.push_back(tmp);
			ttt[0] = ii + 1;
			ttt[1] = CurrV[1];
			ttt[2] = 1;
			LableVec.push_back(ttt);
		}
	}
}
size_t CorrectTrace::findId(VectorVec3d &Curves1, Vec3d &Curves2)
{
	std::vector<double> d;
	for (size_t i = 0; i < Curves1.size(); ++i){
		d.push_back((Curves1[i] - Curves2).norm());
	}
	auto bb = std::min_element(d.begin(), d.end());
	return bb - d.begin();
}

/*correct the connection points*/
void CorrectTrace::RevisedCurvesBifur0Modify(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new)
{
	//[Numx,Numy,Numz]=size(ImageS);
	//	Curve2 = RevisedCurvesFixed(Curve2, ImageS, zeros(1, size(Curve2, 2)), 2, 0.2);
	Vec3d PriorDirc;
	Vec3d ReferenceP;
	Vec3d DircP;
	Curve1_new = Curve1;
	Curve2_new = Curve2;
	size_t Numx = ImageS.x();
	size_t Numy = ImageS.y();
	size_t Numz = ImageS.z();
	std::vector<size_t> LabelMatrix(Curve2_new.size());
	std::for_each(LabelMatrix.begin(), LabelMatrix.end(), [&](size_t& arg){ arg = 0; });
	//time_t mytime = time(0);
	RevisedCurvesFixed(Curve2_new, ImageS, LabelMatrix, 2, 0.2);
	//time_t mytime2 = time(0);
	//std::cout << "Function: RevisedCurvesFixed time= " << mytime2 - mytime << std::endl;
	DVolume PSF3D;
	PSF3D.SetSize(7, 7, 7);
	PSF3D.SetZero();
	for (int i = -3; i <= 3; i++){
		for (int j = -3; j <= 3; j++){
			for (int k = -3; k <= 3; k++){
				PSF3D.operator()(3 + i, 3 + j, 3 + k) = std::exp(-(i *i + j *j + k *k + 0.0) / 6.0);
			}
		}
	}
	//clock_t beg = clock();
	SVolume ReImag;
	FastConv3d(PSF3D, Curve1_new, Curve2_new, Numx, Numy, Numz, ReImag);
	//ReImag.SetSize(Numx, Numy, Numz);
	//ReImag.SetZero();

	//for (int ii = 0; ii < Curve2_new.size(); ii++){
	//	Vec3d datap = Curve2_new[ii];
	//	size_t x = std::min(std::max(std::round(datap[0]), 0.0), (double)Numx);
	//	size_t y = std::min(std::max(std::round(datap[1]), 0.0), (double)Numy);
	//	size_t z = std::min(std::max(std::round(datap[2]), 0.0), (double)Numz);
	//	ReImag.operator()(x - 1, y - 1, z - 1) = 1;
	//	//ReImag.operator()(x, y , z) = 1;
	//}

	//conv3d(ReImag, PSF3D,ReImag);
	//clock_t end = clock();
	//std::cout << "Function: conv3d time= " << end - beg << std::endl;
	//std::cout << ReImag.operator()(82, 68, 45) << std::endl;
	//std::cout << ReImag.operator()(87, 73, 49) << std::endl;
	CurvesCorrforBifurp(Curve2_new, Curve1_new, Curve1_new);

	size_t Index;
	FindSpecialPointInCurve2(Curve2_new, Curve1_new, 4, ReferenceP, Index);
	//std::cout << "ReferenceP= " << ReferenceP << "Index=" << Index<<std::endl;
	if (Index < Curve1_new.size() - 4){
		DircP = Curve1_new[std::min((Curve1_new.size() - 1llu), (Index + 4llu))] - Curve1_new[Index];
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc.setZero();
	//
	//beg = clock();
	for (int i = 0; i < 4; ++i){
		Vec3d ModifiedPoint = RevisedCurveBifurTerminalPoint(Curve1_new[0], ReImag, PriorDirc, ReferenceP);
		Curve1_new[0] = ModifiedPoint;
		//std::cout << "ModifiedPoint= " << ModifiedPoint << std::endl;
		//std::cout << "ReferenceP= " << ReferenceP << std::endl;
		//std::cout << "PriorDirc= " << PriorDirc << std::endl;
		std::vector<size_t> LabelMatrix2(Curve1_new.size());
		std::for_each(LabelMatrix2.begin(), LabelMatrix2.end(), [&](size_t &arg){arg = 0; });
		RevisedCurvesFixed(Curve1_new, ImageS, LabelMatrix2, 2, 2);
		FindSpecialPointInCurve2(Curve2_new, Curve1_new, 4, ReferenceP, Index);
		if (Index < Curve1_new.size() - 4){
			DircP = Curve1_new[std::min((ulong)(Curve1_new.size() - 1), (ulong)(Index + 4))] - Curve1_new[Index];
			PriorDirc = DircP / (DircP.norm() + 0.001);
		}
		else
		{
			PriorDirc = Vec3d::Zero();
		}
	}
	//end = clock();
	//std::cout << "Function: for in RevisedCurvesBifur0Modify time= " << end - beg << std::endl;
}


//********************//
//***2018/01/24****//
void CorrectTrace::RevisedCurvesBifur0Modify(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new, Vec3d &ModifiedPoint)
{
	Vec3d PriorDirc;
	Vec3d ReferenceP;
	Vec3d DircP;
	Curve1_new = Curve1;
	Curve2_new = Curve2;
	size_t Numx = ImageS.x();
	size_t Numy = ImageS.y();
	size_t Numz = ImageS.z();
	std::vector<size_t> LabelMatrix(Curve2_new.size());
	std::for_each(LabelMatrix.begin(), LabelMatrix.end(), [&](size_t& arg){ arg = 0; });
	RevisedCurvesFixed(Curve2_new, ImageS, LabelMatrix, 2, 2);
	VectorVec3d ResampleData;
	ThreeDimdCurveResample(Curve2_new, 1, ResampleData);
	Curve2_new = ResampleData;
	//!!!!!!2/4/22:59
	Vec3d Minxx, Coordinate;
	for (size_t i = 0; i < Curve2_new.size(); i++)
	{
		if (i == 0)
			Minxx = Curve2_new[i];
		else{
			for (int j = 0; j < 3; j++)
			{
				if (Curve2_new[i](j)<Minxx(j))
					Minxx(j) = Curve2_new[i](j);
				if (i == Curve2_new.size() - 1)
					Coordinate(j) = Minxx(j) - 9;
			}
		}
	}
	for (size_t i = 0; i < Curve2_new.size(); i++)
		Curve2_new[i] = Curve2_new[i] - Coordinate;
	for (size_t i = 0; i < Curve1_new.size(); i++)
		Curve1_new[i] = Curve1_new[i] - Coordinate;

	DVolume PSF3D;
	PSF3D.SetSize(7, 7, 7);
	PSF3D.SetZero();
	for (int i = -3; i <= 3; i++){
		for (int j = -3; j <= 3; j++){
			for (int k = -3; k <= 3; k++){
				PSF3D.operator()(3 + i, 3 + j, 3 + k) = std::exp(-(i *i + j *j + k *k + 0.0) / 1.0);
			}
		}
	}

	SVolume ReImag;
	FastConv3d(PSF3D, Curve1_new, Curve2_new, Numx, Numy, Numz, ReImag);
	CurvesCorrforBifurp(Curve2_new, Curve1_new, Curve1_new);

	size_t Index;
	FindSpecialPointInCurve2(Curve2_new, Curve1_new, 4, ReferenceP, Index);
	//std::cout << "ReferenceP= " << ReferenceP << "Index=" << Index<<std::endl;
	if (Index < Curve1_new.size() - 4){
		DircP = Curve1_new[std::min((Curve1_new.size() - 1llu), (Index + 4llu))] - Curve1_new[Index];
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc.setZero();
	//
	//beg = clock();
	for (int i = 0; i < 1; ++i){
		//Vec3d one; one.setOnes();
		//Vec3d tmp = Curve1_new[0] -one;
		ModifiedPoint = RevisedCurveBifurTerminalPoint(Curve1_new[0], ReImag, PriorDirc, ReferenceP);
		Curve1_new[0] = ModifiedPoint;

		for (size_t i = 0; i < Curve1_new.size(); i++)
			Curve1_new[i] = Curve1_new[i] + Coordinate;
		std::vector<size_t> LabelMatrix2(Curve1_new.size());
		std::for_each(LabelMatrix2.begin(), LabelMatrix2.end(), [&](size_t &arg){arg = 0; });
		RevisedCurvesFixed(Curve1_new, ImageS, LabelMatrix2, 2, 2);
	}
	for (size_t i = 0; i < Curve2_new.size(); i++)
		Curve2_new[i] = Curve2_new[i] + Coordinate;
	ModifiedPoint = ModifiedPoint + Coordinate;
}

//********************//
//***2018/01/18****//
void CorrectTrace::RevisedCurvesBifur0Modify_3line(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS, VectorVec3d &Curve2_new, VectorVec3d &Curve1_new)
{
	//[Numx,Numy,Numz]=size(ImageS);
	//	Curve2 = RevisedCurvesFixed(Curve2, ImageS, zeros(1, size(Curve2, 2)), 2, 0.2);
	Vec3d PriorDirc;
	Vec3d ReferenceP;
	Vec3d DircP;
	Curve1_new = Curve1;
	Curve2_new = Curve2;
	size_t Numx = ImageS.x();
	size_t Numy = ImageS.y();
	size_t Numz = ImageS.z();
	std::vector<size_t> LabelMatrix(Curve2_new.size());
	std::for_each(LabelMatrix.begin(), LabelMatrix.end(), [&](size_t& arg){ arg = 0; });
	Vec3d Minxx, Coordinate;
	for (size_t i = 0; i < Curve2_new.size(); i++)
	{
		if (i == 0)
			Minxx = Curve2_new[i];
		else{
			for (int j = 0; j < 3; j++)
			{
				if (Curve2_new[i](j)<Minxx(j))
					Minxx(j) = Curve2_new[i](j);
				if (i == Curve2_new.size() - 1)
					Coordinate(j) = Minxx(j) - 10;
			}
		}
	}
	for (size_t i = 0; i < Curve2_new.size(); i++)
		Curve2_new[i] = Curve2_new[i] - Coordinate;
	for (size_t i = 0; i < Curve1_new.size(); i++)
		Curve1_new[i] = Curve1_new[i] - Coordinate;

	//time_t mytime = time(0);
	RevisedCurvesFixed(Curve2_new, ImageS, LabelMatrix, 2, 0.2);
	//time_t mytime2 = time(0);
	//std::cout << "Function: RevisedCurvesFixed time= " << mytime2 - mytime << std::endl;
	DVolume PSF3D;
	PSF3D.SetSize(7, 7, 7);
	PSF3D.SetZero();
	for (int i = -3; i <= 3; i++){
		for (int j = -3; j <= 3; j++){
			for (int k = -3; k <= 3; k++){
				PSF3D.operator()(3 + i, 3 + j, 3 + k) = std::exp(-(i *i + j *j + k *k + 0.0) / 6.0);
			}
		}
	}
	//clock_t beg = clock();
	SVolume ReImag;
	FastConv3d(PSF3D, Curve1_new, Curve2_new, Numx, Numy, Numz, ReImag);
	//clock_t end = clock();
	VectorVec3d Curve1_new11, Curve1_new10, CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2;
	CurvesCorrforBifurp(Curve2_new, Curve1_new, Curve1_new, Curve1_new11, Curve1_new10);
	//[CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2] = ...
	//distance_check(Curve1, Curve11, Curve10);
	distance_check(Curve1_new, Curve1_new11, Curve1_new10, CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2);
	VectorVec3d ModifiedPoint;
	//main line   ModifiedPoint1
	size_t Index;
	FindSpecialPointInCurve2(Curve2_new, Curve1_new, 4, ReferenceP, Index);
	//std::cout << "ReferenceP= " << ReferenceP << "Index=" << Index<<std::endl;
	if (Index < Curve1_new.size() - 4){
		DircP = Curve1_new[std::min((Curve1_new.size() - 1llu), (Index + 4llu))] - Curve1_new[Index];
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc.setZero();
	ModifiedPoint.push_back(RevisedCurveBifurTerminalPoint(Curve1_new[0], ReImag, PriorDirc, ReferenceP));

	//mian line ModifiedPoint2
	Vec3d ReferenceP2;
	size_t Index2;
	FindSpecialPointInCurve2(CurveBase_bak1, CurveBir_bak1, 1.5, ReferenceP2, Index2);
	if (Index2 < CurveBir_bak1.size() - 4)//matlab___Index<size(CurveBir_bak1,2)-3
	{
		DircP = CurveBir_bak1[std::min((ulong)(CurveBir_bak1.size() - 1), (ulong)(Index2 + 2))] - CurveBir_bak1[Index2];
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc = Vec3d::Zero();
	ModifiedPoint.push_back(RevisedCurveBifurTerminalPoint(CurveBir_bak1[0], ReImag, PriorDirc, ReferenceP2));

	//mian line ModifiedPoint3
	Vec3d ReferenceP3;
	size_t Index3;
	FindSpecialPointInCurve2(CurveBase_bak2, CurveBir_bak2, 1.5, ReferenceP3, Index3);
	if (Index3 < CurveBase_bak2.size() - 4){
		DircP = CurveBase_bak2[std::min((ulong)(CurveBase_bak2.size() - 1), (ulong)(Index3 + 2))] - CurveBase_bak2[Index3];
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc = Vec3d::Zero();
	ModifiedPoint.push_back(RevisedCurveBifurTerminalPoint(CurveBase_bak2[0], ReImag, PriorDirc, ReferenceP3));

	//regress mainline  , determine whether three points are close
	VectorVec3d xxx = ModifiedPoint;
	std::vector<double> pValue;
	weigthvalue(xxx, ReImag, pValue);
	double fenmu = 0;
	for (int i = 0; i < pValue.size(); i++)
	{
		fenmu += pValue[i];
	}
	Vec3d ModifiedPoint00 = (pValue[0] * xxx[0] + pValue[1] * xxx[1] + pValue[2] * xxx[2]) / fenmu;
	Curve1_new[0] = ModifiedPoint00;

	for (size_t i = 0; i < Curve2_new.size(); i++)
		Curve2_new[i] = Curve2_new[i] + Coordinate;
	for (size_t i = 0; i < Curve1_new.size(); i++)
		Curve1_new[i] = Curve1_new[i] + Coordinate;
	//beg = clock();
	//end = clock();
	//std::cout << "Function: for in RevisedCurvesBifur0Modify time= " << end - beg << std::endl;
}

//************************//
//*******2018/1/24******//
void CorrectTrace::RevisedCurvesBifur0Modify_3line(const VectorVec3d &Curve2, const VectorVec3d &Curve1, SVolume &ImageS,
	VectorVec3d &Curve2_new, VectorVec3d &Curve1_new, Vec3d &ModifiedPoint1,
	VectorVec3d &CurveBase_bak1, VectorVec3d &CurveBir_bak1, VectorVec3d &CurveBase_bak2, VectorVec3d &CurveBir_bak2)
{
	Vec3d PriorDirc;
	Vec3d ReferenceP;
	Vec3d DircP;
	Curve1_new = Curve1;
	Curve2_new = Curve2;
	size_t Numx = ImageS.x();
	size_t Numy = ImageS.y();
	size_t Numz = ImageS.z();
	std::vector<size_t> LabelMatrix(Curve2_new.size());
	std::for_each(LabelMatrix.begin(), LabelMatrix.end(), [&](size_t& arg){ arg = 0; });
	RevisedCurvesFixed(Curve2_new, ImageS, LabelMatrix, 1, 0.2);//
	VectorVec3d ResampleData;
	ThreeDimdCurveResample(Curve2_new, 1, ResampleData);
	Curve2_new = ResampleData;
	Vec3d Minxx, Coordinate;
	for (size_t i = 0; i < Curve2_new.size(); ++i)
	{
		if (i == 0)
			Minxx = Curve2_new[i];
		else{
			for (int j = 0; j < Curve2_new[i].size(); ++j)
			{
				if (Curve2_new[i](j)<Minxx(j))
					Minxx(j) = Curve2_new[i](j);
			}
		}
	}
	Vec3d ten(9, 9, 9);
	Coordinate = Minxx - ten;
	for (size_t i = 0; i < Curve2_new.size(); i++)
		Curve2_new[i] = Curve2_new[i] - Coordinate;
	for (size_t i = 0; i < Curve1_new.size(); i++)
		Curve1_new[i] = Curve1_new[i] - Coordinate;

	DVolume PSF3D;
	PSF3D.SetSize(7, 7, 7);
	PSF3D.SetZero();
	for (int i = -3; i <= 3; i++){
		for (int j = -3; j <= 3; j++){
			for (int k = -3; k <= 3; k++){
				PSF3D.operator()(3 + i, 3 + j, 3 + k) = std::exp(-(i *i + j *j + k *k + 0.0) / 1.0);
			}
		}
	}

	SVolume ReImag;
	FastConv3d(PSF3D, Curve1_new, Curve2_new, Numx, Numy, Numz, ReImag);

	VectorVec3d Curve1_new11, Curve1_new10;
	CurvesCorrforBifurp(Curve2_new, Curve1_new, Curve1_new, Curve1_new11, Curve1_new10);

	distance_check(Curve1_new, Curve1_new11, Curve1_new10, CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2);

	//main line   ModifiedPoint1
	size_t Index;
	FindSpecialPointInCurve2(Curve2_new, Curve1_new, 5, ReferenceP, Index);
	//std::cout << "ReferenceP= " << ReferenceP << "Index=" << Index<<std::endl;
	if (Index < Curve1_new.size() - 4){
		DircP = Curve1_new[std::min((Curve1_new.size() - 1llu), (Index + 5llu))] - Curve1_new[Index];//!!!!!!!!!
		PriorDirc = DircP / (DircP.norm() + 0.001);
	}
	else
		PriorDirc.setZero();
	ModifiedPoint1 = RevisedCurveBifurTerminalPoint(Curve1_new[0], ReImag, PriorDirc, ReferenceP);

	//regress mainline  , determine whether three points are close
	Curve1_new[0] = ModifiedPoint1;

	for (size_t i = 0; i < Curve2_new.size(); i++)    Curve2_new[i] = Curve2_new[i] + Coordinate;
	for (size_t i = 0; i < Curve1_new.size(); i++)		Curve1_new[i] = Curve1_new[i] + Coordinate;
	std::vector<size_t> LabelMatrix1(Curve1_new.size());
	std::for_each(LabelMatrix1.begin(), LabelMatrix1.end(), [&](size_t& arg){ arg = 0; });
	RevisedCurvesFixed(Curve1_new, ImageS, LabelMatrix1, 2, 2);
	for (size_t i = 0; i < CurveBase_bak1.size(); i++)		CurveBase_bak1[i] = CurveBase_bak1[i] + Coordinate;
	for (size_t i = 0; i < CurveBir_bak1.size(); i++)		CurveBir_bak1[i] = CurveBir_bak1[i] + Coordinate;
	for (size_t i = 0; i < CurveBase_bak2.size(); i++)		CurveBase_bak2[i] = CurveBase_bak2[i] + Coordinate;
	for (size_t i = 0; i < CurveBir_bak2.size(); i++)		CurveBir_bak2[i] = CurveBir_bak2[i] + Coordinate;
	ModifiedPoint1 = ModifiedPoint1 + Coordinate;
}

void CorrectTrace::RevisedCurvesFixed(VectorVec3d & SingCurveSet, SVolume &ImageS, std::vector<size_t> &LabelMatrix, double ThreShift, double ThreL1)
{
	VectorVec3d CurrCurve = SingCurveSet;
	size_t Numss = CurrCurve.size();
	VectorVec3d NumssMatrix(Numss - 2);
	VectorVec3d NumssMatrix0(Numss - 2);
	VectorVec3d CurrCurve1, d, r, d_new, r_new;
	double ThreShift0 = 0;
	std::for_each(NumssMatrix.begin(), NumssMatrix.end(), [&](Vec3d &arg){arg = Vec3d::Ones()*ThreL1; });

	for (int loop = 1; loop <= 40; ++loop){//20

		if (loop<30){   //if (loop %105 == 0){//10
			ThreShift0 = 1.5 * ThreShift;

			std::for_each(NumssMatrix0.begin(), NumssMatrix0.end(), [&](Vec3d &arg){arg = Vec3d::Ones()*ThreL1*1.5; });
		}
		else
		{
			ThreShift0 = ThreShift;
			NumssMatrix0 = NumssMatrix;
		}
		if (loop > 1){
			//CurvesCorrSubMeanShiftModify000(ImageS, CurrCurve, LabelMatrix, ThreShift0, d, r, CurrCurve1);
			CurvesCorrSubMeanShiftModify0001(ImageS, CurrCurve, LabelMatrix, ThreShift0, d, r, CurrCurve1);
			CurvesCorrSparse(CurrCurve1, r, NumssMatrix0, d_new, r_new);
		}
		else
		{
			InitialDR_matrix(CurrCurve, d, r);
			CurvesCorrSubMeanShiftModify0001(ImageS, CurrCurve, LabelMatrix, ThreShift0, d, r, CurrCurve1);
			CurvesCorrSparse(CurrCurve1, r, NumssMatrix0, d_new, r_new);
		}
		CurrCurve = CurrCurve1;
		d = d_new;
		r = r_new;
	}
	SingCurveSet = CurrCurve;
}
//2018/1/24
void CorrectTrace::CurvesCorrSubMeanShiftModify0001(const SVolume &ImageS, const VectorVec3d &dataPoints,
	const std::vector<size_t> &LabelMatrix, double Lamda, VectorVec3d &d, VectorVec3d &r, VectorVec3d& dataPoints1)
{
	size_t NumP = dataPoints.size();
	dataPoints1 = dataPoints;
	VectorVec3d dataPoints0 = dataPoints1;
	for (int ttt = 0; ttt < 1; ttt++){
		//clock_t beg = clock();
		for (int ii = 0; ii < 3; ++ii){
			for (int i = 1; i < NumP - 1; ++i){
				if (LabelMatrix[i] == 0){
					Vec3d Dircxx = dataPoints1[i] - dataPoints1[i - 1];
					Vec3d Dircyy = dataPoints1[i + 1] - dataPoints1[i];
					Vec3d NDircx = Dircxx / std::max(Dircxx.norm(), 0.001);
					Vec3d NDircy = Dircyy / std::max(Dircyy.norm(), 0.001);
					Vec3d Currdatap;
					if (NDircx.transpose()*NDircy>0.7){
						Vec3d Dircx = NDircx + NDircy;
						Vec3d Dircy, Dircz;
						Dircx = Dircx / Dircx.norm();
						directionc(Dircx, Dircy, Dircz);
						Mat3d SigmaM = (0.1*Dircx*Dircx.transpose() + Dircy*Dircy.transpose() + Dircz*Dircz.transpose()).inverse();
						//Mat3d SigmaM = (  Dircy*Dircy.transpose() + Dircz*Dircz.transpose()).inverse();
						Currdatap = MeanShiftSinglePointModify(dataPoints1[i], SigmaM, ImageS);
					}
					else
					{
						Currdatap = MeanShiftSinglePointModify(dataPoints1[i], Mat3d::Identity(3, 3), ImageS);
					}
					dataPoints0[i] = Currdatap;
				}
			}
		}
		//clock_t end = clock();
		//std::cout << "Function: for in CurvesCorrSubMeanShiftModify000 time= " << end - beg << std::endl;
		VectorVec3d mm;
		Vec3d CurrData;
		for (int i = 0; i < d.size(); i++){
			mm.push_back(d[i] - r[i]);
		}
		for (int ii = 0; ii < 3; ii++){
			CurrData = 2 * dataPoints1[0] + 4 * dataPoints1[2] - dataPoints1[3] + 2 * mm[0] - mm[1];
			CurrData = Lamda*CurrData + dataPoints0[1];
			dataPoints1[1] = CurrData / (1 + 5 * Lamda);
			for (int i = 2; i < NumP - 2; i++){
				if (LabelMatrix[i] == 0){
					CurrData = 4 * dataPoints1[i - 1] + 4 * dataPoints1[i + 1] - dataPoints1[i - 2] - dataPoints1[i + 2];
					//cout << CurrData << endl;
					CurrData = Lamda*(CurrData + 2 * mm[i - 1] - mm[i - 2] - mm[i]) + dataPoints0[i];
					//cout << CurrData << endl;
					dataPoints1[i] = CurrData / (1 + 6 * Lamda);
					//cout << dataPoints1[i] << endl;

				}
			}
			CurrData = 4 * dataPoints1[NumP - 3] + 2 * dataPoints1[NumP - 1] - dataPoints1[NumP - 4] + 2 * mm[NumP - 3] - mm[NumP - 4];
			CurrData = Lamda*CurrData + dataPoints0[NumP - 2];
			dataPoints1[NumP - 2] = CurrData / (1 + 5 * Lamda);
		}
	}
}
void CorrectTrace::CurvesCorrSubMeanShiftModify000(const SVolume &ImageS, const VectorVec3d &dataPoints,
	const std::vector<size_t> &LabelMatrix, double Lamda, VectorVec3d &d, VectorVec3d &r, VectorVec3d& dataPoints1)
{
	size_t NumP = dataPoints.size();
	dataPoints1 = dataPoints;
	VectorVec3d dataPoints0 = dataPoints1;
	for (int ttt = 0; ttt < 1; ttt++){
		//clock_t beg = clock();
		for (int ii = 0; ii < 3; ++ii){
			for (int i = 1; i < NumP - 1; ++i){
				if (LabelMatrix[i] == 0){
					Vec3d Dircxx = dataPoints1[i] - dataPoints1[i - 1];
					Vec3d Dircyy = dataPoints1[i + 1] - dataPoints1[i];
					Vec3d NDircx = Dircxx / std::max(Dircxx.norm(), 0.001);
					Vec3d NDircy = Dircyy / std::max(Dircyy.norm(), 0.001);
					Vec3d Currdatap;
					if (NDircx.transpose()*NDircy>0.7){
						Vec3d Dircx = NDircx + NDircy;
						Vec3d Dircy, Dircz;
						Dircx = Dircx / Dircx.norm();
						directionc(Dircx, Dircy, Dircz);
						Mat3d SigmaM = (100*Dircx*Dircx.transpose() + Dircy*Dircy.transpose() + Dircz*Dircz.transpose()).inverse();
						//Mat3d SigmaM = ( Dircy*Dircy.transpose() + Dircz*Dircz.transpose()).inverse();
						Currdatap = MeanShiftSinglePointModify(dataPoints1[i], SigmaM, ImageS);
					}
					else
					{
						Currdatap = MeanShiftSinglePointModify(dataPoints1[i], Mat3d::Identity(3, 3), ImageS);
					}
					dataPoints0[i] = Currdatap;
				}
			}
		}
		//clock_t end = clock();
		//std::cout << "Function: for in CurvesCorrSubMeanShiftModify000 time= " << end - beg << std::endl;
		VectorVec3d mm;
		Vec3d CurrData;
		for (int i = 0; i < d.size(); i++){
			mm.push_back(d[i] - r[i]);
		}
		for (int ii = 0; ii < 3; ii++){
			CurrData = 2 * dataPoints1[0] + 4 * dataPoints1[2] - dataPoints1[3] + 2 * mm[0] - mm[1];
			CurrData = Lamda*CurrData + dataPoints0[1];
			dataPoints1[1] = CurrData / (1 + 5 * Lamda);
			for (int i = 2; i < NumP - 2; i++){
				if (LabelMatrix[i] == 0){
					CurrData = 4 * dataPoints1[i - 1] + 4 * dataPoints1[i + 1] - dataPoints1[i - 2] - dataPoints1[i + 2];
					//cout << CurrData << endl;
					CurrData = Lamda*(CurrData + 2 * mm[i - 1] - mm[i - 2] - mm[i]) + dataPoints0[i];
					//cout << CurrData << endl;
					dataPoints1[i] = CurrData / (1 + 6 * Lamda);
					//cout << dataPoints1[i] << endl;

				}
			}
			CurrData = 4 * dataPoints1[NumP - 3] + 2 * dataPoints1[NumP - 1] - dataPoints1[NumP - 4] + 2 * mm[NumP - 3] - mm[NumP - 4];
			CurrData = Lamda*CurrData + dataPoints0[NumP - 2];
			dataPoints1[NumP - 2] = CurrData / (1 + 5 * Lamda);
		}
	}
}
void CorrectTrace::directionc(const Vec3d &x1, Vec3d &x2, Vec3d& x3){

	//Vec3d xx = { 1, 0, 0 };
	//Mat3d xt = x1*xx;
	//x2 = x1*xt;
	//x3 = x1*x2;
	x2 = x1;
	//[idxv,idexx]=sort(abs(x1));
	double temp;
	int n = 3;
	int i, j;
	Vec3d idxv = x1.cwiseAbs();//x1 de mo
	Vec3d idexx = { 0, 1, 2 };
	for (j = 0; j <n - 1; j++)
		for (i = 0; i < n - 1 - j; i++)
		{
			if (idxv[i] > idxv[i + 1])
			{
				temp = idxv[i];
				idxv[i] = idxv[i + 1];
				idxv[i + 1] = temp;
				temp = idexx[i];
				idexx[i] = idexx[i + 1];
				idexx[i + 1] = temp;
			}
		}
	double cdd = (idxv[0] * idxv[0] + idxv[1] * idxv[1]) / idxv[2];
	x2[idexx[2]] = -Sign(x1[idexx[2]])*cdd;
	if (cdd != 0){
		x2 = x2 / x2.norm();
	}
	else
	{
		x2(idexx(1)) = 1;
	}
	x3.setZero();
	x3(idexx(0)) = -1;
	MatXd AA = MatXd::Zero(2, 3);
	AA.row(0) = x1.transpose();
	AA.row(1) = x2.transpose();
	MatXd AA1 = MatXd::Zero(2, 2);
	MatXd AA2 = MatXd::Zero(2, 2);
	MatXd AA3 = MatXd::Zero(2, 2);
	AA1.col(0) = AA.col(idexx[0]);
	AA1.col(1) = AA.col(idexx[2]);
	AA2.col(0) = AA.col(idexx[1]);
	AA2.col(1) = AA.col(idexx[2]);
	AA3.col(0) = AA.col(idexx[1]);
	AA3.col(1) = AA.col(idexx[0]);
	//cout << AA1 << endl;
	//cout << AA2 << endl;
	//cout << AA3 << endl;
	//cout << AA1.determinant() << endl;
	//cout << AA2.determinant() << endl;
	//cout << AA3.determinant() << endl;
	x3(idexx(1)) = AA1.determinant() / AA2.determinant();

	x3(idexx(2)) = AA3.determinant() / AA2.determinant();
	x3 = x3 / x3.norm();
}

//SVolume

void CorrectTrace::FastConv3d(const DVolume &fliter, const VectorVec3d &Curve1, const VectorVec3d &Curve2, size_t Numx, size_t Numy, size_t Numz, SVolume &src){
	//VectorVec3d Curve2, Curve1;
	Vec3d cordxyz;
	Vec3d MaxCurve2, MaxCurve1;
	for (int i = 0; i < Curve2.size(); i++){
		if (i == 0)       MaxCurve2 = Curve2[i];
		for (int j = 0; j < 3; j++){
			if (Curve2[i](j)>MaxCurve2(j)) MaxCurve2(j) = Curve2[i](j);
		}
	}
	for (int i = 0; i < Curve1.size(); i++){
		if (i == 0)       MaxCurve1 = Curve1[i];
		for (int j = 0; j < 3; j++){
			if (Curve1[i](j)>MaxCurve1(j)) MaxCurve1(j) = Curve1[i](j);
		}
	}
	for (int i = 0; i < 3; i++)
	{
		cordxyz(i) = std::max(MaxCurve2(i), MaxCurve1(i));
	}
	double Rex = std::round(cordxyz(0) + 10);
	double Rey = std::round(cordxyz(1) + 10);
	double Rez = std::round(cordxyz(2) + 10);
	//double Rex = std::round(std::max(cordxyz()))

	DVolume ReImag;
	VectorVec3d Curve22(Curve2.size() * 2 - 1);
	for (int i = 0; i < Curve22.size(); ++i){
		if (i % 2 == 0){
			Curve22[i] = Curve2[i / 2];
		}
		else
		{
			Curve22[i] = (Curve2[i / 2] + Curve2[i / 2 + 1])*0.5;
		}

	}
	size_t fx = fliter.x();

	int radius = (fx - 1) / 2;
	ReImag.SetSize(Rex, Rey, Rez);
	//ReImag.SetSize(Numx, Numy, Numz);
	ReImag.SetZero();
	src.SetSize(Rex, Rey, Rez);
	//src.SetSize(Numx, Numy, Numz);
	src.SetZero();
	for each(auto &datap in Curve22){
		size_t x = std::min(std::max(std::round(datap[0]), 0.0), (double)Numx - 1);
		size_t y = std::min(std::max(std::round(datap[1]), 0.0), (double)Numy - 1);
		size_t z = std::min(std::max(std::round(datap[2]), 0.0), (double)Numz - 1);
		//ReImag.operator()(x - 1, y - 1, z - 1) = 1;
		if (src.operator()(x, y, z) == 0){
			AddFliter(ReImag, fliter, x, y, z, radius);
			src.operator()(x, y, z) = 1;
		}
	}
	//for each(auto datap in Curve){
	//	size_t x = std::min(std::max(std::round(datap[0]), 0.0), (double)Numx-1);
	//	size_t y = std::min(std::max(std::round(datap[1]), 0.0), (double)Numy-1);
	//	size_t z = std::min(std::max(std::round(datap[2]), 0.0), (double)Numz-1);
	//	//ReImag.operator()(x - 1, y - 1, z - 1) = 1;
	//	if (src.operator()(x, y, z) == 1){
	//		/*AddFliter(ReImag, fliter, x, y, z, radius);*/
	//		src.operator()(x, y, z) = std::round(ReImag.operator()(x, y, z));
	//	}
	//}
	for (size_t i = 0; i < Rex; ++i){
		for (size_t j = 0; j < Rey; ++j){
			for (size_t k = 0; k < Rez; ++k){
				if (ReImag.operator()(i, j, k) < 0.001)continue;
				src.operator()(i, j, k) = std::round(ReImag.operator()(i, j, k));
			}
		}
	}

}
void  CorrectTrace::AddFliter(DVolume & src, const DVolume &fliter, size_t x, size_t y, size_t z, int radius){
	size_t nx = src.x();
	size_t ny = src.y();
	size_t nz = src.z();
	//size_t fx = fliter.x();
	//size_t fy = fliter.y();
	//size_t fz = fliter.z();
	size_t now_x = 0, now_y = 0, now_z = 0;
	for (int i = -radius; i <= radius; i++){
		for (int j = -radius; j <= radius; j++){
			for (int k = -radius; k <= radius; k++){
				now_x = i + x;
				now_y = j + y;
				now_z = k + z;
				if (now_x < 0 || now_y < 0 || now_z < 0)continue;
				if (now_x >= nx || now_y >= ny || now_z >= nz)continue;
				//if (src.operator()(now_x, now_y, now_z) == 0)continue;
				src.operator()(now_x, now_y, now_z) += fliter.operator()(i + radius, j + radius, k + radius);
			}
		}
	}
}
void CorrectTrace::conv3d(const SVolume &src, DVolume &fliter, SVolume &output){
	size_t nx = src.x();
	size_t ny = src.y();
	size_t nz = src.z();
	size_t fx = fliter.x();
	size_t fy = fliter.y();
	size_t fz = fliter.z();
	SVolume mysrc;
	mysrc.QuickCopy(src);
	int radius = 0;
	if (fx == fy&&fy == fz)
	{
		radius = (fx - 1) / 2;
	}
	else
	{
		//error
	}

	output.SetSize(nx, ny, nz);
	output.SetZero();
	for (size_t x = 0; x < nx; x++){
		for (size_t y = 0; y < ny; y++){
			for (size_t z = 0; z < nz; z++){
				output.operator()(x, y, z) = std::round(Getsum(mysrc, fliter, x, y, z, radius));//has
				/*	if (output.operator()(x, y, z) != 0){
				cout << output.operator()(x, y, z) << endl;
				}
				if (mysrc.operator()(x, y, z) != 0){
				cout << mysrc.operator()(x, y, z) << endl;
				}*/
			}
		}
	}
}

//replace
double CorrectTrace::Getsum(SVolume & src, DVolume &fliter, size_t x, size_t y, size_t z, int radius){
	size_t nx = src.x();
	size_t ny = src.y();
	size_t nz = src.z();
	//size_t fx = fliter.x();
	//size_t fy = fliter.y();
	//size_t fz = fliter.z();
	size_t now_x = 0, now_y = 0, now_z = 0;
	double mysum = 0;
	for (int i = -radius; i <= radius; i++){
		for (int j = -radius; j <= radius; j++){
			for (int k = -radius; k <= radius; k++){
				now_x = i + x;
				now_y = j + y;
				now_z = k + z;
				if (now_x<0 || now_y<0 || now_z<0)continue;
				if (now_x >= nx || now_y >= ny || now_z >= nz)continue;
				//if (src.operator()(now_x, now_y, now_z) == 0)continue;
				mysum += (double)src.operator()(now_x, now_y, now_z)*fliter.operator()(i + radius, j + radius, k + radius);
			}
		}
	}
	return mysum;
}
//********************//
//***2018/01/18****//
void CorrectTrace::CurvesCorrforBifurp(const VectorVec3d &Curve1, const VectorVec3d &Curve2, VectorVec3d& Curve3, VectorVec3d& Curve11, VectorVec3d& Curve10){
	Curve3 = Curve2;
	Vec3d datap0 = Curve3[0];
	Vec3d datap1 = Curve3[Curve3.size() - 1];
	size_t Num1 = Curve1.size();
	VectorVec2d MinDistMatrix(Num1);
	double min0, min1;
	size_t Index0, Index1;
	std::for_each(MinDistMatrix.begin(), MinDistMatrix.end(), [&](Vec2d&arg){arg = Vec2d::Zero(); });
	for (size_t i = 0; i < Num1; i++){
		double dst0 = (datap0 - Curve1[i]).norm();
		double dst1 = (datap1 - Curve1[i]).norm();
		MinDistMatrix[i][0] = dst0;
		MinDistMatrix[i][1] = dst1;
		if (i == 0){
			min0 = dst0;
			min1 = dst1;
			Index0 = i;
			Index1 = i;
		}
		else
		{
			if (min0 > dst0){
				min0 = dst0;
				Index0 = i;
			}
			if (min1 > dst1){
				min1 = dst1;
				Index1 = i;
			}
		}
	}
	if (min0 < min1){
		size_t CurrIndex = std::min((size_t)std::max(Index0, (size_t)2), Num1 - 1);//CurrIndex=min(max(Index0,3),Num1);
		for (size_t i = 0; i < CurrIndex + 1; ++i)
			Curve10.push_back(Curve1[i]);
		for (size_t i = 0; i <= Num1 - CurrIndex - 1; ++i)
			Curve11.push_back(Curve1[CurrIndex + i]);
		if (min0 > 0){
			Curve3.insert(Curve3.begin(), Curve11[0]);
		}
	}
	else
	{
		size_t CurrIndex = std::min((size_t)std::max(Index1, (size_t)2), Num1 - 1);
		for (size_t i = 0; i < CurrIndex + 1; ++i)
			Curve10.push_back(Curve1[i]);
		for (size_t i = 0; i <= Num1 - CurrIndex - 1; ++i)
			Curve11.push_back(Curve1[CurrIndex + i]);
		if (min1 > 0){
			Curve3.push_back(Curve1[CurrIndex]);
		}
		Vec3d tmp;
		size_t n = Curve3.size();
		for (int i = 0; i < n / 2; i++){
			tmp = Curve3[i];
			Curve3[i] = Curve3[n - i - 1];
			Curve3[n - i - 1] = tmp;
		}
	}
}

void CorrectTrace::CurvesCorrforBifurp(const VectorVec3d &Curve1, const VectorVec3d &Curve2, VectorVec3d& Curve3){
	Curve3 = Curve2;
	Vec3d datap0 = Curve3[0];
	Vec3d datap1 = Curve3[Curve3.size() - 1];
	size_t Num1 = Curve1.size();
	VectorVec2d MinDistMatrix(Num1);
	double min0, min1;
	size_t Index0, Index1;
	std::for_each(MinDistMatrix.begin(), MinDistMatrix.end(), [&](Vec2d&arg){arg = Vec2d::Zero(); });
	for (size_t i = 0; i < Num1; i++){
		double dst0 = (datap0 - Curve1[i]).norm();
		double dst1 = (datap1 - Curve1[i]).norm();
		MinDistMatrix[i][0] = dst0;
		MinDistMatrix[i][1] = dst1;
		if (i == 0){
			min0 = dst0;
			min1 = dst1;
			Index0 = i;
			Index1 = i;
		}
		else
		{
			if (min0 > dst0){
				min0 = dst0;
				Index0 = i;
			}
			if (min1 > dst1){
				min1 = dst1;
				Index1 = i;
			}
		}
	}
	if (min0 < min1){
		size_t CurrIndex = std::min((size_t)std::max(Index0, (size_t)3), Num1);
		if (min0 > 0){
			Curve3.insert(Curve3.begin(), Curve1[CurrIndex]);
		}
		else
		{
			size_t CurrIndex = std::min((size_t)std::max(Index1, (size_t)3), Num1);
			if (min1 > 0){
				Curve3.push_back(Curve1[CurrIndex]);
			}
			Vec3d tmp;
			size_t n = Curve3.size();
			for (int i = 0; i < n / 2; i++){
				tmp = Curve3[i];
				Curve3[i] = Curve3[n - i - 1];
				Curve3[n - i - 1] = tmp;
			}
		}
	}
}

void CorrectTrace::FindSpecialPointInCurve2(const VectorVec3d &Curve2, const VectorVec3d &Curve1, double Thre, Vec3d& Point, size_t &Index){
	size_t Num1 = Curve1.size();
	size_t Num2 = Curve2.size();
	size_t i;
	for (i = 0; i < Num1; i++){
		Vec3d dataPoint = Curve1[i];
		VectorVec3d PMatrix(Num2);
		double Minv;
		for (int jj = 0; jj < Num2; jj++){
			PMatrix[jj] = dataPoint - Curve2[jj];
			if (jj == 0){
				Minv = std::sqrt(PMatrix[jj].transpose()*PMatrix[jj]);
			}
			else
			{
				if (Minv > std::sqrt(PMatrix[jj].transpose()*PMatrix[jj])){
					Minv = std::sqrt(PMatrix[jj].transpose()*PMatrix[jj]);
				}
			}
		}
		if (Minv > Thre){
			Index = i;
			Point = Curve1[i];
			break;
		}
	}
	if (i == Num1){
		Index = i - 1;
		Point = Curve1[Index];
	}
}

Vec3d CorrectTrace::RevisedCurveBifurTerminalPoint(Vec3d SinglePoint, SVolume &ReImag, Vec3d &PriorDirc, Vec3d &ReferenceP){
	Vec3d x2, x3;
	Vec3d	Currdatap;
	for (int i = 0; i < 20; i++){
		Currdatap = MeanShiftSinglePointModify(SinglePoint, Mat3d::Identity(3, 3), ReImag);
		//Currdatap = MeanShiftSinglePointModify(Currdatap, Mat3d::Identity(3, 3), ReImag);
		if (std::abs(PriorDirc.maxCoeff())>0.001){
			//[x2,x3]=directionc(PriorDirc);
			PriorDirc = PriorDirc / PriorDirc.norm();
			directionc(PriorDirc, x2, x3);
		}
		else
		{
			x2 = Vec3d::Zero();
			x3 = Vec3d::Zero();
		}
		Vec3d    P = ReferenceP - SinglePoint;

		//P = P / std::max(P.norm(), 0.01);
		Vec3d ReferenceP1 = ReferenceP + 2 * P;
		Mat3d SigmaMatrix = x2*x2.transpose() + x3*x3.transpose() + 0.1*PriorDirc*PriorDirc.transpose();
		//Mat3d SigmaMatrix = x2*x2.transpose() + x3*x3.transpose() ;
		//+0.00*PriorDirc*PriorDirc.transpose();
		//Mat3d SigmaMatrix = x2*x2.transpose() + x3*x3.transpose();
		//Currdatap = Currdatap + 0.2*SigmaMatrix*P;
		//(Mat3d::Identity(3, 3) + 0.2*SigmaMatrix).cwiseInverse
		//Mat3d kk = (Mat3d::Identity(3, 3) + 0.2*SigmaMatrix).inverse();
		//std::cout << "kk=" << kk << std::endl;
		//Vec3d kk2 = (Currdatap + 0.2*SigmaMatrix*ReferenceP);
		//std::cout << "kk2=" << kk2 << std::endl;
		Currdatap = (Mat3d::Identity(3, 3) + 0.2*SigmaMatrix).inverse()*(Currdatap + 0.2*SigmaMatrix*ReferenceP1);//!!!
		//std::cout << "Currdatap= " << Currdatap << std::endl;
		SinglePoint = Currdatap;
	}
	return SinglePoint;
}

/*
xiong feng work

use reference!!!
*/
void CorrectTrace::CellCurvesaddedPoints(const CellVec3d& CellCurves, VectorVec3d &LabelVec, CellVec4d&NewCellCurves)
{
	VectorVec3d curCurve; VectorVec4d NewCurve;
	Vec4d NewVec; Vec3d label, tmp3d;
	vector<Vec3d>::iterator Vec3dIt;
	for (size_t i = 0; i < CellCurves.size(); ++i){
		curCurve = CellCurves[i];
		NewCurve.clear();
		for (Vec3dIt = curCurve.begin(); Vec3dIt != curCurve.end(); ++Vec3dIt){
			MakeVec4D(*Vec3dIt, 0, NewVec);
			NewCurve.push_back(NewVec);
		}
		NewCellCurves.push_back(NewCurve);
	}
	VectorVec4d data1;
	size_t Numss = LabelVec.size();
	for (size_t i = 0; i < Numss; ++i)
	{
		Vec3d CurrVec = LabelVec[i];
		VectorVec4d data0 = NewCellCurves[CurrVec(0) - 1];
		VectorVec4d data1 = NewCellCurves[CurrVec(1) - 1];
		if (CurrVec(2) == 0){
			Vec4d Addedpoint4d = data0[0];
			Vec3d Addedpoint;
			Make4DTo3D(Addedpoint4d, Addedpoint);
			AddingPointsTreeModify(data1, Addedpoint);
			NewCellCurves[CurrVec[1] - 1] = data1;
		}
		else
		{
			Vec4d Addedpoint4d = data0[data0.size() - 1];
			Vec3d Addedpoint;
			Make4DTo3D(Addedpoint4d, Addedpoint);
			AddingPointsTreeModify(data1, Addedpoint);
			NewCellCurves[CurrVec[1] - 1] = data1;
		}
	}
	/*
	for (Vec3dIt = LabelVec.begin(); Vec3dIt != LabelVec.end(); ++Vec3dIt){
	label = *Vec3dIt;
	if (label(2) == 0){
	Make4DTo3D(NewCellCurves[label(0) - 1].front(), tmp3d);
	}
	else{
	Make4DTo3D(NewCellCurves[label(0) - 1].back(), tmp3d);
	}
	AddingPointsTreeModify(NewCellCurves[label(1) - 1], tmp3d);
	}*/
}


void CorrectTrace::AddingPointsTreeModify(VectorVec4d &curCurves, Vec3d &AddPoint)
{
	//size_t Numss = curCurves.size();
	//vector<double> MinValues(Numss, 0.0);
	//VectorVec4d Curvess1;
	//for (size_t ii = 0; ii < Numss; ++ii){
	//	Vec3d tmp3d;
	//	Make4DTo3D(curCurves[ii], tmp3d);
	//	MinValues[ii] = (AddPoint - tmp3d).norm();
	//}
	//int MinIndex = distance(MinValues.begin(), min_element(MinValues.begin(), MinValues.end())) + 1;
	//if (MinIndex==0)
	//{
	//	//Curvess1.resize(curCurves.size());
	//	for (size_t i = 0; i < curCurves.size(); i++)
	//	{
	//		Curvess1.push_back(curCurves[i]);
	//	}
	//	Vec4d tmp4d;
	//	MakeVec4D(AddPoint, 1, tmp4d);
	//	Curvess1.insert(Curvess1.begin(), tmp4d);
	//}
	//if (MinIndex==Numss-1)
	//{
	//	for (size_t i = 0; i < curCurves.size(); i++)
	//	{
	//		Curvess1.push_back(curCurves[i]);
	//	}
	//	Vec4d tmp4d;
	//	MakeVec4D(AddPoint, 1, tmp4d);
	//	Curvess1.insert(Curvess1.end(), tmp4d);
	//}
	//if (MinIndex<Numss-1 &&MinIndex >1)
	//{
	//	
	//}

	double normer;
	vector <double> MinValues;
	vector<Vec4d>::iterator Vec4dIt;
	Vec3d subVec, tmp3d;
	for (Vec4dIt = curCurves.begin(); Vec4dIt != curCurves.end(); ++Vec4dIt){
		Make4DTo3D(*Vec4dIt, tmp3d);
		subVec = AddPoint - tmp3d;
		//SubVector3d(AddPoint, tmp3d, subVec);
		normer = subVec.norm();
		MinValues.push_back(normer);
	}
	size_t MinIdex; Vec4d tmp4d;
	MakeVec4D(AddPoint, 1, tmp4d);
	MinIdex = distance(MinValues.begin(), min_element(MinValues.begin(), MinValues.end())) + 1;
	curCurves.insert(curCurves.begin() + MinIdex, tmp4d);

}


void CorrectTrace::RevisedCurvesSetFixed(CellVec4d &CurvesSet, SVolume &ImageS, double Threshift, double ThreL1)
{
	//RevisedCurvesFixed(VectorVec3d & SingCurveSet, SVolume &ImageS, std::vector<size_t> LabelMatrix, double ThreShift, double ThreL1);
	VectorVec4d Curdata, NewCurve; VectorVec3d Curdata1; vector <size_t> Curdata2;  Vec3d tmp3d; Vec4d NewVec, tmp4d;
	vector<Vec4d>::iterator Vec4dIt; vector<Vec3d>::iterator Vec3dIt;
	for (size_t i = 0; i < CurvesSet.size(); ++i){
		Curdata = CurvesSet[i];
		Curdata1.clear();
		Curdata2.clear();
		for (Vec4dIt = Curdata.begin(); Vec4dIt != Curdata.end(); ++Vec4dIt){
			tmp4d = *Vec4dIt;
			Make4DTo3D(*Vec4dIt, tmp3d);
			Curdata1.push_back(tmp3d);
			Curdata2.push_back(tmp4d(3));
		}
		if (Curdata.size() > 10){
			RevisedCurvesFixed(Curdata1, ImageS, Curdata2, Threshift, ThreL1);
			NewCurve.clear();
			for (Vec3dIt = Curdata1.begin(); Vec3dIt != Curdata1.end(); ++Vec3dIt){
				MakeVec4D(*Vec3dIt, 0, NewVec);
				NewCurve.push_back(NewVec);
			}
			CurvesSet[i] = NewCurve;
		}
	}
}

//yi zou
void CorrectTrace::MakeVec4D(const Vec3d &arg, double a4, Vec4d &aim)
{
	aim << arg(0), arg(1), arg(2), a4;
}

void CorrectTrace::Make4DTo3D(const Vec4d &arg, Vec3d &aim)
{
	aim << arg(0), arg(1), arg(2);
}

void CorrectTrace::Make4DTo3Dadv(CellVec4d arg, CellVec3d &aim)
{
	vector<Vec4d>::iterator Vec4dIt; Vec3d tmp3d;
	VectorVec4d Curdata; VectorVec3d Curdata1;
	aim.clear();
	for (size_t i = 0; i < arg.size(); ++i){
		Curdata = arg[i];
		Curdata1.clear();
		for (Vec4dIt = Curdata.begin(); Vec4dIt != Curdata.end(); ++Vec4dIt){
			tmp3d[0] = Vec4dIt[0][0] + paramPack->xMin_;
			tmp3d[1] = Vec4dIt[0][1] + paramPack->yMin_;
			tmp3d[2] = Vec4dIt[0][2] + paramPack->zMin_;
			Curdata1.push_back(tmp3d);
			if (Vec4dIt == Curdata.begin())Curdata1.push_back(tmp3d);
		}
		aim.push_back(Curdata1);
	}
}


//void CorrectTrace::SubVector3d(const Vec3d &arg1, const Vec3d &arg2, Vec3d &aim)
//{
//	aim = arg1 - arg2;
//	//aim << arg1(0) - arg2(0), arg1(1) - arg2(1), arg1(2) - arg2(2);
//}


/*main function of smooth lines*/
void CorrectTrace::weak_signal_analysis(){

	XCellVec3d LocalCurve;
	VectorVec3d LabelVec;
	BifurLocalCurvesExtract(myCellCurves, myConnectSet, LocalCurve, LabelVec);
	size_t Numsss = LabelVec.size();
	VectorVec3d PointsV(Numsss);
	CellVec3d	CellCurves1 = myCellCurves;
	VectorVec3d Curves2, Curves1, Curve2, Curve1;

	for (int i = 0; i < Numsss; i++){
		Curves2 = LocalCurve[i][1];//father
		Curves1 = LocalCurve[i][0];//son
		if (Curves1.size() > 4 && Curves2.size()>2){
			Vec3d M1, M2, M3;
			VectorVec3d CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2;
			//RevisedCurvesBifur0Modify(Curves2, Curves1, src, Curve2, Curve1);
			//RevisedCurvesBifur0Modify_3line(Curves2, Curves1, src, Curve2, Curve1);%2018/1/17
			RevisedCurvesBifur0Modify_3line(Curves2, Curves1, src, Curve2, Curve1, M1,
				CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2);
			if (CurveBir_bak1.size() > 4 && CurveBir_bak2.size() > 4){
				RevisedCurvesBifur0Modify(CurveBase_bak1, CurveBir_bak1, src, CurveBase_bak1, CurveBir_bak1, M2);
				RevisedCurvesBifur0Modify(CurveBase_bak2, CurveBir_bak2, src, CurveBase_bak2, CurveBir_bak2, M3);
				VectorVec3d xxx;
				xxx.push_back(M1);
				xxx.push_back(M2);
				xxx.push_back(M3);
				std::vector<double> pValue;
				weigthvalue(xxx, src, pValue);
				Vec3d ModifiedPoint00 = (pValue[0] * xxx[0] + pValue[1] * xxx[1] + pValue[2] * xxx[2]);
				double fenmu = 0;
				for (int i = 0; i < pValue.size(); i++)
				{
					fenmu += pValue[i];
				}
				Curve1[0] = ModifiedPoint00 / fenmu;
			}
			else
				Curve1[0] = M1;

		}
		else
		{
			Curve1 = Curves1;
			Curve2 = Curves2;
		}
		Vec3d ReferenceP;
		size_t Index;
		VectorVec3d data0 = Curve1;
		PointsV[i] = data0[0];
		Vec3d CurrVec = LabelVec[i];
		VectorVec3d dataCurves = CellCurves1[CurrVec[0] - 1];//dataCurves=CellCurves1{CurrVec(1)};
		if (CurrVec[2] == 1){
			dataCurves[dataCurves.size() - 1] = data0[0];
			//std::reverse(dataCurves.begin(), dataCurves.end());
			/*FindSpecialPointInCurve2(Curve2, dataCurves, 4, ReferenceP,Index);
			dataCurves.erase(dataCurves.begin() + Index, dataCurves.end());
			dataCurves.push_back(data0[0]);*/
		}
		else
		{
			//FindSpecialPointInCurve2(Curve2, dataCurves, 4, ReferenceP, Index);
			//if (Index>0)dataCurves.erase(dataCurves.begin(), dataCurves.begin() + Index);
			/*dataCurves.push_back(data0[0]);*/
			//dataCurves.insert(dataCurves.begin(), data0[0]);
			dataCurves[0] = data0[0];
		}
		CellCurves1[CurrVec[0] - 1] = dataCurves;
	}

	CellVec4d NewCellCurves;
	CellCurvesaddedPoints(CellCurves1, LabelVec, NewCellCurves);
	//mytime = time(0);
	RevisedCurvesSetFixed(NewCellCurves, src, 2 * paramPack->smooth_level, std::min(200*paramPack->smooth_level,(double)6));
	//mytime2 = time(0);
	//std::cout << "Function: RevisedCurvesSetFixed time= " << mytime2 - mytime << std::endl;
	Make4DTo3Dadv(NewCellCurves, FinalCurves);


}

//[CurveBase_bak1, CurveBir_bak1, CurveBase_bak2, CurveBir_bak2] = distance_check(Up_curve, down_curve1, down_curve2)
void CorrectTrace::distance_check(const VectorVec3d &Up_curve, const VectorVec3d &down_curve1, const VectorVec3d &down_curve2,
	VectorVec3d &CurveBase_bak1, VectorVec3d &CurveBir_bak1, VectorVec3d &CurveBase_bak2, VectorVec3d &CurveBir_bak2){
	VectorVec3d Up_curve_new = Up_curve, down_curve1_new = down_curve1, down_curve2_new = down_curve2;
	size_t Num1 = Up_curve.size();
	Vec3d datap00, datap01, datap10, datap11;
	datap00 = down_curve1[0];
	datap01 = down_curve1[down_curve1.size() - 1];
	datap10 = down_curve2[0];
	datap11 = down_curve2[down_curve2.size() - 1];
	VectorVec4d MinDistMatrix(Num1);
	for (size_t i = 0; i < Num1; i++)
	{
		Vec3d currdata = Up_curve[i];
		double dis00 = (currdata - datap00).norm();
		double dis01 = (currdata - datap01).norm();
		double dis10 = (currdata - datap10).norm();
		double dis11 = (currdata - datap11).norm();
		MinDistMatrix[i] << dis00, dis01, dis10, dis11;
	}
	Vec4d idx, zz;
	for (size_t i = 0; i < MinDistMatrix.size(); i++)
	{
		if (i == 0)
		{
			idx = MinDistMatrix[0];
			zz << i, i, i, i;
		}
		else{
			for (size_t j = 0; j < MinDistMatrix[0].size(); j++)
			{
				if (MinDistMatrix[i](j)< idx(j))
				{
					idx(j) = MinDistMatrix[i](j);
					zz(j) = j;
				}
			}
		}
	}
	if (idx(0)<idx(1))
	{
		//down_curve1_new = down_curve1;
		if (zz(0) == 0)
			;//Up_curve_new = Up_curve;
		else{
			//Up_curve_new = Up_curve;
			Vec3d tmp;
			size_t n = Up_curve_new.size();
			for (int i = 0; i < n / 2; i++){
				tmp = Up_curve_new[i];
				Up_curve_new[i] = Up_curve_new[n - i - 1];
				Up_curve_new[n - i - 1] = tmp;
			}
		}
	}
	else
	{
		//down_curve1_new = down_curve1;
		Vec3d tmp;
		size_t n = down_curve1_new.size();
		for (int i = 0; i < n / 2; i++){
			tmp = down_curve1_new[i];
			down_curve1_new[i] = down_curve1_new[n - i - 1];
			down_curve1_new[n - i - 1] = tmp;
		}
		if (zz(1) == 0)
			;// Up_curve_new = Up_curve;
		else
		{
			//Up_curve_new = Up_curve;
			Vec3d tmp;
			size_t n = Up_curve_new.size();
			for (int i = 0; i < n / 2; i++){
				tmp = Up_curve_new[i];
				Up_curve_new[i] = Up_curve_new[n - i - 1];
				Up_curve_new[n - i - 1] = tmp;
			}

		}
	}

	if (idx(2) < idx(3))
		;//down_curve1_new = down_curve1;
	else
	{
		//down_curve2_new = down_curve2;
		Vec3d tmp;
		size_t n = down_curve2_new.size();
		for (int i = 0; i < n / 2; i++){
			tmp = down_curve2_new[i];
			down_curve2_new[i] = down_curve2_new[n - i - 1];
			down_curve2_new[n - i - 1] = tmp;
		}
	}
	CurveBase_bak1 = down_curve1_new;
	Vec3d tmp0;
	size_t n0 = CurveBase_bak1.size();
	for (int i = 0; i < n0 / 2; i++){
		tmp0 = CurveBase_bak1[i];
		CurveBase_bak1[i] = CurveBase_bak1[n0 - i - 1];
		CurveBase_bak1[n0 - i - 1] = tmp0;
	}
	CurveBase_bak1.pop_back();
	//CurveBase_bak1.insert(CurveBase_bak1.end(), Up_curve_new.begin(), Up_curve_new.end());
	for (size_t i = 0; i < Up_curve_new.size(); ++i){
		CurveBase_bak1.push_back(Up_curve_new[i]);
	}
	CurveBir_bak1 = down_curve2_new;


	CurveBase_bak2 = down_curve2_new;
	Vec3d tmp1;
	size_t n1 = CurveBase_bak2.size();
	for (int i = 0; i < n1 / 2; i++){
		tmp1 = CurveBase_bak2[i];
		CurveBase_bak2[i] = CurveBase_bak2[n1 - i - 1];
		CurveBase_bak2[n1 - i - 1] = tmp1;
	}
	CurveBase_bak2.pop_back();
	//CurveBase_bak2.insert(CurveBase_bak2.end(), Up_curve_new.begin(), Up_curve_new.end());
	for (size_t i = 0; i < Up_curve_new.size(); i++){
		CurveBase_bak2.push_back(Up_curve_new[i]);
	}
	CurveBir_bak2 = down_curve1_new;
}

void CorrectTrace::weigthvalue(const VectorVec3d &dataL1, SVolume &L_XX3, std::vector<double> &aa_data){
	int nxss = dataL1.size();
	aa_data.resize(nxss);
	int nx = L_XX3.x();
	int ny = L_XX3.y();
	int nz = L_XX3.z();
	for (int i = 0; i < nxss; i++)
	{
		double x = dataL1[i](0);
		double y = dataL1[i](1);
		double z = dataL1[i](2);
		double dd = 0, ww = 0;
		int idexx = std::max((int)std::min((int)std::round(x), nx), (int)1);
		int idexy = std::max((int)std::min((int)std::round(y), ny), (int)1);
		int idexz = std::max((int)std::min((int)std::round(z), nz), (int)1);
		double w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x + 1), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x - 1), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y - 1), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y + 1), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z - 1), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		idexx = std::max((int)std::min((int)std::round(x), nx), (int)1);
		idexy = std::max((int)std::min((int)std::round(y), ny), (int)1);
		idexz = std::max((int)std::min((int)std::round(z + 1), nz), (int)1);
		w1 = 0.05*std::max((-2)*(pow(x - idexx, 2) + pow(y - idexy, 2) + pow(2 * z - 2 * idexz, 2)), (double)-6);
		dd = dd + exp(w1)*L_XX3.operator()(idexx, idexy, idexz);
		ww = ww + exp(w1);

		aa_data[i] = dd / (ww + 0.0001);
	}
}