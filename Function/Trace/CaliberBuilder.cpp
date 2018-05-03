
#include <Eigen/Eigenvalues>
#include "CaliberBuilder.h"
#include "../IO/imagewriter.h"
#include "../Function/NGUtility.h"

#pragma warning(disable:4996)
CaliberBuilder::CaliberBuilder()
{
    className_ = "CaliberBuilder";
	std::cout << "CaliberBuilder Start!" << std::endl;
}

CaliberBuilder::~CaliberBuilder()
{
}
/*don't use*/
INEURONPROCESSOBJECT_RELEASEDATA_IMPLE(CaliberBuilder)
/*don't use*/
INEURONPROCESSOBJECT_GETOUTPUT_IMPLE(CaliberBuilder, NeuronPopulation)

void CaliberBuilder::SetLinePoints(){
    KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, Currentlines, CurrentIDs,
        paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
    IsForRange = true;
}
/*
set the line to smooth and correct
choose  part of lines in current image
*/
void CaliberBuilder::SetLinePoints(const std::vector<size_t> &curBranchID){
    KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, curBranchID, Currentlines, CurrentIDs,
        paramPack->xMin_, paramPack->xMax_, paramPack->yMin_, paramPack->yMax_, paramPack->zMin_, paramPack->zMax_);
    IsForRange = false;
}
/*
set the line to smooth and correct
local range lines
*/
void CaliberBuilder::SetLinePoints(size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max){
    x_min = x_min > paramPack->xMin_ ? x_min : paramPack->xMin_;
    y_min = y_min > paramPack->yMin_ ? y_min : paramPack->yMin_;
    z_min = z_min > paramPack->zMin_ ? z_min : paramPack->zMin_;
    x_max = x_max < paramPack->xMax_ ? x_max : paramPack->xMax_;
    y_max = y_max < paramPack->yMax_ ? y_max : paramPack->yMax_;
    z_max = z_max < paramPack->zMax_ ? z_max : paramPack->zMax_;
    KPTREE::FindCurrentImagePoints(paramPack->activeTree->m_curveList, Currentlines, CurrentIDs,
        x_min, x_max, y_min, y_max, z_min, z_max);
    IsForRange = true;
}

void CaliberBuilder::SetLinePoints(const std::vector<Line5d>& arg, double xMin, double xMax, double yMin, double yMax, double zMin, double zMax)
{
    Currentlines = arg;
    for (auto &it : Currentlines) {
        for (auto &iter : it) {
            if (iter(0) < xMin) iter(0) = xMin;
            if (iter(0) > xMax) iter(0) = xMax;
            if (iter(1) < yMin) iter(1) = yMin;
            if (iter(1) > yMax) iter(1) = yMax;
            if (iter(2) < zMin) iter(2) = zMin;
            if (iter(2) > zMax) iter(2) = zMax;
        }
    }
}

void CaliberBuilder::SetLinePoints(const std::vector<Line3d>& arg)
{
    Currentlines.clear();
    Line5d tmp5d;
    for (auto &it : arg) {
        tmp5d.clear();
        tmp5d.reserve(it.size());
        for (auto &iter:it) {
            tmp5d.push_back(NGUtility::MakeVec5d(iter));
        }
        Currentlines.emplace_back();
        Currentlines.back().swap(tmp5d);
    }
    /*for (auto &it : Currentlines) {
        for (auto &iter : it) {
        if (iter(0) < xMin) iter(0) = xMin;
        if (iter(0) > xMax) iter(0) = xMax;
        if (iter(1) < yMin) iter(1) = yMin;
        if (iter(1) > yMax) iter(1) = yMax;
        if (iter(2) < zMin) iter(2) = zMin;
        if (iter(2) > zMax) iter(2) = zMax;
        }
        }*/
}

ProcStatPointer CaliberBuilder::Update(){
    //find inner part of target curve and connective curves. get connect information from m_Topo
    if (Currentlines.empty()) {
        MAKEPROCESSSTATUS(resSta, false, className_, "no input lines.");
        return resSta;
    }
    /*if (Currentlines.size() <3) {
        MAKEPROCESSSTATUS(resSta, false, className_, "input trace line is too short. trace is failed");
        return resSta;
        }*/
    if (IsForRange){
        if (!SWCFormatChangeConnect(1.5, paramPack->xMin_, paramPack->yMin_, paramPack->zMin_)){
            MAKEPROCESSSTATUS(resSta, false, className_, "SWCFormatChangeConnect failed.");
            NG_ERROR_MESSAGE("SWCFormatChangeConnect failed.");
            return resSta;
        }
    }
    else{
        if (!FastSWCFormatChange(paramPack->xMin_, paramPack->yMin_, paramPack->zMin_)){
            MAKEPROCESSSTATUS(resSta, false, className_, "FastSWCFormatChange failed.");
            NG_ERROR_MESSAGE("FastSWCFormatChange failed.");
            return resSta;
        }
    }
    NG_CREATE_DYNAMIC_CAST(SVolume, tmpOrig, paramPack->OrigImage);
    if (!tmpOrig || tmpOrig->IsEmpty()){
        MAKEPROCESSSTATUS(resSta, false, className_, "no input image.");
        return resSta;
    }
    //src.QuickCopy(*mysrc);
	if (!caliberPointSet) caliberPointSet = std::make_shared<CellVec3d>();
    NeuritesShapeReconModifyTotal(myCellCurves, myConnectSet, *tmpOrig, *caliberPointSet);
    //caliberPointSet.erase(caliberPointSet.begin() + 1, caliberPointSet.end());
    //map to global space
    for (auto &it : *caliberPointSet) {
        for (auto &pt: it) {
            pt(0) += paramPack->xMin_;
            pt(1) += paramPack->yMin_;
            pt(2) += paramPack->zMin_;
        }
    }

    MAKEPROCESSSTATUS(resSta, true, className_, "");
	std::cout << "CaliberBuilder Completed!" << std::endl;
    return resSta;
}
bool CaliberBuilder::FastSWCFormatChange(double x_min, double y_min, double z_min){
    if (Currentlines.empty()){
        return false;
    }
	std::cout << "Doing FastSWCFormatChange" << std::endl;
    myCellCurves.clear();
    myConnectSet.clear();
    VectorVec3d DataPoints, ResampleData;
    for (size_t i = 0; i < Currentlines.size(); ++i){
        DataPoints.clear();
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
        myCellCurves.push_back(DataPoints);
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
                    if (parentID == CurrentIDs[j][0][0]) {
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
        myConnectSet.push_back(ConVec);
        /*paramPack->activeTree->m_Topo->*/
    }
    return true;
}

/*
change swc format and resample lines
*/
bool CaliberBuilder::SWCFormatChangeConnect(double Thres, double x_min, double y_min, double z_min)
{
    if (Currentlines.empty()){
        return false;
    }
	std::cout << "Doing SWCFormatChangeConnect" << std::endl;
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
        myCellCurves.push_back(DataPoints);
    }
    MatXd DistMatrixStarP, DistMatrixEndP;
    SWCFormatChangeSub1(myCellCurves, DistMatrixStarP, DistMatrixEndP);
    for (int i = 0; i < Currentlines.size(); i++){
        MatXd	CurrVec1 = DistMatrixStarP.col(i);
        MatXd	CurrVec2 = DistMatrixEndP.col(i);
        Eigen::MatrixXf::Index minRow1, minCol1;
        double min1 = CurrVec1.minCoeff(&minRow1, &minCol1);
        Eigen::MatrixXf::Index minRow2, minCol2;
        double min2 = CurrVec2.minCoeff(&minRow2, &minCol2);
        Vec2d ConVec;
        ConVec.setZero();
        if (min1<Thres)
            ConVec[0] = minRow1 + 1;
        if (min2<Thres)
            ConVec[1] = minRow2 + 1;
        myConnectSet.push_back(ConVec);
    }
    return  true;
}

void CaliberBuilder::SWCFormatChangeSub1(const CellVec3d &Curvesets, MatXd &Distances1, MatXd &Distances2){
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
void CaliberBuilder::SWCFormatChangeSub11(const Vec3d &data1, const Vec3d &data2, const VectorVec3d &Currdata0, double &ds1, double &ds2)
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
}
void CaliberBuilder::SampleCylinderDataSetsExtraModify(VectorVec3d &Curves, SVolume & OrigIma, int ScaleSubImage, float Pixelsize, bool IsNew, VectorMatXd &SampleDataset, MatXd &DirecVecs) {
	std::cout << "Doing SampleCylinderDataSetsExtraModify" << std::endl;
    auto Nx = OrigIma.x();
    auto Ny = OrigIma.y();
    auto Nz = OrigIma.z();
    auto NumCurves = Curves.size();
    MatXd DirecVecstmp = DirecVecs;
    if (IsNew){
        CurvesDircs(Curves, DirecVecs);
    }

    for (int i = 0; i < NumCurves - 1; i++){
        Vec3d x1 = DirecVecs.col(i).middleRows(6, 3);
        Vec3d x2 = DirecVecs.col(i).middleRows(9, 3);
        Vec3d y1 = DirecVecs.col(i + 1).middleRows(6, 3);
        Vec3d y2 = DirecVecs.col(i + 1).middleRows(9, 3);
        Vec3d yy1, yy2;
        OrthogonalsystemMatch(x1, x2, y1, y2, yy1, yy2);
        DirecVecs.col(i + 1).middleRows(6, 3) = yy1;
        DirecVecs.col(i + 1).middleRows(9, 3) = yy2;
    }

    for (int i = 0; i < NumCurves; i++){
        MatXd CurrMatrix = MatXd::Zero(ScaleSubImage, ScaleSubImage);
        Vec3d  CenterPoint = DirecVecs.col(i).middleRows(0, 3);
        Vec3d  Xdir = DirecVecs.col(i).middleRows(6, 3);
        Vec3d  Ydir = DirecVecs.col(i).middleRows(9, 3);
        for (int ii = -(ScaleSubImage - 1) / 2; ii <= (ScaleSubImage - 1) / 2; ii++){
            for (int jj = -(ScaleSubImage - 1) / 2; jj <= (ScaleSubImage - 1) / 2; jj++){
                Vec3d  Positions = ii*Pixelsize*Xdir + jj*Pixelsize*Ydir + CenterPoint;
                Positions(0) = std::min(std::max(Positions(0), 0.0), Nx - 1.0);
                Positions(1) = std::min(std::max(Positions(1), 0.0), Ny - 1.0);
                Positions(2) = std::min(std::max(Positions(2), 0.0), Nz - 1.0);
                double v1 = weigthvalue(Positions, OrigIma);
                CurrMatrix(ii + (ScaleSubImage - 1) / 2, jj + (ScaleSubImage - 1) / 2) = v1;
            }
        }
        SampleDataset.push_back(CurrMatrix);
        //cout << "CurrMatrix=" << endl << CurrMatrix << endl;
    }
    if (!IsNew){
        DirecVecs = DirecVecstmp;
    }
}
void CaliberBuilder::SBGCSTotalTT(const MatXd& ImageF, MatXd &ImagPrio, double Lamda0, double Threv, double LamdaPrio, MatXd & USeg, MatXd &U, Vec2d & c1_c2){
    double Threv0 = 0.5*(ImageF.maxCoeff() + ImageF.minCoeff());
    U = ImageF;
    for (int xx = 0; xx < ImageF.rows(); xx++){
        for (int yy = 0; yy < ImageF.cols(); yy++){
            U.operator()(xx, yy) = ImageF.operator()(xx, yy)>Threv0;
        }
    }
    //cout << U;
    MatXd Ux, Uy;
    FusedImageD(U, Ux, Uy);
    MatXd dx, dy, bx, by;
    dx.setZero(Ux.rows(), Ux.cols());
    dy.setZero(Uy.rows(), Uy.cols());
    bx = dx;
    by = dy;
    MatXd R;
    SBGCSsub1(ImageF, U, Threv, R, c1_c2);
    //cout << ImageF << endl << endl;
    //cout << R;
    if (c1_c2[0] - c1_c2[1]>100){
        Threv = 0.3;
        LamdaPrio = 0.1;
    }
    else if (c1_c2[0] - c1_c2[1]>50)
    {
        Threv = 0.3;
        LamdaPrio = 2;
    }
    else
    {
        Threv = 0.5;
        LamdaPrio = 4;
    }
    USeg = U;
    MatXd X, LamdaVectorX, LamdaVectorY, Y;
    for (int i = 0; i < 10; i++)
    {
        SBGCSsub1(ImageF, U, Threv, R, c1_c2);
        double tmp = 10 / (R.cwiseAbs()).maxCoeff();
        R /= tmp;
        //R = 10 * R.array() / (R.cwiseAbs()).maxCoeff();
        for (int xx = 0; xx < U.rows(); xx++){
            for (int yy = 0; yy < U.cols(); yy++){
                USeg.operator()(xx, yy) = U.operator()(xx, yy)>Threv;
            }
        }
        if (USeg.sum() < 11)break;
        //cout << "R=" << endl << R << endl;
        double Objectvalue = 0;
        for (int jj = 0; jj < 5;++ jj){
            //cout << "U=" << endl << U << endl << "Objectvalue=" << Objectvalue<<endl;
            NueriteShapeIterativeOP(U, R, ImagPrio, Lamda0, Threv, LamdaPrio, dx, dy, bx, by, Objectvalue);
        }
        //cout << U<<endl;
        FusedImageD(U, Ux, Uy);
        //MatXd X, LamdaVector, Y;
        //X = Ux + bx;
        X = Ux; X += bx;
        //LamdaVector = 1 / Lamda0*MatXd::Ones(Ux.rows(), Ux.cols());
        LamdaVectorX.setOnes(Ux.rows(), Ux.cols());
        LamdaVectorX *= 1.0 / Lamda0;
        SoftThresholdingOperation(X, LamdaVectorX, dx);
        //MatXd X, LamdaVector, Y;
        //X = Uy + by;
        X = Uy; X += by;
        //LamdaVector = 1 / Lamda0*MatXd::Ones(Uy.rows(), Uy.cols());
        LamdaVectorY.setOnes(Uy.rows(), Uy.cols());
        LamdaVectorY *= 1.0 / Lamda0;
        SoftThresholdingOperation(X, LamdaVectorY, dy);
        //bx = bx + (Ux - dx);
        bx += Ux; bx -= dx;
        //by = by + (Uy - dy);
        by += Uy;
        by -= dy;
    }
    //cout <<"U="<< endl<<U << endl;
    for (int xx = 0; xx < U.rows(); xx++){
        for (int yy = 0; yy < U.cols(); yy++){
            USeg.operator()(xx, yy) = U.operator()(xx, yy)>Threv;
        }
    }
}

void CaliberBuilder::FusedImageD(MatXd& RecImage, MatXd& Imag0, MatXd& Imag1){
    auto Nx = RecImage.rows();
    auto Ny = RecImage.cols();
    Imag0.setZero(Nx, Ny - 1);
    Imag1.setZero(Nx - 1, Ny);
    for (int ii = 0; ii < Nx; ii++){
        Imag0.row(ii) = RecImage.row(ii).tail(Ny - 1) - RecImage.row(ii).head(Ny - 1);
    }
    for (int ii = 0; ii < Ny; ii++){
        Imag1.col(ii) = RecImage.col(ii).tail(Nx - 1) - RecImage.col(ii).head(Nx - 1);
    }
}
void CaliberBuilder::SBGCSsub1(const MatXd &ImageF, MatXd& U, double Thres, MatXd & R, Vec2d & c1_c2){
    MatXd ss = U;
    for (int xx = 0; xx < ImageF.rows(); xx++){
        for (int yy = 0; yy < ImageF.cols(); yy++){
            ss.operator()(xx, yy) = U.operator()(xx, yy)>Thres;
        }
    }
    //MatXd sb = MatXd::Ones(ss.rows(), ss.cols()) - ss;
    MatXd sb; sb.setOnes(ss.rows(), ss.cols()); sb -= ss;
    c1_c2[0] = (ss.array()*ImageF.array()).sum() / std::max(ss.sum(), 1.0);
    c1_c2[1] = (sb.array()*ImageF.array()).sum() / std::max(sb.sum(), 1.0);
    R = (ImageF.array() - c1_c2[0]).square() - (ImageF.array() - c1_c2[1]).square();
}

void CaliberBuilder::NueriteShapeIterativeOP(MatXd &SegImage, MatXd &SegImageInf, MatXd &SegImageprio, double Lamda, double Mu, double Lamda0, MatXd &Ua, MatXd & Ub, MatXd & Va, MatXd &Vb, double &Objectvalue){
    //MatXd Imag0, Imag1, GradientImage2, GradientImage3, GradientImage4, GradientImage;

    FusedImageD(SegImage, NSIOPImag0, NSIOPImag1);
    FusedImageDTranspose(NSIOPImag0, NSIOPImag1, NSIOPGradientImage2);
    FusedImageDTranspose(Ua, Ub, NSIOPGradientImage3);
    FusedImageDTranspose(Va, Vb, NSIOPGradientImage4);
    //NSIOPGradientImage = NSIOPGradientImage2 - NSIOPGradientImage3 + NSIOPGradientImage4 + Mu / Lamda*SegImageInf + Lamda0 / Lamda*(SegImage - SegImageprio);
    NSIOPGradientImage = NSIOPGradientImage2; NSIOPGradientImage -= NSIOPGradientImage3; NSIOPGradientImage += NSIOPGradientImage4;
    NSIOPtmp = SegImageInf; NSIOPtmp *= Mu / Lamda; NSIOPGradientImage += NSIOPtmp;
    NSIOPtmp = SegImage; NSIOPtmp -= SegImageprio; NSIOPtmp *= Lamda0 / Lamda; NSIOPGradientImage += NSIOPtmp;

    double tmp = 10.0 / (std::abs(NSIOPGradientImage.maxCoeff()) + 0.01);
    NSIOPGradientImage *= tmp;
    //NSIOPGradientImage = 10 * NSIOPGradientImage / (std::abs(NSIOPGradientImage.maxCoeff()) + 0.01);
    double tt1 = NueriteShapeObjectFun(SegImage, SegImageInf, SegImageprio, Va, Vb, Ua, Ub, Lamda, Mu, Lamda0);
    //MatXd JJk;
    double tt2 = 10000;
    for (int i = 0; i < 10; i++){
        //JJk = SegImage - (1.0 / std::pow(2, i + 1))*NSIOPGradientImage;
        NSIOPJJk = NSIOPGradientImage; NSIOPJJk *= -(1.0 / std::pow(2, i + 1)); NSIOPJJk += SegImage;
        for (int xx = 0; xx < NSIOPJJk.rows(); xx++){
            for (int yy = 0; yy < NSIOPJJk.cols(); yy++){
                NSIOPJJk(xx, yy) = std::min(std::max(NSIOPJJk(xx, yy), 0.0), 1.0);
            }
        }
        tt2 = NueriteShapeObjectFun(SegImage, SegImageInf, SegImageprio, Va, Vb, Ua, Ub, Lamda, Mu, Lamda0);
        if (tt2 < tt1)	break;
    }
    SegImage.swap(NSIOPJJk);
    //SegImage = NSIOPJJk;s
    Objectvalue = tt2;
    //cout << GradientImage << endl;
}

void CaliberBuilder::FusedImageDTranspose(MatXd&Imag0, MatXd&Imag1, MatXd&GradientImage2){
    auto Nx0 = Imag0.rows();
    auto Ny0 = Imag0.cols();
    auto Nx1 = Imag1.rows();
    auto Ny1 = Imag1.cols();
    //MatXd Imag00 = MatXd::Zero(Nx0, Ny0 + 1);
    FDTImag00.setZero(Nx0, Ny0 + 1);
    //MatXd Imag11 = MatXd::Zero(Nx1 + 1, Ny1);
    FDTImag11.setZero(Nx0, Ny0 + 1);
    for (int i = 0; i < Nx0; i++){
        auto &ss = Imag0.row(i);
        FDTImag00(i, 0) = -ss(0);
        FDTImag00.row(i).middleCols(1, Ny0 - 1) = ss.head(Ny0 - 1) - ss.tail(Ny0 - 1);
        FDTImag00(i, Ny0) = ss(Ny0 - 1);
    }
    for (int i = 0; i < Ny1; i++){
        auto &ss = Imag1.col(i);
        FDTImag11(0, i) = -ss(0);
        FDTImag11.col(i).middleRows(1, Nx1 - 1) = ss.head(Nx1 - 1) - ss.tail(Nx1 - 1);
        FDTImag11(Nx1, i) = ss(Nx1 - 1);
    }
    GradientImage2 = FDTImag00;
    GradientImage2 += FDTImag11;
}
double  CaliberBuilder::NueriteShapeObjectFun(MatXd &SegImage, MatXd & SegImageInf, MatXd &SegImagePrio, MatXd &Va, MatXd &Vb, MatXd &Ua, MatXd & Ub, double Lamda, double Mu, double Lamda0){

    //MatXd SS = SegImage.array()*SegImageInf.array();
    NSOFSS = SegImage;
    //NSOFSS = NSOFSS.array() * SegImageInf.array();
    auto nx = NSOFSS.rows(); auto ny = NSOFSS.cols();
    for (auto i = 0; i < nx; ++i) {
        for (auto j = 0; j < ny; ++j) {
            NSOFSS(i, j) *= SegImageInf(i, j);
        }
    }
    double ss = 0.5*Lamda*NSOFSS.sum();
    //NSOFSS = (SegImage - SegImagePrio);
    NSOFSS = SegImage;
    NSOFSS -= SegImagePrio;
    //NSOFSS = NSOFSS.array()*NSOFSS.array();
    nx = NSOFSS.rows(); ny = NSOFSS.cols();
    for (auto i = 0; i < nx; ++i) {
        for (auto j = 0; j < ny; ++j) {
            NSOFSS(i, j) = NSOFSS(i, j) * NSOFSS(i, j);
        }
    }
    ss += 0.5*Lamda0*NSOFSS.sum();
    //MatXd Imag0, Imag1;
    FusedImageD(SegImage, NSOFImag0, NSOFImag1);
    //NSOFSS = (NSOFImag0 - Ua + Va).array().square();
    NSOFSS = NSOFImag0; NSOFSS -= Ua; NSOFSS += Va;
    //NSOFSS = NSOFSS.array().square();
    nx = NSOFSS.rows(); ny = NSOFSS.cols();
    for (auto i = 0; i < nx; ++i) {
        for (auto j = 0; j < ny; ++j) {
            NSOFSS(i, j) = NSOFSS(i, j) * NSOFSS(i, j);
        }
    }
    ss +=  Mu*NSOFSS.sum();
    //NSOFSS = (NSOFImag1 - Ub + Vb).array().square();
    NSOFSS = NSOFImag1; NSOFSS -= Ub; NSOFSS += Vb;
    NSOFSS = NSOFSS.array().square();
    ss += Mu*NSOFSS.sum();
    return ss;
}

void CaliberBuilder::SoftThresholdingOperation(MatXd &X, MatXd & LamdaVector, MatXd &Y){
    auto Nxx = X.rows();
    auto Nyy = X.cols();
    Y.setZero(Nxx, Nyy);
    for (int xx = 0; xx < Nxx; xx++){
        for (int yy = 0; yy < Nyy; yy++){
            auto Currdata = X(xx, yy);
            auto CurrThreshhold = LamdaVector(xx, yy);
            if (std::abs(Currdata) > CurrThreshhold){
                Y(xx, yy) = NGUtility::sign(Currdata)*(std::abs(Currdata) - CurrThreshhold);
            }
        }
    }
}
void CaliberBuilder::NeuritesShapeReconModifyTotal(CellVec3d &curveSet, VectorVec2d& connectSet, SVolume &origImag, CellVec3d &shapePointSet){

	std::cout << "Doing NeuritesShapeReconModifyTotal" << std::endl;
    CellVec3d curveSet0 = curveSet;
    auto Numss = curveSet0.size();
    VectorVecXd LabelVecSet0;
    for (int ii = 0; ii < Numss; ++ii){
        Vec2d  ConMatrix = connectSet[ii];
		//curveSet0[ii] = curveSet[ii];
		VectorVec3d &datap = curveSet0[ii];
        if (ConMatrix(0) != 0){
            datap.erase(datap.begin());
        }
        if (ConMatrix(1) != 0){
            datap.erase(datap.end() - 1);
        }
        //curveSet0[ii] = datap;
        LabelVecSet0.push_back(VecXd::Zero(datap.size()));
    }
	//force each curves to 3 node
	//for (auto &it : curveSet0){
	//	if (it.size() ==  1) {
	//		VectorVec3d datapp(2*it.size());///////////  3 ??
	//		for (int n = 0; n < it.size(); ++n){
	//			datapp[2 * n] = it[n];
	//			datapp[2 * n + 1] = it[n];
	//		}
	//		std::swap(it, datapp);
	//	}
	//	else if (it.size() == 2){
	//		VectorVec3d datapp;
	//		datapp.resize(2 * it.size() - 1);
	//		for (int n = 0; n < it.size() - 1; ++n){
	//			datapp[2 * n] = it[n];
	//			datapp[2 * n + 1] = (it[n] + it[n + 1])*0.5;
	//		}
	//		std::swap(it, datapp);
	//	}
	//}
    auto Numx = origImag.x();
    auto Numy = origImag.y();
    auto Numz = origImag.z();
    SVolume ImageLabel;
    ImageLabel.SetSize(Numx, Numy, Numz);
    ImageLabel.SetZero();
    for (int ii = 0; ii < Numss; ++ii){
        VectorVec3d datap = curveSet0[ii];
        auto nump = datap.size();
		VectorVec3d datapp;
		datapp.resize(2 * nump - 1);
		for (int n = 0; n < nump - 1; ++n){
			datapp[2 * n] = datap[n];
			datapp[2 * n + 1] = (datap[n] + datap[n + 1])*0.5;
		}
		for (int kk = 0; kk < datapp.size(); ++kk){
            Vec3d p;
            p[0] = std::round(datapp[kk][0]);
            p[1] = std::round(datapp[kk][1]);
            p[2] = std::round(datapp[kk][2]);
            auto x1 = std::min(std::max(p[0], 0.0), Numx - 1.0);
            auto y1 = std::min(std::max(p[1], 0.0), Numy - 1.0);
            auto z1 = std::min(std::max(p[2], 0.0), Numz - 1.0);
            ImageLabel(x1, y1, z1) = ii + 1;
        }
    }
    CellVec3d ShapePointSet1, ShapePointSet2;
    VectorVecXd LabelVecSet1, LabelVecSet2;
    NeuritesShapeReconModify(curveSet0, origImag, ImageLabel, 21, LabelVecSet0, true, shapePointSet, LabelVecSet1);
    NeuritesShapeReconModify(curveSet0, origImag, ImageLabel, 11, LabelVecSet1, false, ShapePointSet2, LabelVecSet2);
    assert(ShapePointSet2.size() == shapePointSet.size());
    for (int ii = 0; ii < shapePointSet.size(); ++ii) {
        std::move(ShapePointSet2[ii].begin(), ShapePointSet2[ii].end(), std::back_inserter(shapePointSet[ii]));
    }
}
void CaliberBuilder::NeuritesShapeReconModify(CellVec3d & CurveSet, SVolume & OrigImage, SVolume & ImagLabel, int SampleImgsize, VectorVecXd &LabelSet, bool IsNew, CellVec3d & shapePointSet, VectorVecXd &LabelVecSet){
	std::cout << "Doing NeuritesShapeReconModify" << std::endl;
    auto curveSetNum = CurveSet.size();
    //Nums = 1;
    for (auto jj = 0; jj < curveSetNum; ++jj) {
        VectorMatXd SampleDataset;
        VectorMatXd SegmentDataSet;
        VectorVec2d MatrixV;
        //MatXd USeg, U;
        //Vec2d c1_c2;
        //MatXd MMkk;
        MatXd DirecVecs, DirecVecsNew;
        MatXd ImagPrio;
        //VectorMatXd ImagPrioList;
        CellVec7d dataPointSet;
        VectorVec3d dataSet;
        VecXd Labelvec2;
        auto Curves = CurveSet[jj];
        auto Labelvec = LabelSet[jj];
        VectorVec3d CurvesNew;
        for (int nn = 0; nn < Labelvec.size(); nn++){
            if (Labelvec[nn] == 0){
                CurvesNew.push_back(Curves[nn]);
            }
        }
        if (CurvesNew.size() > 0)
        {
            if (!IsNew)
            {
                CurvesDircs(Curves, DirecVecsNew);
                DirecVecs = MatXd(12, CurvesNew.size());
                int kk = 0;
                for (int nn = 0; nn < Labelvec.size(); nn++){
                    if (Labelvec[nn] == 0){
                        DirecVecs.col(kk) = DirecVecsNew.col(nn);
                        kk++;
                    }
                }
            }
            SampleCylinderDataSetsExtraModify(CurvesNew, OrigImage, SampleImgsize, 0.5, IsNew, SampleDataset, DirecVecs);
            //cout << "SampleDataset[0]=" << endl << SampleDataset[0]<< endl;
            //double *RR = new double[SampleDataset.size()];
            //std::vector<double> RR(SampleDataset.size());
            auto beg = clock();
            MatrixV.resize(SampleDataset.size());
            SegmentDataSet.resize(SampleDataset.size());
            //ImagPrioList.resize(SampleDataset.size());
            //for (auto k = 0; k < SampleDataset.size(); ++k)
            //    ImagPrioList[k].setZero(SampleDataset[k].rows(), SampleDataset[k].cols());
            for (int ii = 0; ii < SampleDataset.size(); ++ii){
                const MatXd &MMkk = SampleDataset[ii];
                MatXd USeg, U;
                Vec2d c1_c2;
                //SBGCSTotalTT(MMkk, ImagPrioList[ii], 1, 0.5, 2, USeg, U, c1_c2);
                SBGCSTotalTT(MMkk, ImagPrio.setZero(MMkk.rows(), MMkk.cols()), 1, 0.5, 2, USeg, U, c1_c2);
                SegmentDataSet[ii].swap(USeg);
                //SegmentDataSet.push_back(USeg);
                //SegmentDataSet1.push_back(U);
                //MatrixV.push_back(c1_c2);
                MatrixV[ii].swap(c1_c2);
                //RR[ii] = std::sqrt(USeg.sum()*0.25 / 3.14);
                //if (RR[ii]>5) printf("R>5\n");
            }
            auto fin = clock();
            //printf("%d ms eclipsed in SBGCSTotalTT. \n", int(fin - beg));
            //beg = clock();
            SegmentpointModify(SegmentDataSet, DirecVecs, OrigImage, MatrixV, dataPointSet);
            //fin = clock();
            //printf("%d ms eclipsed in SegmentpointModify. \n", int(fin - beg));
            //beg = clock();
            SegmentPointSceening(dataPointSet, ImagLabel, OrigImage, jj + 1, dataSet, Labelvec2);
            //fin = clock();
            //printf("%d ms eclipsed in SegmentPointSceening. \n", int(fin - beg));
            //ShapePointSet.push_back(dataSet);
            shapePointSet.emplace_back();
            shapePointSet.back().swap(dataSet);
            LabelVecSet.push_back(Labelvec2);
        }
        else  {
            dataSet.clear();
            shapePointSet.push_back(dataSet);
        }
        //std::cout << "line[" << jj << "]  size= " << shapePointSet.back().size() << "   " << 1.0*(jj + 1) / curveSetNum << " is complete !" << std::endl;
    }
}
void CaliberBuilder::principald(VectorVec3d &dataL, Vec3d&x1){
    auto nx = dataL.size();
    Vec3d nn = Vec3d::Zero();
    for (int i = 0; i < nx - 1; i++){
        nn = nn + (dataL[i + 1] - dataL[i]).norm()*(dataL[i + 1] - dataL[i]);
    }
    x1 = nn / nn.norm();
}
void CaliberBuilder::OrthogonalsystemMatch(Vec3d&x1, Vec3d&x2, Vec3d&y1, Vec3d&y2, Vec3d&yy1, Vec3d&yy2){
    Vec3d y1_new, y2_new;
    if (std::abs(y1.transpose()*x1)>0.1){
        y1_new = NGUtility::sign(y1.transpose()*x1)*y1;
    }
    else
    {
        y1_new = y1;
    }
    if (std::abs(y2.transpose()*x2)>0.1){
        y2_new = NGUtility::sign(y2.transpose()*x2)*y2;
    }
    else
    {
        y2_new = y2;
    }
    double s11 = x1.transpose()*y1_new;
    double s12 = x2.transpose()*y2_new;
    double s1 = s11 + s12;
    double s21 = x1.transpose()*y2_new;
    double s22 = y1_new.transpose()*x2;
    double s2 = s21 - s22;
    if ((s1*s1 + s2*s2) > 0.000001){
        double t = std::sqrt(1 / (s1*s1 + s2*s2));
        double a = s1*t;
        double b = s2*t;
        yy1 = a*y1_new + b*y2_new;
        yy1 = NGUtility::sign(x1.transpose()*yy1)*yy1;
        yy2 = -b*y1_new + a*y2_new;
        yy2 = NGUtility::sign(x2.transpose()*yy2)*yy2;
    }
    else{
        yy1 = y1_new;
        yy2 = y2_new;
    }

}
double CaliberBuilder::weigthvalue(Vec3d &dataL1, SVolume & L_XX3){
    auto nx = L_XX3.x();
    auto ny = L_XX3.y();
    auto nz = L_XX3.z();
    auto x = dataL1(0);
    auto y = dataL1(1);
    auto z = dataL1(2);
    //cout << "(" << std::round(x) << " " << std::round(y) << " " << std::round(z) << ")=" << L_XX3(std::round(x), std::round(y), std::round(z)) << endl;
    double dd = 0, ww = 0, idexx, idexy, idexz, w1;

    idexx = std::max(std::min(std::round(x), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);

    idexx = std::max(std::min(std::round(x + 1), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    idexx = std::max(std::min(std::round(x - 1), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    idexx = std::max(std::min(std::round(x), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y + 1), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    idexx = std::max(std::min(std::round(x), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y - 1), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    idexx = std::max(std::min(std::round(x), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z + 1), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    idexx = std::max(std::min(std::round(x), nx - 1.0), 0.0);
    idexy = std::max(std::min(std::round(y), ny - 1.0), 0.0);
    idexz = std::max(std::min(std::round(z - 1), nz - 1.0), 0.0);
    w1 = 0.05*std::max(-2 * ((x - idexx)*(x - idexx) + (y - idexy) * (y - idexy) + (2 * z - 2 * idexz)* (2 * z - 2 * idexz)), -6.0);
    dd = dd + std::exp(w1)*L_XX3(idexx, idexy, idexz);
    ww = ww + std::exp(w1);
    return dd / (ww + 0.0001);
}
void CaliberBuilder::SegmentpointModify(VectorMatXd &Segmentdata, MatXd &DircVecs, SVolume & OrigImag, VectorVec2d&ThresMatrix, CellVec7d &dataPointSet){
	std::cout << "Doing SegmentpointModify" << std::endl;
    auto Imsizex = OrigImag.x();
    auto Imsizey = OrigImag.y();
    auto Imsizez = OrigImag.z();
    auto Numz = Segmentdata.size();
    auto Numx = Segmentdata[0].rows();
    auto Numy = Segmentdata[0].cols();
    //SVolume ImageLabel;
    VectorVec7d dataPoint;
    //ImageLabel.SetSize(Imsizex, Imsizey, Imsizez);
    //ImageLabel.SetZero();
    commonLabel_.SetSize(Imsizex, Imsizey, Imsizez);
    commonLabel_.SetZero();
    auto cx = (Numx - 1) / 2;
    auto cy = (Numy - 1) / 2;
    for (int i = 0; i < Numz; i++){
        auto MatrixLabel = Segmentdata[i];
        //cout << "MatrixLabel="<<endl<< MatrixLabel<<endl;
        LabelImage(MatrixLabel);
        //ImageLabel.SetZero();
        commonLabel_.SetZero();
        //Vec2d Threv=1.0 / 3.0 * (ThresMatrix[std::max(i - 1, 0)] + ThresMatrix[i] + ThresMatrix[std::min(i + 1, Numz - 1)]);
        //cout << "MatrixLabel=" << endl << MatrixLabel << endl;
        for (int ii = 0; ii < Numx; ii++){
            for (int jj = 0; jj < Numy; jj++){
                if (MatrixLabel(ii, jj) > 0){
                    Vec3d p0 = DircVecs.col(i).middleRows(0, 3);
                    Vec3d px = DircVecs.col(i).middleRows(6, 3);
                    Vec3d py = DircVecs.col(i).middleRows(9, 3);
                    auto ii0 = 0.5*(ii - cx);
                    auto jj0 = 0.5*(jj - cy);
                    Vec3d pp0 = p0 + ii0*px + jj0*py;
                    pp0(0) = std::max(pp0(0), 0.0);
                    pp0(1) = std::max(pp0(1), 0.0);
                    pp0(2) = std::max(pp0(2), 0.0);
                    pp0(0) = std::min(std::round(pp0(0)), Imsizex - 1.0);
                    pp0(1) = std::min(std::round(pp0(1)), Imsizey - 1.0);
                    pp0(2) = std::min(std::round(pp0(2)), Imsizez - 1.0);
                    //if (ImageLabel(pp0(0), pp0(1), pp0(2)) == 0){
                    if (commonLabel_(pp0(0), pp0(1), pp0(2)) == 0){
                        Vec7d Point;
                        Point.middleRows(0, 3) = pp0;
                        Point(3) = OrigImag(pp0(0), pp0(1), pp0(2));
                        Point.middleRows(4, 2) = 1.0 / 3.0 * (ThresMatrix[std::max(i - 1, 0)] + ThresMatrix[i] + ThresMatrix[std::min(i + 1.0, Numz - 1.0)]);
                        Point(6) = MatrixLabel(ii, jj);
                        dataPoint.push_back(Point);
                        //ImageLabel(pp0(0), pp0(1), pp0(2)) = 1;
                        commonLabel_(pp0(0), pp0(1), pp0(2)) = 1;
                    }
                }
            }
        }
        dataPointSet.push_back(dataPoint);
        dataPoint.clear();
    }
}
void CaliberBuilder::LabelImage(MatXd&LabelImag){
    MatXd BinaryImage = LabelImag;
    auto Numx = BinaryImage.rows();
    auto Numy = BinaryImage.cols();
    LabelImag.setZero(Numx, Numy);
    VectorVec3d PointSet;
    VectorVec2d Datapoints, Seedp;
    int  ss = 0;
    for (int i = 0; i < Numx; i++){
        for (int j = 0; j < Numy; j++){
            Seedp.clear();
            Datapoints.clear();
            if (BinaryImage(i, j) == 1){
                Vec2d seed;
                seed(0) = i;
                seed(1) = j;
                Seedp.push_back(seed);
                LabelImageExtra(BinaryImage, Seedp, Datapoints);
                ss = ss + 1;
                for (int n = 0; n < Datapoints.size(); n++){
                    Vec3d tmp;
                    tmp.head(2) = Datapoints[n];
                    tmp(2) = ss;
                    PointSet.push_back(tmp);
                }
            }
        }
    }
    for (int ii = 0; ii < PointSet.size(); ii++){
        LabelImag(PointSet[ii](0), PointSet[ii](1)) = PointSet[ii](2);
    }
}
void CaliberBuilder::LabelImageExtra(MatXd& BinaryImage, VectorVec2d& Seedp, VectorVec2d&Datapoints){
    //auto Numx = BinaryImage.rows();
    //auto Numy = BinaryImage.cols();
    VectorVec2d datapoints;
    VectorVec2d Seedp_tmp = Seedp;
    while (true)
    {
        LabelImageExtra_other(BinaryImage, Seedp_tmp, datapoints);
        if (datapoints.size() > 0){
            Datapoints.insert(Datapoints.end(), datapoints.begin(), datapoints.end());
            Seedp_tmp = datapoints;
        }
        else
        {
            break;
        }
    }
}
void CaliberBuilder::LabelImageExtra_other(MatXd& BinaryImage, VectorVec2d& Seedp, VectorVec2d&Datapoints){
    auto Numx = BinaryImage.rows();
    auto Numy = BinaryImage.cols();
    auto Num1 = Seedp.size();
    Datapoints.clear();
    for (int i = 0; i < Num1; i++){
        auto datap = Seedp[i];
        for (int Indexx = std::max(datap(0) - 1, 0.0); Indexx <= std::min(datap(0) + 1, Numx - 1.0); Indexx++){
            for (int Indeyy = std::max(datap(1) - 1, 0.0); Indeyy <= std::min(datap(1) + 1, Numy - 1.0); Indeyy++){
                if (BinaryImage(Indexx, Indeyy) == 1){
                    Vec2d Currdata;
                    Currdata(0) = Indexx;
                    Currdata(1) = Indeyy;
                    Datapoints.push_back(Currdata);
                    BinaryImage(Indexx, Indeyy) = 0;
                }
            }
        }
    }
}

void CaliberBuilder::SegmentPointSceening(CellVec7d&dataPointSet, SVolume &ImagLabel, SVolume & OrigImage, int CurveLabel, VectorVec3d & dataSet, VecXd& Labelvec){
	
    auto Numss = dataPointSet.size();
    Labelvec = VecXd::Zero(Numss);
    VectorVec3d dataPointS, dataS;
    for (int i = 0; i < Numss; i++){
        int kkss = 0;
        auto Currdata = dataPointSet[i];
        VectorVec7d CurrdataNew;
        if (Currdata.size() > 0)
        {
            kkss = 0;
            int Num0 = Currdata[Currdata.size() - 1](6);
            for (int ii = 0; ii < Num0; ii++)
            {
                for (int n = 0; n < Currdata.size(); n++)
                {
                    if (Currdata[n](6) == ii + 1){
                        CurrdataNew.push_back(Currdata[n]);
                        Vec3d tmp = Currdata[n].middleRows(0, 3);
                        dataPointS.push_back(tmp);
                    }
                }
                int sst = SegmentPointScreeningSub1(dataPointS, ImagLabel, CurveLabel);
                if (sst == 1){
                    FilleddataPointModify(CurrdataNew, OrigImage, dataS);
                    dataSet.insert(dataSet.end(), dataS.begin(), dataS.end());
                    kkss = 1;
                }
                dataS.clear();
                dataPointS.clear();
                CurrdataNew.clear();
            }
            Labelvec(i) = kkss;
        }
    }
}
int  CaliberBuilder::SegmentPointScreeningSub1(VectorVec3d &dataPointS, SVolume & ImagLabel, int CurveLabel){
	
    auto Numss = dataPointS.size();
    auto Numx = ImagLabel.x();
    auto Numy = ImagLabel.y();
    auto Numz = ImagLabel.z();
    int sst = 0;
    for (int i = 0; i < Numss; i++){
        Vec3d p = dataPointS[i].middleRows(0, 3);
        for (int ii = std::min(std::max(p(0) - 1, 0.0), Numx - 1.0); ii <= std::min(std::max(p(0) + 1, 0.0), Numx - 1.0); ii++)
        {
            for (int jj = std::min(std::max(p(1) - 1, 0.0), Numy - 1.0); jj <= std::min(std::max(p(1) + 1, 0.0), Numy - 1.0); jj++)
            {
                for (int ij = std::min(std::max(p(2) - 1, 0.0), Numz - 1.0); ij <= std::min(std::max(p(2) + 1, 0.0), Numz - 1.0); ij++)
                {
                    if (ImagLabel(ii, jj, ij) == CurveLabel){
                        sst = 1;
                        break;
                    }

                }
                if (sst == 1)break;
            }
            if (sst == 1)break;
        }
        if (sst == 1)break;
    }
    return sst;
}
void CaliberBuilder::FilleddataPointModify(VectorVec7d &dataPoint, SVolume & OrigImag, VectorVec3d &dataPointS){
    auto Numx = OrigImag.x();
    auto Numy = OrigImag.y();
    auto Numz = OrigImag.z();
    //SVolume ImagLabel;
    //ImagLabel.SetSize(Numx, Numy, Numz);
    //ImagLabel.SetZero();
    commonLabel_.SetSize(Numx, Numy, Numz);
    commonLabel_.SetZero();
    auto NumPoint = dataPoint.size();

    for (int ii = 0; ii < NumPoint; ii++){
        //ImagLabel(dataPoint[ii](0), dataPoint[ii](1), dataPoint[ii](2)) = dataPoint[ii](3);
        commonLabel_(dataPoint[ii](0), dataPoint[ii](1), dataPoint[ii](2)) = dataPoint[ii](3);
    }
    for (int ii = 0; ii < NumPoint; ii++){

        for (int ik = std::min(std::max(dataPoint[ii](0) - 2, 0.0), Numx - 1.0); ik <= std::min(std::max(dataPoint[ii](0) + 2, 0.0), Numx - 1.0); ik++)
        {
            for (int jk = std::min(std::max(dataPoint[ii](1) - 2, 0.0), Numy - 1.0); jk <= std::min(std::max(dataPoint[ii](1) + 2, 0.0), Numy - 1.0); jk++)
            {
                for (int mk = std::min(std::max(dataPoint[ii](2) - 2, 0.0), Numz - 1.0); mk <= std::min(std::max(dataPoint[ii](2) + 2, 0.0), Numz - 1.0); mk++)
                {
                    double Thre = std::max(dataPoint[ii](3), 0.5*(dataPoint[ii](4) + dataPoint[ii](5)));
                    Thre = std::min(1.1*Thre, dataPoint[ii](4));
                    //cout << "OrigImag=" << OrigImag(ik, jk, mk) << "  ImagLabel=  " << ImagLabel(ik, jk, mk)<< endl;
                    //if (ImagLabel(ik, jk, mk) == 0 && OrigImag(ik, jk, mk) > Thre){
                    if (commonLabel_(ik, jk, mk) == 0 && OrigImag(ik, jk, mk) > Thre){
                        Vec3d Point;
                        Point(0) = ik;
                        Point(1) = jk;
                        Point(2) = mk;
                        dataPointS.push_back(Point);
                        //ImagLabel(ik, jk, mk) = 1;
                        commonLabel_(ik, jk, mk) = 1;
                    }
                }
            }
        }
    }
}
void CaliberBuilder::CurvesDircs(VectorVec3d &Curves, MatXd&DirecVecs){
	std::cout << "Doing CurvesDircs" << std::endl;
    auto NumCurves = Curves.size();
    DirecVecs = MatXd::Zero(12, NumCurves);
    VectorVec3d tmp;
    Vec3d x1, x2, x3;
    for (int i = 1; i < NumCurves - 1; i++){
        tmp.insert(tmp.begin(), Curves.begin() + std::max(i - 2, 0), Curves.begin() + std::min(i + 3.0, NumCurves+0.0));
        principald(tmp, x1);
        directionc(x1, x2, x3);
        DirecVecs.col(i).middleRows(0, 3) = Curves[i];
        DirecVecs.col(i).middleRows(3, 3) = x1;
        DirecVecs.col(i).middleRows(6, 3) = x2;
        DirecVecs.col(i).middleRows(9, 3) = x3;
        tmp.clear();
    }
    DirecVecs.col(0) = DirecVecs.col(1);
    DirecVecs.col(NumCurves - 1) = DirecVecs.col(NumCurves - 2);
}
//int CaliberBuilder::Sign(double n)
//{
//	if (0 == n)
//		return 0;
//	else if (n > 0)
//		return 1;
//	else
//		return -1;
//}
void CaliberBuilder::directionc(const Vec3d &x1, Vec3d &x2, Vec3d& x3){

    //Vec3d xx = { 1, 0, 0 };
    //Mat3d xt = x1*xx;
    //x2 = x1*xt;
    //x3 = x1*x2;
    x2 = x1;
    //[idxv,idexx]=sort(abs(x1));
    double temp;
    int n = 3;
    int i, j;
    Vec3d idxv = x1.cwiseAbs();
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
    x2[idexx[2]] = -NGUtility::sign(x1[idexx[2]])*cdd;
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