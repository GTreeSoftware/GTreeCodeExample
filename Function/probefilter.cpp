/*
 * Copyright (c)2013-2015  Zhou Hang, Shaoqun Zeng, Tingwei Quan
 * Britton Chance Center for Biomedical Photonics, Huazhong University of Science and Technology
 * All rights reserved.
 chenyijun zhushi 
 */
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <omp.h>
#include <Eigen/SVD>
#ifdef _WIN32
#include <ctime>
#define M_E        2.71828182845904523536
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#else
#include <sys/time.h>
#endif
#include "probefilter.h"
#include "../ngtypes/volume.h"
#include "../ngtypes/soma.h"
#include "volumealgo.h"
#include "contourutil.h"
#include "Function/Trace/traceutil.h"
#include "Function/NGUtility.h"

//struct Vec4d_3th_less{
//    bool operator() (const Vec4d& lhs, const Vec4d& rhs){
//		return lhs(3) < rhs(3);
//	}
//};
//
//struct Vec4d_3th_great{
//    bool operator() (const Vec4d& lhs, const Vec4d& rhs){
//		return lhs(3) > rhs(3);
//	}
//};
//
//struct Vec3i_less{
//	bool operator() (const Vec3i& lhs, const Vec3i &rhs){
//		if (lhs(0) != rhs(0))  return lhs(0) < rhs(0);
//		else if (lhs(1) != rhs(1))  return lhs(1) < rhs(1);
//		else if (lhs(2) != rhs(2))  return lhs(2) < rhs(2);
//		return false;
//	}
//};


ProbeFilter::ProbeFilter()
{
    className_ = std::string("ProbFilter");
	m_Source = std::shared_ptr<Soma>(new Soma(this));
    
    threadNum = 4;
    k_theta = 40;
    k_phi = 20;
	slice_ = 0.5;
    minRadius_ = 4;
    minConnectNum_ = 100;
    //-------2014-5-15--------//
    maxConnectNum_ = 4000;
    maxBoundingBoxSize_ = 125000;//50*50*50
}

ProbeFilter::~ProbeFilter()
{

}

ProcStatPointer ProbeFilter::Update()
{
    if (m_Input->GetIdentifyName() != std::string("Volume") || m_BinImg->GetIdentifyName() != std::string("Volume")){
        MAKEPROCESSSTATUS(resSta, false, className_, "no input data.");
        return resSta;
    }
	std::shared_ptr<const SVolume> origImg = std::dynamic_pointer_cast<const SVolume>(m_Input);
	std::shared_ptr<const CVolume> binImg = std::dynamic_pointer_cast<const CVolume>(m_BinImg);
	std::shared_ptr<Soma> m_Soma = std::dynamic_pointer_cast<Soma>(m_Source);

    int width	= origImg->x();
    int height = origImg->y();
    int frame	= origImg->z();

    xVoxel_ = origImg->XResolution();//
    yVoxel_ = origImg->YResolution();//
    zVoxel_ = origImg->ZResolution();//

    minConnectNum_ = NGUtility::Round(4.0 * minRadius_ * minRadius_ * minRadius_);//2015-5-7//Round(double(50) * 2.0 / (xVoxel_*yVoxel_*zVoxel_));//TODO
    maxConnectNum_ = NGUtility::Round(double(4000) * 2.0 / (xVoxel_*yVoxel_*zVoxel_));

    maxBoundingBoxSize_ = NGUtility::Round(double(maxBoundingBoxSize_) * 2.0 / (xVoxel_*yVoxel_*zVoxel_));
	std::shared_ptr<CVolume> erobinImg = std::tr1::shared_ptr<CVolume>(new CVolume());
    erobinImg->SetSize(width, height, frame);

	std::shared_ptr<CVolume> tmpBinImg = std::tr1::shared_ptr<CVolume>(new CVolume());
    tmpBinImg->QuickCopy(*binImg);


    VectorVec3i	 erodedPtSet_;
    std::vector<int> connectPtSetNum_;
    VectorVec3i	 connectPtSet_;

    VectorVec3i tmpBinPtSet(*m_BinPtSet);//copy
    VectorVec3i tmpConnectPtSet;
    std::vector<int> tmpConnectNum;
    const int erodeProcessNum = 50;

#ifdef _WIN32
    clock_t beg = clock();
#else
    timeval start1, end1;
    gettimeofday(&start1, 0);
#endif

    /*50 erode operations, extract all connected domain*/
    for (int i = 1; i < erodeProcessNum; ++i){//101
        printf("%d erode process\n", i);

        Corrode(*tmpBinImg, tmpBinPtSet, double(i), *erobinImg, erodedPtSet_);
        GetConnectedDomain(*erobinImg, erodedPtSet_,
            *tmpBinImg, tmpBinPtSet, tmpConnectPtSet, tmpConnectNum);
        //
        VectorVec3i::size_type sz1 = connectPtSetNum_.size();
        VectorVec3i::size_type sz2 = connectPtSet_.size();
        connectPtSetNum_.resize(sz1 + tmpConnectNum.size());
        connectPtSet_.resize(sz2 + tmpConnectPtSet.size());
        std::copy(tmpConnectNum.begin(), tmpConnectNum.end(), connectPtSetNum_.begin() + sz1);
        std::copy(tmpConnectPtSet.begin(), tmpConnectPtSet.end(), connectPtSet_.begin() + sz2);
        printf("  connect-area new : %d\n", (int)tmpConnectNum.size() );
        printf("  connect-area sum : %d\n", (int)connectPtSetNum_.size() );

        if (tmpBinPtSet.empty())  break;
    }
    /*localize soma*/
    if (!connectPtSetNum_.empty()){
        SomaLocalize(*origImg, connectPtSet_, connectPtSetNum_, *m_Soma);
    }

    //multiply the z resolution
    //xy resolution are both turn into 1 um
    /*Cell tmpCell;
    for(size_t i = 0 ; i < m_Soma->size(); ++i){
        tmpCell = m_Soma->GetCell(i);
        tmpCell.x *= xVoxel_;
        tmpCell.y *= yVoxel_;
        tmpCell.z *= zVoxel_;
        m_Soma->GetCell(i) = tmpCell;
    }*/
    //move out to add offset.

    /*clear*/
    tmpBinPtSet.clear();
    tmpConnectPtSet.clear();
    tmpConnectNum.clear();
    //Clear();

    ///=================================
    printf("Finding Soma Completed.\n");

    ///=============time=============///
    double time;
#ifdef _WIN32
    clock_t end = clock();
    time = double(end - beg);
#else
    gettimeofday(&end1, 0);
    time = 1000000*(end1.tv_sec-start1.tv_sec)+end1.tv_usec-start1.tv_usec;
    time /= 1000;
#endif
//    double time=(double)(end-start)/CLOCKS_PER_SEC;
    printf("Recognize Soma used %lf ms.\n", time);
    printf("there are %d soma.\n", (int)m_Soma->size() );
    MAKEPROCESSSTATUS(resSta, true, className_, "");
    return resSta;
}

ConstIDataPointer ProbeFilter::GetOutput()
{
	if(!m_Source) m_Source = std::tr1::shared_ptr<Soma>(new Soma(this));
    return m_Source;
}

IDataPointer ProbeFilter::ReleaseData()
{
    m_Source->ReleaseProcessObject();
    IDataPointer tData(m_Source);
    m_Source.reset();
    return tData;
}

void ProbeFilter::SetThreadNum(int t)
{
    threadNum = t;
}

void ProbeFilter::SetInputBin(ConstIDataPointer p)
{
    m_BinImg = p;
}

void ProbeFilter::SetInputBinPtSet(BinPtSetPointer p)
{
    m_BinPtSet = p;
}

void ProbeFilter::SetConnectDomainRange(int m1, int m2)
{
    if(m2 <= m1) return;
    minConnectNum_ = m1;
    maxConnectNum_ = m2;
}

void ProbeFilter::SetMinRadius(double arg)
{
    minRadius_ = arg;
}

double ProbeFilter::GetMinRadius() const
{
    return minRadius_;
}

void ProbeFilter::ExtractLocalDomain(const Vec3d &initSoma, const SVolume &v, SVolume &locOrigImg, Vec3d &locSoma)
{
    int vx = v.x();
    int vy = v.y();
    int vz = v.z();
    int minX = std::max((int)initSoma(0) - 30, 0);
    int maxX = std::min((int)initSoma(0) + 30, vx - 1);
    int minY = std::max((int)initSoma(1) - 30, 0);
    int maxY = std::min((int)initSoma(1) + 30, vy - 1);
    int minZ = std::max((int)initSoma(2) - 15, 0);//30
    int maxZ = std::min((int)initSoma(2) + 15, vz - 1);
    locSoma = Vec3d(initSoma(0) - minX, initSoma(1) - minY, initSoma(2) - minZ);
    int subx = maxX - minX + 1;
    int suby = maxY - minY + 1;
    int subz = maxZ - minZ + 1;
    locOrigImg.SetSize(subx, suby, subz);
    const SVolume &pptr = v;//&v.Image();
    for (int i = minX; i <= maxX; ++i){
        for (int j = minY; j <= maxY; ++j){
            for (int ij = minZ; ij <= maxZ; ++ij)
                locOrigImg(i - minX, j - minY, ij - minZ) = pptr(i, j, ij);
        }
    }//for
}

void ProbeFilter::GetSubDomain(const Vec3i &initPoint, VectorVec3i &connectPtSet, CVolume &eroImage)
{
    int i = 0;
    connectPtSet.clear();//
    VectorVec3i tempPos1;		//point set
    VectorVec3i tempPos2;		//tmp data
    tempPos2.push_back(initPoint);			//
    ///----------------initialize center point -----------------
    eroImage( initPoint(0) , initPoint(1) ,  initPoint(2) ) = 0;
    connectPtSet.push_back(initPoint);//
    tempPos1.push_back(initPoint);//for region growing
    ///----------------start ---------------------
    while(!tempPos2.empty() && i < 1280 ){
        ++i;
        //////////--- region growing  ---//////////////////
        SubDomainExpand( tempPos1, tempPos2, eroImage);
        //resize
        VectorVec3i::size_type sz = connectPtSet.size();
        connectPtSet.resize(sz + tempPos2.size());
        std::copy(tempPos2.begin(), tempPos2.end(), connectPtSet.begin() + sz);
        tempPos1 = tempPos2;
    }
}

void ProbeFilter::SubDomainExpand(const VectorVec3i &initPoints, VectorVec3i &connectPtSet, CVolume &eroImage)
{
    int i,j,ij;
    int n, t_Sum;;
    int xMin,xMax;
    int yMin,yMax;
    int zMin,zMax;
    int nx	= eroImage.x();
    int ny	= eroImage.y();
    int nz	= eroImage.z();
    int nxx	= initPoints.size();

    Vec3i pt;
    pt.setZero();
    connectPtSet.clear();//must clear
    ///-------------------------------------
    for (n = 0 ; n < nxx; ++n){
        const Vec3i& it = initPoints[n];
        int id1 = it(0);
        int id2 = it(1);
        int id3 = it(2);

        if (id1 > 1 && id1 < nx-2 &&id2 > 1 && id2 < ny-2 && id3 > 1 && id3 < nz-2){
            xMin = id1 -1; xMax = id1 +1;
            yMin = id2 -1; yMax = id2 +1;
            zMin = id3 -1; zMax = id3 +1;

            //cal connective domain signal
            t_Sum = 0;
            for (i = xMin; i <= xMax; ++i)
                for (ij = zMin; ij <= zMax; ++ij)
                    for (j = yMin; j <= yMax; ++j){
                        t_Sum += eroImage(i , j , ij );
                    }
            t_Sum /= 255;

            //extract
            if (t_Sum > 4){
                for (i = xMin; i <= xMax; ++i)
                    for (j = yMin; j <= yMax; ++j)
                        for (ij = zMin; ij <= zMax; ++ij)
                            if (eroImage(i ,j , ij ) > 100){
                                eroImage(i ,j , ij ) = 0;
                                pt << i, j, ij;
                                connectPtSet.push_back(pt);
                            }
                        }
        }
    }//for
}

bool ProbeFilter::IsCxDomainBoundingBoxTooLarge(const VectorVec3i &locPtSet)
{
    int xMin, xMax, yMin, yMax, zMin, zMax;
    CalculateCxDomainBoundingBox(locPtSet, xMin, xMax, yMin, yMax, zMin, zMax);
    //------------2014-5-15--------------------//
    if ( (xMax - xMin + 1) * (yMax - yMin + 1) * (zMax - zMin + 1) < maxBoundingBoxSize_ )
        return true;
    else return false;
}

void ProbeFilter::CalculateCxDomainBoundingBox(const VectorVec3i &locPtSet, int &xMin, int &xMax, int &yMin, int &yMax, int &zMin, int &zMax)
{
    xMin = yMin = zMin = 100000;
    xMax = yMax = zMax =  0;

    for(int i = 0; i < (int)locPtSet.size(); ++i){
        const Vec3i& it = locPtSet[i];
        if (it(0) < xMin) xMin = it(0);
        if (it(0) > xMax) xMax = it(0);
        if (it(1) < yMin) yMin = it(1);
        if (it(1) > yMax) yMax = it(1);
        if (it(2) < zMin) zMin = it(2);
        if (it(2) > zMax) zMax = it(2);
    }
}



bool ProbeFilter::Corrode(const CVolume &binImage, const VectorVec3i &binPtSet,
                          const double eroIntensity,
                          CVolume &eroBinImg, VectorVec3i &eroPtSet)
{
    int xMin,xMax;
    int yMin,yMax;
    int zMin,zMax;
    int binDomainSum;

    //int nXX = ero_pt_set.size();
    int nx	= binImage.x();
    int ny	= binImage.y();
    int nz	= binImage.z();

    eroBinImg.SetZero();
    eroPtSet.clear();
///-------------------erode ---------------------------
    //VectorVec3i::const_iterator itend = binPtSet.end();
    //for (VectorVec3i::const_iterator it = binPtSet.begin(); it != itend; ++it)
    for(size_t i = 0; i < binPtSet.size(); ++i){
        const Vec3i& it = binPtSet[i];
        binDomainSum = 0;
        Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, it(0), it(1), it(2),
            1, 1, 1,
            0, nx-1, 0, ny-1, 0, nz-1);
        ///-----------------neighborhood--------------------------------
        GetAreaSum(binImage, xMin, xMax, yMin, yMax, zMin, zMax, binDomainSum);
        binDomainSum /= 255;
        ///--------------------judge---------------------------
        if (binDomainSum > int(27.0f*( 1.0f/3.0f + 0.01f * eroIntensity))) {
            eroBinImg( it(0), it(1) , it(2)) = 255;
            eroPtSet.push_back(it);
        }
    }
    if (eroPtSet.empty()) return false;
    return true;
}

void ProbeFilter::GetConnectedDomain(const CVolume &eroImage, const VectorVec3i &eroPtSet,
                                     CVolume &remainEroImg, VectorVec3i &remainBinPtSet,
                                     VectorVec3i &connectPtSet, std::vector<int> &connectNum)
{
   int n;
   int xMin,xMax;
   int yMin,yMax;
   int zMin,zMax;
   int ax = eroImage.x();
   int ay = eroImage.y();
   int az = eroImage.z();
   int nxx = int(eroPtSet.size());
   int sum = 0;

   remainBinPtSet.clear();
   remainEroImg.SetZero();
   CVolume tmpEroImg;
   //Quickstd::copyVoxel(eroImage, tmpEroImg;);
   tmpEroImg.QuickCopy(eroImage);
   VectorVec3i curExtractPtSet;
   connectPtSet.clear();
   connectNum.clear();
    ///--------------start---------------
   Vec3i pt;
   for (n = 0; n < nxx; ++n){
       pt(0) = eroPtSet[n](0);
       pt(1) = eroPtSet[n](1);
       pt(2) = eroPtSet[n](2);
       /////////////////////-- get seeds --/////////////////////////////
       if (pt(0) > 1 && pt(0) < ax-1 && pt(1) > 1 && pt(1) < ay-1 && pt(2) > 1 && pt(2) < az-1){
           sum = 0;
           Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, pt(0), pt(1), pt(2),
                       1,1,1);
           GetAreaSum(tmpEroImg, xMin, xMax, yMin, yMax, zMin, zMax, sum);
           sum /= 255;
       }
       else  sum = 0;

       ///----------------extract---------------------
       curExtractPtSet.clear();//warning
       if (sum > 20){
           /*extract area*/
           GetSubDomain( pt , curExtractPtSet, tmpEroImg);
           /*whether connective domain, less than given threshold, limited bounding box*/
           if (minConnectNum_ < (int)curExtractPtSet.size() && maxConnectNum_ > (int)curExtractPtSet.size()
               && IsCxDomainBoundingBoxTooLarge(curExtractPtSet) ){
               connectNum.push_back(curExtractPtSet.size());
               VectorVec3i::size_type sz = connectPtSet.size();
               connectPtSet.resize(sz + curExtractPtSet.size());
               std::copy(curExtractPtSet.begin(), curExtractPtSet.end(), connectPtSet.begin() + sz);
           }
           /*update next erode image*/
           else{
               VectorVec3i::size_type sz = remainBinPtSet.size();
               remainBinPtSet.resize(sz + curExtractPtSet.size());
               std::copy(curExtractPtSet.begin(), curExtractPtSet.end(), remainBinPtSet.begin() + sz);
               int num_cur_extract_pt_set = curExtractPtSet.size();
               for (VectorVec3i::size_type ii = 0; ii < (VectorVec3i::size_type)num_cur_extract_pt_set; ++ii){
                   remainEroImg( curExtractPtSet[ii](0), curExtractPtSet[ii](1),
                       curExtractPtSet[ii](2) ) = 255;
               }

           }
       }
   }
}

void ProbeFilter::SomaLocalize(const SVolume &v, const VectorVec3i &connectPtSet, const std::vector<int> &connectNum,
                               Soma &resultSoma)
{
    ///---------------------initialize ---------------------------------
    int i, j;//openMP use int, not size_t
    resultSoma.clear();
#ifdef _WIN32
	std::tr1::shared_ptr<const CVolume> binImg = std::tr1::dynamic_pointer_cast<const CVolume>(m_BinImg);
#else
	std::shared_ptr<const CVolume> binImg = std::dynamic_pointer_cast<const CVolume>(m_BinImg);
#endif
    
    /*the initial position of N-th connective domain*/
    std::vector<int> conPtSetPosition(connectNum.size()+1);
    conPtSetPosition[0] = 0;
    std::partial_sum(connectNum.begin(), connectNum.end(), conPtSetPosition.begin()+1);
    int conPtSetSum = int(conPtSetPosition.size());
    int finalSomaIndex = 0;
    omp_set_num_threads(threadNum);
#pragma omp parallel for private(i,j) schedule(dynamic)
    for ( i = 0; i < conPtSetSum -1; ++i ){//TODO
        //--------test----2014-6-3---//
        //printf("%d\n",i);
        VectorVec3i innerPtSet;
        //Vec3d initSomaForward;
        //////////////////////////////////
        CVolume	locBinImg;	//
        SVolume	locOrigImg;	//local original image

        /*get singla connective domain*/
        VectorVec3i curConPtSet;//
        curConPtSet.assign(connectPtSet.begin() + conPtSetPosition[i],
            connectPtSet.begin() + conPtSetPosition[i+1] );

        VectorVec4d locSeed;//private data
        VectorVec3d globalSeed;//private data

        /*get candidate seeds*/
        CreateSeeds(curConPtSet, v, locBinImg, locOrigImg, locSeed, globalSeed);
////---draw connect domain------//
//        for(size_t kk = 0; kk < curConPtSet.size();++kk){
//            globalSeed.push_back(Vec3d(curConPtSet[kk](0),curConPtSet[kk](1),curConPtSet[kk](2)));
//            locSeed.push_back(Vec4d(curConPtSet[kk](0),curConPtSet[kk](1),curConPtSet[kk](2),2.0));
//        }

        VectorVec4d::size_type seedSum = locSeed.size();
        if (seedSum > 0 && locSeed[0](0) > -1){
            VectorVec3d locPtSet;
            for (j = 0; j < (int)seedSum; ++j){
                locPtSet.push_back(Vec3d(locSeed[j](0), locSeed[j](1), locSeed[j](2)));
            }

            VectorVec3d resultLocPtSet;//new position of seeds
            std::vector<double> radius;//the radii of seeds

            _SomaRecognition(locBinImg, locPtSet, locOrigImg, resultLocPtSet, radius);

			locBinImg.SetZero();//save locBinImg for shape reflection
            locOrigImg.SetSize(0,0,0);

            VectorVec3d finPos;
            for(VectorVec4d::size_type n = 0; n < locSeed.size(); ++n){
                finPos.push_back(Vec3d(resultLocPtSet[n](0) - locSeed[n](0),
                                        resultLocPtSet[n](1) - locSeed[n](1),
                                        resultLocPtSet[n](2) - locSeed[n](2) ));//the offset distance
            }
            seedSum = radius.size();

			//screening little somas
			Vec3d tmpSet = globalSeed[0] - locSeed[0].block(0,0,3,1);
            Vec3i offSet(NGUtility::Round(tmpSet(0)), NGUtility::Round(tmpSet(1)), NGUtility::Round(tmpSet(2)));
			VectorVec4d gloSomaStack;
			for(size_t ij = 0; ij < locSeed.size(); ++ij){
                if(radius[ij] >= minRadius_ ){
                    gloSomaStack.push_back(Vec4d(
                        globalSeed[ij](0) + finPos[ij](0),
                        globalSeed[ij](1) + finPos[ij](1),
                        globalSeed[ij](2) + finPos[ij](2),
                        radius[ij]));
                }
			}
            if (gloSomaStack.empty()) continue;
            /**** 
                the following is Matlab function SomaShapeFeedback
            *****/
			std::vector<int> flagYN(gloSomaStack.size(), 1);
            //contour , surface , volume
            //double meanUx(0.0);
            //double majAxis(0.0);
            //double minAxis(0.0);
            //double sur(0.0);
            //double volume(0.0);
            //std::vector<double> somaMeanRadius;//(seedSum, 0.0);//平均半径
            //std::vector<double> somaSurface;//(seedSum, 0.0);
            //std::vector<double> somaVolume;//(seedSum, 0.0);
            //std::vector<double> somaMajAxis;//(seedSum, 0.0);
            //std::vector<double> somaMinAxis;//(seedSum, 0.0);
//            int s = 0;
            double minDiameter = minRadius_ / slice_;// used to identify columnar
            std::vector<VectorVec3i> interPtSets;
			std::vector<std::vector<std::vector<double> > > maxRayLengths;
			std::vector<std::vector<VectorVec3i> > localConnectPtSets;
			std::vector<double> backThrevs;

            for(size_t somaIndex = 0; somaIndex < gloSomaStack.size(); ++somaIndex){
                Vec3d initSoma = gloSomaStack[somaIndex].block(0,0,3,1);
                double backThrev;
                VectorVec3i curInterPtSet;
                std::vector<VectorVec3i> localConnectPtSet;
                MatXd rayLength;
                std::vector<std::vector<double> > maxRayLength;
                SomaShape(initSoma, v, *binImg, curInterPtSet, rayLength, maxRayLength, localConnectPtSet, backThrev);
				//cal soma info //2015-5-7
				//MatXd newRayLength;
				//CorrectRayMatrix(rayLength, k_theta, k_phi, newRayLength);
				//meanUx = newRayLength.mean();
				//VectorVec3d contourPtSet;
				//VectorVec3d contourPatchPtSet;
				//MatXd resultRayLength;
				//GetContour(initSoma, newRayLength, k_theta, k_phi, contourPtSet);//轮廓
				//GetAxis(newRayLength, majAxis, minAxis);//长短轴
				////const int k_theta = 40, k_phi = 20;
				//ContourUtil::CalAreaVolume(contourPtSet, k_theta, k_phi,
				//	contourPatchPtSet, sur, volume);//
				//somaMeanRadius.push_back(meanUx);
				//somaSurface.push_back(sur);
				//somaVolume.push_back(volume);
				//somaMajAxis.push_back(majAxis);
				//somaMinAxis.push_back(minAxis);
                //CalSomaInfo(initSoma, smoothRay, majAxis, minAxis, sur, volume);
                interPtSets.push_back(VectorVec3i());
                interPtSets.back().swap(curInterPtSet);
				maxRayLengths.push_back(std::vector<std::vector<double> >());
				maxRayLengths.back().swap(maxRayLength);
				localConnectPtSets.push_back(std::vector<VectorVec3i>());
				localConnectPtSets.back().swap(localConnectPtSet);
				backThrevs.push_back(backThrev);
			}

			for(size_t somaIndex = 0; somaIndex < gloSomaStack.size(); ++somaIndex){
				//misty area
                std::vector<double> tmpSortRayLen;
                for (size_t ii = 0; ii < maxRayLengths[somaIndex].size(); ++ii) {
                    for (size_t jj = 0; jj < maxRayLengths[somaIndex][ii].size(); ++jj) {
                        tmpSortRayLen.push_back(maxRayLengths[somaIndex][ii][jj]);
                    }
                }
                std::sort(tmpSortRayLen.begin(), tmpSortRayLen.end());
                tmpSortRayLen.erase(tmpSortRayLen.begin() + NGUtility::Round(0.6 * tmpSortRayLen.size()), tmpSortRayLen.end());
                // get mean value.
                double tmpSortRayLenMean =
                    std::accumulate(tmpSortRayLen.begin(), tmpSortRayLen.end(), 0.0)
                    / double(tmpSortRayLen.size());
                //standard variation---------------------------------------
                double tmpSortRayLenMeanSqrtSum(0);
                int numTmpSortRayLen = tmpSortRayLen.size();
                for (int ii = 0; ii< numTmpSortRayLen; ++ii)
                {
                    tmpSortRayLenMeanSqrtSum +=
                        (tmpSortRayLen[ii] - tmpSortRayLenMean) * (tmpSortRayLen[ii] - tmpSortRayLenMean);
                }
                double tmpSortRayLenStd = std::sqrt(tmpSortRayLenMeanSqrtSum / (tmpSortRayLen.size() - 1));

                //
                if (tmpSortRayLenMean - tmpSortRayLenStd < minDiameter) {//XXX 
                    flagYN[somaIndex] = 0;//Misty area
                } else{

                    //int xLower, xUpper, yLower, yUpper, zLower, zUpper;
                    //CalculateCxDomainBoundingBox(curConPtSet, xLower, xUpper, yLower, yUpper, zLower, zUpper);
                    /*double rangeX = xUpper - xLower;
                    double rangeY = yUpper - yLower;
                    double rangeZ = zUpper - zLower;
                    double maxXYZ = std::max( rangeZ, std::max(rangeX, rangeY));*/
                    //if (maxXYZ > minRadius_ * 20 && std::accumulate(flagYN.begin(), flagYN.end(), 0) > 0) {
					//if (maxXYZ > 60.0 / (xVoxel_ * yVoxel_ * zVoxel_) && std::accumulate(flagYN.begin(), flagYN.end(), 0) > 0) {
					size_t sumSomaNodes = 0;
					for (size_t kk = 0; kk < interPtSets.size(); ++kk) {
						sumSomaNodes += interPtSets.size();
					}
					if (sumSomaNodes < curConPtSet.size() * 0.9 && std::accumulate(flagYN.begin(), flagYN.end(), 0) > 0) {
                        size_t maxConPtSetNum = 0;
                        size_t idMax = 10000000;
                        //just find the largest connect domain
                        for (size_t ii = 0; ii < localConnectPtSets[somaIndex].size(); ++ii) {
                            if(localConnectPtSets[somaIndex][ii].size() > maxConPtSetNum ) {
                                maxConPtSetNum = localConnectPtSets[somaIndex][ii].size();
                                idMax = ii;
                            }
                        }
                        //the neighborhood area is large enough.
                        double curRadius = gloSomaStack[somaIndex](3);
                        if (double(maxConPtSetNum) > 2 * curRadius * curRadius  ) {
                            //soma center
                            Vec3i curCenter(0,0,0);
                            /*for (size_t ii = 0; ii < curInterPtSet.size(); ++ii) {
                                curCenter += curInterPtSet[ii];
                            }
                            curCenter /= curInterPtSet.size();*/
                            for (size_t ii = 0; ii < interPtSets[somaIndex].size(); ++ii) {
                                curCenter += interPtSets[somaIndex][ii];
                            }
                            curCenter /= interPtSets[somaIndex].size();
                            //connect domain center
                            Vec3i curMaxConCenter(0,0,0);
                            VectorVec3i& maxConPtSet = localConnectPtSets[somaIndex][idMax];
                            for (size_t ii = 0; ii < maxConPtSet.size(); ++ii) {
                                curMaxConCenter += maxConPtSet[ii];
                            }
                            curMaxConCenter /= maxConPtSet.size();
                            Vec3d dir;
                            dir << curMaxConCenter(0) - curCenter(0), curMaxConCenter(1) - curCenter(1), curMaxConCenter(2) - curCenter(2);
                            if (dir.norm() < 0.001) dir << 0,0,0;
                            else dir.normalize();
                            //forward test
                            int forwardStep=0;
                            int backwardStep=0;
                            Vec3i curPt(0,0,0);
                            for (size_t ii = 0; ii <= 6 * curRadius; ++ii) {
                                curPt << curCenter(0) + ii * dir(0), curCenter(1) + ii * dir(1), curCenter(2) + ii * dir(2);
                                curPt(0) = std::max(0, std::min(v.x() - 1, curPt(0)));
                                curPt(1) = std::max(0, std::min(v.y() - 1, curPt(1)));
                                curPt(2) = std::max(0, std::min(v.z() - 1, curPt(2)));
                                if( v(curPt(0), curPt(1), curPt(2)) < backThrevs[somaIndex] ) break;
                                else ++forwardStep;
                            }
                            //backward test
                            for (size_t ii = 0; ii <= 6 * curRadius; ++ii) {
                                curPt << curCenter(0) - ii * dir(0), curCenter(1) - ii * dir(1), curCenter(2) - ii * dir(2);
                                curPt(0) = std::max(0, std::min(v.x() - 1, curPt(0)));
                                curPt(1) = std::max(0, std::min(v.y() - 1, curPt(1)));
                                curPt(2) = std::max(0, std::min(v.z() - 1, curPt(2)));
                                if( v(curPt(0), curPt(1), curPt(2)) < backThrevs[somaIndex] ) break;
                                else ++backwardStep;
                            }
                            //if( std::max(backwardStep, forwardStep) > 4 * curRadius) flagYN[somaIndex] = 1;//TODO
							if(  backwardStep + forwardStep  > 8 * curRadius) flagYN[somaIndex] = 0;//TODO
                        }
                    }
                }
            }

            //remove redundant soma
            if (std::accumulate(flagYN.begin(), flagYN.end(), 0) > 1) {
                for(size_t ii = 0; ii < gloSomaStack.size() - 1; ++ii){
                    for(size_t jj = ii+1; jj < gloSomaStack.size() ; ++jj){
                        if( 1 == flagYN[ii] && 1 == flagYN[jj] ){
                            double overlapRate = GetOverlapRate(interPtSets[ii], interPtSets[jj]);
                            if (overlapRate > 0.8) {
                                flagYN[jj] = 0;
                            }
                        }
                    }
                }
            }
            //save true soma
            VectorVec4d finalSoma;
            std::vector<double> finalSomaMeanRadius;//(seedSum, 0.0);
			std::vector<double> finalSomaSurface;//(seedSum, 0.0);
			std::vector<double> finalSomaVolume;//(seedSum, 0.0);
			std::vector<double> finalSomaMajAxis;//(seedSum, 0.0);
			std::vector<double> finalSomaMinAxis;//(seedSum, 0.0);
            for (size_t somaIndex = 0; somaIndex < flagYN.size(); ++somaIndex) {
				if( 1 == flagYN[somaIndex]) {//
					finalSoma.push_back(gloSomaStack[somaIndex]);
					//2015-5-7
					finalSomaMeanRadius.push_back(0);
					finalSomaSurface.push_back(0);
					finalSomaVolume.push_back(0);
					finalSomaMajAxis.push_back(0);
					finalSomaMinAxis.push_back(0);
					/*finalSomaMeanRadius.push_back(somaMeanRadius[somaIndex]);
					finalSomaSurface.push_back(somaSurface[somaIndex]);
					finalSomaVolume.push_back(somaVolume[somaIndex]);
					finalSomaMajAxis.push_back(somaMajAxis[somaIndex]);
					finalSomaMinAxis.push_back(somaMinAxis[somaIndex]);*/
				}
                //if( 1 ) finalSoma.push_back(gloSomaStack[somaIndex]);
            }
            
#pragma omp critical
            {
                //2015-1-5
				if(!finalSoma.empty()){
                    
                    for (size_t k = 0; k < finalSoma.size(); ++k){
                    		resultSoma.push_back(	Cell(finalSomaIndex++,
                    			finalSoma[k](0),finalSoma[k](1),finalSoma[k](2), finalSoma[k](3), 
                    			0, 0,                 
                    			finalSomaMeanRadius[k], finalSomaSurface[k], finalSomaVolume[k])//Cell
                    			);
                    }
					printf("Thread %d : find %lu somas.\n", omp_get_thread_num(), finalSoma.size());//线程号
				}
            }
        }
    }
}

void ProbeFilter::CreateSeeds(const VectorVec3i &curConPtSet, const SVolume &origImage,
                              CVolume &dstBinImg, SVolume &dstOrigImg, VectorVec4d &locSeed, VectorVec3d &globalSeed)
{
    locSeed.clear();
    globalSeed.clear();

    int i,j,ij;
    int xMin(100000), xMax(0);
    int yMin(100000), yMax(0);
    int zMin(100000), zMax(0);
    int nx,ny,nz;
    int nxx,nyy,nzz,nXY;
    int maxNum = 0;

    VectorVec4d tmpGlobalNode;

    int *tmpBinImg;
    int *tmpOrigImg;
    ///--------------------bounding box------------------------
    CalculateCxDomainBoundingBox(curConPtSet, xMin, xMax, yMin, yMax, zMin, zMax);

    nx = xMax - xMin + 1;
    ny = yMax - yMin + 1;
    nz = zMax - zMin + 1;
    ///============================
    nxx = nx + 6;//padarray 3 pixels
    nyy = ny + 6;
    nzz = nz + 6;
    nXY = nxx * nyy;
    ///----------------------------------------------------------------------
    dstBinImg.SetSize(nxx,nyy,nzz);
    dstOrigImg.SetSize(nxx,nyy,nzz);
    ///-----------------------extract local signal
    for(VectorVec3i::size_type n = 0; n < curConPtSet.size(); ++n){
        const Vec3i& it = curConPtSet[n];
        dstBinImg(it(0)-xMin+3 , it(1) - yMin +3 , it(2) - zMin + 3 ) = 1;//warning it is 1 rather than 255
        dstOrigImg(it(0)-xMin+3 , it(1) - yMin +3 , it(2) - zMin + 3 ) =
            origImage(it(0) , it(1) , it(2) );
    }

    ///=====================================
    //tmpBinImg.SetSize(nxx, nyy, nzz);
    //tmpOrigImg.SetSize(nxx, nyy, nzz);
    tmpBinImg = new int[nxx * nyy * nzz];
    memset(tmpBinImg, 0, sizeof(int) * nxx * nyy * nzz);
    tmpOrigImg = new int[nxx * nyy * nzz];
    memset(tmpOrigImg, 0, sizeof(int) * nxx * nyy * nzz);
    ///============special convolution===============//
    int xLen = CalcConv3dRadius(xVoxel_);
    int yLen = CalcConv3dRadius(yVoxel_);
    int zLen = CalcConv3dRadius(zVoxel_);
    //int xLen = Round(3.0 / xVoxel_ );
    //int yLen = Round(3.0 / yVoxel_ );
    //int zLen = Round(3.0 / zVoxel_ );
    Conv3d(tmpBinImg,dstBinImg, xLen, yLen,zLen);
    Conv3d(tmpOrigImg,dstOrigImg, xLen, yLen,zLen);

    /*get seeds , which are the peak of neighborhood*/
    maxNum = 0;
    int seedDistance =  2;//???
    //int seedDistance = 4;
    //--------------2014-5-15----------------------//
    int binValue = 0.3*double( (2*xLen+1) * (2*yLen+1) * (2*zLen+1));//70
    for (j = seedDistance ; j <= nyy - seedDistance - 1; ++j){//
        for (i = seedDistance ; i <= nxx - seedDistance - 1; ++i){
            for (ij = seedDistance ; ij <= nzz - seedDistance - 1; ++ij){
                if (tmpBinImg[i + j *nxx + ij * nXY ] > binValue){
                    bool isMax = true;
                    int tmp_Value= tmpOrigImg[i + j * nxx +  ij * nXY];
                    /*whether peak value?*/
                    for(int ii = i - seedDistance; ii <= i + seedDistance; ++ii){
                        for(int jj = j - seedDistance; jj <= j +seedDistance; ++jj){
                            for(int iijj = ij - seedDistance; iijj <= ij + seedDistance; ++iijj){
                                if(tmp_Value < tmpOrigImg[ii + jj * nxx +  iijj * nXY]){
                                    isMax = false;
                                    break;
                                }
                            }
                            if(!isMax) break;
                        }
                        if(!isMax) break;
                    }
                    ///complete
                    if (isMax){
                        ++maxNum;
                        locSeed.push_back(Vec4d(i,j,ij,tmp_Value));
                    }
                }
            }
        }
    }

    delete[] tmpBinImg;
    delete[] tmpOrigImg;

    if (maxNum < 1){
        locSeed.clear();
        locSeed.push_back(Vec4d(-1,-1,-1,-1));
    }
    /**/
    else if (maxNum > 1){
        VectorVec4d::size_type sz = tmpGlobalNode.size();
        tmpGlobalNode.resize(sz + locSeed.size());
        std::copy(locSeed.begin(), locSeed.end(), tmpGlobalNode.begin() + sz);
        int maxIndex = 0;
        int maxValue = 1;
        double seedLeastDistance = 4.0 * std::pow(2.0 / (xVoxel_ * yVoxel_ * zVoxel_), 1.0/3.0);
        while(maxValue > 0){
            /*VectorVec4d::const_iterator it = std::max_element(tmpGlobalNode.begin(), tmpGlobalNode.end(), [](Vec4d& lhs, Vec4d& rhs){
                    return lhs(3) < rhs(3);
            });*/
			VectorVec4d::const_iterator it = std::max_element(tmpGlobalNode.begin(), tmpGlobalNode.end(), Vec4d_3th_less());
			
            maxValue = (*it)(3);
            maxIndex = int(it - tmpGlobalNode.begin());//warning

            if (maxValue > 0){
                tmpGlobalNode[maxIndex](3) = -1;

                for (i = 0; i < maxNum; ++i){	//remove near seeds
                    if (i != maxIndex){
                        Vec3d tempPt1;
                        tempPt1 << tmpGlobalNode[maxIndex](0), tmpGlobalNode[maxIndex](1), tmpGlobalNode[maxIndex](2);
                        Vec3d tempPt2;
                        tempPt2 << tmpGlobalNode[i](0), tmpGlobalNode[i](1), tmpGlobalNode[i](2);
                        ///==============get euclidean distance=============
                        double normal = (tempPt1 - tempPt2).norm();
                        if (normal < seedLeastDistance) tmpGlobalNode[i](3) = 0;
                    }
                }
            }
        }

        VectorVec4d tmpData;
        for(VectorVec4d::size_type n = 0; n < tmpGlobalNode.size(); ++n){
            if(tmpGlobalNode[n](3) == -1)
                tmpData.push_back(tmpGlobalNode[n]);
        }
        std::swap(locSeed, tmpData);
    }

    /*sort ascend*/
    /*std::sort(locSeed.begin(), locSeed.end(), [](const Vec4d& lhs, const Vec4d& rhs){
        return lhs(3) > rhs(3);
    });*/
	std::sort(locSeed.begin(), locSeed.end(), Vec4d_3th_great());

    /*seeds num*/
    int min_10 = std::min((int)locSeed.size(), 100);
    tmpGlobalNode.clear();
    VectorVec4d tmpLocPos;
    for(int j = 0; j < min_10; ++j)
        tmpLocPos.push_back(locSeed[j]);
    //tmp_locPos.assign(loc_seed.begin(), loc_seed.begin() + min_10);
    locSeed.clear();
    locSeed.assign(tmpLocPos.begin(), tmpLocPos.end() );
    tmpLocPos.clear();

    /*padarray 3 pixels*/
    Vec3d Offset(xMin - 3, yMin -3, zMin -3);
    for(int i = 0; i < (int)locSeed.size(); ++i)
        globalSeed.push_back(Vec3d(locSeed[i](0) + Offset(0), locSeed[i](1) + Offset(1), locSeed[i](2) + Offset(2)));
}

void ProbeFilter::RayBurstShape(const Vec3d &initSoma, const SVolume &v, MatXd &resultRayLength, DVolume &smoothRay)
{
    resultRayLength.setZero();
    smoothRay.SetSize(0,0,0);

    //double slice_ is global varient , of which the value is 0.5   //double slice_ = 0.5;//(double)(minLen) / 82.0;
    const int blocksum = 41;//41;

    std::vector<double> lineSegment;
    for(int i = 0; i < blocksum; ++i){
        lineSegment.push_back(double(i) * slice_);
    }
    //generate_n(back_inserter(lineSegment), blocksum, GenArray<double>(0.0, slice_));//41个

    const int Theta = 40;
    const int Phi = 20;

    SVolume locOrigImg;
    Vec3d locSoma;

    ExtractLocalDomain(initSoma, v, locOrigImg, locSoma);

    DVolume sphereRayWet;//(lineSegment.size(), Theta, Phi);
    sphereRayWet.SetSize(lineSegment.size(), Theta, Phi);
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;

    //Volumn SubVol;//(locOrigImg.x(), locOrigImg.y(), locOrigImg.z());
    //int subx = locOrigImg.x();
    //int suby = locOrigImg.y();
    //int subz = locOrigImg.z();

    double segmentDense(0.0), segmentWet(0.0);
    double x,y,z;

    int numLineSegment;
    for (int j = 1; j <= Phi; ++j){
        for (int i = 1; i <= Theta; ++i){
            numLineSegment = lineSegment.size();
            for (int k = 0; k < numLineSegment; ++k){
                x = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::cos(a * (double)i * PI_180) + locSoma(0);
                y = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::sin(a * (double)i * PI_180) + locSoma(1);
                z = lineSegment[k] * std::cos(b * (double)j * PI_180) + locSoma(2);
                segmentDense = segmentWet = 0.0;
                ContourUtil::CalculateSphereOneNode(locOrigImg, 1.0, x, y, z, segmentDense, segmentWet);
                sphereRayWet(k, i-1, j-1) = segmentDense / (segmentWet + 0.0001);
            }
        }
    }//for
    //column storage
    int lenR0_1 = int(lineSegment.size()) - 1;
    std::vector<double> outerShell;
    for (int j = 0; j < Phi; ++j){
        for (int i = 0; i < Theta; ++i){
            outerShell.push_back(sphereRayWet(lenR0_1, i, j));
        }
    }

    std::vector<double> boundaryBack;
    ContourUtil::GetBoundaryBack(outerShell, 4, boundaryBack);

    //three_vs = mean(Lssx1)+3.5*std(Lssx1)
    double constrictionThrev(0.0);
    double constrictionThrevMean =
        std::accumulate(boundaryBack.begin(), boundaryBack.end(), 0.0) / boundaryBack.size();
    //get standard deviation---------------------------------------
    double constrictionThrevMeanSqrtSum(0);
    int numBoundaryBack = boundaryBack.size();
    for (int i = 0; i< numBoundaryBack; ++i)
    {
        constrictionThrevMeanSqrtSum +=
            (boundaryBack[i] - constrictionThrevMean) * (boundaryBack[i] - constrictionThrevMean);
    }
    constrictionThrev = 3.5 * std::sqrt(constrictionThrevMeanSqrtSum / (boundaryBack.size() - 1));
    constrictionThrev += constrictionThrevMean;
    //boundary value threshold
    for (int i = 0; i < Phi; ++i){
        for (int j = 0; j < Theta; ++j)
            sphereRayWet(lenR0_1, j, i) = boundaryBack[i * Theta + j];
    }

    std::vector<std::vector<double> > rayLimit;
    ContourUtil::GetRayLimit(sphereRayWet, constrictionThrev, rayLimit);
    ContourUtil::GetGradientVectorFlow(sphereRayWet, smoothRay);

    resultRayLength = 3.0 * MatXd::Ones(Theta + 2, Phi + 2);//ray length
    std::vector<double> lineSegLength;//segment length
    //generate_n(back_inserter(lineSegLength), blocksum, GenArray<double>(1.0, 1.0));//
    for(int i = 0; i < blocksum; ++i){//TODO:2014-5-21 here is a big bug!!!
        lineSegLength.push_back(double(i + 1));
    }

    double curRayLength(0);
    std::vector<double> curSmoothRay( smoothRay.x() );
    std::vector<double> distWet( smoothRay.x() );

    double reduceSmooth[2];
    reduceSmooth[0] = 0.9;
    reduceSmooth[1] = 1.0 - reduceSmooth[0];

    int repeat = 50;//!!!!!!
    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                curRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    curSmoothRay[ij] = smoothRay(ij, i - 1, j - 1);
                    distWet[ij] = curSmoothRay[ij] * std::exp( -0.05 *std::pow( lineSegLength[ij] - curRayLength, 2.0));
                }

                resultRayLength(i,j) += reduceSmooth[0] * ( ( (inner_product(lineSegLength.begin(), lineSegLength.end(),
                    distWet.begin(), 0.0)) /
                    (accumulate(distWet.begin(), distWet.end(), 0.0) ) ) - resultRayLength(i,j) )
                    - reduceSmooth[1] * (4.0 * resultRayLength(i,j)
                    - resultRayLength(i-1,j)
                    - resultRayLength(i+1,j)
                    - resultRayLength(i,j-1)
                    - resultRayLength(i,j+1));

                resultRayLength(i,j) = std::min(resultRayLength(i,j), rayLimit[i - 1][ j - 1]);
            }
            resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for
}

void ProbeFilter::_SomaRecognition(const CVolume &locBinImg, const VectorVec3d &locPtSet, const SVolume &locOrigImg,
                                   VectorVec3d &resultLocPtSet, std::vector<double> &radius)
{
    for(VectorVec3d::size_type i = 0; i < locPtSet.size(); ++i){
        resultLocPtSet.push_back(Vec3d(locPtSet[i](0), locPtSet[i](1), locPtSet[i](2)));
    }
    int i, mm;
    VectorVec3d::size_type locPtSetNum = resultLocPtSet.size();
    double maxAbsR;
    double rThre = 2.0f;
    radius.clear();
    radius.resize(locPtSetNum);
    std::vector<double> absR( radius.size() );
    /*initial radius*/
    std::vector<double> L1Wet(locPtSetNum);
    std::fill(radius.begin(), radius.end(), rThre);
    /*wet value*/
    std::fill(L1Wet.begin(), L1Wet.end(), 0.02f);//warning!! parameter
    //optimalize procedure
    for (i = 2; i <= 100; ++i){
        std::vector<double> newRadius;
        _CalculateL1Minimum(locBinImg, resultLocPtSet, L1Wet, radius, newRadius);//update radii
        for (mm = 1; mm < 3; ++mm){//
            VectorVec3d temp;
            _PostionUpdate(locOrigImg,resultLocPtSet, radius, temp);//update seeds position
            resultLocPtSet = temp;
        }
        std::vector<double> tmpR(radius.size());
        int numRadius = radius.size();
        for (int j = 0; j < numRadius; ++j){
            tmpR[j] = radius[j] - newRadius[j];//
        }
        /*increment of radii*/
        std::vector<double>::const_iterator it_MAX = std::max_element(radius.begin(),radius.end());
        double normalization = std::sqrt(std::inner_product(tmpR.begin(), tmpR.end(), tmpR.begin(), 0.0));//normalize
        /*convergence condition*/
        if (normalization < 0.005f  && *it_MAX > rThre)//
            break;
        else
            std::copy(newRadius.begin(), newRadius.end(), radius.begin());
    }

    /*modify wet value*/
    //w=.025*max(abs(r))./(abs(r)+0.1);
    int numRadius = radius.size();
    for (int j = 0; j < numRadius; ++j){
        absR[j] = std::abs(radius[j]);
    }

    std::vector<double>::const_iterator abs_it_MAX = std::max_element(absR.begin(),absR.end());
    maxAbsR = *abs_it_MAX;
    numRadius = radius.size();
    for (int j = 0; j < numRadius; ++j){
        L1Wet[j] = std::min(0.025 * maxAbsR / (absR[j] + 0.1), 1.0);
    }

    ///--------------------repeat-------------------------------------------
    for (i = 2; i <= 100; ++i){
        std::vector<double> newRadius;
        _CalculateL1Minimum(locBinImg, resultLocPtSet, L1Wet, radius, newRadius);

        for (mm = 1; mm <= 2; ++mm){
            VectorVec3d temp;
            _PostionUpdate(locOrigImg,resultLocPtSet, radius, temp);
            resultLocPtSet = temp;
        }

        std::vector<double> tmpR(radius.size());
        numRadius = radius.size();
        for (int j = 0; j < numRadius; ++j){
            tmpR[j] = radius[j] - newRadius[j];
        }

        std::vector<double>::const_iterator it_MAX = std::max_element(radius.begin(),radius.end());
        double normalization = std::sqrt(std::inner_product(tmpR.begin(), tmpR.end(), tmpR.begin(), 0.0));//normalize

        if (normalization < 0.005f ) //&& *it_MAX > r_Thre
            break;
        else std::copy(newRadius.begin(), newRadius.end(), radius.begin());
    }
    //w=.025*max(abs(r))./(abs(r)+0.1);
    numRadius = radius.size();
    for (int j = 0; j < numRadius; ++j){
        absR[j] = std::abs(radius[j]);
    }
    std::vector<double>::const_iterator tmp_it_MAX = std::max_element(absR.begin(),absR.end());
    maxAbsR = (*abs_it_MAX);
    numRadius = radius.size();
    for (int j = 0; j< numRadius;++j){
        L1Wet[j] = std::min(0.025 * maxAbsR / (absR[j] + 0.01), 1.0);
    }
    ///----------------------repeat-----------------------------------------
    for (i = 2; i <= 100; ++i){
        std::vector<double> newRadius;
        _CalculateL1Minimum(locBinImg, resultLocPtSet, L1Wet, radius, newRadius);
        for (mm = 1; mm <= 2; ++mm){
            VectorVec3d temp;
            _PostionUpdate(locOrigImg,resultLocPtSet, radius, temp);
            resultLocPtSet = temp;
        }
        std::vector<double> tmpR(radius.size());
        numRadius = radius.size();
        for (int j = 0; j < numRadius; ++j){
            tmpR[j] = radius[j] - newRadius[j];
        }
        std::vector<double>::const_iterator it_MAX = std::max_element(radius.begin(),radius.end());
        double normalization = std::sqrt(std::inner_product(tmpR.begin(), tmpR.end(), tmpR.begin(), 0.0));//normalize
        if (normalization < 0.005f ) //&& *it_MAX > r_Thre
            break;
        else std::copy(newRadius.begin(), newRadius.end(), radius.begin());
    }
    ///----------end-------------------------------------------------------
}

void ProbeFilter::_PostionUpdate(const SVolume &locOrigImg, const VectorVec3d &locPtSet, const std::vector<double> &radius,
                                 VectorVec3d &resultLocPtSet)
{
    VectorVec3d::size_type nx = locPtSet.size();
    resultLocPtSet.clear();
    for (VectorVec3d::size_type i = 0; i < nx; ++i){
        Vec3d newPtSet;
        _CalculateNewPosition(locOrigImg, locPtSet[i], radius[i], newPtSet);
        resultLocPtSet.push_back(newPtSet);
    }
}

void ProbeFilter::_CalculateNewPosition(const SVolume &locOrigImg, const Vec3d &locPtSet, const double radius, Vec3d &resultLocPos)
{
    int i, j, ij;
    int nx = locOrigImg.x();
    int ny = locOrigImg.y();
    int nz = locOrigImg.z();
    //int nxy = nx * ny;
    int n = (int)radius;//int n = (int)round<double>(radius - 0.5f);
    double locImgIndense = 0.0;
    double locImgWet[3] = {0.0f,0.0f,0.0f};
    Vec3d roundPtSet(int(locPtSet(0) + 0.5), int(locPtSet(1) + 0.5), int(locPtSet(2)+0.5));
    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, roundPtSet(0), roundPtSet(1), roundPtSet(2),
        n, n, n, 0, nx-1, 0, ny-1, 0, nz-1);
    ///---------------------------------------------------------------------------///
    for (j = yMin; j <= yMax; ++j){
        for (i = xMin; i <= xMax; ++i){
            for (ij = zMin; ij <= zMax; ++ij) {
                double rDist[] = {xVoxel_ * ((double)i-locPtSet(0)),
                    yVoxel_ * ((double)j - locPtSet(1)),
                    zVoxel_ * ((double)ij-locPtSet(2))};

                if (std::sqrt(std::inner_product(rDist,rDist+3, rDist, 0.0)) < radius){//seed radius
                    double temp   = (double)locOrigImg(i , j , ij );
                    locImgIndense	+= temp;
                    locImgWet[0]	+= temp * rDist[0];
                    locImgWet[1]	+= temp * rDist[1];
                    locImgWet[2]	+= temp * rDist[2];
                }
            }
        }
    }
    ///--------------------------------------------------------------------------------------///
    double temp = locImgIndense + 0.01f; //avoid 0
    resultLocPos(0) = locPtSet(0) + locImgWet[0] / temp;
    resultLocPos(1) = locPtSet(1) + locImgWet[1] / temp;
    resultLocPos(2) = locPtSet(2) + locImgWet[2] / temp;
}

void ProbeFilter::_CalculateL1Minimum(const CVolume &locBinImg, const VectorVec3d &locPtSet, const std::vector<double> &L1Wet,
                                      const std::vector<double> &srcRadius, std::vector<double> &resultRadius)
{
    resultRadius = srcRadius;

    int i,j;
    double gAbsMax;
    int	nss;
    const int N = 18;

    std::vector<double> rGrad;
    _CalculateGrad(locBinImg, locPtSet, L1Wet, resultRadius, rGrad);//partial derivative
    //g=g./(max(abs(g))+0.01);
    std::vector<double> absG( rGrad.size() );
    int numRGrad = rGrad.size();
    for (i = 0; i < numRGrad; ++i){
        absG[i] = std::abs(rGrad[i]);
    }

    std::vector<double>::const_iterator it_abs_g_MAX = std::max_element(absG.begin(), absG.end());
    gAbsMax = *it_abs_g_MAX + 0.01;//avoid 0

    numRGrad = rGrad.size();
    for (i = 0; i < numRGrad;++i){
        rGrad[i] /= gAbsMax;//normalize
    }
    nss = rGrad.size();
    for (i = 0; i < nss;++i){
        if (rGrad[i] > 0.0 && resultRadius[i] < 0.01)
            rGrad[i] = 0.0;
    }

    std::vector<double> calL1Value(N);
    double	curL1Value=0.0;

    /*L1 minimum sphere function*/
    _CalculateL1Value(locBinImg, locPtSet, L1Wet, resultRadius, curL1Value);

    for (i = 0 ; i <N; ++i){//
        std::vector<double> tmpR1(resultRadius.size());
        int numResultRadius = resultRadius.size();
        for ( j = 0; j < numResultRadius; ++j ){
            tmpR1[j] = std::max(0.0, resultRadius[j] - std::pow(0.5,(int)i+1) * rGrad[j]);//g is negtive
        }

        _CalculateL1Value(locBinImg, locPtSet, L1Wet, tmpR1, calL1Value[i]);//get sphere function
    }

    int numCalL1Value = calL1Value.size();
    for (i = 0; i < numCalL1Value; ++i){
        calL1Value[i] = calL1Value[i] - curL1Value;
    }
    std::vector<double>::const_iterator it_MIN = std::min_element(calL1Value.begin(), calL1Value.end());//minimize
    double v = *it_MIN;
    int idx = int(it_MIN - calL1Value.begin());

    if (v < 0.0f){
        int numResultRadius = resultRadius.size();
        for (i = 0; i < numResultRadius;++i){
            resultRadius[i] = std::max(0.0, resultRadius[i] - std::pow(0.5, (int)idx+1) * rGrad[i]);
        }
    }
}

void ProbeFilter::_CalculateGrad(const CVolume &locBinImg, const VectorVec3d &locPtSet,
                                 const std::vector<double> &rWet, const std::vector<double> &radius, std::vector<double> &rGrad)
{
    rGrad.clear();
    VectorVec3d::size_type seednum = locPtSet.size();
    double curL1Value(0.0);
    double calL1Value(0.0);
    _CalculateL1Value(locBinImg, locPtSet, rWet, radius, curL1Value);

    double accu = 0.5;//_bMerge? 0.5f : 1.0f;

    for (VectorVec3d::size_type i = 0; i < seednum; ++i){
        std::vector<double> increRadius(radius);
        increRadius[i] += accu; //0.5f
        _CalculateL1Value(locBinImg, locPtSet, rWet, increRadius, calL1Value);
        rGrad.push_back((calL1Value - curL1Value) / accu ); //
    }
}

void ProbeFilter::_CalculateL1Value(const CVolume &locBinImg, const VectorVec3d &locPtSet, const std::vector<double> &L1Wet,
                                    const std::vector<double> &radius, double &L1Value)
{
    size_t np = locPtSet.size();
    std::vector<double> H;
    /*H.assign(IIm.head(), IIm.tail());*/
    for (int ij = 0; ij < locBinImg.z(); ++ij){
        for (int j = 0; j < locBinImg.y(); ++j){
            for ( int i = 0; i < locBinImg.x(); ++i){
                H.push_back((double)locBinImg(i,j,ij));
            }
        }
    }

    int vx = locBinImg.x();
    int vy = locBinImg.y();
    int vz = locBinImg.z();
    std::vector<double> M;
    for (size_t i= 0; i < np; ++i){
        //M = calculate_M( pars( : , i ) , nx , ny , nz , r(i) );
        _CalculateSphereWet(locPtSet[i], vx, vy , vz , radius[i], M);//sphere function
        size_t num_M = M.size();
        for (size_t j = 0; j < num_M; ++j){
            H[j] = H[j] - M[j];
        }
    }
    L1Value = 0.0f;
    L1Value = std::inner_product(H.begin(),H.end(),H.begin(),0.0);
    L1Value = pow(L1Value, 1.0/3.0);
    L1Value += std::inner_product(L1Wet.begin(),L1Wet.end(),radius.begin(),0.0);
}

void ProbeFilter::_CalculateSphereWet(const Vec3d &seed, const int vx, const int vy, const int vz,
                                      const double radius, std::vector<double> &sphereWet)
{
    //cal sphere function
    int i, j, ij;
    sphereWet.clear();
    sphereWet.resize(vx * vy * vz);
    //TODO:2014-5-17
    int wetRadius = int(std::sqrt(3.0) * radius + 0.5) + 1;//threshold value is 3.0

    Vec3i rSeed(int(seed(0) + 0.5), int(seed(1) + 0.5), int(seed(2) + 0.5));//int

    int xMin, xMax, yMin, yMax, zMin, zMax;
    Get3DRegion(xMin, xMax, yMin, yMax, zMin, zMax, rSeed(0), rSeed(1), rSeed(2),
        wetRadius, wetRadius, wetRadius, 0, vx-1, 0, vy-1, 0, vz-1);

    int xy = vx * vy;
    for (j = yMin; j <= yMax; ++j){
        for (i = xMin; i <= xMax; ++i){
            for (ij = zMin; ij <= zMax; ++ij){
                double rDist[] = {xVoxel_ * ((double)i-seed(0)),
                                  yVoxel_ * ((double)j - seed(1)),
                                  zVoxel_ * ((double)ij-seed(2))};
                double nr = std::sqrt(std::inner_product(rDist,rDist+3, rDist, 0.0));
                sphereWet[i + j * vx + ij * xy] = (nr < radius) ? 1.0f : std::pow(M_E, - std::pow(nr-radius,2)/0.1f);//exp
            }
        }
    }
}

void ProbeFilter::CalSomaInfo(const Vec3d &initSoma, const DVolume &smoothRay,
                              double &majAxis, double &minAxis, double &somaSurface, double &somaVolume)
{
    VectorVec3d contourPtSet;
    VectorVec3d contourPatchPtSet;
    MatXd resultRayLength;
    GetContour(initSoma, smoothRay, contourPtSet, resultRayLength);//contour
    GetAxis(resultRayLength, majAxis, minAxis);//long short axis
    //const int k_theta = 40, k_phi = 20;
    ContourUtil::CalAreaVolume(contourPtSet, k_theta, k_phi,
        contourPatchPtSet, somaSurface, somaVolume);//surface and volume
}

void ProbeFilter::GetContour(const Vec3d &initSoma, const DVolume &smoothRay,
                             VectorVec3d &contourPtSet, MatXd &resultRayLength)
{
     //double slice_ is global varient , of which the value is 0.5  //double slice = 0.5;//(double)(minLen) / 82.0;
    const int blocksum = 31;//41;
    const int Theta = 40;
    const int Phi = 20;
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;

    resultRayLength = 3.0 * MatXd::Ones(Theta + 2, Phi + 2);//ray length
    std::vector<double> lineSegLength;
    //generate_n(back_inserter(line_seg_length), blocksum, GenArray<double>(1.0, 1.0));
    for(int i = 0; i < blocksum; ++i){
        lineSegLength.push_back(double(i + 1));
    }

    double curRayLength(0);
    std::vector<double> curSmoothRay(smoothRay.x());
    std::vector<double> distWet(smoothRay.x());

    double reduceSmooth[2];
    reduceSmooth[0] = 0.9;
    reduceSmooth[1] = 1.0 - reduceSmooth[0];

    int repeat = 100;

    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                curRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    curSmoothRay[ij] = smoothRay(ij, i - 1, j - 1);
                    distWet[ij] = curSmoothRay[ij] * std::exp( -0.05 *std::pow( lineSegLength[ij] - curRayLength, 2.0));
                }

                resultRayLength(i,j) += reduceSmooth[0] * ( ( (inner_product(lineSegLength.begin(), lineSegLength.end(), distWet.begin(), 0.0)) /
                    (accumulate(distWet.begin(), distWet.end(), 0.0) ) ) - resultRayLength(i,j) )
                    - reduceSmooth[1] * (4.0 * resultRayLength(i,j)
                    - resultRayLength(i-1,j)
                    - resultRayLength(i+1,j)
                    - resultRayLength(i,j-1)
                    - resultRayLength(i,j+1));

            }
            resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for

    double xx(0.0), yy(0.0), zz(0.0);
    Vec3d tmp;
    //double k(0.0);
    VectorVec3d tmpContour;

    for (int j = 2; j <= Phi  +1; ++j){
        for (int i = 2; i <= Theta + 1;++i){
            if (resultRayLength(i-1, j - 1)  > 0.0){
                double k =  slice_ * resultRayLength(i - 1, j - 1);//	slice_ * meanUx;
                xx = (k * std::sin(b * (double)(j-1) * PI_180) * std::cos(a * (double)(i-1) * PI_180) + initSoma(0));
                yy = (k * std::sin(b * (double)(j-1) * PI_180) * std::sin(a * (double)(i-1) * PI_180) + initSoma(1));
                zz = (k * std::cos(b * (double)(j-1) * PI_180) + initSoma(2));
                tmp(0) = xx;
                tmp(1) = yy;
                tmp(2) = zz;
                tmpContour.push_back(tmp);
            }
        }
    }

    if(tmpContour.empty())
        return;
    Vec3d sum(0.0,0.0,0.0);
    for (VectorVec3d::const_iterator it = tmpContour.begin() + 1; it != tmpContour.begin() + 41; ++it){
        sum(0) += (*it)(0);
        sum(1) += (*it)(1);
        sum(2) += (*it)(2);
    }
    sum(0) /= 40.0f; sum(1) /= 40.0f; sum(2) /= 40.0f;
    contourPtSet.push_back(sum);

    std::copy(tmpContour.begin(), tmpContour.end(), std::back_inserter(contourPtSet));
    contourPtSet.erase(contourPtSet.end() - 39, contourPtSet.end());
}

void ProbeFilter::GetContour(const Vec3d &initSoma, const MatXd &resultRayLength,const int Theta, const int Phi,
                    VectorVec3d &contourPtSet)
{
	//const int Theta = 40;
    //const int Phi = 20;
	const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;
	double xx(0.0), yy(0.0), zz(0.0);
    Vec3d tmp;
    //double k(0.0);
    VectorVec3d tmpContour;

    for (int j = 2; j <= Phi  +1; ++j){
        for (int i = 2; i <= Theta + 1;++i){
            if (resultRayLength(i-1, j - 1)  > 0.0){
                double k =  slice_ * resultRayLength(i - 1, j - 1);//	slice * meanUx;
                xx = (k * std::sin(b * (double)(j-1) * PI_180) * std::cos(a * (double)(i-1) * PI_180) + initSoma(0));
                yy = (k * std::sin(b * (double)(j-1) * PI_180) * std::sin(a * (double)(i-1) * PI_180) + initSoma(1));
                zz = (k * std::cos(b * (double)(j-1) * PI_180) + initSoma(2));
                tmp(0) = xx;
                tmp(1) = yy;
                tmp(2) = zz;
                tmpContour.push_back(tmp);
            }
        }
    }

    if(tmpContour.empty())
        return;
    Vec3d sum(0.0,0.0,0.0);
    for (VectorVec3d::const_iterator it = tmpContour.begin() + 1; it != tmpContour.begin() + 41; ++it){
        sum(0) += (*it)(0);
        sum(1) += (*it)(1);
        sum(2) += (*it)(2);
    }
    sum(0) /= 40.0f; sum(1) /= 40.0f; sum(2) /= 40.0f;
    contourPtSet.push_back(sum);

    std::copy(tmpContour.begin(), tmpContour.end(), std::back_inserter(contourPtSet));
    contourPtSet.erase(contourPtSet.end() - 39, contourPtSet.end());
}

void ProbeFilter::GetAxis(const MatXd &resultRayLength, double &majAxis, double &minAxis)
{
    MatXd blockUx = resultRayLength.block(1,1, k_theta, k_phi);
    int nx = int(blockUx.rows());
    int ny = int(blockUx.cols());
    std::vector<double> R;
    for (int i = 0; i < nx / 2; ++i){
        for ( int j = 0; j < ny - 1; ++j ){
            R.push_back( blockUx(i,j) + blockUx(i + nx / 2, ny - j - 2) );
        }
    }
    std::sort(R.begin(), R.end());
    //
    minAxis = 0.2 * std::accumulate(R.begin(), R.begin() + 8, 0.0) / 8.0;
    majAxis = 0.2 * std::accumulate(R.rbegin(), R.rbegin() + 8, 0.0) / 8.0;
}

int ProbeFilter::CalcConv3dRadius(const double resolution)
{
    if(resolution >0.5 && resolution <= 1.5){
        return 3;
    }else if(resolution > 1.5 && resolution <= 2.5){
        return 2;
    }else if(resolution > 2.5){
        return 1;
    }
    return 0;
}

void ProbeFilter::GetInnerPtSet(const MatXd &resultRayLength, int Theta, int Phi, double slice, double xVox, double yVox, double zVox,
                                VectorVec3i& innerPtSet)
{
    Vec3i tmp;
    double xx, yy, zz;
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;
    for (int i = 2; i <= Theta + 1;++i){
        for (int j = 2; j <= Phi  +1; ++j){
            if (resultRayLength(i-1, j - 1)  > 0.0){
                for (double k = 1.0; k <= slice * resultRayLength(i - 1, j - 1); ++k){//warning!!
                    xx = k * std::sin(b * (double)(j-1) * PI_180) * std::cos(a * (double)(i-1) * PI_180);
                    yy = k * std::sin(b * (double)(j-1) * PI_180) * std::sin(a * (double)(i-1) * PI_180);
                    zz = k * std::cos(b * (double)(j-1) * PI_180);

                    tmp(0) = NGUtility::Round(xx * xVox);
                    tmp(1) = NGUtility::Round(yy * yVox);
                    tmp(2) = NGUtility::Round(zz * zVox);
                    innerPtSet.push_back(tmp);
                }
            }
        }
    }
    std::sort(innerPtSet.begin(), innerPtSet.end(), Vec3i_less());
    innerPtSet.erase(std::unique(innerPtSet.begin(), innerPtSet.end()), innerPtSet.end());
}

void ProbeFilter::CalcInnerRegionPCA(const VectorVec3i &innerPtSet, Vec3d& firDir, double &lambda1, double &lambda2, double &lambda3)
{
    MatXd regionMatrix(3, innerPtSet.size());
    for(size_t i = 0; i < innerPtSet.size(); ++i){
        regionMatrix(0,i) = innerPtSet[i](0);
        regionMatrix(1,i) = innerPtSet[i](1);
        regionMatrix(2,i) = innerPtSet[i](2);
    }
    Mat3d PCAMatrix = regionMatrix * regionMatrix.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(PCAMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV );
    firDir = svd.matrixU().col(0);
    Vec3d e= svd.singularValues();
    lambda1 = e(0);
    lambda2 = e(1);
    lambda3 = e(2);
}

int ProbeFilter::JudgeRegionShapeUsingEigenValue(double lambda1, double lambda2, double lambda3)
{
    if(lambda1 / std::sqrt(lambda2 * lambda2 + lambda3 * lambda3) > 2.0){
        return 0;
    }
//    if(lambda1 / std::sqrt(std::abs(lambda2 * lambda3)) > 3){
//        return 0;
//    }
//    if(lambda1 / lambda3 < 5)
//        return 0;
    return 1;
}

void ProbeFilter::CorrectRayMatrix(const MatXd &rayLen, int Theta, int Phi, MatXd &newRayLen)
{
    Vec3i tmp;
    newRayLen = rayLen;
    double xx, yy, zz;
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;
    double len;
    for (int i = 2; i <= Theta + 1;++i){
        for (int j = 2; j <= Phi  +1; ++j){
            if (rayLen(i-1, j - 1)  > 0.0){
                double k = rayLen(i - 1, j - 1);
                xx = k * std::sin(b * (double)(j-1) * PI_180) * std::cos(a * (double)(i-1) * PI_180);
                yy = k * std::sin(b * (double)(j-1) * PI_180) * std::sin(a * (double)(i-1) * PI_180);
                zz = k * std::cos(b * (double)(j-1) * PI_180);
                xx *= xVoxel_;
                yy *= yVoxel_;
                zz *= zVoxel_;
                len = std::sqrt(xx * xx + yy * yy + zz * zz);
                newRayLen(i-1, j - 1) = len;
            }
        }
    }
}



bool ProbeFilter::IsInVolume(const Vec3d &pos, int lx, int ly, int lz, int hx, int hy, int hz)
{
    int sx = NGUtility::Round(pos(0));
    int sy = NGUtility::Round(pos(1));
    int sz = NGUtility::Round(pos(2));
    if(sx < lx || sx > hx || sy < ly || sy > hy || sz < lz || sz > hz){
        return false;
    }
    return true;
}

void ProbeFilter::MergeSomasByShape(const SVolume &v, StackVec4d& gloSomaStack, 
		const Vec3i& offSet, MatXd& newUx, CVolume &locBinImg )
{
	//initial
	int sx = locBinImg.x();
	int sy = locBinImg.y();
	int sz = locBinImg.z();
	Vec3d tmpSet(offSet(0), offSet(1), offSet(2));
	StackVec4d tmpGloSomaStack;
	VectorVec4d mergePt;
	MatXd Ux;
	DVolume smoothRay;
	Vec4d popPt;
	Vec3i locPti;
	bool firstUxUsed = false;
	//loop
	do{
		Vec4d head = gloSomaStack.top();
		gloSomaStack.pop();
		mergePt.clear();
		if(firstUxUsed){//need to cal newUx again
			RayBurstShape(head.block(0,0,3,1), v, Ux, smoothRay);
			CorrectRayMatrix(Ux, k_theta, k_phi, newUx);// introduce voxel size
		}
		if(!firstUxUsed)
			firstUxUsed = true;
		locBinImg.SetZero();
		FillShape(head.block(0,0,3,1) - tmpSet, newUx, locBinImg);
		//search merge somas
		while(!gloSomaStack.empty()){
			popPt = gloSomaStack.top();
			gloSomaStack.pop();
			//limit index
            locPti << NGUtility::Round(popPt(0)) - offSet(0), NGUtility::Round(popPt(1)) - offSet(1), NGUtility::Round(popPt(2)) - offSet(2);
			locPti(0) = std::max(std::min(sx -1, locPti(0)), 0);
			locPti(1) = std::max(std::min(sy -1, locPti(1)), 0);
			locPti(2) = std::max(std::min(sz -1, locPti(2)), 0);
			if(locBinImg(locPti(0), locPti(1), locPti(2)) > 0)// I do not sure 255 or 1
				mergePt.push_back(popPt);
			else
				tmpGloSomaStack.push(popPt);
		}
		for(size_t i = 0; i < mergePt.size(); ++i)
			head += mergePt[i];
		head.block(0,0,3,1) /= double(mergePt.size()) + 1.0;//do not divide radius
		tmpGloSomaStack.push(head);//remained soma
#ifdef _WIN32
		//gloSomaStack.swap(tmpGloSomaStack);//tmpGloSomaStack is empty
		std::swap(gloSomaStack, tmpGloSomaStack);
#else
		gloSomaStack.swap(tmpGloSomaStack);//tmpGloSomaStack is empty
#endif
		
	}while( !mergePt.empty());
}

void ProbeFilter::FillShape(const Vec3d& locSoma, MatXd& newUx, CVolume& locBinImg)
{
	//
	int sx = locBinImg.x();
	int sy = locBinImg.y();
	int sz = locBinImg.z();
    Vec3i locSomai(NGUtility::Round(locSoma(0)), NGUtility::Round(locSoma(1)), NGUtility::Round(locSoma(2)));
	//get inner point set
	Vec3i tmp;
	VectorVec3i innerPtSet;
    double xx, yy, zz;
    const double a = 360.0 / (double)k_theta;
    const double b = 180.0 / (double)k_phi;
    const double PI_180 = M_PI / 180.0;
	/*FILE* fp = fopen("D:/7_Ux.txt","w");
	for(int i = 0; i < newUx.rows(); ++i){
		for(int j = 0; j < newUx.cols(); ++ j)
			fprintf(fp, "%lf ", newUx(i,j));
		fprintf(fp,"\n");
	}
	fclose(fp);*/
    for (int i = 2; i <= k_theta + 1;++i){
        for (int j = 2; j <= k_phi  +1; ++j){
            if (newUx(i-1, j - 1)  > 0.0){
                for (double k = 1.0; k <= slice_ * newUx(i - 1, j - 1); ++k){//warning!!
                    xx = k * std::sin(b * (double)(j-1) * PI_180) * std::cos(a * (double)(i-1) * PI_180);
                    yy = k * std::sin(b * (double)(j-1) * PI_180) * std::sin(a * (double)(i-1) * PI_180);
                    zz = k * std::cos(b * (double)(j-1) * PI_180);
					//
                    tmp(0) = NGUtility::Round(xx);//Round(xx * xVoxel_);
                    tmp(1) = NGUtility::Round(yy);//Round(yy * yVoxel_);
                    tmp(2) = NGUtility::Round(zz);//Round(zz * zVoxel_);
                    innerPtSet.push_back(tmp);
					for(int ii = tmp(0) - 1 ; ii <= tmp(0) + 1; ++ii){
						for(int jj = tmp(1) - 1 ; jj <= tmp(1) + 1; ++jj){
							for(int kk = tmp(2) - 1 ; kk <= tmp(2) + 1; ++kk){
								innerPtSet.push_back(Vec3i(ii,jj,kk));
							}
						}
					}
                }
            }
        }
    }
    std::sort(innerPtSet.begin(), innerPtSet.end(), Vec3i_less());
	innerPtSet.erase(std::unique(innerPtSet.begin(), innerPtSet.end()), innerPtSet.end());
	/*FILE* fp = fopen("D:/7_crop.txt","w");
	for(size_t i = 0; i < innerPtSet.size(); ++i){
		fprintf(fp, "%d %d %d\n",innerPtSet[i](0),innerPtSet[i](1),innerPtSet[i](2));
	}
	fclose(fp);*/
	//fill locBinImg
	//FILE* fp = fopen("D:/7_crop.txt","w");
	Vec3i tmpPos;
	for(size_t i = 0; i < innerPtSet.size(); ++i){
		tmpPos << std::max(std::min(innerPtSet[i](0) + locSomai(0), sx - 1), 0),
			std::max(std::min(innerPtSet[i](1) + locSomai(1), sy -1 ), 0),
			std::max(std::min(innerPtSet[i](2) + locSomai(2), sz -1 ), 0);
		locBinImg( tmpPos(0), tmpPos(1), tmpPos(2) ) = 255;
		//fprintf(fp, "%d %d %d\n",tmpPos(0),tmpPos(1),tmpPos(2));
	}
	//fclose(fp);
}

void ProbeFilter::SomaShape( const Vec3d& initSoma, const SVolume& origImg, const CVolume& binImg, 
                            VectorVec3i& interPtSet, MatXd& resultRayLength, std::vector<std::vector<double> >& rayLimit,
                            std::vector<VectorVec3i> &connectPtSet, double& constrictionThrev)
{
     //double slice_ is global varient , of which the value is 0.5  //double slice = 0.5;//(double)(minLen) / 82.0;
    const int blocksum = 41;

    std::vector<double> lineSegment;
    for(int i = 0; i < blocksum; ++i){
        lineSegment.push_back(double(i) * slice_);
    }

    const int Theta = 40;
    const int Phi = 20;

    const SVolume& locOrigImg = origImg;
    const Vec3d& locSoma = initSoma;

    DVolume sphereRayWet;//(lineSegment.size(), Theta, Phi);
    sphereRayWet.SetSize(lineSegment.size(), Theta, Phi);
    const double a = 360.0 / (double)Theta;
    const double b = 180.0 / (double)Phi;
    const double PI_180 = M_PI / 180.0;

    double segmentDense(0.0), segmentWet(0.0);
    double x,y,z;

    //int subx = locOrigImg.x();
    //int suby = locOrigImg.y();
    //int subz = locOrigImg.z();

    int numLineSegment;
    numLineSegment = lineSegment.size();
    for (int k = 0; k < numLineSegment; ++k){
    //for (int k = numLineSegment - 1; k < numLineSegment; ++k){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                x = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::cos(a * (double)i * PI_180) + locSoma(0);
                y = lineSegment[k] * std::sin(b * (double)j * PI_180) * std::sin(a * (double)i * PI_180) + locSoma(1);
                z = lineSegment[k] * std::cos(b * (double)j * PI_180) + locSoma(2);
                segmentDense = segmentWet = 0.0;
                //ContourUtil::CalculateSphereOneNode(locOrigImg,  x, y, z, 1.0, segmentDense, segmentWet);
                ContourUtil::CalculateSphereOneNode(locOrigImg, 1.0, x, y, z, segmentDense, segmentWet);
                sphereRayWet(k, i-1, j-1) = segmentDense / (segmentWet + 0.0001);
            }
        }
    }//for
    //
    int lenR0_1 = int(lineSegment.size()) - 1;
    std::vector<double> outerShell;
    for (int j = 0; j < Phi; ++j){
        for (int i = 0; i < Theta; ++i){
            outerShell.push_back(sphereRayWet(lenR0_1, i, j));
        }
    }

    std::vector<double> boundaryBack;
    ContourUtil::GetBoundaryBack(outerShell, 4, boundaryBack);

    //three_vs = mean(Lssx1)+3.5*std(Lssx1)
    constrictionThrev = 0.0;
    double constrictionThrevMean =
        std::accumulate(boundaryBack.begin(), boundaryBack.end(), 0.0)
        / double(boundaryBack.size());
    //---------------------------------------
    /*double constrictionThrevMeanSqrtSum(0);
    size_t numBoundaryBack = boundaryBack.size();
    for (size_t i = 0; i< numBoundaryBack; ++i)
    {
        constrictionThrevMeanSqrtSum +=
            (boundaryBack[i] - constrictionThrevMean) * (boundaryBack[i] - constrictionThrevMean);
    }
    constrictionThrev = 3.5 * std::sqrt(constrictionThrevMeanSqrtSum / (boundaryBack.size() - 1));
	constrictionThrev += constrictionThrevMean;*/
    constrictionThrev = 1.2 * constrictionThrevMean;//2015-4-27

    for (int i = 0; i < Phi; ++i){
        for (int j = 0; j < Theta; ++j)
            sphereRayWet(lenR0_1, j, i) = boundaryBack[i * Theta + j];
    }

    Volume<double> smoothRay;//2015-5-7
    //std::vector<std::vector<double> > rayLimit;
    ContourUtil::GetRayLimit(sphereRayWet, constrictionThrev, rayLimit);
    TraceUtil::GetGradientVectorFlowForTrace(sphereRayWet, smoothRay);

    //MatXd resultRayLength = 3.0 * MatXd::Ones(Theta + 2, Phi + 2);
    std::vector<double> lineSegLength;
    //generate_n(back_inserter(lineSegLength), blocksum, GenArray<double>(1.0, 1.0));
    for(int i = 0; i < blocksum; ++i){
        lineSegLength.push_back(double(i + 1));
    }

    double curRayLength(0);
    std::vector<double> curSmoothRay( smoothRay.x() );
    std::vector<double> distWet( smoothRay.x() );

    double reduceSmooth[2];
    reduceSmooth[0] = 0.9;
    reduceSmooth[1] = 1.0 - reduceSmooth[0];

    resultRayLength=6.0*MatXd::Ones(Theta+2, Phi+2);// = rayLimit;
    /*for(size_t i = 0; i < rayLimit.size();++i){
    for(size_t j = 0; j < rayLimit[0].size(); ++j){
    resultRayLength(i+1,j+1)=rayLimit[i][j];
    }
    }*/

    int repeat = 100;//
    for (int jj = 0; jj < repeat;++jj){
        for (int i = 1; i <= Theta; ++i){
            for (int j = 1; j <= Phi; ++j){
                curRayLength = resultRayLength(i, j);
                for (int ij = 0; ij < smoothRay.x(); ++ij){
                    curSmoothRay[ij] = smoothRay(ij, i - 1, j - 1);
                    distWet[ij] = curSmoothRay[ij] * std::exp( -0.05 *std::pow( lineSegLength[ij] - curRayLength, 2.0));
                }

                resultRayLength(i,j) += reduceSmooth[0] * ( ( (std::inner_product(lineSegLength.begin(), lineSegLength.end(),
                    distWet.begin(), 0.0)) /
                    (std::accumulate(distWet.begin(), distWet.end(), 0.0) + 0.001 ) ) - resultRayLength(i,j) )
                    - reduceSmooth[1] * (4.0 * resultRayLength(i,j)
                    - resultRayLength(i-1,j)
                    - resultRayLength(i+1,j)
                    - resultRayLength(i,j-1)
                    - resultRayLength(i,j+1));

                resultRayLength(i,j) = std::min(resultRayLength(i,j), rayLimit[i - 1][ j - 1]);
            }
            resultRayLength = (resultRayLength.array() > 0.0).select(resultRayLength, 0.0);
            resultRayLength.row(0) = resultRayLength.row(2);
            resultRayLength.row(Theta + 1) = resultRayLength.row(Theta - 1);
            resultRayLength.col(0) = resultRayLength.col(2);
            resultRayLength.col(Phi + 1) = resultRayLength.col(Phi - 1);
        }
    }//for

    //bounding box
    int xBoundMin = 1000000;
    int yBoundMin = 1000000;
    int zBoundMin = 1000000;
    int xBoundMax = 0;
    int yBoundMax = 0;
    int zBoundMax = 0;
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    //get shape inner point set
    interPtSet.reserve(Theta * Phi * lineSegment.size());
    double xx(0.0), yy(0.0), zz(0.0);
    Vec3i tmp;
    for(int j = 1; j < Phi + 1; ++j){
        for(int i = 1; i < Theta + 1; ++i){
            if(resultRayLength(i, j) > 1.0){
                for(double k = 1.0; k <= slice_ * resultRayLength(i,j); k += slice_){
                    xx = (k * std::sin(b * (double)(j) * PI_180) * std::cos(a * (double)(i) * PI_180) + initSoma(0));
                    yy = (k * std::sin(b * (double)(j) * PI_180) * std::sin(a * (double)(i) * PI_180) + initSoma(1));
                    zz = (k * std::cos(b * (double)(j) * PI_180) + initSoma(2));
                    tmp(0) = int(NGUtility::Round(xx));
                    tmp(1) = int(NGUtility::Round(yy));
                    tmp(2) = int(NGUtility::Round(zz));
                    // far away from boundary.
                    if ( tmp(0) >= nx-1 || tmp(0) < 1 || tmp(1) >= ny -1 || tmp(1) < 1 || tmp(2) >= nz -1 || tmp(2) < 1   ) {
                        continue;
                    }
                    for(int l = std::max(tmp(0) - 1, 0); l < std::min(tmp(0) +2, nx); ++ l)
                        for(int m = std::max(tmp(1) - 1, 0); m < std::min(tmp(1) +2, ny); ++m)
                            for(int n = std::max(tmp(2) - 1, 0); n < std::min(tmp(2) +2, nz); ++n){
                                interPtSet.push_back(Vec3i(l,m,n));
                                if ( xBoundMin > l) xBoundMin = l;
                                if ( yBoundMin > m) yBoundMin = m;
                                if ( zBoundMin > n) zBoundMin = n;
                                if ( xBoundMax < l) xBoundMax = l;
                                if ( yBoundMax < m) yBoundMax = m;
                                if ( zBoundMax < n) zBoundMax = n;
                            }
                }
            }
        }
    }

    //masaka, empty interPtSet
    if (interPtSet.empty()) {
        connectPtSet.clear();
        return;
    }

    std::sort(interPtSet.begin(), interPtSet.end(), Vec3i_less());
    interPtSet.erase(std::unique(interPtSet.begin(), interPtSet.end()), interPtSet.end());
    //extract local binary image
    //int xOffset = (xBoundMax - xBoundMin)/2 + 5;
    //int yOffset = (yBoundMax - yBoundMin)/2 + 5;
    //int zOffset = (zBoundMax - zBoundMin)/2 + 5;
    int xOffset = std::max(xBoundMax - int(initSoma(0)),  int(initSoma(0)) - xBoundMin ) + 5;
    int yOffset = std::max(yBoundMax - int(initSoma(1)),  int(initSoma(1)) - yBoundMin ) + 5;
    int zOffset = std::max(zBoundMax - int(initSoma(2)),  int(initSoma(2)) - zBoundMin ) + 5;
    int xMin = std::max(0, NGUtility::Round(initSoma(0) - xOffset));
    int yMin = std::max(0, NGUtility::Round(initSoma(1) - yOffset));
    int zMin = std::max(0, NGUtility::Round(initSoma(2) - zOffset));
    int xMax = std::min(nx - 1, NGUtility::Round(initSoma(0) + xOffset));
    int yMax = std::min(ny - 1, NGUtility::Round(initSoma(1) + yOffset));
    int zMax = std::min(nz - 1, NGUtility::Round(initSoma(2) + zOffset));
    int globalOffsetX = xMin;
    int globalOffsetY = yMin;
    int globalOffsetZ = zMin;
    //extract local binary image
    CVolume locBinImg;
    locBinImg.SetSize(xMax - xMin + 1, yMax - yMin + 1, zMax - zMin + 1);
    for (int i = xMin; i <= xMax; ++i ) 
        for (int j = yMin; j <= yMax; ++j) 
            for (int ij = zMin; ij <= zMax; ++ij) 
                locBinImg(i - xMin, j - yMin, ij - zMin) = binImg(i,j,ij);
    //extract local binary image point set, just for test
    /*VectorVec3i locBinImgPtSet;
    for (int i = 0; i < locBinImg.x(); ++i)
    for (int j = 0; j < locBinImg.y(); ++j)
    for(int ij = 0; ij < locBinImg.z(); ++ij)
    if( 0 < locBinImg(i,j,ij) )
    locBinImgPtSet.push_back(Vec3i(i,j,ij));*/

    //remove soma shape
    VectorVec3i locInterPtSet;
    locInterPtSet.reserve(interPtSet.size());
    Vec3i globalOffset(globalOffsetX, globalOffsetY, globalOffsetZ);
    for(size_t i = 0; i < interPtSet.size(); ++i){
        locInterPtSet.push_back(interPtSet[i] - globalOffset);
        locBinImg(locInterPtSet[i](0), locInterPtSet[i](1), locInterPtSet[i](2)) = 0;
    }

    //get shell
    VectorVec3i shellPtSet;
    int shellThick = 2;
    int ss = 0;
    for (size_t i = 0; i < interPtSet.size(); ++i) {
        int idx1 = std::max(locInterPtSet[i](0) - shellThick, 0);
        int idx2 = std::min(locInterPtSet[i](0) + shellThick, locBinImg.x() - 1);
        int idy1 = std::max(locInterPtSet[i](1) - shellThick, 0);
        int idy2 = std::min(locInterPtSet[i](1) + shellThick, locBinImg.y() - 1);
        int idz1 = std::max(locInterPtSet[i](2) - shellThick, 0);
        int idz2 = std::min(locInterPtSet[i](2) + shellThick, locBinImg.z() - 1);        
        ss = 0;
        for(int l = idx1; l<= idx2; ++l)
            for(int m = idy1; m <= idy2; ++m)
                for(int n = idz1; n <= idz2; ++n)
                    ss += locBinImg(l,m,n);
        if (ss > 0) {
            for(int l = idx1; l<= idx2; ++l)
                for(int m = idy1; m <= idy2; ++m)
                    for(int n = idz1; n <= idz2; ++n){
                        if( locBinImg(l,m,n) == 255){
                            shellPtSet.push_back(Vec3i(l,m,n) + globalOffset);
                            locBinImg(l,m,n) = 0;
                        }
                    }
        }
    }
    std::sort(shellPtSet.begin(), shellPtSet.end(), Vec3i_less());
    shellPtSet.erase(std::unique(shellPtSet.begin(), shellPtSet.end()), shellPtSet.end());
    //get connectArea point set.
    NeighborhoodPtSet(shellPtSet, connectPtSet);
}

void ProbeFilter::NeighborhoodPtSet( const VectorVec3i& ptSet, std::vector<VectorVec3i>& connectPtSet )
{
    connectPtSet.clear();
    if(ptSet.empty()) return;
    VecXi flag;
    flag.resize( ptSet.size());
    flag.setOnes();
    size_t idxx;
    VectorVec3i tmpConnectPtSet;
    while (flag.sum() > 0) {
        tmpConnectPtSet.clear();
        //get initial point
        for(int i = 0; i < flag.rows(); ++i){
            if( 1 == flag(i) ){
                idxx = i;
                break;
            }
        }
        tmpConnectPtSet.push_back(ptSet[idxx]);
        flag(idxx) = 0;
        //scooping area
        VectorVec3i newPts = tmpConnectPtSet;
        while(!newPts.empty()){
            VectorVec3i addPts;
            for (size_t i = 0; i < newPts.size(); ++i) {
                for (size_t j =0 ; j < ptSet.size(); ++j) {
                    if (flag(j) == 1 && std::abs(ptSet[j](0) - newPts[i](0)) < 2 && std::abs(ptSet[j](1) - newPts[i](1)) < 2
                        && std::abs(ptSet[j](2) - newPts[i](2)) < 2) {
                            addPts.push_back(ptSet[j]);
                            flag(j) = 0;
                    }
                }
            }
            newPts = addPts;
            std::copy(addPts.begin(), addPts.end(), std::back_inserter(tmpConnectPtSet));
        }
        connectPtSet.push_back(VectorVec3i());
        connectPtSet.back().swap(tmpConnectPtSet);
    }
}

double ProbeFilter::GetOverlapRate( const VectorVec3i& ptSet1 , const VectorVec3i& ptSet2 )
{
    size_t nx1 = ptSet1.size();
    size_t nx2 = ptSet2.size();
    int xBoundMin = 10000;
    int yBoundMin = 10000;
    int zBoundMin = 10000;
    int xBoundMax = 0;
    int yBoundMax = 0;
    int zBoundMax = 0;
    //get bounding box
    for (size_t i = 0; i < nx1; ++i) {
        //min
        if (ptSet1[i](0) < xBoundMin) xBoundMin = ptSet1[i](0); 
        if (ptSet1[i](1) < yBoundMin) yBoundMin = ptSet1[i](1);
        if (ptSet1[i](2) < zBoundMin) zBoundMin = ptSet1[i](2);
        //max
        if (ptSet1[i](0) > xBoundMax) xBoundMax = ptSet1[i](0); 
        if (ptSet1[i](1) > yBoundMax) yBoundMax = ptSet1[i](1);
        if (ptSet1[i](2) > zBoundMax) zBoundMax = ptSet1[i](2);
    }
    //offset
    int xGlobalOffSet = xBoundMin ;
    int yGlobalOffSet = yBoundMin ;
    int zGlobalOffSet = zBoundMin ;
    //create local binary image
    CVolume binImg;
    binImg.SetSize(xBoundMax - xBoundMin + 1, yBoundMax - yBoundMin + 1, zBoundMax - zBoundMin + 1);
    binImg.SetZero();
    //fill binary image
    for (size_t i = 0; i < nx1; ++i) {
        binImg(ptSet1[i](0) - xGlobalOffSet, ptSet1[i](1) - yGlobalOffSet, ptSet1[i](2) - zGlobalOffSet) = 1;
    }
    int overlapPtNum(0);
    for (size_t i = 0; i < nx2; ++i) {
        if( ptSet2[i](0) < xBoundMin || ptSet2[i](1) < yBoundMin || ptSet2[i](2) < zBoundMin || 
            ptSet2[i](0) > xBoundMax || ptSet2[i](1) > yBoundMax || ptSet2[i](2) > zBoundMax )
            continue;
        if( 1 == binImg(ptSet2[i](0) - xGlobalOffSet, ptSet2[i](1) - yGlobalOffSet, ptSet2[i](2) - zGlobalOffSet) )
            ++ overlapPtNum;
    }
    return double(overlapPtNum) / double(std::min(nx1, nx2));
}

bool ProbeFilter::WriteSwcTest( const char* path, const VectorVec3i& ptSet)
{
    FILE* fp = fopen(path, "w");
    int index = 1;
    for (size_t i = 0; i < ptSet.size(); ++i) {
        fprintf(fp, "%d 1 %d %d %d 1 -1\n", index ++, ptSet[i](0), ptSet[i](1), ptSet[i](2));
    }
    fclose(fp);
    return true;
}
