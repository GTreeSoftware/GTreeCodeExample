
#include "contourutil.h"
//#include "Core/algo_stl-inl.h"
#include "../ngtypes/volume.h"
//#include "common.h"
#include <deque>
#include <algorithm>
#include <numeric>
#include <Eigen/Geometry>
#include "Function/volumealgo.h"
#include "NGUtility.h"
using namespace std;

/**/
void ContourUtil::CalAreaVolume(const VectorVec3d &dataPP, const int Theta, const int Phi,
                     VectorVec3d &Patch, double &area, double &cellVol)
{
    double ratio(2.0);//
    int count(0);
    area = 0.0f;
    cellVol = 0.0f;
    Patch.clear();
    Vec3d center(dataPP[0](0), dataPP[0](1), dataPP[0](2)*ratio);//
    Vec3d A;
    Vec3d B;
    Vec3d C;
    Vec3d AB;
    Vec3d AC;
    Vec3d OA;
    Vec3d OB;
    Vec3d OC;

    /**/
    for (int i=0; i < Theta - 1; ++i)
    {
        Patch.push_back(dataPP[1]);//A
        Patch.push_back(dataPP[i + 2]);//B
        Patch.push_back(dataPP[i + 3]);//C
        /*ABC */
        A(0) = dataPP[1](0); A(1) = dataPP[1](1); A(2) = dataPP[1](2)*ratio;
        B(0) = dataPP[i + 2](0); B(1) = dataPP[i + 2](1); B(2) = dataPP[i + 2](2)*ratio;
        C(0) = dataPP[i + 3](0); C(1) = dataPP[i + 3](1); C(2) = dataPP[i + 3](2)*ratio;
        AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
        ++count;
    }
    //i = 40
    Patch.push_back(dataPP[1]);
    Patch.push_back(dataPP[Theta + 1]);
    Patch.push_back(dataPP[2]);
    ///*ABC */
    A(0) = dataPP[1](0); A(1) = dataPP[1](1); A(2) = dataPP[1](2)*ratio;
    B(0) = dataPP[Theta + 1](0); B(1) = dataPP[Theta + 1](1); B(2) = dataPP[Theta + 1](2)*ratio;
    C(0) = dataPP[2](0); C(1) = dataPP[2](1); C(2) = dataPP[2](2)*ratio;
    AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
    ++count;

    /**/
    for (int i = 0; i < Phi - 2; ++i)//
    {
        int j;
        for (j = 2 + i * Theta; j < 1 + (i + 1) * Theta; ++j)//
        {
            /**/
            Patch.push_back(dataPP[j]);
            Patch.push_back(dataPP[j + Theta]);
            Patch.push_back(dataPP[j + Theta + 1]);//anti clock wise
            /*ABC */
            A(0) = dataPP[j](0); A(1) = dataPP[j](1); A(2) = dataPP[j](2)*ratio;
            B(0) = dataPP[j + Theta](0); B(1) = dataPP[j + Theta](1); B(2) = dataPP[j + Theta](2)*ratio;
            C(0) = dataPP[j + Theta + 1](0); C(1) = dataPP[j + Theta + 1](1); C(2) = dataPP[j + Theta + 1](2)*ratio;
            AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
            ++count;

            /**/
            Patch.push_back(dataPP[j + 1]);
            Patch.push_back(dataPP[j]);
            Patch.push_back(dataPP[j + Theta + 1]);
            /*ABC */
            A(0) = dataPP[j+1](0); A(1) = dataPP[j+1](1); A(2) = dataPP[j+1](2)*ratio;
            B(0) = dataPP[j](0); B(1) = dataPP[j](1); B(2) = dataPP[j](2)*ratio;
            C(0) = dataPP[j + Theta + 1](0); C(1) = dataPP[j + Theta + 1](1); C(2) = dataPP[j + Theta + 1](2)*ratio;
            AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
            ++count;
        }

        /*j = 41*/
        Patch.push_back(dataPP[j]);//41
        Patch.push_back(dataPP[j + Theta]);//81
        Patch.push_back(dataPP[j + 1]);//42
        /*ABC */
        A(0) = dataPP[j](0); A(1) = dataPP[j](1); A(2) = dataPP[j](2)*ratio;
        B(0) = dataPP[j + Theta](0); B(1) = dataPP[j + Theta](1); B(2) = dataPP[j + Theta](2)*ratio;
        C(0) = dataPP[j + 1](0); C(1) = dataPP[j + 1](1); C(2) = dataPP[j + 1](2)*ratio;
        AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
        ++count;

        Patch.push_back(dataPP[j+1 - Theta]);//2
        Patch.push_back(dataPP[j]);//41
        Patch.push_back(dataPP[j + 1]);//42
        /*ABC */
        A(0) = dataPP[j+1 - Theta](0); A(1) = dataPP[j+1 - Theta](1); A(2) = dataPP[j+1 - Theta](2)*ratio;
        B(0) = dataPP[j](0); B(1) = dataPP[j](1); B(2) = dataPP[j](2)*ratio;
        C(0) = dataPP[j + 1](0); C(1) = dataPP[j + 1](1); C(2) = dataPP[j + 1](2)*ratio;
        AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
        ++count;
    }

    int i=0;
    for (; i < Theta-1 ; ++i)
    {
        Patch.push_back(dataPP[dataPP.size() - 1]);
        Patch.push_back(dataPP[dataPP.size() - i - 2]);
        Patch.push_back(dataPP[dataPP.size() - i - 3]);

        A(0) = dataPP[dataPP.size() - 1](0); A(1) = dataPP[dataPP.size() - 1](1); A(2) = dataPP[dataPP.size() - 1](2)*ratio;
        B(0) = dataPP[dataPP.size() - i - 2](0); B(1) = dataPP[dataPP.size() - i - 2](1); B(2) = dataPP[dataPP.size() - i - 2](2)*ratio;
        C(0) = dataPP[dataPP.size() - i - 3](0); C(1) = dataPP[dataPP.size() - i - 3](1); C(2) = dataPP[dataPP.size() - i - 3](2)*ratio;
        AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
        ++count;
    }
    //i = 40
    Patch.push_back(dataPP[dataPP.size() - 1]);
    Patch.push_back(dataPP[dataPP.size() - i - 2]);
    Patch.push_back(dataPP[dataPP.size() - 2]);
    /*求ABC */
    A(0) = dataPP[dataPP.size() - 1](0); A(1) = dataPP[dataPP.size() - 1](1); A(2) = dataPP[dataPP.size() - 1](2)*ratio;
    B(0) = dataPP[dataPP.size() - i - 2](0); B(1) = dataPP[dataPP.size() - i - 2](1); B(2) = dataPP[dataPP.size() - i - 2](2)*ratio;
    C(0) = dataPP[dataPP.size()  - 2](0); C(1) = dataPP[dataPP.size() - 2](1); C(2) = dataPP[dataPP.size() - 2](2)*ratio;
    AddAreaVolume(center, A, B, C, AB, AC, OA, OB, OC, area, cellVol);//
    ++count;

    printf("%d\n", count);
}

/**/
void ContourUtil::AddAreaVolume(const Vec3d &center, const Vec3d &A, const Vec3d &B, const Vec3d &C,
                    Vec3d &AB, Vec3d &AC, Vec3d &OA, Vec3d &OB, Vec3d &OC,
                    double &area, double &cellVol)
{
    //
    /*Area*/
    AB = B - A;//AB
    AC = C - A;//AB
    area += 0.5f * (std::abs)((AB.cross(AC)).norm());
    /*Volume*/
    OA = A - center;
    OB = B - center;
    OC = C - center;
    cellVol += 1.0f / 6.0f * (std::abs)((OA.cross(OB)).dot(OC));
}


/**/
void ContourUtil::GetRayLimit(const DVolume &sphereRayWet,
                 const double constrictionThrev,
                 std::vector<std::vector<double> > &rayLimit)
{
    int nx = sphereRayWet.x();
    int ny = sphereRayWet.y();
    int nz = sphereRayWet.z();
    vector<double> sphereRay;
    vector<double> tmpUzz;
    sphereRay.clear();
    int arra(0);

    for (int i = 0; i < ny; ++i){
        tmpUzz.clear();
        for (int j = 0; j < nz; ++j){
            for (int ij = 0; ij < nx; ++ij)
                sphereRay.push_back(sphereRayWet(ij, i, j));

            CalculateOneRayLimit(sphereRay, constrictionThrev, arra);
            tmpUzz.push_back(arra);
            sphereRay.clear();
        }
        rayLimit.push_back(tmpUzz);
    }
}
/**/
void ContourUtil::CalculateOneRayLimit(const vector<double> &ray,
                          const double constrictionThrev,
                          int &oneRayLimit)
{
    size_t nx = ray.size();
    oneRayLimit = 0;

    for (size_t i = 0; i < nx; ++i){
        if (ray[i] < constrictionThrev) break;
        else  ++oneRayLimit;
    }
}

void ContourUtil::CalculateSphereOneNode(const CVolume &locOrigImg,

                            const double value, 
                            const double x, 
                            const double y, 
                            const double z, 
                            double &segmentDense,
                            double& segmentWet)
{
    int subx = locOrigImg.x();
    int suby = locOrigImg.y();
    int subz = locOrigImg.z();
    double rayNodeDist(0.0);

    double idexx = double(std::max(std::min(NGUtility::Round(x), subx - 1), 0));
    double idexy = double(std::max(std::min(NGUtility::Round(y), suby - 1), 0));
    double idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist= value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x + 1.0), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x - 1.0), subx - 1), 0));//-1.0 + 0.5
    idexy = double(std::max(std::min(NGUtility::Round(y), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y - 1.0), suby - 1), 0));//-1.0 + 0.5
    idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y + 1.0), suby - 1), 0));//+1.0 + 0.5
    idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z - 1.0), subz - 1), 0));//-1.0 + 0.5

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z + 1.0), subz - 1), 0));//+1.0 + 0.5

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);
}

//2015-6-8
void ContourUtil::CalculateSphereOneNode(const SVolume &locOrigImg,

                            const double value,
                            const double x,
                            const double y,
                            const double z,
                            double &segmentDense,
                            double& segmentWet)
{
    int subx = locOrigImg.x();
    int suby = locOrigImg.y();
    int subz = locOrigImg.z();
    double rayNodeDist(0.0);

    double idexx = double(std::max(std::min(NGUtility::Round(x), subx - 1), 0));
    double idexy = double(std::max(std::min(NGUtility::Round(y), suby - 1), 0));
    double idexz = double(std::max(std::min(NGUtility::Round(z), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist= value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x +1.0), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z ), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x - 1.0), subx - 1), 0));//-1.0 + 0.5
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z ), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y - 1.0), suby - 1), 0));//-1.0 + 0.5
    idexz = double(std::max(std::min(NGUtility::Round(z ), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y + 1.0), suby - 1), 0));//+1.0 + 0.5
    idexz = double(std::max(std::min(NGUtility::Round(z ), subz - 1), 0));

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z - 1.0), subz - 1), 0));//-1.0 + 0.5

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);

    idexx = double(std::max(std::min(NGUtility::Round(x ), subx - 1), 0));
    idexy = double(std::max(std::min(NGUtility::Round(y ), suby - 1), 0));
    idexz = double(std::max(std::min(NGUtility::Round(z + 1.0), subz - 1), 0));//+1.0 + 0.5

    rayNodeDist = -2.0 * ( (x - idexx) * (x - idexx) + (y - idexy) * (y - idexy) +
        4.0 * (z - idexz) * (z - idexz) );
    rayNodeDist=value * std::max(rayNodeDist,-6.0);
    segmentDense += std::exp(rayNodeDist) * (double)locOrigImg(int(idexx), int(idexy), int(idexz));
    segmentWet += std::exp(rayNodeDist);
}

void ContourUtil::GetGradientVectorFlow(const DVolume &sphereRayWet, DVolume &smoothRay)
{
    smoothRay.SetSize(sphereRayWet.x(), sphereRayWet.y(), sphereRayWet.z());

    std::vector<double> smoothOneRayWet;

    for (int i = 0; i < sphereRayWet.y(); ++i)//1
    {
        for (int j = 0; j < sphereRayWet.z(); ++j)//2
        {
            std::vector<double> tmpSphere;
            smoothOneRayWet.clear();
            for (int ij = 0; ij < sphereRayWet.x(); ++ij)//
                tmpSphere.push_back((double)sphereRayWet(ij,i,j));

            SmoothGradientCurves(tmpSphere, smoothOneRayWet);

            int numSmoothOneRayWet =int( smoothOneRayWet.size());
            for (int ij = 0; ij < numSmoothOneRayWet; ++ij)
                smoothRay(ij,i,j) = smoothOneRayWet[ij];
        }
    }
}



void ContourUtil::SmoothGradientCurves(const std::vector<double> &oneRayWet, std::vector<double> &smoothOneRayWet)
{
    smoothOneRayWet.clear();

    std::deque<double> tmpDiff;//原始梯度
    //diff_vector<double>(oneRayWet, tmp_diff);
    std::adjacent_difference(oneRayWet.begin(), oneRayWet.end(), std::back_inserter(tmpDiff));
    tmpDiff.pop_back();
    std::vector<double> diffRay;

    //NG_Util::Interpl_2_Mean(tmpDiff, diffRay);//插值
    for (std::vector<double>::size_type i = 0; i < tmpDiff.size() - 1; ++i)
    {
        diffRay.push_back(tmpDiff[i]);
        diffRay.push_back((tmpDiff[i+1] + tmpDiff[i]) / 2.0);
    }
    diffRay.push_back(tmpDiff[tmpDiff.size() - 1]);
    int diffRayLength = int(diffRay.size());

    std::vector<double> yy = diffRay;
    std::vector<double> yy1 = yy;

    /*平滑梯度*///只是循环100次
    for (int ddk = 1; ddk < 101; ++ddk){
        for (int i =1; i < diffRayLength - 1; ++i){
            yy1[i] = ((std::abs)(yy[i]) * yy[i] + 2.0 * (yy[i-1]+yy[i+1])) / ((std::abs)(yy[i]) + 4.0);
            yy = yy1;
        }
    }

    std::vector<size_t> tts;
    //generate_n(back_inserter(tts), (int)(yy.size() / 2) + 1,GenArray<int>(0, 2));
    for(std::vector<double>::size_type i = 0; i < yy.size(); i+=2){
        tts.push_back(i);
    }

    size_t num_tts = tts.size();
    for (size_t i = 0; i < num_tts; ++i){
        smoothOneRayWet.push_back(yy[tts[i]]);
    }

    size_t numSmoothOneRayWet = smoothOneRayWet.size();
    for (size_t i = 0; i < numSmoothOneRayWet; ++i){
        smoothOneRayWet[i] = smoothOneRayWet[i] > 0.0? 0.0 : (-smoothOneRayWet[i]);
    }
}

void ContourUtil::GetBoundaryBack(const std::vector<double> &outerShell,
                     const double threv, 
                     std::vector<double> &boundaryBack)
{
    int nxx = int(outerShell.size());
    boundaryBack.clear();
    std::vector<double> adjustShell;
    double tmp;
    double third = 1.0/3.0;
    for (int i = 0; i < nxx; ++i){
        //2014-4-19
        tmp = outerShell[i] < 400.0 ? outerShell[i] : 400.0;//400.0
        adjustShell.push_back(tmp);
    }

    for (int i = 0; i < 100; ++i){
        for (int j = 1; j < nxx - 1; ++j)
            adjustShell[j] = third * (adjustShell[j - 1] + adjustShell[j] + adjustShell[j + 1]);

        adjustShell[nxx - 1] = 0.5 * (adjustShell[nxx - 1] + adjustShell[nxx - 2]);
        adjustShell[0] = 0.5 * (adjustShell[0] + adjustShell[1]);
    }
    //boundaryBack.resize(outerShell.size());
    boundaryBack = outerShell;
    for (int i = 0; i < nxx; ++i){
        if (std::abs(outerShell[i]-adjustShell[i]) > threv * std::sqrt(adjustShell[i]))
            boundaryBack[i] = adjustShell[i];
    }
}
