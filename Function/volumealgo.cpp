#include "volumealgo.h"
#include "../ngtypes/volume.h"
#include "NGUtility.h"
#pragma warning(disable:4996)


//void GetAreaSum(const Volume &dst, const int xMin, const int xMax, const int yMin, const int yMax, const int zMin, const int zMax, int &value)



void Get3DRegion(int &xMin, int &xMax, int &yMin, int &yMax, int &zMin, int &zMax,
                 const int xCenter, const int yCenter, const int zCenter,
                 const int xOffset, const int yOffset, const int zOffset)
{
    xMin = xCenter - xOffset;
    xMax = xCenter + xOffset;
    yMin = yCenter - yOffset;
    yMax = yCenter + yOffset;
    zMin = zCenter - zOffset;
    zMax = zCenter + zOffset;
}


void Get3DRegion(int &xMin, int &xMax, int &yMin, int &yMax, int &zMin, int &zMax,
                 const int xCenter, const int yCenter, const int zCenter,
                 const int xOffset, const int yOffset, const int zOffset,
                 const int xlower, const int xupper, const int ylower, const int yupper, const int zlower, const int zupper)
{
    xMin = std::max(xlower, xCenter - xOffset);
    xMax = std::min(xupper, xCenter + xOffset);
    yMin = std::max(ylower, yCenter - yOffset);
    yMax = std::min(yupper, yCenter + yOffset);
    zMin = std::max(zlower, zCenter - zOffset);
    zMax = std::min(zupper, zCenter + zOffset);
}

//2015-6-18
void ExtractLocalDomain(const Vec3d &initPoint, const Volume<unsigned short> &origImg,
                        Volume<unsigned short> &locOrigImg, Vec3d &locPoint)
{
    //Volumn vol(v.Width(), v.Height(), v.Depth());
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    int minX = std::max(NGUtility::Round(initPoint(0) - 30.0 ), 0);
    int maxX = std::min(NGUtility::Round(initPoint(0) + 30.0), nx - 1);
    int minY = std::max(NGUtility::Round(initPoint(1) - 30.0), 0);
    int maxY = std::min(NGUtility::Round(initPoint(1) + 30.0), ny - 1);
    int minZ = std::max(NGUtility::Round(initPoint(2) - 30.0), 0);//2015-6-18
    int maxZ = std::min(NGUtility::Round(initPoint(2) + 30.0), nz - 1);

    locPoint = Vec3d(initPoint(0) - minX, initPoint(1) - minY, initPoint(2) - minZ);

    //Volumn SubVol(maxX - minX + 1, maxY - minY + 1, maxZ - minZ + 1);
    int sx = maxX - minX + 1;
    int sy = maxY - minY + 1;
    int sz = maxZ - minZ + 1;

    locOrigImg.SetSize(sx, sy, sz);
    for (int i = minX; i <= maxX; ++i){
        for (int j = minY; j <= maxY; ++j){
            for (int ij = minZ; ij <= maxZ; ++ij)
                locOrigImg(i - minX, j - minY, ij - minZ) = origImg(i, j, ij);
        }
    }//for
}

//2015-6-18
void ExtractLocalDomain(const Vec3d &initPoint, const Volume<unsigned char> &origImg,
                        Volume<unsigned char> &locOrigImg, Vec3d &locPoint)
{
    //Volumn vol(v.Width(), v.Height(), v.Depth());
    int nx = origImg.x();
    int ny = origImg.y();
    int nz = origImg.z();
    int minX = std::max(NGUtility::Round(initPoint(0) - 30.0), 0);
    int maxX = std::min(NGUtility::Round(initPoint(0) + 30.0), nx - 1);
    int minY = std::max(NGUtility::Round(initPoint(1) - 30.0), 0);
    int maxY = std::min(NGUtility::Round(initPoint(1) + 30.0), ny - 1);
    int minZ = std::max(NGUtility::Round(initPoint(2) - 30.0), 0);//2015-6-18
    int maxZ = std::min(NGUtility::Round(initPoint(2) + 30.0), nz - 1);

    locPoint = Vec3d(initPoint(0) - minX, initPoint(1) - minY, initPoint(2) - minZ);

    //Volumn SubVol(maxX - minX + 1, maxY - minY + 1, maxZ - minZ + 1);
    int sx = maxX - minX + 1;
    int sy = maxY - minY + 1;
    int sz = maxZ - minZ + 1;

    locOrigImg.SetSize(sx, sy, sz);
    for (int i = minX; i <= maxX; ++i){
        for (int j = minY; j <= maxY; ++j){
            for (int ij = minZ; ij <= maxZ; ++ij)
                locOrigImg(i - minX, j - minY, ij - minZ) = origImg(i, j, ij);
        }
    }//for
}

//2015-6-8
int AbsDiff(const int a, const int b)
{
    return std::abs(a - b);
}

//2015-10-9
void LimitRange( int xMin, int xMax, int yMin, int yMax, int zMin, int zMax, int &dstXMin, int &dstXMax, int &dstYMin, int &dstYMax, int &dstZMin, int &dstZMax )
{
	dstXMin = std::max(xMin, dstXMin);
	dstXMax = std::min(xMax, dstXMax);
	dstYMin = std::max(yMin, dstYMin);
	dstYMax = std::min(yMax, dstYMax);
	dstZMin = std::max(zMin, dstZMin);
	dstZMax = std::min(zMax, dstZMax);
}


void GenerateSerialVector( double start, double finish, double step, std::vector<double>& dst)
{
    dst.clear();
    for (double i = start; i <= finish; i+= step) {
        dst.push_back(i);
    }
}


